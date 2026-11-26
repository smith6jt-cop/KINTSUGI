"""The remote (streamable HTTP) transport: authentication, host checks, confinement.

The app under test is the one ``kintsugi mcp start --transport streamable-http``
serves (:func:`kintsugi.mcp.server.create_http_app`), driven in-process through
Starlette's test client: bearer-token access, the full OAuth 2.1 flow claude.ai
uses (dynamic registration, /login with the secret, PKCE token exchange),
DNS-rebinding protection and the confinement of ``project_path`` to the served
project.
"""

from __future__ import annotations

import base64
import hashlib
import json
import secrets
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

from kintsugi.mcp._sdk import MCP_AVAILABLE

pytestmark = pytest.mark.skipif(not MCP_AVAILABLE, reason="mcp>=2.0 not installed")

SECRET = "s" * 40
CLAUDE_CALLBACK = "https://claude.ai/api/mcp/auth_callback"
JSON_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json, text/event-stream",
}
INIT = {
    "jsonrpc": "2.0",
    "id": 1,
    "method": "initialize",
    "params": {
        "protocolVersion": "2025-06-18",
        "capabilities": {},
        "clientInfo": {"name": "pytest", "version": "0"},
    },
}


def _project(path: Path) -> Path:
    from kintsugi.mcp.session import PROJECT_MARKER

    path.mkdir(parents=True, exist_ok=True)
    (path / PROJECT_MARKER).write_text(json.dumps({"config": {"name": path.name}}))
    return path


@pytest.fixture
def project(tmp_path):
    from kintsugi.mcp.session import set_project_root

    root = _project(tmp_path / "served")
    set_project_root(root)
    return root


@pytest.fixture
def client(project):
    from starlette.testclient import TestClient

    from kintsugi.mcp.server import HttpConfig, create_http_app

    app, _ = create_http_app(HttpConfig(secret=SECRET, host="127.0.0.1", port=8000))
    with TestClient(app, base_url="http://127.0.0.1:8000") as c:
        yield c


def _rpc_payloads(response) -> list[dict]:
    """JSON-RPC messages from a JSON or SSE response body."""
    if response.headers.get("content-type", "").startswith("application/json"):
        body = response.json()
        return body if isinstance(body, list) else [body]
    out = []
    for line in response.text.splitlines():
        if line.startswith("data:") and line[5:].strip():
            out.append(json.loads(line[5:]))
    return out


def _session(client, token: str) -> dict[str, str]:
    headers = {**JSON_HEADERS, "Authorization": f"Bearer {token}"}
    r = client.post("/mcp", json=INIT, headers=headers)
    assert r.status_code == 200, r.text
    session_id = r.headers.get("mcp-session-id")
    if session_id:
        headers["mcp-session-id"] = session_id
    proto = _rpc_payloads(r)[0]["result"]["protocolVersion"]
    headers["mcp-protocol-version"] = proto
    client.post(
        "/mcp", json={"jsonrpc": "2.0", "method": "notifications/initialized"}, headers=headers
    )
    return headers


def _call(client, headers, method: str, params: dict | None = None, rid: int = 2) -> dict:
    r = client.post(
        "/mcp",
        json={"jsonrpc": "2.0", "id": rid, "method": method, "params": params or {}},
        headers=headers,
    )
    assert r.status_code == 200, r.text
    return _rpc_payloads(r)[-1]


# ---------------------------------------------------------------------------
# Bearer token
# ---------------------------------------------------------------------------


def test_requests_without_a_token_are_refused(client):
    r = client.post("/mcp", json=INIT, headers=JSON_HEADERS)
    assert r.status_code == 401
    assert "resource_metadata" in r.headers.get("www-authenticate", "")


def test_a_wrong_token_is_refused(client):
    r = client.post("/mcp", json=INIT, headers={**JSON_HEADERS, "Authorization": "Bearer nope"})
    assert r.status_code == 401


def test_the_static_token_lists_every_tool(client):
    from kintsugi.mcp.tool_specs import TOOL_SPECS

    headers = _session(client, SECRET)
    result = _call(client, headers, "tools/list")["result"]
    assert {t["name"] for t in result["tools"]} == {s["name"] for s in TOOL_SPECS}


def test_an_unknown_host_header_is_refused(client):
    headers = {**JSON_HEADERS, "Authorization": f"Bearer {SECRET}", "Host": "evil.example"}
    r = client.post("/mcp", json=INIT, headers=headers)
    assert r.status_code == 421


def test_health_endpoint_needs_no_token(client):
    assert client.get("/healthz").json() == {"status": "ok"}


# ---------------------------------------------------------------------------
# Confinement
# ---------------------------------------------------------------------------


def test_project_path_outside_the_served_project_is_refused(client, tmp_path):
    other = _project(tmp_path / "other")
    headers = _session(client, SECRET)
    reply = _call(
        client,
        headers,
        "tools/call",
        {"name": "list_channels", "arguments": {"project_path": str(other)}},
    )
    result = reply["result"]
    assert result["isError"] is True
    assert "outside the project this server serves" in result["content"][0]["text"]


def test_project_path_inside_the_served_project_is_allowed(client, project):
    headers = _session(client, SECRET)
    reply = _call(
        client,
        headers,
        "tools/call",
        {"name": "list_channels", "arguments": {"project_path": str(project)}},
    )
    assert "outside the project" not in json.dumps(reply)


def test_check_project_path_rejects_symlink_escape(project, tmp_path):
    from kintsugi.mcp.session import (
        ProjectOutsideRoot,
        check_project_path,
        resolve_project_path,
        set_confined,
    )

    outside = _project(tmp_path / "outside")
    link = project / "link"
    link.symlink_to(outside, target_is_directory=True)
    check_project_path(link)  # not confined: anything goes
    set_confined(True)
    with pytest.raises(ProjectOutsideRoot):
        check_project_path(link)
    with pytest.raises(ProjectOutsideRoot):
        resolve_project_path(str(project / ".." / "outside"))
    assert resolve_project_path(None) == str(project.resolve())


# ---------------------------------------------------------------------------
# OAuth (claude.ai custom connectors)
# ---------------------------------------------------------------------------


def _pkce() -> tuple[str, str]:
    verifier = secrets.token_urlsafe(48)
    challenge = base64.urlsafe_b64encode(hashlib.sha256(verifier.encode()).digest())
    return verifier, challenge.rstrip(b"=").decode()


def _register(client, redirect: str = CLAUDE_CALLBACK):
    return client.post(
        "/register",
        json={
            "redirect_uris": [redirect],
            "client_name": "pytest connector",
            "grant_types": ["authorization_code", "refresh_token"],
            "response_types": ["code"],
            "token_endpoint_auth_method": "none",
        },
    )


def _authorize(client, client_id: str, challenge: str) -> str:
    r = client.get(
        "/authorize",
        params={
            "response_type": "code",
            "client_id": client_id,
            "redirect_uri": CLAUDE_CALLBACK,
            "code_challenge": challenge,
            "code_challenge_method": "S256",
            "state": "xyz",
            "scope": "kintsugi",
            "resource": "http://127.0.0.1:8000/mcp",
        },
        follow_redirects=False,
    )
    assert r.status_code in (302, 307), r.text
    login = urlparse(r.headers["location"])
    assert login.path == "/login"
    return parse_qs(login.query)["request"][0]


def test_oauth_metadata_is_published(client):
    meta = client.get("/.well-known/oauth-authorization-server").json()
    assert meta["issuer"].rstrip("/") == "http://127.0.0.1:8000"
    assert "S256" in meta["code_challenge_methods_supported"]
    resource = client.get("/.well-known/oauth-protected-resource/mcp").json()
    assert resource["resource"] == "http://127.0.0.1:8000/mcp"


def test_registration_refuses_unlisted_redirects(client):
    r = _register(client, "https://attacker.example/callback")
    assert r.status_code == 400
    assert r.json()["error"] == "invalid_redirect_uri"


def test_full_oauth_flow_yields_a_working_token(client):
    reg = _register(client)
    assert reg.status_code == 201, reg.text
    client_id = reg.json()["client_id"]
    verifier, challenge = _pkce()
    request_id = _authorize(client, client_id, challenge)

    page = client.get("/login", params={"request": request_id})
    assert page.status_code == 200 and "KINTSUGI_MCP_TOKEN" in page.text

    bad = client.post("/login", data={"request": request_id, "secret": "wrong"})
    assert bad.status_code == 401

    ok = client.post(
        "/login", data={"request": request_id, "secret": SECRET}, follow_redirects=False
    )
    assert ok.status_code == 302
    target = urlparse(ok.headers["location"])
    assert f"{target.scheme}://{target.netloc}{target.path}" == CLAUDE_CALLBACK
    query = parse_qs(target.query)
    assert query["state"] == ["xyz"]

    token = client.post(
        "/token",
        data={
            "grant_type": "authorization_code",
            "code": query["code"][0],
            "redirect_uri": CLAUDE_CALLBACK,
            "client_id": client_id,
            "code_verifier": verifier,
            "resource": "http://127.0.0.1:8000/mcp",
        },
    )
    assert token.status_code == 200, token.text
    access = token.json()["access_token"]

    headers = _session(client, access)
    assert _call(client, headers, "tools/list")["result"]["tools"]

    replay = client.post(
        "/token",
        data={
            "grant_type": "authorization_code",
            "code": query["code"][0],
            "redirect_uri": CLAUDE_CALLBACK,
            "client_id": client_id,
            "code_verifier": verifier,
        },
    )
    assert replay.status_code == 400


def test_login_locks_out_after_repeated_failures(client):
    from kintsugi.mcp.remote_auth import MAX_FAILURES

    client_id = _register(client).json()["client_id"]
    _, challenge = _pkce()
    request_id = _authorize(client, client_id, challenge)
    for _ in range(MAX_FAILURES):
        client.post("/login", data={"request": request_id, "secret": "wrong"})
    locked = client.post("/login", data={"request": request_id, "secret": SECRET})
    assert locked.status_code == 401
    assert "Too many failed attempts" in locked.text


def test_unknown_login_request_is_rejected(client):
    assert client.get("/login", params={"request": "missing"}).status_code == 400


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


def test_short_or_missing_secret_is_refused():
    from kintsugi.mcp.remote_auth import SecretError, validate_secret

    with pytest.raises(SecretError):
        validate_secret(None)
    with pytest.raises(SecretError):
        validate_secret("short")
    assert validate_secret(SECRET) == SECRET


def test_redirect_allowlist():
    from kintsugi.mcp.remote_auth import redirect_allowed

    assert redirect_allowed(CLAUDE_CALLBACK)
    assert redirect_allowed("http://localhost:33418/callback")
    assert redirect_allowed("http://127.0.0.1:9/cb")
    assert not redirect_allowed("https://localhost.evil.example/cb")
    assert not redirect_allowed("https://claude.ai/other")
    assert redirect_allowed("https://x.example/cb", ("https://x.example/cb",))


@pytest.mark.parametrize(
    ("args", "env", "message"),
    [
        ([], {"KINTSUGI_MCP_TOKEN": SECRET}, "requires --project-dir"),
        (["--project-dir", "{p}"], {}, "KINTSUGI_MCP_TOKEN is not set"),
        (
            ["--project-dir", "{p}", "--host", "0.0.0.0"],
            {"KINTSUGI_MCP_TOKEN": SECRET},
            "--public-url",
        ),
        (
            ["--project-dir", "{p}", "--public-url", "http://example.org"],
            {"KINTSUGI_MCP_TOKEN": SECRET},
            "must use https",
        ),
    ],
)
def test_cli_refuses_unsafe_http_configurations(tmp_path, monkeypatch, args, env, message):
    from click.testing import CliRunner

    from kintsugi.cli import main

    monkeypatch.delenv("KINTSUGI_MCP_TOKEN", raising=False)
    monkeypatch.delenv("KINTSUGI_MCP_PUBLIC_URL", raising=False)
    proj = _project(tmp_path / "p")
    argv = ["mcp", "start", "--transport", "streamable-http"]
    argv += [a.format(p=proj) for a in args]
    result = CliRunner().invoke(main, argv, env=env)
    assert result.exit_code == 2
    assert message in " ".join(result.output.split())
