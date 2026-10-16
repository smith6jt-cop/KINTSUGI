"""The agentic eval's provider adapters (evals/mcp/providers.py), without network access."""

from __future__ import annotations

import io
import json
import sys
import threading
import urllib.error
from email.message import Message
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pytest

EVALS = Path(__file__).resolve().parents[1] / "evals" / "mcp"
sys.path.insert(0, str(EVALS))
try:
    import providers
finally:
    sys.path.remove(str(EVALS))

KEY = "sk-test-0123456789"
TOOLS = [
    {
        "name": "denoise",
        "description": "Remove noise.",
        "input_schema": {
            "type": "object",
            "title": "denoiseArguments",
            "properties": {
                "channel": {"type": "string", "title": "Channel", "description": "Image."},
                "method": {
                    "type": "string",
                    "enum": ["median", "uniform"],
                    "default": "median",
                    "title": "Method",
                },
                "filter_size": {
                    "anyOf": [{"type": "integer", "minimum": 1, "maximum": 15}, {"type": "null"}],
                    "default": None,
                    "description": "Kernel size.",
                    "title": "Filter Size",
                },
                "title": {"type": "string", "format": "uri", "title": "Title"},
            },
            "required": ["channel"],
            "additionalProperties": False,
        },
    }
]


class _Response(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


@pytest.fixture
def captured(monkeypatch):
    """Replace urlopen; record requests and answer with the queued responses."""
    calls: list[dict] = []
    queue: list = []

    def fake_urlopen(request, timeout=None):
        calls.append(
            {
                "url": request.full_url,
                "headers": {k.lower(): v for k, v in request.header_items()},
                "body": json.loads(request.data),
                "timeout": timeout,
            }
        )
        item = queue.pop(0)
        if isinstance(item, Exception):
            raise item
        return _Response(json.dumps(item).encode())

    monkeypatch.setattr(providers.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(providers.time, "sleep", lambda s: None)
    for name in (*providers.API_KEY_ENV.values(), *providers.BASE_URL_ENV.values()):
        monkeypatch.delenv(name, raising=False)
    return calls, queue


def _http_error(code: int, body: dict | str, retry_after: str | None = None):
    headers = Message()
    if retry_after is not None:
        headers["retry-after"] = retry_after
    text = body if isinstance(body, str) else json.dumps(body)
    return urllib.error.HTTPError("u", code, "err", headers, io.BytesIO(text.encode()))


def test_anthropic_request_and_tool_choice(captured):
    calls, queue = captured
    queue.append(
        {"content": [{"type": "text", "text": "ok"}, {"type": "tool_use", "name": "denoise"}]}
    )
    chosen = providers.select_tool("anthropic", "m-1", "clean CD3", TOOLS, api_key=KEY)
    assert chosen == "denoise"
    (call,) = calls
    assert call["url"] == "https://api.anthropic.com/v1/messages"
    assert call["headers"]["x-api-key"] == KEY
    assert call["headers"]["anthropic-version"] == "2023-06-01"
    assert call["headers"]["content-type"] == "application/json"
    body = call["body"]
    assert body["model"] == "m-1" and body["max_tokens"] == providers.MAX_TOKENS
    assert body["tool_choice"] == {"type": "auto"}
    assert body["messages"] == [{"role": "user", "content": "clean CD3"}]
    # The Messages API receives the advertised schema unchanged.
    assert body["tools"] == [
        {
            "name": "denoise",
            "description": "Remove noise.",
            "input_schema": TOOLS[0]["input_schema"],
        }
    ]


def test_openai_request_and_tool_choice(captured, monkeypatch):
    calls, queue = captured
    monkeypatch.setenv("OPENAI_API_KEY", KEY)
    queue.append({"choices": [{"message": {"tool_calls": [{"function": {"name": "denoise"}}]}}]})
    assert providers.select_tool("openai", "m-2", "clean CD3", TOOLS) == "denoise"
    (call,) = calls
    assert call["url"] == "https://api.openai.com/v1/chat/completions"
    assert call["headers"]["authorization"] == f"Bearer {KEY}"
    body = call["body"]
    assert body["model"] == "m-2" and body["tool_choice"] == "auto"
    assert "max_tokens" not in body and "max_completion_tokens" not in body
    (tool,) = body["tools"]
    assert tool["type"] == "function"
    assert tool["function"]["name"] == "denoise"
    assert tool["function"]["parameters"] == providers.sanitize_schema(TOOLS[0]["input_schema"])


@pytest.mark.parametrize(
    ("provider", "response"),
    [
        ("anthropic", {"content": [{"type": "text", "text": "Which channel?"}]}),
        ("openai", {"choices": [{"message": {"content": "Which channel?", "tool_calls": None}}]}),
    ],
)
def test_a_text_answer_selects_no_tool(captured, provider, response):
    _, queue = captured
    queue.append(response)
    assert providers.select_tool(provider, "m", "p", TOOLS, api_key=KEY) is None


def test_sanitize_schema():
    clean = providers.sanitize_schema(TOOLS[0]["input_schema"])
    text = json.dumps(clean)
    assert '"default"' not in text and '"format"' not in text and "anyOf" not in text
    assert "title" not in clean and "title" not in clean["properties"]["channel"]
    # A parameter called "title" is a name, not the keyword: it stays.
    assert clean["properties"]["title"] == {"type": "string"}
    assert clean["properties"]["filter_size"] == {
        "type": "integer",
        "minimum": 1,
        "maximum": 15,
        "description": "Kernel size.",
    }
    assert clean["properties"]["method"]["enum"] == ["median", "uniform"]
    assert clean["required"] == ["channel"] and clean["additionalProperties"] is False
    # The input is not modified.
    assert TOOLS[0]["input_schema"]["properties"]["filter_size"]["anyOf"]


def test_sanitize_keeps_a_real_union():
    schema = {"anyOf": [{"type": "string"}, {"type": "integer"}, {"type": "null"}]}
    assert providers.sanitize_schema(schema) == schema


def test_sanitize_every_advertised_schema():
    pytest.importorskip("mcp")
    from kintsugi.mcp._sdk import MCP_AVAILABLE

    if not MCP_AVAILABLE:
        pytest.skip("mcp 2.x not installed")
    sys.path.insert(0, str(EVALS))
    try:
        import harness
    finally:
        sys.path.remove(str(EVALS))

    def keywords(node, found):
        if isinstance(node, dict):
            for key, value in node.items():
                if key == "properties":
                    for sub in value.values():
                        keywords(sub, found)
                else:
                    found.add(key)
                    keywords(value, found)
        elif isinstance(node, list):
            for item in node:
                keywords(item, found)
        return found

    for tool in harness.load_tools():
        found = keywords(providers.sanitize_schema(tool["input_schema"]), set())
        assert not found & {"title", "default", "format"}, tool["name"]
        assert "anyOf" not in found, tool["name"]  # every union in the surface is T | None


def test_base_url_and_local_server_without_key(captured, monkeypatch):
    calls, queue = captured
    monkeypatch.setenv("OPENAI_BASE_URL", "http://localhost:11434/v1/")
    queue.append({"choices": [{"message": {"tool_calls": [{"function": {"name": "denoise"}}]}}]})
    assert providers.select_tool("openai", "llama", "p", TOOLS) == "denoise"
    assert calls[0]["url"] == "http://localhost:11434/v1/chat/completions"
    assert "authorization" not in calls[0]["headers"]
    assert providers.endpoint_host("openai") == "localhost:11434"
    assert providers.endpoint_host("openai", "https://u:p@example.org/v1") == "example.org"


@pytest.mark.parametrize("provider", ["anthropic", "openai"])
def test_missing_key_for_a_hosted_api(captured, provider):
    with pytest.raises(providers.ProviderError, match=providers.API_KEY_ENV[provider]):
        providers.select_tool(provider, "m", "p", TOOLS)
    assert captured[0] == []


def test_http_error_message_never_contains_the_key(captured):
    _, queue = captured
    queue.append(_http_error(401, {"error": {"message": f"invalid x-api-key {KEY}"}}))
    with pytest.raises(providers.ProviderError) as info:
        providers.select_tool("anthropic", "m", "p", TOOLS, api_key=KEY)
    assert "HTTP 401" in str(info.value) and "invalid x-api-key" in str(info.value)
    assert KEY not in str(info.value)


def test_transient_errors_are_retried(captured, monkeypatch):
    calls, queue = captured
    waits: list[float] = []
    monkeypatch.setattr(providers.time, "sleep", waits.append)
    queue.extend(
        [
            _http_error(429, {"error": {"message": "slow down"}}, retry_after="3"),
            _http_error(529, "overloaded"),
            {"content": [{"type": "tool_use", "name": "denoise"}]},
        ]
    )
    assert providers.select_tool("anthropic", "m", "p", TOOLS, api_key=KEY) == "denoise"
    assert len(calls) == 3 and waits == [3.0, 4.0]

    queue.extend([_http_error(503, "down")] * (providers.RETRIES + 1))
    with pytest.raises(providers.ProviderError, match="HTTP 503"):
        providers.select_tool("anthropic", "m", "p", TOOLS, api_key=KEY)


def test_unreachable_endpoint(captured):
    _, queue = captured
    queue.append(urllib.error.URLError("Connection refused"))
    with pytest.raises(providers.ProviderError, match="Could not reach"):
        providers.select_tool("openai", "m", "p", TOOLS, base_url="http://127.0.0.1:9/v1")


def test_resolve_provider():
    assert providers.resolve_provider("openai", {"ANTHROPIC_API_KEY": "k"}) == "openai"
    assert providers.resolve_provider(None, {"ANTHROPIC_API_KEY": "k"}) == "anthropic"
    assert providers.resolve_provider(None, {"OPENAI_API_KEY": "k"}) == "openai"
    assert providers.resolve_provider(None, {"OPENAI_BASE_URL": "http://h/v1"}) == "openai"
    with pytest.raises(ValueError, match="ANTHROPIC_API_KEY or OPENAI_API_KEY"):
        providers.resolve_provider(None, {"ANTHROPIC_API_KEY": ""})
    with pytest.raises(ValueError, match="pass --provider"):
        providers.resolve_provider(None, {"ANTHROPIC_API_KEY": "a", "OPENAI_API_KEY": "b"})
    with pytest.raises(ValueError, match="Unknown provider"):
        providers.resolve_provider("acme", {})


def test_harness_requires_a_model_for_agentic(capsys):
    sys.path.insert(0, str(EVALS))
    try:
        import harness
    finally:
        sys.path.remove(str(EVALS))
    with pytest.raises(SystemExit) as info:
        harness.main(["--agentic", "--provider", "openai"])
    assert info.value.code == 2
    assert "--agentic requires --model" in capsys.readouterr().err


def test_real_http_round_trip(monkeypatch):
    """Both adapters against a local HTTP server: real sockets, real urllib."""
    for name in (*providers.API_KEY_ENV.values(), *providers.BASE_URL_ENV.values()):
        monkeypatch.delenv(name, raising=False)
    seen: list[tuple[str, dict]] = []

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):  # noqa: N802 (http.server API)
            body = json.loads(self.rfile.read(int(self.headers["content-length"])))
            seen.append((self.path, body))
            if self.path.endswith("/messages"):
                answer = {"content": [{"type": "tool_use", "name": body["tools"][0]["name"]}]}
            else:
                name = body["tools"][0]["function"]["name"]
                answer = {"choices": [{"message": {"tool_calls": [{"function": {"name": name}}]}}]}
            data = json.dumps(answer).encode()
            self.send_response(200)
            self.send_header("content-type", "application/json")
            self.send_header("content-length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def log_message(self, *args):
            pass

    server = HTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = f"http://127.0.0.1:{server.server_port}"
    try:
        assert (
            providers.select_tool("anthropic", "m", "p", TOOLS, base_url=base, api_key=KEY)
            == "denoise"
        )
        assert providers.select_tool("openai", "m", "p", TOOLS, base_url=f"{base}/v1") == "denoise"
    finally:
        server.shutdown()
        server.server_close()
    assert [path for path, _ in seen] == ["/v1/messages", "/v1/chat/completions"]
