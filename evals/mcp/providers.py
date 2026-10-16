"""Model providers for the agentic tool-selection eval (standard library only).

Each adapter sends one chat request that offers the live MCP tool definitions to a
model and returns the name of the first tool the model calls (``None`` when it
answers in text). Two wire formats cover the common platforms:

``anthropic``
    The Anthropic Messages API (``POST {base}/v1/messages``).
``openai``
    The OpenAI Chat Completions API (``POST {base}/chat/completions``), which
    many platforms and local servers also implement: set ``OPENAI_BASE_URL`` (or
    ``--base-url``) to point it at one, for example Ollama at
    ``http://localhost:11434/v1`` or a vLLM server.

No vendor SDK is needed: a request is one JSON POST through :mod:`urllib`.
"""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.request
from collections.abc import Mapping
from typing import Any
from urllib.parse import urlsplit

PROVIDERS = ("anthropic", "openai")

API_KEY_ENV = {"anthropic": "ANTHROPIC_API_KEY", "openai": "OPENAI_API_KEY"}
BASE_URL_ENV = {"anthropic": "ANTHROPIC_BASE_URL", "openai": "OPENAI_BASE_URL"}
DEFAULT_BASE_URL = {
    "anthropic": "https://api.anthropic.com",
    "openai": "https://api.openai.com/v1",
}

ANTHROPIC_VERSION = "2023-06-01"
MAX_TOKENS = 512

# Statuses worth retrying: rate limits, transient server errors, overload.
RETRY_STATUSES = frozenset({408, 409, 429, 500, 502, 503, 504, 529})
RETRIES = 2
MAX_RETRY_WAIT_S = 30.0

# JSON Schema keywords that carry no constraint and that some OpenAI-compatible
# servers reject or mishandle.
_DROP_KEYWORDS = frozenset({"title", "default", "format"})


class ProviderError(RuntimeError):
    """A provider request failed (configuration, network or HTTP error)."""


def resolve_base_url(provider: str, explicit: str | None = None) -> str:
    """The endpoint base: ``explicit``, else the provider's environment variable, else default."""
    _check_provider(provider)
    url = explicit or os.environ.get(BASE_URL_ENV[provider]) or DEFAULT_BASE_URL[provider]
    return url.rstrip("/")


def endpoint_host(provider: str, explicit: str | None = None) -> str:
    """Host (and port) of the endpoint, without credentials or path, for reports."""
    parts = urlsplit(resolve_base_url(provider, explicit))
    host = parts.hostname or ""
    return f"{host}:{parts.port}" if parts.port else host


def resolve_provider(explicit: str | None, env: Mapping[str, str] | None = None) -> str:
    """The provider to use: ``explicit``, else the only one configured in ``env``.

    A provider counts as configured when its API key variable is set; ``openai``
    also counts when ``OPENAI_BASE_URL`` is set (a local server may need no key).
    """
    if explicit:
        _check_provider(explicit)
        return explicit
    env = os.environ if env is None else env
    configured = [p for p in PROVIDERS if env.get(API_KEY_ENV[p])]
    if env.get(BASE_URL_ENV["openai"]) and "openai" not in configured:
        configured.append("openai")
    if len(configured) == 1:
        return configured[0]
    if not configured:
        raise ValueError(
            "No model provider configured: set "
            + " or ".join(API_KEY_ENV[p] for p in PROVIDERS)
            + " (for a local OpenAI-compatible server, set OPENAI_BASE_URL), "
            "or pass --provider."
        )
    raise ValueError(
        f"Several providers are configured ({', '.join(configured)}): pass --provider."
    )


def sanitize_schema(schema: Any) -> Any:
    """A JSON Schema without annotation-only keywords, for OpenAI-compatible servers.

    Drops ``title``, ``default`` and ``format``, and turns an optional parameter's
    ``anyOf: [T, {"type": "null"}]`` into ``T``. Constraints (``enum``, bounds,
    ``required``, ``items``, ``additionalProperties``) are kept.
    """
    if isinstance(schema, list):
        return [sanitize_schema(item) for item in schema]
    if not isinstance(schema, dict):
        return schema

    out: dict[str, Any] = {}
    for key, value in schema.items():
        if key in _DROP_KEYWORDS:
            continue
        if key in ("properties", "$defs", "definitions") and isinstance(value, dict):
            # Keys here are parameter names, not keywords: a parameter may be called "title".
            out[key] = {name: sanitize_schema(sub) for name, sub in value.items()}
        else:
            out[key] = sanitize_schema(value)

    options = out.get("anyOf")
    if isinstance(options, list):
        non_null = [o for o in options if not (isinstance(o, dict) and o.get("type") == "null")]
        if len(non_null) == 1 and len(non_null) < len(options) and isinstance(non_null[0], dict):
            del out["anyOf"]
            out = {**non_null[0], **out}
    return out


def select_tool(
    provider: str,
    model: str,
    prompt: str,
    tools: list[dict[str, Any]],
    *,
    base_url: str | None = None,
    api_key: str | None = None,
    timeout: float = 60,
) -> str | None:
    """Offer ``tools`` to ``model`` with ``prompt``; return the first tool it calls.

    ``tools`` are the advertised MCP definitions (``name``, ``description``,
    ``input_schema``). Raises :class:`ProviderError` on a configuration, network
    or HTTP error; the API key never appears in its message.
    """
    _check_provider(provider)
    root = resolve_base_url(provider, base_url)
    key = api_key or os.environ.get(API_KEY_ENV[provider]) or ""
    if not key and (provider == "anthropic" or root == DEFAULT_BASE_URL["openai"]):
        raise ProviderError(f"{provider}: set {API_KEY_ENV[provider]}.")

    if provider == "anthropic":
        url = f"{root}/v1/messages"
        headers = {"x-api-key": key, "anthropic-version": ANTHROPIC_VERSION}
        body: dict[str, Any] = {
            "model": model,
            "max_tokens": MAX_TOKENS,
            "tools": [
                {
                    "name": t["name"],
                    "description": t["description"],
                    "input_schema": t["input_schema"],
                }
                for t in tools
            ],
            "tool_choice": {"type": "auto"},
            "messages": [{"role": "user", "content": prompt}],
        }
        data = _post(url, headers, body, key, timeout)
        for block in data.get("content") or []:
            if isinstance(block, dict) and block.get("type") == "tool_use":
                return block.get("name")
        return None

    url = f"{root}/chat/completions"
    headers = {"Authorization": f"Bearer {key}"} if key else {}
    # No token limit: its field name differs between compatible servers, and the
    # answer is a single tool call.
    body = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "tools": [
            {
                "type": "function",
                "function": {
                    "name": t["name"],
                    "description": t["description"],
                    "parameters": sanitize_schema(t["input_schema"]),
                },
            }
            for t in tools
        ],
        "tool_choice": "auto",
    }
    data = _post(url, headers, body, key, timeout)
    choices = data.get("choices") or [{}]
    calls = (choices[0].get("message") or {}).get("tool_calls") or []
    if calls:
        return (calls[0].get("function") or {}).get("name")
    return None


def _check_provider(provider: str) -> None:
    if provider not in PROVIDERS:
        raise ValueError(f"Unknown provider {provider!r}; choose from {', '.join(PROVIDERS)}.")


def _post(
    url: str, headers: dict[str, str], body: dict[str, Any], key: str, timeout: float
) -> dict[str, Any]:
    payload = json.dumps(body).encode()
    request_headers = {"content-type": "application/json", **headers}
    for attempt in range(RETRIES + 1):
        request = urllib.request.Request(url, data=payload, headers=request_headers, method="POST")
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode() or "{}")
        except urllib.error.HTTPError as exc:
            detail = _redact(_error_detail(exc), key)
            if exc.code in RETRY_STATUSES and attempt < RETRIES:
                time.sleep(_retry_wait(exc, attempt))
                continue
            raise ProviderError(f"HTTP {exc.code} from {url}: {detail}") from None
        except urllib.error.URLError as exc:
            raise ProviderError(f"Could not reach {url}: {_redact(str(exc.reason), key)}") from None
        except json.JSONDecodeError:
            raise ProviderError(f"{url} did not return JSON.") from None
    raise AssertionError("unreachable")  # pragma: no cover


def _error_detail(exc: urllib.error.HTTPError) -> str:
    try:
        text = exc.read().decode(errors="replace")
    except Exception:
        return exc.reason or ""
    try:
        error = json.loads(text).get("error")
    except (ValueError, AttributeError):
        return text[:300]
    if isinstance(error, dict) and error.get("message"):
        return str(error["message"])[:300]
    return str(error or text)[:300]


def _retry_wait(exc: urllib.error.HTTPError, attempt: int) -> float:
    try:
        wait = float(exc.headers.get("retry-after", ""))
    except (TypeError, ValueError):
        wait = 2.0 * 2**attempt
    return min(max(wait, 0.0), MAX_RETRY_WAIT_S)


def _redact(text: str, key: str) -> str:
    return text.replace(key, "<redacted>") if key else text
