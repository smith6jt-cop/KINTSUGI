"""Authentication for the remote (streamable HTTP) KINTSUGI MCP server.

The stdio server needs no authentication: the MCP client starts it as a
subprocess with the user's own permissions. Over HTTP anyone who can reach the
port can talk to it, so ``kintsugi mcp start --transport streamable-http``
requires a secret, ``KINTSUGI_MCP_TOKEN``, and accepts it two ways:

* **Bearer token** — clients that let the user set a header (Claude Code:
  ``claude mcp add --transport http ... --header "Authorization: Bearer $TOKEN"``)
  send the secret itself.
* **OAuth 2.1** — clients that only speak OAuth (claude.ai custom connectors)
  register dynamically, are sent to a ``/login`` page where the user types the
  same secret, and receive short-lived access tokens bound to this server
  (PKCE, the token endpoint and the metadata documents are the MCP SDK's).

This is a *single-user* authorization server: there is one secret and one
principal. Registration is limited to redirect URIs on an allowlist (claude.ai's
callback and loopback addresses by default), failed logins are throttled, and
all state is in memory, so a restart signs every client out.

Nothing here runs unless the HTTP transport is selected; importing the module
requires the ``mcp`` extra.
"""

from __future__ import annotations

import hmac
import html
import secrets
import time
from dataclasses import dataclass, field
from urllib.parse import urlparse

from mcp.server.auth.provider import (
    AccessToken,
    AuthorizationCode,
    AuthorizationParams,
    AuthorizeError,
    RefreshToken,
    RegistrationError,
    TokenError,
    construct_redirect_uri,
)
from mcp.server.auth.settings import AuthSettings, ClientRegistrationOptions, RevocationOptions
from mcp.shared.auth import OAuthClientInformationFull, OAuthToken

TOKEN_ENV = "KINTSUGI_MCP_TOKEN"
MIN_SECRET_LENGTH = 32

#: claude.ai's OAuth callback for custom connectors.
CLAUDE_AI_CALLBACK = "https://claude.ai/api/mcp/auth_callback"
DEFAULT_REDIRECT_URIS: tuple[str, ...] = (CLAUDE_AI_CALLBACK,)
LOOPBACK_HOSTS = frozenset({"localhost", "127.0.0.1", "::1"})

SCOPE = "kintsugi"
SUBJECT = "kintsugi-owner"
STATIC_CLIENT_ID = "kintsugi-static-token"

LOGIN_PATH = "/login"
CODE_TTL_S = 300
REQUEST_TTL_S = 600
ACCESS_TTL_S = 3600
REFRESH_TTL_S = 30 * 24 * 3600
MAX_FAILURES = 5
LOCKOUT_S = 300

# Bounds on state that unauthenticated callers can create (registration and
# authorization requests need no secret, so they must not grow without limit).
MAX_CLIENTS = 64
MAX_PENDING = 32
CLIENT_IDLE_TTL_S = 24 * 3600
RATE_WINDOW_S = 60
MAX_UNAUTH_PER_WINDOW = 30
MAX_CONSUMED_REFRESH = 1024


class SecretError(ValueError):
    """``KINTSUGI_MCP_TOKEN`` is missing or too short for a network listener."""


def validate_secret(secret: str | None) -> str:
    """Return ``secret`` or raise :class:`SecretError` if it is unusable."""
    if not secret:
        raise SecretError(
            f"{TOKEN_ENV} is not set. Generate one with: "
            "python -c 'import secrets; print(secrets.token_urlsafe(32))'"
        )
    if len(secret) < MIN_SECRET_LENGTH:
        raise SecretError(f"{TOKEN_ENV} must be at least {MIN_SECRET_LENGTH} characters long.")
    return secret


def redirect_allowed(uri: str, extra: tuple[str, ...] = ()) -> bool:
    """True for an exact allowlisted URI or an ``http`` loopback URI (any port/path).

    Loopback redirects are how desktop and CLI clients (Claude Code) receive the
    code; RFC 8252 section 7.3 allows any port for them.
    """
    if uri in DEFAULT_REDIRECT_URIS or uri in extra:
        return True
    parsed = urlparse(uri)
    return parsed.scheme == "http" and (parsed.hostname or "") in LOOPBACK_HOSTS


@dataclass
class _PendingLogin:
    client_id: str
    params: AuthorizationParams
    expires_at: float


@dataclass
class _Throttle:
    failures: int = 0
    locked_until: float = 0.0

    def locked(self, now: float) -> bool:
        return now < self.locked_until

    def fail(self, now: float) -> None:
        self.failures += 1
        if self.failures >= MAX_FAILURES:
            self.locked_until = now + LOCKOUT_S
            self.failures = 0

    def reset(self) -> None:
        self.failures = 0
        self.locked_until = 0.0


@dataclass
class _RateLimit:
    """Sliding-window limit on unauthenticated registration/authorization calls."""

    limit: int = MAX_UNAUTH_PER_WINDOW
    window: float = RATE_WINDOW_S
    events: list[float] = field(default_factory=list)

    def allow(self, now: float) -> bool:
        self.events = [t for t in self.events if now - t < self.window]
        if len(self.events) >= self.limit:
            return False
        self.events.append(now)
        return True


def _evict_oldest(store: dict, limit: int) -> None:
    """Drop the oldest entries (dicts keep insertion order) until below ``limit``."""
    while len(store) >= limit:
        del store[next(iter(store))]


@dataclass
class SingleUserOAuthProvider:
    """In-memory OAuth authorization server for one secret-holding user.

    Implements the MCP SDK's ``OAuthAuthorizationServerProvider`` protocol;
    :meth:`complete_login` is the second half of the authorization step, called
    by the ``/login`` route once the user has typed the secret.
    """

    secret: str
    public_url: str
    resource_url: str
    extra_redirect_uris: tuple[str, ...] = ()
    clients: dict[str, OAuthClientInformationFull] = field(default_factory=dict)
    pending: dict[str, _PendingLogin] = field(default_factory=dict)
    codes: dict[str, AuthorizationCode] = field(default_factory=dict)
    access_tokens: dict[str, AccessToken] = field(default_factory=dict)
    refresh_tokens: dict[str, RefreshToken] = field(default_factory=dict)
    throttle: _Throttle = field(default_factory=_Throttle)
    rate: _RateLimit = field(default_factory=_RateLimit)
    #: client_id -> last time it registered or authorized (idle expiry).
    client_seen: dict[str, float] = field(default_factory=dict)
    #: token (access or refresh) -> token family; one family per authorization.
    family_of: dict[str, str] = field(default_factory=dict)
    #: refresh tokens already exchanged -> (family, expiry); replay revokes the family.
    consumed: dict[str, tuple[str, float]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        validate_secret(self.secret)
        self.public_url = self.public_url.rstrip("/")

    # -- secret ------------------------------------------------------------

    def secret_matches(self, candidate: str) -> bool:
        """Constant-time comparison with the configured secret."""
        return hmac.compare_digest(candidate.encode(), self.secret.encode())

    # -- clients -----------------------------------------------------------

    async def get_client(self, client_id: str) -> OAuthClientInformationFull | None:
        self._expire(time.time())
        return self.clients.get(client_id)

    async def register_client(self, client_info: OAuthClientInformationFull) -> None:
        uris = [str(u) for u in (client_info.redirect_uris or [])]
        if not uris:
            raise RegistrationError("invalid_redirect_uri", "redirect_uris is required")
        refused = [u for u in uris if not redirect_allowed(u, self.extra_redirect_uris)]
        if refused:
            raise RegistrationError(
                "invalid_redirect_uri",
                f"redirect URI not allowed by this server: {', '.join(refused)}",
            )
        if client_info.client_id is None:
            raise RegistrationError("invalid_client_metadata", "client_id was not assigned")
        now = time.time()
        if not self.rate.allow(now):
            raise RegistrationError(
                "invalid_client_metadata", "too many requests; try again in a minute"
            )
        self._expire(now)
        _evict_oldest(self.clients, MAX_CLIENTS)
        self.clients[client_info.client_id] = client_info
        self.client_seen[client_info.client_id] = now

    # -- authorization -----------------------------------------------------

    async def authorize(
        self, client: OAuthClientInformationFull, params: AuthorizationParams
    ) -> str:
        if client.client_id is None:
            raise AuthorizeError("unauthorized_client", "unknown client")
        now = time.time()
        if not self.rate.allow(now):
            raise AuthorizeError("temporarily_unavailable", "too many requests")
        self._expire(now)
        self.client_seen[client.client_id] = now
        _evict_oldest(self.pending, MAX_PENDING)
        request_id = secrets.token_urlsafe(32)
        self.pending[request_id] = _PendingLogin(
            client_id=client.client_id, params=params, expires_at=now + REQUEST_TTL_S
        )
        return f"{self.public_url}{LOGIN_PATH}?request={request_id}"

    def pending_login(self, request_id: str) -> _PendingLogin | None:
        """The live authorization request ``request_id``, or ``None``."""
        self._expire(time.time())
        return self.pending.get(request_id)

    def complete_login(self, request_id: str, candidate: str) -> str:
        """Check the typed secret and return the client redirect carrying the code.

        Raises ``PermissionError`` for a wrong secret or while locked out, and
        ``KeyError`` for an unknown or expired request.
        """
        now = time.time()
        if self.throttle.locked(now):
            raise PermissionError("Too many failed attempts; try again in a few minutes.")
        pending = self.pending_login(request_id)
        if pending is None:
            raise KeyError("This sign-in request has expired; start again from your client.")
        if not self.secret_matches(candidate):
            self.throttle.fail(now)
            raise PermissionError("Incorrect secret.")
        self.throttle.reset()
        del self.pending[request_id]

        params = pending.params
        code = secrets.token_urlsafe(32)
        self.codes[code] = AuthorizationCode(
            code=code,
            scopes=params.scopes or [SCOPE],
            expires_at=now + CODE_TTL_S,
            client_id=pending.client_id,
            code_challenge=params.code_challenge,
            redirect_uri=params.redirect_uri,
            redirect_uri_provided_explicitly=params.redirect_uri_provided_explicitly,
            resource=params.resource or self.resource_url,
            subject=SUBJECT,
        )
        return construct_redirect_uri(str(params.redirect_uri), code=code, state=params.state)

    async def load_authorization_code(
        self, client: OAuthClientInformationFull, authorization_code: str
    ) -> AuthorizationCode | None:
        code = self.codes.get(authorization_code)
        if code is None or code.client_id != client.client_id:
            return None
        if code.expires_at < time.time():
            self.codes.pop(authorization_code, None)
            return None
        return code

    async def exchange_authorization_code(
        self, client: OAuthClientInformationFull, authorization_code: AuthorizationCode
    ) -> OAuthToken:
        if self.codes.pop(authorization_code.code, None) is None:
            raise TokenError("invalid_grant", "authorization code already used")
        self.client_seen[authorization_code.client_id] = time.time()
        return self._issue(
            authorization_code.client_id,
            authorization_code.scopes,
            authorization_code.resource,
            family=secrets.token_urlsafe(16),
        )

    # -- refresh -----------------------------------------------------------

    async def load_refresh_token(
        self, client: OAuthClientInformationFull, refresh_token: str
    ) -> RefreshToken | None:
        replay = self.consumed.get(refresh_token)
        if replay is not None:
            # A rotated-out refresh token was presented again: it was stolen or
            # leaked, so revoke every token descended from the same authorization.
            self._revoke_family(replay[0])
            return None
        token = self.refresh_tokens.get(refresh_token)
        if token is None or token.client_id != client.client_id:
            return None
        if token.expires_at is not None and token.expires_at < time.time():
            self.refresh_tokens.pop(refresh_token, None)
            return None
        return token

    async def exchange_refresh_token(
        self,
        client: OAuthClientInformationFull,
        refresh_token: RefreshToken,
        scopes: list[str],
    ) -> OAuthToken:
        # Rotate: the presented refresh token is single-use; remember it so a
        # later replay revokes the whole family.
        if self.refresh_tokens.pop(refresh_token.token, None) is None:
            raise TokenError("invalid_grant", "refresh token already used")
        family = self.family_of.pop(refresh_token.token, None) or secrets.token_urlsafe(16)
        _evict_oldest(self.consumed, MAX_CONSUMED_REFRESH)
        self.consumed[refresh_token.token] = (family, time.time() + REFRESH_TTL_S)
        return self._issue(
            refresh_token.client_id,
            scopes or refresh_token.scopes,
            refresh_token.resource,
            family=family,
        )

    # -- access tokens -----------------------------------------------------

    async def load_access_token(self, token: str) -> AccessToken | None:
        if self.secret_matches(token):
            return AccessToken(
                token=token,
                client_id=STATIC_CLIENT_ID,
                scopes=[SCOPE],
                resource=self.resource_url,
                subject=SUBJECT,
            )
        access = self.access_tokens.get(token)
        if access is None:
            return None
        if access.expires_at is not None and access.expires_at < time.time():
            self.access_tokens.pop(token, None)
            return None
        return access

    async def revoke_token(self, token: AccessToken | RefreshToken) -> None:
        self.access_tokens.pop(token.token, None)
        self.refresh_tokens.pop(token.token, None)
        self.family_of.pop(token.token, None)

    # -- helpers -----------------------------------------------------------

    def _issue(
        self, client_id: str, scopes: list[str], resource: str | None, *, family: str
    ) -> OAuthToken:
        now = int(time.time())
        access = secrets.token_urlsafe(32)
        refresh = secrets.token_urlsafe(32)
        resource = resource or self.resource_url
        self.access_tokens[access] = AccessToken(
            token=access,
            client_id=client_id,
            scopes=scopes,
            expires_at=now + ACCESS_TTL_S,
            resource=resource,
            subject=SUBJECT,
        )
        self.refresh_tokens[refresh] = RefreshToken(
            token=refresh,
            client_id=client_id,
            scopes=scopes,
            expires_at=now + REFRESH_TTL_S,
            resource=resource,
            subject=SUBJECT,
        )
        self.family_of[access] = family
        self.family_of[refresh] = family
        return OAuthToken(
            access_token=access,
            expires_in=ACCESS_TTL_S,
            scope=" ".join(scopes),
            refresh_token=refresh,
        )

    def _revoke_family(self, family: str) -> None:
        for token in [t for t, f in self.family_of.items() if f == family]:
            self.access_tokens.pop(token, None)
            self.refresh_tokens.pop(token, None)
            del self.family_of[token]

    def _expire(self, now: float) -> None:
        for key in [k for k, v in self.pending.items() if v.expires_at < now]:
            del self.pending[key]
        for key in [k for k, (_, exp) in self.consumed.items() if exp < now]:
            del self.consumed[key]
        live = {t.client_id for t in self.refresh_tokens.values()}
        live |= {p.client_id for p in self.pending.values()}
        for cid in [
            c
            for c in self.clients
            if c not in live and now - self.client_seen.get(c, now) > CLIENT_IDLE_TTL_S
        ]:
            del self.clients[cid]
            self.client_seen.pop(cid, None)


def auth_settings(public_url: str, resource_url: str) -> AuthSettings:
    """``AuthSettings`` for :class:`SingleUserOAuthProvider` served at ``public_url``."""
    # model_validate (not the constructor) so the URLs keep their exact string
    # form: OAuth issuer comparison is exact and AnyHttpUrl("https://h") would
    # gain a trailing slash.
    return AuthSettings.model_validate(
        {
            "issuer_url": public_url.rstrip("/"),
            "resource_server_url": resource_url,
            "client_registration_options": ClientRegistrationOptions(
                enabled=True, valid_scopes=[SCOPE], default_scopes=[SCOPE]
            ),
            "revocation_options": RevocationOptions(enabled=True),
            "required_scopes": [SCOPE],
            "validate_token_resource": True,
        }
    )


# ---------------------------------------------------------------------------
# /login page
# ---------------------------------------------------------------------------

_PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>KINTSUGI MCP sign-in</title>
<style>
body{{font-family:system-ui,sans-serif;max-width:28rem;margin:4rem auto;padding:0 1rem;
color:#1b1b1b;background:#fafafa}}
input,button{{font:inherit;padding:.5rem;width:100%;box-sizing:border-box;margin-top:.5rem}}
.err{{color:#a00}}
</style></head><body>
<h1>KINTSUGI MCP server</h1>
<p>Client <code>{client}</code> asks to use this server's image-processing tools.
Enter the server secret (<code>KINTSUGI_MCP_TOKEN</code>) to allow it.</p>
{error}
<form method="post" action="{action}">
<input type="hidden" name="request" value="{request}">
<label>Secret<input type="password" name="secret" autocomplete="current-password" required
autofocus></label>
<button type="submit">Allow</button>
</form></body></html>"""


def render_login(client_id: str, request_id: str, error: str | None = None) -> str:
    """The sign-in form for the pending request ``request_id``."""
    return _PAGE.format(
        client=html.escape(client_id),
        request=html.escape(request_id, quote=True),
        action=LOGIN_PATH,
        error=f'<p class="err">{html.escape(error)}</p>' if error else "",
    )


def register_login_route(server, provider: SingleUserOAuthProvider) -> None:  # noqa: ANN001
    """Add ``GET/POST /login`` (and ``GET /healthz``) to ``server``."""
    from starlette.requests import Request
    from starlette.responses import HTMLResponse, JSONResponse, RedirectResponse, Response

    headers = {
        "Cache-Control": "no-store",
        "X-Frame-Options": "DENY",
        "Content-Security-Policy": "default-src 'none'; style-src 'unsafe-inline'; "
        "form-action 'self' https:; frame-ancestors 'none'",
        "Referrer-Policy": "no-referrer",
    }

    @server.custom_route(LOGIN_PATH, methods=["GET", "POST"], include_in_schema=False)
    async def login(request: Request) -> Response:
        if request.method == "GET":
            request_id = request.query_params.get("request", "")
            pending = provider.pending_login(request_id)
            if pending is None:
                return HTMLResponse(
                    "Sign-in request expired; start again from your client.",
                    status_code=400,
                    headers=headers,
                )
            return HTMLResponse(render_login(pending.client_id, request_id), headers=headers)

        form = await request.form()
        request_id = str(form.get("request", ""))
        candidate = str(form.get("secret", ""))
        try:
            target = provider.complete_login(request_id, candidate)
        except KeyError as exc:
            return HTMLResponse(html.escape(str(exc.args[0])), status_code=400, headers=headers)
        except PermissionError as exc:
            pending = provider.pending_login(request_id)
            client_id = pending.client_id if pending else "unknown"
            return HTMLResponse(
                render_login(client_id, request_id, str(exc)), status_code=401, headers=headers
            )
        return RedirectResponse(target, status_code=302, headers=headers)

    @server.custom_route("/healthz", methods=["GET"], include_in_schema=False)
    async def healthz(request: Request) -> Response:  # noqa: ARG001
        return JSONResponse({"status": "ok"})
