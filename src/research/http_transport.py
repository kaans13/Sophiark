"""Small, bounded JSON transport for optional Research Context providers.

Importing this module performs no network activity.  The requests session is
created lazily on the first explicit provider call, and every response is
bounded before JSON decoding so an external service cannot consume unbounded
memory in the application process.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from types import MappingProxyType
from typing import Any, Mapping, Protocol

from .config import ProviderTimeout


DEFAULT_MAX_RESPONSE_BYTES = 4 * 1024 * 1024


@dataclass(frozen=True, slots=True)
class JsonHttpResponse:
    """Transport-neutral subset of an HTTP JSON response."""

    status_code: int
    data: Any
    headers: Mapping[str, str]
    url: str

    def __post_init__(self) -> None:
        if isinstance(self.status_code, bool) or not isinstance(self.status_code, int):
            raise TypeError("status_code must be an integer")
        object.__setattr__(
            self,
            "headers",
            MappingProxyType({str(key).casefold(): str(value) for key, value in self.headers.items()}),
        )
        object.__setattr__(self, "url", str(self.url).strip())

    def header(self, name: str) -> str | None:
        return self.headers.get(str(name).casefold())


class JsonTransport(Protocol):
    """Injectable transport contract used by all concrete providers."""

    def get_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None,
        timeout: ProviderTimeout,
        headers: Mapping[str, str] | None = None,
    ) -> JsonHttpResponse:
        ...


class RequestsJsonTransport:
    """Lazy requests-based implementation with finite connect/read limits."""

    def __init__(
        self,
        *,
        user_agent: str = "SOPHIARK-Research-Context/1",
        max_response_bytes: int = DEFAULT_MAX_RESPONSE_BYTES,
    ) -> None:
        if isinstance(max_response_bytes, bool) or not isinstance(max_response_bytes, int):
            raise TypeError("max_response_bytes must be an integer")
        if max_response_bytes <= 0:
            raise ValueError("max_response_bytes must be positive")
        self._user_agent = str(user_agent).strip() or "SOPHIARK-Research-Context/1"
        self._max_response_bytes = max_response_bytes
        self._session: Any = None

    @property
    def max_response_bytes(self) -> int:
        return self._max_response_bytes

    def _get_session(self) -> Any:
        if self._session is not None:
            return self._session
        import requests

        # A fresh session per concurrent call avoids sharing mutable requests
        # session state across ProviderExecutor worker threads.
        return requests.Session()

    def get_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None,
        timeout: ProviderTimeout,
        headers: Mapping[str, str] | None = None,
    ) -> JsonHttpResponse:
        if not isinstance(timeout, ProviderTimeout):
            raise TypeError("timeout must be ProviderTimeout")
        request_headers = {
            "Accept": "application/json",
            "User-Agent": self._user_agent,
            **{str(key): str(value) for key, value in (headers or {}).items()},
        }
        session = self._get_session()
        owns_session = self._session is None
        try:
            response = session.get(
                str(url),
                params=dict(params or {}),
                headers=request_headers,
                timeout=(timeout.connect_seconds, timeout.read_seconds),
                stream=True,
            )
        except Exception as exc:
            # Keep requests optional at import time and translate its exception
            # hierarchy to the runtime's transport-neutral boundaries.
            try:
                import requests

                if isinstance(exc, requests.exceptions.Timeout):
                    raise TimeoutError(str(exc) or "External provider timed out") from exc
                if isinstance(exc, requests.exceptions.ConnectionError):
                    raise ConnectionError(str(exc) or "External provider is offline") from exc
            except ImportError:
                pass
            raise
        finally:
            if owns_session and "response" not in locals():
                session.close()

        try:
            content_length = response.headers.get("Content-Length")
            if content_length is not None and int(content_length) > self._max_response_bytes:
                raise ValueError(
                    f"Provider response exceeds {self._max_response_bytes} byte safety limit"
                )
            content = bytearray()
            for chunk in response.iter_content(chunk_size=64 * 1024):
                if not chunk:
                    continue
                content.extend(chunk)
                if len(content) > self._max_response_bytes:
                    raise ValueError(
                        f"Provider response exceeds {self._max_response_bytes} byte safety limit"
                    )
            if content:
                try:
                    encoding = response.encoding or "utf-8"
                    data = json.loads(bytes(content).decode(encoding))
                except (UnicodeDecodeError, ValueError, json.JSONDecodeError) as exc:
                    raise ValueError("Provider returned invalid JSON") from exc
            else:
                data = None
            return JsonHttpResponse(
                status_code=int(response.status_code),
                data=data,
                headers=dict(response.headers),
                url=str(response.url),
            )
        finally:
            response.close()
            if owns_session:
                session.close()

    def close(self) -> None:
        if self._session is not None:
            self._session.close()
            self._session = None


__all__ = [
    "DEFAULT_MAX_RESPONSE_BYTES",
    "JsonHttpResponse",
    "JsonTransport",
    "RequestsJsonTransport",
]
