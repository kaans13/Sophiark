"""Safety tests for the lazy, bounded Phase 7 HTTP transport."""

from __future__ import annotations

import subprocess
import sys

import pytest
import requests

from src.research.config import ProviderTimeout
from src.research.http_transport import RequestsJsonTransport


TIMEOUT = ProviderTimeout(0.1, 0.2, 0.3)


class FakeResponse:
    def __init__(self, chunks, *, content_length=None) -> None:
        self.status_code = 200
        self.headers = {} if content_length is None else {"Content-Length": str(content_length)}
        self.url = "https://official.example/record"
        self.encoding = "utf-8"
        self._chunks = list(chunks)
        self.closed = False

    def iter_content(self, chunk_size):
        assert chunk_size == 64 * 1024
        yield from self._chunks

    def close(self):
        self.closed = True


class FakeSession:
    def __init__(self, response=None, error=None) -> None:
        self.response = response
        self.error = error
        self.kwargs = None

    def get(self, url, **kwargs):
        self.kwargs = kwargs
        if self.error is not None:
            raise self.error
        return self.response

    def close(self):
        pass


def test_transport_construction_is_lazy_and_import_does_not_load_requests() -> None:
    script = """
import sys
import src.research.http_transport as module
transport = module.RequestsJsonTransport()
assert transport._session is None
assert 'requests' not in sys.modules
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout


def test_transport_streams_with_component_timeouts_and_decodes_bounded_json() -> None:
    response = FakeResponse([b'{"ok":', b'true}'])
    session = FakeSession(response)
    transport = RequestsJsonTransport(max_response_bytes=100)
    transport._session = session

    result = transport.get_json("https://official.example/record", params={"limit": 3}, timeout=TIMEOUT)

    assert result.data == {"ok": True}
    assert session.kwargs["timeout"] == (0.1, 0.2)
    assert session.kwargs["stream"] is True
    assert session.kwargs["params"] == {"limit": 3}
    assert response.closed is True


@pytest.mark.parametrize(
    "response",
    (
        FakeResponse([b"{}"], content_length=101),
        FakeResponse([b"x" * 101]),
    ),
)
def test_transport_rejects_oversized_response_and_closes_it(response: FakeResponse) -> None:
    transport = RequestsJsonTransport(max_response_bytes=100)
    transport._session = FakeSession(response)
    with pytest.raises(ValueError, match="safety limit"):
        transport.get_json("https://official.example/record", params=None, timeout=TIMEOUT)
    assert response.closed is True


@pytest.mark.parametrize(
    ("error", "expected"),
    (
        (requests.exceptions.Timeout("slow"), TimeoutError),
        (requests.exceptions.ConnectionError("offline"), ConnectionError),
    ),
)
def test_requests_failures_are_translated_to_transport_neutral_errors(error, expected) -> None:
    transport = RequestsJsonTransport()
    transport._session = FakeSession(error=error)
    with pytest.raises(expected):
        transport.get_json("https://official.example/record", params=None, timeout=TIMEOUT)

