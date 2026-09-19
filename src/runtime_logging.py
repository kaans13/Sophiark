"""Cross-platform runtime logging helpers with no import side effects."""

from __future__ import annotations

import sys
from typing import TextIO


def ensure_utf8_stream(stream: TextIO | None) -> None:
    reconfigure = getattr(stream, "reconfigure", None)
    if callable(reconfigure):
        try:
            reconfigure(encoding="utf-8", errors="replace")
        except (OSError, ValueError):
            # Embedded hosts may expose a stream that cannot be reconfigured.
            pass


def ensure_utf8_standard_streams() -> None:
    ensure_utf8_stream(sys.stdout)
    ensure_utf8_stream(sys.stderr)
