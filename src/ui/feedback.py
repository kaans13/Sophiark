"""User-safe status messaging; technical details belong in logs, not the UI."""

from __future__ import annotations

import logging
import streamlit as st

log = logging.getLogger(__name__)


def user_error(message: str, *, exception: Exception | None = None) -> None:
    if exception is not None:
        log.error(message, exc_info=(type(exception), exception, exception.__traceback__))
    st.error(message)


def user_warning(message: str, *, exception: Exception | None = None) -> None:
    if exception is not None:
        log.warning("%s: %s", message, exception)
    st.warning(message)
