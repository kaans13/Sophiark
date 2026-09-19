"""Version-aware reuse of immutable product resources.

This module caches I/O and target-independent preparation only. Scientific
working graphs are always returned as detached copies to their callers.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from threading import Event, Lock
from types import MappingProxyType
from typing import Callable, Mapping

import pandas as pd

from .capabilities import CapabilityReport, detect_capabilities
from .datasets import DataRegistry
from .paths import ProjectPaths


_IGRAPH_BUILD_LOCK = Lock()


def _reset_igraph_rng(seed: int = 42) -> None:
    """Honor the frozen engines' declared seed on igraph versions without a seed kwarg."""
    try:
        import igraph as ig
        import random

        setter = getattr(ig, "set_random_number_generator", None)
        if callable(setter):
            setter(random.Random(seed))
    except (ImportError, RuntimeError):
        # Engine code retains its existing compatibility fallback.
        return


def file_identity(path: str | Path) -> tuple[str, int | None, int | None]:
    resolved = Path(path).resolve()
    try:
        stat = resolved.stat()
    except OSError:
        return str(resolved), None, None
    return str(resolved), stat.st_size, stat.st_mtime_ns


@lru_cache(maxsize=8)
def _mapping_rows(identity: tuple[str, int | None, int | None]) -> tuple[tuple[str, str], ...]:
    path = Path(identity[0])
    if not path.is_file():
        return ()
    frame = pd.read_csv(path, usecols=lambda column: column in {"gene", "Symbol"}).dropna()
    return tuple(zip(frame["gene"].astype(str), frame["Symbol"].astype(str)))


def symbol_mapping(path: str | Path) -> dict[str, str]:
    """Return a fresh mapping backed by version-aware cached immutable rows."""
    return dict(_mapping_rows(file_identity(path)))


@dataclass(frozen=True, slots=True)
class ProductDataSnapshot:
    identity: tuple[tuple[str, int | None, int | None], ...]
    versions: Mapping[str, str]
    fingerprints: Mapping[str, str]
    capabilities: Mapping[str, str]
    capability_report: CapabilityReport


@lru_cache(maxsize=4)
def _registry(
    root: str, catalog_identity: tuple[str, int | None, int | None],
) -> DataRegistry:
    return DataRegistry.from_json(catalog_identity[0], root=root)


@lru_cache(maxsize=4)
def _snapshot(
    root: str,
    catalog_identity: tuple[str, int | None, int | None],
    resource_identity: tuple[tuple[str, int | None, int | None], ...],
) -> ProductDataSnapshot:
    registry = _registry(root, catalog_identity)
    health = registry.health_report(fingerprint=True)
    capabilities = detect_capabilities(registry)
    return ProductDataSnapshot(
        identity=(catalog_identity, *resource_identity),
        versions=MappingProxyType({name: item.version for name, item in health.items()}),
        fingerprints=MappingProxyType({
            name: item.fingerprint for name, item in health.items()
            if item.fingerprint is not None
        }),
        capabilities=MappingProxyType({
            name.value: status.value for name, status in capabilities.statuses.items()
        }),
        capability_report=capabilities,
    )


def product_data_snapshot(paths: ProjectPaths | None = None) -> ProductDataSnapshot:
    active = paths or ProjectPaths.discover()
    catalog_identity = file_identity(active.config / "datasets.json")
    registry = _registry(str(active.root), catalog_identity)
    identities = tuple(
        file_identity(path) if (path := registry.resolve(registry.get(name))) is not None
        else (f"missing:{name}", None, None)
        for name in registry.names()
    )
    return _snapshot(str(active.root), catalog_identity, identities)


class DetachedPreparedContextCache:
    """Small LRU for immutable base contexts; consumers receive detached state."""

    def __init__(self, maxsize: int = 2) -> None:
        self.maxsize = maxsize
        self._items: OrderedDict[tuple[object, ...], tuple[object, pd.DataFrame]] = OrderedDict()
        self._lock = Lock()
        self._building: dict[tuple[object, ...], Event] = {}

    def get_or_build(
        self,
        key: tuple[object, ...],
        builder: Callable[[], tuple[object, pd.DataFrame]],
    ) -> tuple[object, pd.DataFrame]:
        while True:
            with self._lock:
                cached = self._items.get(key)
                if cached is not None:
                    self._items.move_to_end(key)
                    break
                ready = self._building.get(key)
                if ready is None:
                    ready = Event()
                    self._building[key] = ready
                    is_builder = True
                else:
                    is_builder = False
            if not is_builder:
                ready.wait()
                continue
            try:
                # python-igraph exposes one process-wide RNG.  Two tissue builds
                # must not interleave between resetting that RNG and Louvain.
                with _IGRAPH_BUILD_LOCK:
                    _reset_igraph_rng()
                    built_graph, built_scores = builder()
                cached = (built_graph.copy(), built_scores.copy(deep=True))
                with self._lock:
                    self._items[key] = cached
                    self._items.move_to_end(key)
                    while len(self._items) > self.maxsize:
                        self._items.popitem(last=False)
                break
            finally:
                with self._lock:
                    self._building.pop(key, None)
                    ready.set()
        graph, scores = cached
        return graph.copy(), scores.copy(deep=True)

    def clear(self) -> None:
        with self._lock:
            self._items.clear()


def clear_runtime_resource_caches() -> None:
    _mapping_rows.cache_clear()
    _snapshot.cache_clear()
    _registry.cache_clear()
