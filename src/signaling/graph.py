"""Indexed directed graph and bounded signed traversal for Layer 2."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field

from .models import PathSummary, SignalingDataset, SignalingEdge, SignStatus


def path_sign(edges: tuple[SignalingEdge, ...]) -> str:
    """Return a deterministic net sign; uncertainty is never guessed away."""

    product = 1
    for edge in edges:
        if edge.sign_status in {SignStatus.UNSIGNED, SignStatus.CONFLICTING}:
            return "uncertain"
        product *= edge.canonical_sign
    return "activation" if product > 0 else "inhibition"


def _path_category(edges: tuple[SignalingEdge, ...]) -> str:
    if any(edge.sign_status is SignStatus.CONFLICTING for edge in edges):
        return "uncertain"
    if any(edge.sign_status is SignStatus.UNSIGNED for edge in edges):
        return "unsigned"
    return path_sign(edges)


@dataclass(slots=True)
class DirectedSignalingGraph:
    dataset: SignalingDataset
    nodes: dict = field(init=False, repr=False)
    outgoing: dict = field(init=False, repr=False)
    incoming: dict = field(init=False, repr=False)
    aliases: dict = field(init=False, repr=False)
    _distance_cache: dict = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.nodes = {node.node_id: node for node in self.dataset.nodes}
        self.outgoing: dict[str, tuple[SignalingEdge, ...]] = {}
        self.incoming: dict[str, tuple[SignalingEdge, ...]] = {}
        outgoing: dict[str, list[SignalingEdge]] = {}
        incoming: dict[str, list[SignalingEdge]] = {}
        aliases: dict[str, set[str]] = {}
        for node in self.dataset.nodes:
            for value in (node.node_id, node.layer1_id, node.symbol, node.source_identifier):
                if value:
                    aliases.setdefault(str(value).strip().casefold(), set()).add(node.node_id)
        for edge in self.dataset.edges:
            outgoing.setdefault(edge.source, []).append(edge)
            incoming.setdefault(edge.target, []).append(edge)
        self.outgoing = {key: tuple(sorted(value, key=lambda edge: edge.target)) for key, value in outgoing.items()}
        self.incoming = {key: tuple(sorted(value, key=lambda edge: edge.source)) for key, value in incoming.items()}
        self.aliases = {key: tuple(sorted(value)) for key, value in aliases.items()}
        self._distance_cache: dict[tuple[str, bool, int], dict[str, int]] = {}

    def resolve_node(self, value: object) -> str | None:
        candidates = self.aliases.get(str(value or "").strip().casefold(), ())
        return candidates[0] if len(candidates) == 1 else None

    def successors(self, node: object) -> tuple[str, ...]:
        resolved = self.resolve_node(node)
        if resolved is None:
            return ()
        return tuple(edge.target for edge in self.outgoing.get(resolved, ()))

    def predecessors(self, node: object) -> tuple[str, ...]:
        resolved = self.resolve_node(node)
        if resolved is None:
            return ()
        return tuple(edge.source for edge in self.incoming.get(resolved, ()))

    def bounded_distances(self, source: object, *, max_depth: int = 4, reverse: bool = False) -> dict[str, int]:
        if max_depth < 0:
            raise ValueError("max_depth cannot be negative")
        resolved = self.resolve_node(source)
        if resolved is None:
            return {}
        key = (resolved, reverse, max_depth)
        if key in self._distance_cache:
            return dict(self._distance_cache[key])
        distances = {resolved: 0}
        queue = deque([resolved])
        adjacency = self.incoming if reverse else self.outgoing
        while queue:
            node = queue.popleft()
            depth = distances[node]
            if depth >= max_depth:
                continue
            for edge in adjacency.get(node, ()):
                neighbor = edge.source if reverse else edge.target
                if neighbor not in distances:
                    distances[neighbor] = depth + 1
                    queue.append(neighbor)
        distances.pop(resolved, None)
        self._distance_cache[key] = dict(distances)
        return distances

    def bounded_paths(
        self, source: object, target: object, *, max_depth: int = 4, max_paths: int = 2048,
    ) -> tuple[tuple[SignalingEdge, ...], ...]:
        """Enumerate only bounded simple paths with a hard result cap."""

        if max_depth < 1 or max_paths < 1:
            return ()
        start, finish = self.resolve_node(source), self.resolve_node(target)
        if start is None or finish is None or start == finish:
            return ()
        # Reverse-distance pruning prevents exploring branches that cannot
        # reach the requested endpoint within the remaining depth budget.
        can_reach_finish = self.bounded_distances(finish, max_depth=max_depth, reverse=True)
        can_reach_finish[finish] = 0
        paths: list[tuple[SignalingEdge, ...]] = []
        queue = deque([(start, (), frozenset({start}))])
        while queue and len(paths) < max_paths:
            node, edges, visited = queue.popleft()
            if len(edges) >= max_depth:
                continue
            for edge in self.outgoing.get(node, ()):
                if edge.target in visited:
                    continue
                candidate = (*edges, edge)
                if edge.target == finish:
                    paths.append(candidate)
                    if len(paths) >= max_paths:
                        break
                else:
                    remaining = max_depth - len(candidate)
                    if can_reach_finish.get(edge.target, max_depth + 1) <= remaining:
                        queue.append((edge.target, candidate, visited | {edge.target}))
        return tuple(paths)

    def summarize_paths(self, source: object, target: object, *, max_depth: int = 4) -> PathSummary:
        paths = self.bounded_paths(source, target, max_depth=max_depth)
        if not paths:
            return PathSummary(None, (), 0, 0, 0, 0, 0, "no_path")
        signs = tuple(_path_category(path) for path in paths)
        positive, negative = signs.count("activation"), signs.count("inhibition")
        unsigned = signs.count("unsigned")
        uncertain = signs.count("uncertain")
        shortest = min(paths, key=lambda item: (len(item), tuple(edge.target for edge in item)))
        example = (shortest[0].source, *(edge.target for edge in shortest))
        if uncertain or unsigned or (positive and negative):
            net = "uncertain"
        elif positive:
            net = "activation"
        elif negative:
            net = "inhibition"
        else:
            net = "uncertain"
        return PathSummary(
            shortest_distance=len(shortest), example_path=tuple(example), path_count=len(paths),
            positive_paths=positive, negative_paths=negative, unsigned_paths=unsigned,
            uncertain_paths=uncertain, net_sign=net,
        )

    def relation(self, source: object, target: object, *, max_depth: int = 4) -> tuple[str, PathSummary]:
        start, finish = self.resolve_node(source), self.resolve_node(target)
        if start is None or finish is None:
            return "unresolved", PathSummary(None, (), 0, 0, 0, 0, 0, "no_path")
        # One cached traversal per perturbation target rejects unrelated pairs
        # before bounded path enumeration. This is the critical dense-graph guard.
        downstream_reachable = finish in self.bounded_distances(start, max_depth=max_depth)
        upstream_reachable = finish in self.bounded_distances(start, max_depth=max_depth, reverse=True)
        empty = PathSummary(None, (), 0, 0, 0, 0, 0, "no_path")
        if not downstream_reachable and not upstream_reachable:
            return "no directed path", empty
        downstream = self.summarize_paths(start, finish, max_depth=max_depth) if downstream_reachable else empty
        upstream = self.summarize_paths(finish, start, max_depth=max_depth) if upstream_reachable else empty
        if downstream_reachable and upstream_reachable:
            return "bidirectionally reachable", downstream
        if downstream_reachable:
            if self.nodes[start].layer1_id is None or self.nodes[finish].layer1_id is None:
                return "signaling-only relationship", downstream
            relation = "direct downstream" if downstream.shortest_distance == 1 else "indirect downstream"
            return relation, downstream
        if upstream_reachable:
            return "upstream", upstream
        return "no directed path", empty

    def path_evidence(self, path: tuple[str, ...]) -> tuple[tuple[str, ...], tuple[str, ...]]:
        resources: set[str] = set()
        references: set[str] = set()
        for source, target in zip(path, path[1:]):
            for edge in self.outgoing.get(source, ()):
                if edge.target == target:
                    resources.update(edge.resources)
                    references.update(edge.references)
        return tuple(sorted(resources)), tuple(sorted(references))
