"""Cross-layer Phase 6 guards for the zero-external-call result path."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

from src.research.provider_cache import CacheFirstProviderClient
from src.research.provider_runtime import ProviderExecutor
from src.research.query_planner import QueryPlanner
from tests.test_research_service import _service, _snapshot


def test_offline_bundle_and_query_plan_never_enter_l3_provider_execution() -> None:
    with patch.object(
        ProviderExecutor,
        "execute",
        side_effect=AssertionError("L3 provider execution is forbidden on result path"),
    ), patch.object(
        CacheFirstProviderClient,
        "execute",
        side_effect=AssertionError("cache/provider client is explicit-only"),
    ):
        bundle = _service().build_offline_bundle(_snapshot())
        plan = QueryPlanner().plan(
            bundle,
            expected_snapshot_id=bundle.snapshot_id,
        )

    assert bundle.healthy
    assert plan.snapshot_id == bundle.snapshot_id
    assert plan.provider_request_count <= bundle.config.quick_context_budget.max_provider_requests


def test_application_facade_import_does_not_load_external_runtime_or_cache() -> None:
    project_root = Path(__file__).parents[1]
    script = """
import sys
import src.research.facade
forbidden = {
    'src.research.provider_runtime',
    'src.research.provider_cache',
    'src.research.query_planner',
}
loaded = sorted(forbidden.intersection(sys.modules))
if loaded:
    raise SystemExit('unexpected critical-path imports: ' + ','.join(loaded))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=project_root,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout

