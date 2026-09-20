from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_inner_divergence_norm_scope import (
    ScopeAuditError,
    audit,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs" / "kokuno_inner_divergence_norm_scope.json"


def _mutated_contract(tmp_path: Path, mutate) -> Path:
    payload = json.loads(CONTRACT.read_text(encoding="utf-8"))
    mutate(payload)
    path = tmp_path / "mutated_scope.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def test_scope_audit_accepts_frozen_contract() -> None:
    report = audit()
    assert report["ok"] is True
    assert report["same_numeric_threshold"] is True
    assert report["metrics_equivalent"] is False
    assert report["scientific_promotions"] == []


@pytest.mark.parametrize(
    "mutate",
    [
        lambda payload: payload["metric_scope"].__setitem__("metrics_equivalent", True),
        lambda payload: payload["metric_scope"].__setitem__(
            "upstream_pass_can_discharge_cr001_divergence_L2", True
        ),
        lambda payload: payload["metric_scope"].__setitem__(
            "upstream_pass_can_discharge_complete_candidate_divergence_acceptance", True
        ),
        lambda payload: payload["metric_scope"].__setitem__(
            "upstream_pass_can_discharge_pde_validation", True
        ),
        lambda payload: payload["canonical_cr001_snapshot"].__setitem__(
            "divergence_L2", 2.0e-5
        ),
        lambda payload: payload["canonical_cr001_snapshot"].__setitem__(
            "pde_residual_L2", 2.0e-3
        ),
        lambda payload: payload["status_separation"].__setitem__(
            "pde_validated", True
        ),
        lambda payload: payload["status_separation"].__setitem__(
            "openai_field_identified", True
        ),
        lambda payload: payload["status_separation"].__setitem__(
            "callable_velocity_delivery_blocked_by_pde_pending", True
        ),
        lambda payload: payload["classification"].__setitem__(
            "public_source_fact", ["the 384-sample RMS is an OpenAI/source validation norm"]
        ),
    ],
)
def test_scope_audit_fails_closed_on_claim_or_constraint_mutation(
    tmp_path: Path, mutate
) -> None:
    path = _mutated_contract(tmp_path, mutate)
    with pytest.raises(ScopeAuditError):
        audit(path)
