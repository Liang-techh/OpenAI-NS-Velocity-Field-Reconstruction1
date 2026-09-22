from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_relative_swirl_delta_target_scope import (
    ScopeAuditError,
    audit_repository,
    audit_upstream_source,
    autonomous_delta_nonuniqueness_witness,
    validate_contract_payload,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs/kokuno_relative_swirl_delta_target_scope.json"
UPSTREAM = (
    ROOT
    / "src/openai_ns_reconstruction/kokuno_current_relative_swirl_angular_correction.py"
)


def _payload() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_exact_repository_scope_audit_is_fail_closed_and_non_promoting() -> None:
    receipt = audit_repository(ROOT)
    assert receipt["exact_dependency_head"] == "6ddf35b9dc726c720033ee44161e29efc77d755e"
    assert receipt["exact_dependency_source_blob"] == "37e6b0814fbad14084a8c8c4815122111252ac9b"
    assert receipt["canonical_constraints_blob"] == "6c559e42895a606e2ef025ade4cb448966d75814"
    assert receipt["upstream_scope"] == {
        "delta_target_materialized": True,
        "absolute_base_materialized": False,
        "total_target_materialized": False,
        "explicit_base_target_required": True,
        "explicit_base_target_eta_required": True,
    }
    assert receipt["scientific_state_promoted"] is False
    assert receipt["candidate_bytes_changed_by_this_audit"] is False
    assert receipt["cr001_threshold_changed_by_this_audit"] is False


def test_autonomous_nonuniqueness_witness_same_delta_has_distinct_totals() -> None:
    witness = autonomous_delta_nonuniqueness_witness()
    assert witness["role"] == "autonomous_mechanics_only"
    assert witness["delta"] != 0.0
    assert witness["base_a"] != witness["base_b"]
    assert witness["total_a"] != witness["total_b"]
    assert witness["same_delta_distinct_totals"] is True
    assert witness["source_evidence"] is False
    assert witness["openai_numerical_evidence"] is False
    assert witness["pde_evidence"] is False


def test_contract_rejects_delta_to_absolute_or_total_target_promotion() -> None:
    payload = _payload()
    for key in (
        "imported_base_absolute_target_materialized",
        "current_absolute_entering_I_materialized",
        "current_total_relative_swirl_target_materialized",
        "current_cartesian_relative_swirl_composed",
        "kokuno_velocity_export_ready",
        "visual_correspondence_verified",
        "heldout_complete_ns_residual_assessed",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        mutated = copy.deepcopy(payload)
        mutated["materialization_state"][key] = True
        with pytest.raises(ScopeAuditError):
            validate_contract_payload(mutated)


def test_contract_rejects_evidence_firewall_relaxation() -> None:
    payload = _payload()
    for key in payload["evidence_transfer_rules"]:
        mutated = copy.deepcopy(payload)
        mutated["evidence_transfer_rules"][key] = False
        with pytest.raises(ScopeAuditError):
            validate_contract_payload(mutated)


def test_contract_rejects_cr001_threshold_or_force_drift() -> None:
    payload = _payload()

    mutated = copy.deepcopy(payload)
    mutated["canonical_cr001"]["momentum_max"] = 0.002
    with pytest.raises(ScopeAuditError):
        validate_contract_payload(mutated)

    mutated = copy.deepcopy(payload)
    mutated["canonical_cr001"]["divergence_L2"] = 2e-5
    with pytest.raises(ScopeAuditError):
        validate_contract_payload(mutated)

    mutated = copy.deepcopy(payload)
    mutated["canonical_cr001"]["forcing_mode"] = "residual_defined"
    with pytest.raises(ScopeAuditError):
        validate_contract_payload(mutated)


def test_contract_requires_four_distinct_nonempty_provenance_classes() -> None:
    payload = _payload()
    mutated = copy.deepcopy(payload)
    mutated["provenance_classes"].pop("pending_unknown")
    with pytest.raises(ScopeAuditError):
        validate_contract_payload(mutated)

    mutated = copy.deepcopy(payload)
    mutated["provenance_classes"]["public_source_fact"] = []
    with pytest.raises(ScopeAuditError):
        validate_contract_payload(mutated)


def test_upstream_truth_boundary_rejects_total_target_promotion() -> None:
    source = UPSTREAM.read_text(encoding="utf-8")
    mutated = source.replace(
        '"current_relative_swirl_total_target_materialized": False,',
        '"current_relative_swirl_total_target_materialized": True,',
        1,
    )
    assert mutated != source
    with pytest.raises(ScopeAuditError):
        audit_upstream_source(mutated)


def test_upstream_helper_cannot_default_missing_base_to_zero() -> None:
    source = UPSTREAM.read_text(encoding="utf-8")
    mutated = source.replace(
        "base_target: Any,\n        base_target_eta: Any,",
        "base_target: Any = 0.0,\n        base_target_eta: Any = 0.0,",
        1,
    )
    assert mutated != source
    with pytest.raises(ScopeAuditError):
        audit_upstream_source(mutated)


def test_upstream_helper_keeps_explicit_base_provenance_warning() -> None:
    source = UPSTREAM.read_text(encoding="utf-8")
    mutated = source.replace(
        "There is intentionally no default for either base quantity.",
        "The base quantity may be inferred automatically.",
        1,
    )
    assert mutated != source
    with pytest.raises(ScopeAuditError):
        audit_upstream_source(mutated)


def test_mechanics_witness_refuses_degenerate_inputs() -> None:
    with pytest.raises(ValueError):
        autonomous_delta_nonuniqueness_witness(delta=0.0)
    with pytest.raises(ValueError):
        autonomous_delta_nonuniqueness_witness(base_a=1.0, base_b=1.0)
