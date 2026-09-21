from __future__ import annotations

import copy
import json
from pathlib import Path

from openai_ns_reconstruction.audit_kokuno_current_cartesian_partial_velocity_api_scope import (
    CONSTRAINTS_PATH,
    PROJECT_STATUS_PATH,
    SCOPE_PATH,
    audit_scope,
)


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_current_partial_cartesian_velocity_scope_is_consistent() -> None:
    assert audit_scope() == []


def test_callable_does_not_promote_global_api_totality() -> None:
    scope = _load(SCOPE_PATH)
    scope["delivery_contract"]["callable_implies_global_api_totality"] = True
    errors = audit_scope(scope=scope)
    assert "delivery_forbidden_promotion_callable_implies_global_api_totality" in errors


def test_save_load_does_not_promote_global_api_totality() -> None:
    scope = _load(SCOPE_PATH)
    scope["delivery_contract"]["save_load_implies_global_api_totality"] = True
    errors = audit_scope(scope=scope)
    assert "delivery_forbidden_promotion_save_load_implies_global_api_totality" in errors


def test_semantic_identity_does_not_promote_global_api_totality() -> None:
    scope = _load(SCOPE_PATH)
    scope["delivery_contract"]["semantic_identity_implies_global_api_totality"] = True
    errors = audit_scope(scope=scope)
    assert "delivery_forbidden_promotion_semantic_identity_implies_global_api_totality" in errors


def test_partial_kokuno_provider_cannot_be_marked_export_ready() -> None:
    scope = _load(SCOPE_PATH)
    scope["delivery_contract"]["kokuno_velocity_export_ready"] = True
    errors = audit_scope(scope=scope)
    assert "delivery_forbidden_promotion_kokuno_velocity_export_ready" in errors


def test_beyond_xh_cannot_be_promoted_without_materialization() -> None:
    scope = _load(SCOPE_PATH)
    scope["truth_boundary"]["velocity_beyond_xh_materialized"] = True
    errors = audit_scope(scope=scope)
    assert "truth_boundary_promotion_velocity_beyond_xh_materialized" in errors


def test_mechanics_witness_cannot_be_laundered_into_public_source_fact() -> None:
    scope = _load(SCOPE_PATH)
    scope["mechanics_witness"]["classification"] = "public_source_fact"
    errors = audit_scope(scope=scope)
    assert "mechanics_witness_source_laundering" in errors


def test_pde_pending_must_not_become_callable_delivery_blocker() -> None:
    scope = _load(SCOPE_PATH)
    scope["delivery_contract"]["pde_pending_blocks_callable_velocity_delivery"] = True
    errors = audit_scope(scope=scope)
    assert "delivery_forbidden_promotion_pde_pending_blocks_callable_velocity_delivery" in errors


def test_canonical_eq45_delivery_remains_independent() -> None:
    status = _load(PROJECT_STATUS_PATH)
    status["states"]["velocity_export_ready"] = False
    errors = audit_scope(project_status=status)
    assert "canonical_eq45_velocity_export_ready_lost" in errors


def test_cr001_momentum_threshold_relaxation_is_rejected() -> None:
    constraints = _load(CONSTRAINTS_PATH)
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    errors = audit_scope(constraints=constraints)
    assert "cr001_momentum_gate_drift" in errors


def test_cr001_divergence_threshold_relaxation_is_rejected() -> None:
    constraints = _load(CONSTRAINTS_PATH)
    constraints["validation"]["thresholds"]["divergence_L2"] = 1.0e-4
    errors = audit_scope(constraints=constraints)
    assert "cr001_divergence_gate_drift" in errors


def test_residual_defined_free_force_firewall_is_required() -> None:
    constraints = _load(CONSTRAINTS_PATH)
    constraints["forcing"]["restriction"] = "Only a,c may be fitted."
    errors = audit_scope(constraints=constraints)
    assert "cr001_free_force_firewall_lost" in errors


def test_candidate_collapse_firewall_is_required() -> None:
    constraints = _load(CONSTRAINTS_PATH)
    constraints["nontriviality"]["enforcement"] = "normalize initial velocity parameters"
    errors = audit_scope(constraints=constraints)
    assert "cr001_collapse_firewall_lost" in errors


def test_visual_or_pde_promotion_is_rejected() -> None:
    for key in ("visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        scope = _load(SCOPE_PATH)
        scope["truth_boundary"][key] = True
        errors = audit_scope(scope=scope)
        assert f"truth_boundary_promotion_{key}" in errors


def test_four_way_provenance_partition_must_remain_explicit() -> None:
    scope = _load(SCOPE_PATH)
    scope["provenance"]["other"] = ["ambiguous"]
    errors = audit_scope(scope=scope)
    assert "four_way_provenance_partition_drift" in errors
