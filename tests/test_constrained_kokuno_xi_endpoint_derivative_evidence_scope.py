from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_xi_endpoint_derivative_evidence_scope import (
    SCOPE_PATH,
    audit_repository,
    audit_scope,
    mechanics_witness,
)


ROOT = Path(__file__).resolve().parents[1]


def _scope() -> dict:
    return json.loads((ROOT / SCOPE_PATH).read_text(encoding="utf-8"))


def test_registered_scope_audits_fail_closed() -> None:
    receipt = audit_repository(ROOT)
    assert receipt["interior_independent_derivative_audit_registered"] is True
    assert receipt["Xi_public_handoff_replay_registered"] is True
    assert receipt["Xi_independent_endpoint_derivative_verified"] is False
    assert receipt["actual_candidate_Xi_derivative_handoff_verified"] is False
    assert receipt["kokuno_velocity_export_ready"] is False
    assert receipt["pde_validated"] is False
    assert receipt["canonical_eq45_velocity_delivery_remains_independent"] is True
    assert len(receipt["receipt_sha256"]) == 64


def test_mechanics_witness_separates_replay_from_endpoint_accuracy() -> None:
    replay = mechanics_witness(_scope())
    assert replay["production_handoff_replay_abs_error"] == pytest.approx(0.0, abs=1e-15)
    assert replay["independent_three_point_backward_derivative"] == pytest.approx(2.0, abs=1e-12)
    assert replay["production_vs_independent_abs_error"] == pytest.approx(0.1, abs=1e-12)
    assert replay["replay_consistency_implies_endpoint_accuracy"] is False


def test_rejects_endpoint_replay_promoted_to_independent_validation() -> None:
    scope = copy.deepcopy(_scope())
    scope["governed_distinctions"][
        "Xi_handoff_derivative_replay_is_independent_endpoint_derivative_validation"
    ] = True
    with pytest.raises(ValueError, match="unsupported Xi endpoint evidence promotion"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_independent_endpoint_derivative_promotion() -> None:
    scope = copy.deepcopy(_scope())
    scope["governed_distinctions"]["Xi_radial_derivatives_independently_verified"] = True
    with pytest.raises(ValueError, match="unsupported Xi endpoint evidence promotion"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_interior_audit_laundered_into_endpoint_coverage() -> None:
    scope = copy.deepcopy(_scope())
    scope["registered_upstream_evidence"]["agent4_independent_audit"]["independent_sample_x_range"] = [
        100.35,
        110.0,
    ]
    with pytest.raises(ValueError, match="A4 independent X sampling range drift"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_unregistered_independent_Xi_operator() -> None:
    scope = copy.deepcopy(_scope())
    scope["registered_upstream_evidence"]["agent4_independent_audit"][
        "independent_endpoint_derivative_operator_applied_at_Xi"
    ] = True
    with pytest.raises(ValueError, match="unregistered independent Xi endpoint operator"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_mechanics_witness_as_public_source_fact() -> None:
    scope = copy.deepcopy(_scope())
    scope["provenance_classes"]["public_source_fact"] = [
        "The shared-bug mechanics_witness proves the Kokuno Xi derivative is wrong."
    ]
    with pytest.raises(ValueError, match="public-source fact"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_dropping_pending_endpoint_safe_validation() -> None:
    scope = copy.deepcopy(_scope())
    scope["provenance_classes"]["pending_unknown"] = [
        item
        for item in scope["provenance_classes"]["pending_unknown"]
        if "endpoint-safe" not in item
    ]
    with pytest.raises(ValueError, match="pending independent Xi endpoint"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_cr001_gate_relaxation() -> None:
    scope = copy.deepcopy(_scope())
    scope["canonical_constraints"]["momentum_max_gate"] = 1.0e-2
    with pytest.raises(ValueError, match="CR001 constraint drift: momentum_max_gate"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_free_residual_defined_force() -> None:
    scope = copy.deepcopy(_scope())
    scope["canonical_constraints"]["residual_defined_pointwise_force_forbidden"] = False
    with pytest.raises(ValueError, match="free residual-defined forcing"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_delivery_or_pde_state_coupling() -> None:
    scope = copy.deepcopy(_scope())
    scope["delivery_and_scientific_state"]["pde_validated"] = True
    with pytest.raises(ValueError, match="delivery/PDE/visual/exactness state coupling"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_actual_candidate_handoff_promotion() -> None:
    scope = copy.deepcopy(_scope())
    scope["governed_distinctions"]["actual_candidate_Xi_derivative_handoff_verified"] = True
    with pytest.raises(ValueError, match="unsupported Xi endpoint evidence promotion"):
        audit_scope(scope, repository_root=ROOT)
