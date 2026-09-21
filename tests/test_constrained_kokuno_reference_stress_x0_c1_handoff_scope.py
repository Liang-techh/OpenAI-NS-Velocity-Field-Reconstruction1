from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_reference_stress_x0_c1_handoff_scope import (
    SCOPE_PATH,
    audit_repository,
    audit_scope,
    mechanics_witness,
)


ROOT = Path(__file__).resolve().parents[1]


def _scope() -> dict:
    return json.loads((ROOT / SCOPE_PATH).read_text(encoding="utf-8"))


def test_scope_audits_exact_parent_constraints_and_fail_closed_state() -> None:
    receipt = audit_repository(ROOT)
    assert receipt["audited_parent_head"] == "8bb1139f2b2d9824fd3bbac693ffa4bf2e71dea2"
    assert receipt["audited_parent_source_blob"] == "733f76ab40bd48d916301ea9ed2d4c8f9e1eff09"
    assert receipt["canonical_constraints_blob"] == "6c559e42895a606e2ef025ade4cb448966d75814"
    assert receipt["x0_value_handoff_registered"] is True
    assert receipt["right_branch_derivative_audit_registered"] is True
    assert receipt["x0_c1_reference_stress_handoff_verified"] is False
    assert receipt["source_stress_free_derivative_compatibility_verified"] is False
    assert receipt["source_prepared_reference_stress_identity_verified"] is False
    assert receipt["pde_validated"] is False
    assert receipt["canonical_eq45_velocity_delivery_remains_independent"] is True
    assert len(receipt["receipt_sha256"]) == 64


def test_mechanics_witness_only_proves_value_match_does_not_imply_slope_match() -> None:
    witness = mechanics_witness(_scope())
    assert witness["equal_value"] is True
    assert witness["different_derivative"] is True
    assert witness["left_value_at_x0"] == pytest.approx(witness["right_value_at_x0"], abs=1e-15)
    assert witness["value_jump"] == pytest.approx(0.0, abs=1e-15)
    assert witness["left_derivative"] == pytest.approx(-0.2)
    assert witness["right_derivative"] == pytest.approx(0.4)
    assert witness["derivative_jump"] == pytest.approx(0.6)


def test_rejects_value_handoff_laundered_into_c1_splice() -> None:
    scope = _scope()
    scope["governed_distinctions"]["x0_value_equality_implies_c1_splice"] = True
    with pytest.raises(ValueError, match="C1/source-stress promotion"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_right_branch_audit_laundered_into_left_right_derivative_match() -> None:
    scope = _scope()
    scope["governed_distinctions"][
        "right_branch_derivative_consistency_implies_left_right_derivative_match"
    ] = True
    with pytest.raises(ValueError, match="C1/source-stress promotion"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_source_derivative_or_fixed_kappa_promotion() -> None:
    for key in (
        "source_stress_free_derivative_compatibility_verified",
        "source_prepared_reference_stress_identity_verified",
        "autonomous_reference_stress_authorized_as_actual_fixed_kappa_stress",
        "autonomous_reference_stress_authorized_as_agent3_correction_target",
    ):
        scope = _scope()
        scope["governed_distinctions"][key] = True
        with pytest.raises(ValueError, match="C1/source-stress promotion"):
            audit_scope(scope, repository_root=ROOT)


def test_rejects_dropping_pending_left_derivative_compatibility() -> None:
    scope = _scope()
    scope["provenance_classes"]["pending_unknown"] = [
        item
        for item in scope["provenance_classes"]["pending_unknown"]
        if "left radial derivative" not in item
    ]
    with pytest.raises(ValueError, match="pending left/right derivative compatibility"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_mechanics_witness_laundered_into_public_source_fact() -> None:
    scope = _scope()
    scope["provenance_classes"]["public_source_fact"].append(
        "The equal-value/different-slope mechanics witness is a public source fact."
    )
    with pytest.raises(ValueError, match="public-source fact laundering"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_cr001_gate_relaxation_and_shortcuts() -> None:
    scope = _scope()
    scope["canonical_constraints"]["momentum_l2_gate"] = 2.0e-3
    with pytest.raises(ValueError, match="CR001 constraint drift"):
        audit_scope(scope, repository_root=ROOT)

    scope = _scope()
    scope["canonical_constraints"]["divergence_max_gate"] = 2.0e-5
    with pytest.raises(ValueError, match="CR001 constraint drift"):
        audit_scope(scope, repository_root=ROOT)

    scope = _scope()
    scope["canonical_constraints"]["residual_defined_pointwise_force_forbidden"] = False
    with pytest.raises(ValueError, match="pointwise forcing"):
        audit_scope(scope, repository_root=ROOT)

    scope = _scope()
    scope["canonical_constraints"]["collapsed_candidate_forbidden"] = False
    with pytest.raises(ValueError, match="candidate collapse"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_delivery_pde_visual_or_exactness_coupling() -> None:
    for key, value in (
        ("canonical_eq45_velocity_delivery_remains_independent", False),
        ("kokuno_velocity_export_ready", True),
        ("visual_correspondence_verified", True),
        ("pde_validated", True),
        ("paper_exact", True),
        ("openai_field_identified", True),
    ):
        scope = _scope()
        scope["delivery_and_scientific_state"][key] = value
        with pytest.raises(ValueError, match="state coupling or promotion"):
            audit_scope(scope, repository_root=ROOT)


def test_rejects_upstream_x0_audit_scope_rewrite() -> None:
    scope = _scope()
    scope["upstream_pins"]["agent4_reference_stress_audit"]["exact_x0_check"] = (
        "C1 derivative compatibility proved"
    )
    with pytest.raises(ValueError, match="X0 audit scope drift"):
        audit_scope(scope, repository_root=ROOT)
