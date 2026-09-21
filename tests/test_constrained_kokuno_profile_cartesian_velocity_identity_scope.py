from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_profile_cartesian_velocity_identity_scope import (
    SCOPE_PATH,
    audit_repository,
    audit_scope,
    mechanics_witness,
)


ROOT = Path(__file__).resolve().parents[1]


def _scope() -> dict:
    return json.loads((ROOT / SCOPE_PATH).read_text(encoding="utf-8"))


def test_scope_audits_exact_parent_constraints_and_fail_closed_delivery_state() -> None:
    receipt = audit_repository(ROOT)
    assert receipt["audited_parent_head"] == "744a4fbc642b861f61da91de713dd06adc357f09"
    assert receipt["audited_parent_source_blob"] == "502e64e97903853a16025815eabe31bf8b2111d0"
    assert receipt["canonical_constraints_blob"] == "6c559e42895a606e2ef025ade4cb448966d75814"
    assert receipt["source_coordinate_profile_executable"] is True
    assert receipt["profile_save_load_semantic_identity_registered"] is True
    assert receipt["profile_radial_derivative_audit_registered"] is True
    assert receipt["kokuno_cartesian_velocity_identity_materialized"] is False
    assert receipt["kokuno_velocity_export_ready"] is False
    assert receipt["pde_validated"] is False
    assert receipt["canonical_eq45_velocity_delivery_remains_independent"] is True
    assert len(receipt["receipt_sha256"]) == 64


def test_mechanics_witness_only_proves_profile_packet_needs_a_bound_component_map() -> None:
    witness = mechanics_witness(_scope())
    assert witness["same_profile_packet"] is True
    assert witness["different_cartesian_velocity"] is True
    assert witness["cartesian_a"] == pytest.approx([0.2, 0.3, 1.0], abs=1e-15)
    assert witness["cartesian_b"] == pytest.approx([0.2, 1.0, 0.3], abs=1e-15)
    assert witness["cartesian_difference_l2"] > 0.0


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("profile_identity_implies_cartesian_velocity_identity", True),
        ("profile_semantic_sha_is_velocity_candidate_sha", True),
        ("profile_radial_derivatives_imply_cartesian_velocity_jacobian", True),
        ("profile_values_alone_fix_unique_cartesian_component_map", True),
        ("kokuno_cartesian_coordinate_time_pullback_materialized", True),
        ("kokuno_cartesian_component_map_bound_to_profile_artifact", True),
        ("kokuno_axis_regular_cartesian_velocity_materialized", True),
        ("kokuno_global_outer_join_materialized", True),
        ("kokuno_velocity_provider_save_load_identity_materialized", True),
        ("kokuno_complete_candidate_api_ready", True),
        ("kokuno_velocity_export_ready", True),
    ],
)
def test_rejects_profile_to_cartesian_promotion_without_materialized_contract(
    key: str, value: bool
) -> None:
    scope = _scope()
    scope["governed_distinctions"][key] = value
    with pytest.raises(ValueError, match="profile-to-Cartesian promotion"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_dropping_physical_pullback_or_velocity_save_load_requirements() -> None:
    scope = _scope()
    scope["provenance_classes"]["pending_unknown"] = [
        item
        for item in scope["provenance_classes"]["pending_unknown"]
        if "physical inputs (x,y,z,t)" not in item
    ]
    with pytest.raises(ValueError, match="physical coordinate/time pullback"):
        audit_scope(scope, repository_root=ROOT)

    scope = _scope()
    scope["provenance_classes"]["pending_unknown"] = [
        item
        for item in scope["provenance_classes"]["pending_unknown"]
        if "actual Kokuno Cartesian velocity provider" not in item
    ]
    with pytest.raises(ValueError, match="Cartesian velocity save/load identity"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_delivery_requirement_deletion_or_rewrite() -> None:
    scope = _scope()
    scope["cartesian_velocity_delivery_requirements"].pop(2)
    with pytest.raises(ValueError, match="delivery requirements drift"):
        audit_scope(scope, repository_root=ROOT)

    scope = _scope()
    scope["cartesian_velocity_delivery_requirements"][0] = (
        "F/U/E profile configuration alone is the public Cartesian velocity API."
    )
    with pytest.raises(ValueError, match="delivery requirements drift"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_mechanics_witness_laundered_into_public_source_fact() -> None:
    scope = _scope()
    scope["provenance_classes"]["public_source_fact"].append(
        "The two-component-map mechanics counterexample is a public source fact."
    )
    with pytest.raises(ValueError, match="public-source fact laundering"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_profile_audit_relabelled_as_complete_ns_residual() -> None:
    scope = _scope()
    scope["upstream_pins"]["agent4_profile_audit"]["source_coordinate_profile_audit_only"] = False
    with pytest.raises(ValueError, match="Agent-4 audit boundary drift"):
        audit_scope(scope, repository_root=ROOT)

    scope = _scope()
    scope["upstream_pins"]["agent4_profile_audit"]["complete_ns_residual_audit"] = True
    with pytest.raises(ValueError, match="Agent-4 audit boundary drift"):
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


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("canonical_eq45_velocity_delivery_remains_independent", False),
        ("canonical_eq45_velocity_export_ready", False),
        ("kokuno_velocity_export_ready", True),
        ("visual_correspondence_verified", True),
        ("pde_validated", True),
        ("paper_exact", True),
        ("openai_field_identified", True),
        ("blowup_proved", True),
    ],
)
def test_rejects_delivery_pde_visual_or_exactness_coupling(key: str, value: bool) -> None:
    scope = _scope()
    scope["delivery_and_scientific_state"][key] = value
    with pytest.raises(ValueError, match="state coupling or promotion"):
        audit_scope(scope, repository_root=ROOT)


def test_rejects_profile_source_exactness_laundering() -> None:
    scope = _scope()
    scope["upstream_pins"]["agent1_fixed_kappa_profile"]["source_exact_fixed_kappa_continuation"] = True
    with pytest.raises(ValueError, match="source-exactness boundary drift"):
        audit_scope(scope, repository_root=ROOT)
