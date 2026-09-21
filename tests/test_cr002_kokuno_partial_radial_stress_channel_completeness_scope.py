from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.audit_kokuno_partial_radial_stress_channel_completeness_scope import (
    audit_scope,
    enforce_scope,
    load_scope,
    validate_scope,
)


def _scope() -> dict:
    return load_scope()


def _assert_rejected(payload: dict, token: str) -> None:
    errors = validate_scope(payload)
    assert errors
    assert any(token in item for item in errors)


def test_registered_scope_passes_and_reports_missing_radial_channel() -> None:
    scope = _scope()
    enforce_scope(scope)
    receipt = audit_scope()
    assert receipt["status"] == "pass"
    assert receipt["materialized_stress_channels"] == ["theta", "axial"]
    assert receipt["missing_radial_stress_channel"] is True
    assert receipt["full_three_component_stress_materialized"] is False
    assert receipt["complete_ns_defect_materialized"] is False
    assert receipt["cartesian_correction_velocity_materialized"] is False
    assert receipt["kokuno_velocity_export_ready"] is False
    assert receipt["canonical_eq45_velocity_export_ready"] is True
    assert receipt["pde_validated"] is False


def test_radial_mean_cannot_be_laundered_into_inverted_radial_stress() -> None:
    scope = _scope()
    scope["upstream_agent3"]["radial_cylindrical_stress_inverted"] = True
    _assert_rejected(scope, "radial_stress_inversion_laundered")

    scope = _scope()
    scope["governed_distinctions"]["recorded_radial_mean_is_inverted_radial_stress"] = True
    _assert_rejected(scope, "recorded_radial_mean_is_inverted_radial_stress_laundered")


def test_two_channels_cannot_be_promoted_to_full_three_component_stress() -> None:
    scope = _scope()
    scope["governed_distinctions"][
        "theta_axial_two_channel_stress_is_full_three_component_cylindrical_stress"
    ] = True
    _assert_rejected(scope, "full_three_component_cylindrical_stress_laundered")

    scope = _scope()
    scope["readiness"]["kokuno_full_three_component_cylindrical_stress_materialized"] = True
    _assert_rejected(scope, "kokuno_full_three_component_cylindrical_stress_materialized")


def test_a4_two_channel_audit_cannot_be_promoted_to_radial_or_full_stress_validation() -> None:
    scope = _scope()
    scope["upstream_agent4"]["radial_cylindrical_stress_independently_audited"] = True
    _assert_rejected(scope, "a4_radial_stress_audit_laundered")

    scope = _scope()
    scope["governed_distinctions"][
        "theta_axial_a4_audit_is_full_three_component_stress_validation"
    ] = True
    _assert_rejected(scope, "full_three_component_stress_validation_laundered")


@pytest.mark.parametrize(
    ("key", "token"),
    [
        ("kokuno_complete_ns_defect_materialized", "complete_ns_defect"),
        ("kokuno_complete_correction_target_materialized", "complete_correction_target"),
        ("kokuno_cartesian_correction_velocity_materialized", "cartesian_correction_velocity"),
        ("kokuno_velocity_export_ready", "kokuno_velocity_export_ready"),
        ("kokuno_visual_correspondence_verified", "visual_correspondence"),
        ("pde_validated", "pde_validated"),
        ("paper_exact", "paper_exact"),
        ("openai_field_identified", "openai_field_identified"),
    ],
)
def test_scoped_stress_cannot_promote_scientific_or_delivery_states(
    key: str, token: str
) -> None:
    scope = _scope()
    scope["readiness"][key] = True
    _assert_rejected(scope, token)


def test_complete_defect_correction_and_pde_implications_remain_false() -> None:
    keys = (
        "scoped_two_channel_stress_is_complete_ns_defect",
        "scoped_two_channel_stress_is_authorized_correction_target",
        "scoped_two_channel_stress_is_cartesian_correction_velocity",
        "scoped_two_channel_stress_implies_pde_validation",
    )
    for key in keys:
        scope = _scope()
        scope["governed_distinctions"][key] = True
        _assert_rejected(scope, key)


def test_pending_missing_radial_completion_requirement_cannot_be_dropped() -> None:
    scope = _scope()
    scope["pending_requirements"].remove(
        "materialize_or_mathematically_eliminate_the_missing_radial_cylindrical_stress_channel"
    )
    _assert_rejected(scope, "pending_requirements_dropped")


def test_mechanics_witness_is_projection_only_and_not_source_data() -> None:
    scope = _scope()
    witness = scope["mechanics_witness"]
    assert witness["vector_a_radial_theta_axial"][1:] == witness[
        "vector_b_radial_theta_axial"
    ][1:]
    assert witness["vector_a_radial_theta_axial"][0] != witness[
        "vector_b_radial_theta_axial"
    ][0]
    assert witness["not_kokuno_or_openai_data"] is True

    mutated = copy.deepcopy(scope)
    mutated["mechanics_witness"]["classification"] = "public_source_fact"
    _assert_rejected(mutated, "mechanics_witness_source_laundered")


def test_cr001_threshold_force_and_collapse_contract_is_immutable() -> None:
    scope = _scope()
    scope["cr001"]["momentum_max_gate"] = 2.0e-3
    _assert_rejected(scope, "cr001_contract_drift")

    scope = _scope()
    scope["cr001"]["residual_defined_pointwise_free_force_allowed"] = True
    _assert_rejected(scope, "cr001_contract_drift")

    scope = _scope()
    scope["cr001"]["collapsed_candidate_allowed"] = True
    _assert_rejected(scope, "cr001_contract_drift")

    scope = _scope()
    scope["cr001"]["post_hoc_threshold_relaxation_allowed"] = True
    _assert_rejected(scope, "cr001_contract_drift")


def test_kokuno_missing_channel_cannot_downgrade_independent_eq45_delivery() -> None:
    scope = _scope()
    scope["canonical_eq45_delivery"]["velocity_export_ready"] = False
    _assert_rejected(scope, "eq45_delivery_incorrectly_downgraded")


def test_parent_and_upstream_identities_are_pinned() -> None:
    scope = _scope()
    scope["parent_a5"]["head"] = "0" * 40
    _assert_rejected(scope, "parent_a5_identity_drift")

    scope = _scope()
    scope["upstream_agent3"]["source_blob"] = "0" * 40
    _assert_rejected(scope, "a3_source_blob_drift")

    scope = _scope()
    scope["upstream_agent4"]["source_blob"] = "0" * 40
    _assert_rejected(scope, "a4_source_blob_drift")
