from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_axial_cap_poloidal_representation_governance import (
    audit_axial_cap_poloidal_scope,
    load_constraints,
    load_project_status,
    load_scope,
)


def test_current_axial_cap_scope_passes():
    report = audit_axial_cap_poloidal_scope()
    assert report["status"] == "governance_pass"
    assert report["mode_name"] == "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL"
    assert report["upstream_pr"] == 287
    assert report["canonical_thresholds_unchanged"] is True
    assert report["diagnostic_amplitude_selected"] is False
    assert report["nonlinear_reach_is_capacity_only"] is True
    assert report["axial_cap_capacity_status"] == "integrated_prematerialization"
    assert report["live_next_integration_task_preserved"] is True
    assert report["velocity_export_ready"] is True
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False


def test_rejects_source_laundering_and_bound_promotion():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["source_classification"]["axial_cap_band_and_envelope"] = "public_source_fact"
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(scope=mutated)

    for key in (
        "diagnostic_coefficients_select_value",
        "diagnostic_coefficients_define_materialization_bound",
        "inherited_profile_limit_is_materialization_bound",
    ):
        mutated = deepcopy(scope)
        mutated["capacity_evidence_only"][key] = True
        with pytest.raises(ValueError):
            audit_axial_cap_poloidal_scope(scope=mutated)


def test_rejects_nonlinear_morphology_claim_laundering():
    scope = load_scope()
    for key in (
        "zero_plus_minus_span_means_zero_capacity",
        "reach_change_selects_coefficient_sign",
        "reach_change_establishes_linear_sensitivity",
        "reach_change_establishes_openai_tip_geometry",
        "reach_change_establishes_visual_correspondence",
    ):
        mutated = deepcopy(scope)
        mutated["nonlinear_morphology_semantics"][key] = True
        with pytest.raises(ValueError):
            audit_axial_cap_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["capacity_evidence_only"]["morphology_capacity_is_pde_improvement"] = True
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(scope=mutated)


def test_rejects_live_route_drift():
    scope = load_scope()

    mutated = deepcopy(scope)
    mutated["integration_routing_scope"]["live_integrated_mode"] = "COMPACT_C4_ODD_Z_POLOIDAL"
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["integration_routing_scope"]["axial_cap_capacity_is_integrated_prematerialization_evidence"] = False
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["integration_routing_scope"]["additional_capacity_growth_allowed_before_live_materialization_screen"] = True
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["integration_routing_scope"]["promotion_requires"].remove("new_candidate_sha256")
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(scope=mutated)


def test_rejects_missing_future_identity_or_fresh_revalidation():
    scope = load_scope()
    for key in (
        "explicit_autonomous_coefficient_bound_required",
        "new_representation_family_identity_required",
        "new_candidate_sha_required",
        "fresh_full_per_component_momentum_divergence_required",
        "model_selection_data_may_not_be_reused_as_acceptance_data",
    ):
        mutated = deepcopy(scope)
        mutated["future_nonzero_materialization"][key] = False
        with pytest.raises(ValueError):
            audit_axial_cap_poloidal_scope(scope=mutated)


def test_rejects_cr001_and_project_status_drift():
    constraints = load_constraints()
    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(constraints=mutated_constraints)

    mutated_constraints = deepcopy(constraints)
    mutated_constraints["forcing"]["parameters"]["c"] = [0.0, 20.0]
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(constraints=mutated_constraints)

    status = load_project_status()
    mutated_status = deepcopy(status)
    mutated_status["next_integration_task"] = "materialize AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL now"
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["active_scientific_route"] = "materialize_compact_poloidal_child"
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["latest_integrated_axial_cap_capacity"]["production_coefficient_value_selected"] = True
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["integration_policy"]["avoid"].remove(
        "additional_capacity_basis_growth_before_a_current_integrated_direction_is_materialized_and_screened"
    )
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["states"]["paper_exact"] = True
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_scope(project_status=mutated_status)
