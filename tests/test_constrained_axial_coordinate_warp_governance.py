from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_axial_coordinate_warp_governance import (
    audit_axial_coordinate_warp_governance,
    load_constraints,
    load_project_status,
    load_scope,
)


def test_current_coordinate_warp_governance_passes():
    report = audit_axial_coordinate_warp_governance()
    assert report["status"] == "governance_pass"
    assert report["upstream_pr"] == 367
    assert report["upstream_status"] == "open_unintegrated_experimental_evidence"
    assert report["beta_0p40_is_diagnostic_only"] is True
    assert report["warp_is_navier_stokes_invariance"] is False
    assert report["parent_pde_evidence_transfer_allowed"] is False
    assert report["live_route_preserved"] is True
    assert report["canonical_thresholds_unchanged"] is True
    assert report["velocity_export_ready"] is True
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False


def test_rejects_source_laundering_or_hidden_map_recovery():
    scope = load_scope()

    mutated = deepcopy(scope)
    mutated["source_classification"]["beta_0p40_diagnostic_crossing"] = "public_source_fact"
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(scope=mutated)

    mutated = deepcopy(scope)
    mutated["source_classification"]["openai_hidden_coordinate_or_aspect_map"] = "autonomous_design"
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(scope=mutated)


def test_rejects_capacity_to_selection_or_visual_claim_promotion():
    scope = load_scope()
    for key in (
        "beta_0p40_selects_production_value",
        "screen_beta_grid_defines_production_bound",
        "positive_jacobian_establishes_openai_coordinate_map",
        "target_free_morphology_is_visual_correspondence",
        "lower_collar_or_radial_cost_is_candidate_superiority",
    ):
        mutated = deepcopy(scope)
        mutated["screen_evidence_only"][key] = True
        with pytest.raises(ValueError):
            audit_axial_coordinate_warp_governance(scope=mutated)


def test_rejects_divergence_energy_or_parent_pde_inheritance_laundering():
    scope = load_scope()
    for key in (
        "piola_identity_is_registered_full_candidate_divergence_acceptance",
        "fd_divergence_is_full_pde_acceptance",
        "common_scale_is_bookkeeping_only",
        "parent_pressure_force_or_pde_evidence_transfer_allowed",
    ):
        mutated = deepcopy(scope)
        mutated["screen_evidence_only"][key] = True
        with pytest.raises(ValueError):
            audit_axial_coordinate_warp_governance(scope=mutated)

    mutated = deepcopy(scope)
    mutated["representation"]["navier_stokes_invariant_under_warp"] = True
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(scope=mutated)


def test_rejects_derivative_or_materialization_shortcuts():
    scope = load_scope()

    mutated = deepcopy(scope)
    mutated["derivative_semantics"]["future_time_dependent_beta_may_inherit_static_warp_pde_evidence"] = True
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(scope=mutated)

    for key in (
        "explicit_autonomous_beta_bound_required",
        "explicit_nonzero_beta_value_required",
        "explicit_normalization_rule_required",
        "new_representation_family_identity_required",
        "new_candidate_sha_required",
        "rebuild_or_refit_restricted_pressure_force_contract",
        "fresh_full_per_component_momentum_divergence_required",
        "unused_acceptance_data_required",
    ):
        mutated = deepcopy(scope)
        mutated["future_warp_materialization"][key] = False
        with pytest.raises(ValueError):
            audit_axial_coordinate_warp_governance(scope=mutated)


def test_rejects_silent_replacement_of_live_axial_cap_route():
    scope = load_scope()
    status = load_project_status()

    mutated = deepcopy(scope)
    mutated["integration_routing_scope"]["coordinate_warp_may_silently_replace_live_route"] = True
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(scope=mutated)

    mutated_status = deepcopy(status)
    mutated_status["active_scientific_route"] = "materialize_axial_coordinate_warp_child"
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["latest_integrated_axial_cap_capacity"]["candidate_sha_created"] = True
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["next_integration_task"] = "switch directly to coordinate warp"
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(project_status=mutated_status)


def test_rejects_cr001_threshold_force_or_truth_state_drift():
    constraints = load_constraints()
    status = load_project_status()

    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(constraints=mutated_constraints)

    mutated_constraints = deepcopy(constraints)
    mutated_constraints["forcing"]["parameters"]["a"] = [0.0, 20.0]
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(constraints=mutated_constraints)

    mutated_status = deepcopy(status)
    mutated_status["states"]["visual_correspondence_verified"] = True
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["states"]["pde_validated"] = True
    with pytest.raises(ValueError):
        audit_axial_coordinate_warp_governance(project_status=mutated_status)
