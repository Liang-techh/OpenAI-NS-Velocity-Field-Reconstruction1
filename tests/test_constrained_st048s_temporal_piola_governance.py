from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st048s_temporal_piola_governance import (
    audit_st048s_temporal_piola_governance,
    beta_gamma,
    beta_gamma_prime,
    load_constraints,
    load_project_status,
    load_scope,
)


def test_current_temporal_piola_governance_passes():
    report = audit_st048s_temporal_piola_governance()
    assert report["status"] == "governance_pass"
    assert report["upstream_pr"] == 418
    assert report["gamma_0p025_is_diagnostic_only"] is True
    assert report["beta_prime_endpoint_abs_at_gamma_0p025"] == 0.2
    assert report["u_t_map_motion_terms_required"] is True
    assert report["parent_pde_evidence_transfer_allowed"] is False
    assert report["live_route_preserved"] is True
    assert report["canonical_thresholds_unchanged"] is True
    assert report["velocity_export_ready"] is True
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False


def test_beta_schedule_and_time_derivative_are_explicit():
    assert [beta_gamma(t, 0.025) for t in (0.25, 0.5, 0.75)] == [0.1, 0.075, 0.1]
    assert [beta_gamma_prime(t, 0.025) for t in (0.25, 0.5, 0.75)] == [-0.2, 0.0, 0.2]


def test_rejects_omitted_map_motion_or_parent_pde_inheritance():
    scope = load_scope()
    for key, bad_value in (
        ("omitting_map_motion_terms_allowed", True),
        ("parent_static_warp_u_t_receipt_transfer_allowed", True),
        ("parent_pressure_force_momentum_receipt_transfer_allowed", True),
        ("H_t_contains_beta_prime", False),
        ("H_zt_contains_beta_prime", False),
    ):
        mutated = deepcopy(scope)
        mutated["time_derivative_contract"][key] = bad_value
        with pytest.raises(ValueError):
            audit_st048s_temporal_piola_governance(scope=mutated)


def test_rejects_source_laundering_and_production_selection():
    scope = load_scope()

    mutated = deepcopy(scope)
    mutated["source_classification"]["gamma_0p025_clean_crossing"] = "public_source_fact"
    with pytest.raises(ValueError):
        audit_st048s_temporal_piola_governance(scope=mutated)

    mutated = deepcopy(scope)
    mutated["source_classification"]["openai_hidden_time_dependent_coordinate_map"] = "autonomous_design"
    with pytest.raises(ValueError):
        audit_st048s_temporal_piola_governance(scope=mutated)

    for key in (
        "selects_production_gamma",
        "screen_grid_defines_production_gamma_bound",
        "morphology_crossing_is_visual_correspondence",
        "rank_two_sensitivity_is_pde_leverage_proof",
        "energy_support_core_divergence_preflights_are_full_pde_acceptance",
    ):
        mutated = deepcopy(scope)
        mutated["screen_evidence_only"][key] = True
        with pytest.raises(ValueError):
            audit_st048s_temporal_piola_governance(scope=mutated)


def test_rejects_materialization_shortcuts_or_live_route_replacement():
    scope = load_scope()
    status = load_project_status()

    for key in (
        "explicit_autonomous_beta_and_gamma_bounds_required",
        "new_representation_identity_required",
        "new_candidate_sha_required",
        "analytic_or_independently_checked_u_t_chain_rule_required",
        "rebuild_or_refit_restricted_pressure_force_required",
        "fresh_full_momentum_divergence_required",
        "unused_acceptance_data_required",
    ):
        mutated = deepcopy(scope)
        mutated["future_materialization"][key] = False
        with pytest.raises(ValueError):
            audit_st048s_temporal_piola_governance(scope=mutated)

    mutated = deepcopy(scope)
    mutated["integration_routing_scope"]["may_silently_replace_live_route"] = True
    with pytest.raises(ValueError):
        audit_st048s_temporal_piola_governance(scope=mutated)

    mutated_status = deepcopy(status)
    mutated_status["active_scientific_route"] = "materialize_st048s_temporal_piola_child"
    with pytest.raises(ValueError):
        audit_st048s_temporal_piola_governance(project_status=mutated_status)


def test_rejects_cr001_threshold_force_or_truth_state_drift():
    constraints = load_constraints()
    status = load_project_status()

    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError):
        audit_st048s_temporal_piola_governance(constraints=mutated_constraints)

    mutated_constraints = deepcopy(constraints)
    mutated_constraints["forcing"]["parameters"]["a"] = [0.0, 20.0]
    with pytest.raises(ValueError):
        audit_st048s_temporal_piola_governance(constraints=mutated_constraints)

    mutated_status = deepcopy(status)
    mutated_status["states"]["visual_correspondence_verified"] = True
    with pytest.raises(ValueError):
        audit_st048s_temporal_piola_governance(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["states"]["pde_validated"] = True
    with pytest.raises(ValueError):
        audit_st048s_temporal_piola_governance(project_status=mutated_status)
