from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st048s_swirl_gain_governance import (
    audit_st048s_swirl_gain_governance,
    load_constraints,
    load_project_status,
    load_scope,
)


def test_current_swirl_gain_governance_passes():
    report = audit_st048s_swirl_gain_governance()
    assert report["status"] == "governance_pass"
    assert report["upstream_pr"] == 428
    assert report["kappa_0p05_is_diagnostic_only"] is True
    assert report["constant_kappa_adds_no_kappa_prime_term"] is True
    assert report["temporal_piola_chain_terms_still_required"] is True
    assert report["axis_regular_extension_required"] is True
    assert report["parent_pde_evidence_transfer_allowed"] is False
    assert report["live_route_preserved"] is True
    assert report["canonical_thresholds_unchanged"] is True
    assert report["velocity_export_ready"] is True
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False


def test_rejects_derivative_axis_or_pde_inheritance_shortcuts():
    scope = load_scope()
    for key, bad_value in (
        ("temporal_piola_map_motion_terms_still_required", False),
        ("future_time_dependent_kappa_requires_kappa_prime_u_theta_term", False),
        ("future_time_dependent_scale_requires_scale_prime_u_term", False),
        ("cylindrical_to_cartesian_export_must_use_regular_axis_extension", False),
        ("near_axis_finite_cartesian_velocity_required", False),
        ("fresh_numerical_divergence_recheck_required_after_materialization", False),
        ("parent_pressure_force_momentum_receipt_transfer_allowed", True),
        ("parent_pde_validation_transfer_allowed", True),
    ):
        mutated = deepcopy(scope)
        mutated["derivative_and_axis_contract"][key] = bad_value
        with pytest.raises(ValueError):
            audit_st048s_swirl_gain_governance(scope=mutated)


def test_rejects_source_laundering_and_diagnostic_kappa_promotion():
    scope = load_scope()

    mutated = deepcopy(scope)
    mutated["source_classification"]["kappa_0p05_clean_winding_crossing"] = "public_source_fact"
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(scope=mutated)

    mutated = deepcopy(scope)
    mutated["source_classification"]["openai_hidden_swirl_gain_or_coefficient"] = "autonomous_design"
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(scope=mutated)

    for key in (
        "selects_production_kappa",
        "screen_grid_defines_production_kappa_bound",
        "winding_proxy_crossing_is_visual_correspondence",
        "rank_three_sensitivity_is_pde_leverage_proof",
        "divergence_support_energy_core_preflights_are_full_pde_acceptance",
    ):
        mutated = deepcopy(scope)
        mutated["screen_evidence_only"][key] = True
        with pytest.raises(ValueError):
            audit_st048s_swirl_gain_governance(scope=mutated)


def test_rejects_lineage_identity_or_live_route_shortcuts():
    scope = load_scope()
    status = load_project_status()

    mutated = deepcopy(scope)
    mutated["lineage_and_routing_scope"]["branch_tip_is_scientific_identity"] = True
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(scope=mutated)

    mutated = deepcopy(scope)
    mutated["lineage_and_routing_scope"]["exact_parent_commit_representation_identity_and_candidate_sha_required_before_promotion"] = False
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(scope=mutated)

    mutated = deepcopy(scope)
    mutated["lineage_and_routing_scope"]["downstream_green_ci_promotes_parent_or_child"] = True
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(scope=mutated)

    mutated = deepcopy(scope)
    mutated["lineage_and_routing_scope"]["may_silently_replace_canonical_velocity_or_live_route"] = True
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(scope=mutated)

    mutated_status = deepcopy(status)
    mutated_status["active_scientific_route"] = "materialize_st048s_temporal_piola_swirl_gain_child"
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(project_status=mutated_status)


def test_rejects_materialization_without_new_identity_or_fresh_validation():
    scope = load_scope()
    for key in (
        "explicit_autonomous_beta_gamma_kappa_bounds_required",
        "new_representation_identity_required",
        "new_candidate_sha_required",
        "save_load_replay_required",
        "analytic_or_independently_checked_temporal_chain_rule_required",
        "axis_regular_cartesian_export_check_required",
        "rebuild_or_refit_restricted_pressure_force_required",
        "fresh_full_momentum_divergence_required",
        "unused_acceptance_data_required",
    ):
        mutated = deepcopy(scope)
        mutated["future_materialization"][key] = False
        with pytest.raises(ValueError):
            audit_st048s_swirl_gain_governance(scope=mutated)


def test_rejects_cr001_threshold_force_or_truth_state_drift():
    constraints = load_constraints()
    status = load_project_status()

    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(constraints=mutated_constraints)

    mutated_constraints = deepcopy(constraints)
    mutated_constraints["forcing"]["parameters"]["c"] = [0.0, 20.0]
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(constraints=mutated_constraints)

    mutated_status = deepcopy(status)
    mutated_status["states"]["visual_correspondence_verified"] = True
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["states"]["pde_validated"] = True
    with pytest.raises(ValueError):
        audit_st048s_swirl_gain_governance(project_status=mutated_status)
