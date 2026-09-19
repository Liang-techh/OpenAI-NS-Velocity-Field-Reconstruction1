from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction import constrained_st052m_redistribution_evidence_governance as gov


def _audit_mutated(path: tuple[str, ...], value):
    contract = deepcopy(gov.load_contract())
    cursor = contract
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    return gov.audit_contract(
        contract,
        gov.load_constraints(),
        gov.load_parent_governance(),
        gov.load_project_status(),
        gov.load_velocity_delivery_contract(),
    )


def test_current_repository_passes_evidence_composition_audit():
    result = gov.audit_current_repository()
    assert result == {
        "task_id": gov.TASK_ID,
        "status": "pass",
        "parent_candidate": "ST052-M",
        "transformed_child_materialized": False,
        "parent_pde_receipt_transfers": False,
        "positive_material_path_transfer": True,
        "materialization_priority_supported": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "canonical_velocity_unchanged": True,
    }


def test_transform_is_axis_regular_sign_preserving_and_support_preserving_by_contract():
    c = gov.load_contract()
    transform = c["frozen_transform"]
    rep = c["representation_audit"]
    assert transform["identity_near_axis_below_r_0p30"] is True
    assert transform["multiplier_global_lower_bound_from_profile_range"] > 0.0
    assert transform["multiplier_global_upper_bound_from_profile_range"] == pytest.approx(1.05)
    assert rep["near_axis_parent_swirl_order_is_unchanged"] is True
    assert rep["finite_positive_multiplier_preserves_parent_swirl_sign"] is True
    assert rep["pointwise_multiplier_preserves_parent_zero_support"] is True
    assert rep["axisymmetric_swirl_only_multiplier_adds_no_theta_divergence_term"] is True
    assert rep["analytic_divergence_structure_preservation_is_full_pde_preservation"] is False


def test_positive_material_path_receipt_is_kept_separate_from_outer_radius_tradeoff():
    paths = gov.load_contract()["child_material_path_evidence"]
    assert paths["mean_absolute_turns_relative_change"] > 0.0
    assert paths["maximum_absolute_turns_relative_change"] > 0.0
    assert paths["radial_contraction_magnitude_relative_change"] > 0.0
    assert paths["mean_pair_axial_separation_change_relative"] > 0.0
    assert paths["winding_relative_change_by_seed_radius"]["0.6"] > 0.0
    assert paths["winding_relative_change_by_seed_radius"]["0.9"] > 0.0
    assert paths["winding_relative_change_by_seed_radius"]["1.2"] < 0.0
    assert paths["visual_correspondence_verified"] is False
    assert paths["pde_evidence"] is False


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("parent_pde_evidence", "may_transfer_to_transformed_child"), True),
        (("representation_audit", "pressure_can_be_inherited_without_rebuild"), True),
        (("representation_audit", "restricted_force_fit_can_be_inherited_without_rebuild"), True),
        (("representation_audit", "momentum_residual_can_be_inherited_from_parent"), True),
        (("representation_audit", "reference_energy_normalization_is_navier_stokes_invariance"), True),
        (("child_eulerian_evidence", "full_registered_divergence_acceptance_replayed"), True),
        (("child_eulerian_evidence", "full_momentum_replayed"), True),
        (("child_eulerian_evidence", "pde_validated"), True),
        (("child_material_path_evidence", "visual_correspondence_verified"), True),
        (("child_material_path_evidence", "pde_evidence"), True),
        (("child_delivery_identity", "materialized_versioned_candidate_exists"), True),
        (("child_delivery_identity", "candidate_sha256_assigned"), True),
        (("child_delivery_identity", "save_load_contract_exists"), True),
        (("child_delivery_identity", "unified_velocity_api_registered"), True),
        (("child_delivery_identity", "production_candidate_selected"), True),
        (("evidence_composition_rules", "parent_pde_plus_child_eulerian_equals_child_pde"), True),
        (("evidence_composition_rules", "parent_pde_plus_child_material_paths_equals_child_pde"), True),
        (("evidence_composition_rules", "child_eulerian_plus_child_material_paths_equals_visual_correspondence"), True),
        (("evidence_composition_rules", "clean_fd_divergence_plus_parent_momentum_equals_child_full_acceptance"), True),
        (("evidence_composition_rules", "positive_material_path_transfer_selects_production_candidate"), True),
        (("independent_truth_states", "transformed_child_velocity_export_ready"), True),
        (("independent_truth_states", "visual_correspondence_verified"), True),
        (("independent_truth_states", "pde_validated"), True),
        (("independent_truth_states", "paper_exact"), True),
        (("independent_truth_states", "openai_field_identified"), True),
    ],
)
def test_forbidden_promotions_fail_closed(path, value):
    with pytest.raises(ValueError):
        _audit_mutated(path, value)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("frozen_transform", "alpha"), 2.4),
        (("frozen_transform", "kappa"), 0.025),
        (("frozen_transform", "exact_eulerian_normalization_receipt"), 1.0),
        (("snapshot", "parent_pr_head"), "0" * 40),
        (("snapshot", "eulerian_transform_head"), "0" * 40),
        (("snapshot", "material_path_head"), "0" * 40),
        (("parent_identity", "raw_candidate_sha256"), "0" * 64),
        (("parent_identity", "modifier_sha256"), "0" * 64),
        (("canonical_cr001", "validation_seed"), 9175291),
        (("canonical_cr001", "pde_residual_max"), 0.01),
        (("canonical_cr001", "divergence_max"), 1e-4),
        (("canonical_delivery", "candidate_family"), "ST052-M"),
        (("canonical_delivery", "velocity_api"), "diagnostic:velocity"),
    ],
)
def test_identity_threshold_and_lineage_drift_fail_closed(path, value):
    with pytest.raises(ValueError):
        _audit_mutated(path, value)


def test_parent_pde_receipt_remains_metric_qualified_and_noncanonical():
    c = gov.load_contract()
    pde = c["parent_pde_evidence"]
    assert pde["sampled_max_improved"] is True
    assert pde["fixed_time_volume_L2_improved"] is False
    assert pde["formal_1e3_momentum_gates_pass"] is False
    assert pde["registered_CR001_seed_914027_replayed"] is False
    assert pde["may_transfer_to_transformed_child"] is False


def test_materialization_route_requires_new_identity_and_fresh_child_validation():
    req = gov.load_contract()["future_materialization_requirements"]
    assert req["assign_new_candidate_id_and_sha256"] is True
    assert req["create_reproducible_save_load_recipe"] is True
    assert req["register_explicit_velocity_x_y_z_t_api_or_adapter"] is True
    assert req["rebuild_or_refit_compatible_pressure"] is True
    assert req["rebuild_or_refit_only_the_existing_restricted_two_parameter_force_family"] is True
    assert req["do_not_use_free_or_residual_defined_force"] is True
    assert req["run_fresh_full_momentum_and_divergence_validation_on_exact_materialized_child"] is True
    assert req["run_registered_CR001_acceptance_contract_before_formal_pde_promotion"] is True
