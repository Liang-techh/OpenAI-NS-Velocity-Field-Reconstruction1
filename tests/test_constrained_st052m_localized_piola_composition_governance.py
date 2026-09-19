from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction import constrained_st052m_localized_piola_composition_governance as gov


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


def test_current_repository_passes_localized_piola_composition_audit():
    result = gov.audit_current_repository()
    assert result == {
        "task_id": gov.TASK_ID,
        "status": "pass",
        "taper_child_materialized": False,
        "central_coordinate_map_identity": True,
        "central_final_velocity_identity": False,
        "taper_material_paths_integrated": False,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "canonical_velocity_unchanged": True,
    }


def test_coordinate_map_identity_is_not_final_velocity_identity_after_scale():
    c = gov.load_contract()
    taper = c["frozen_taper"]
    rep = c["representation_audit"]
    sampled = c["sampled_taper_evidence"]
    assert taper["map_identity_for_abs_z_over_2_le_0p50"] is True
    assert taper["taper_energy_scale_reported"] == pytest.approx(0.9929528816)
    assert taper["taper_energy_scale_reported"] != pytest.approx(1.0)
    assert rep["coordinate_map_identity_implies_final_velocity_identity_after_common_scale"] is False
    assert rep["central_final_velocity_is_positive_common_scale_times_parent"] is True
    assert rep["central_velocity_amplitude_is_exactly_unchanged"] is False
    assert sampled["central_coordinate_map_identity_error"] == 0.0
    assert sampled["central_angular_rate_retention_r_0p6"] == pytest.approx(
        taper["taper_energy_scale_reported"], abs=2e-6
    )


def test_callable_piola_structure_does_not_promote_pde_or_support_identity():
    c = gov.load_contract()
    rep = c["representation_audit"]
    assert rep["transform_is_applied_to_callable_parent_not_trilinear_sample_grid"] is True
    assert rep["continuous_piola_identity_conditionally_preserves_divergence"] is True
    assert rep["common_constant_scale_preserves_divergence_identity"] is True
    assert rep["a_ge_1_and_z_unchanged_preserve_parent_zero_extension_outside_registered_support"] is True
    assert rep["nonzero_support_region_is_proved_identical_to_parent"] is False
    assert rep["piola_divergence_identity_is_navier_stokes_invariance"] is False
    assert rep["common_energy_scale_is_navier_stokes_invariance"] is False
    assert rep["pressure_can_be_inherited_without_rebuild"] is False
    assert rep["restricted_force_fit_can_be_inherited_without_rebuild"] is False
    assert rep["momentum_residual_can_be_inherited_from_parent"] is False


def test_parent_material_paths_and_pde_receipts_do_not_transfer_through_taper():
    c = gov.load_contract()
    own = c["evidence_ownership_and_noninheritance"]
    assert gov.load_parent_governance()["child_material_path_evidence"][
        "positive_material_path_transfer_under_frozen_protocol"
    ] is True
    assert own["redistributed_parent_material_paths_belong_to_taper_child"] is False
    assert own["redistributed_parent_positive_material_path_transfer_is_taper_material_path_evidence"] is False
    assert own["taper_material_paths_integrated"] is False
    assert own["st052m_parent_pde_receipt_belongs_to_taper_child"] is False


def test_order64_and_custom_fd_screen_are_not_cr001_acceptance():
    c = gov.load_contract()
    sampled = c["sampled_taper_evidence"]
    cr = c["canonical_cr001"]
    assert sampled["energy_quadrature_order"] == 64
    assert cr["energy_quadrature_orders_per_axis"] == [24, 48, 96]
    assert sampled["registered_CR001_energy_quadrature_ladder_replayed"] is False
    assert sampled["cartesian_fd_divergence_max"] < cr["divergence_max"]
    assert sampled["formal_CR001_divergence_acceptance_replayed"] is False
    assert sampled["formal_CR001_full_momentum_replayed"] is False
    assert sampled["pde_validated"] is False


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("representation_audit", "coordinate_map_identity_implies_final_velocity_identity_after_common_scale"), True),
        (("representation_audit", "central_velocity_amplitude_is_exactly_unchanged"), True),
        (("representation_audit", "piola_divergence_identity_is_navier_stokes_invariance"), True),
        (("representation_audit", "common_energy_scale_is_navier_stokes_invariance"), True),
        (("representation_audit", "pressure_can_be_inherited_without_rebuild"), True),
        (("representation_audit", "restricted_force_fit_can_be_inherited_without_rebuild"), True),
        (("representation_audit", "momentum_residual_can_be_inherited_from_parent"), True),
        (("sampled_taper_evidence", "registered_CR001_energy_quadrature_ladder_replayed"), True),
        (("sampled_taper_evidence", "formal_CR001_divergence_acceptance_replayed"), True),
        (("sampled_taper_evidence", "formal_CR001_full_momentum_replayed"), True),
        (("sampled_taper_evidence", "visual_correspondence_verified"), True),
        (("sampled_taper_evidence", "pde_validated"), True),
        (("evidence_ownership_and_noninheritance", "redistributed_parent_material_paths_belong_to_taper_child"), True),
        (("evidence_ownership_and_noninheritance", "redistributed_parent_positive_material_path_transfer_is_taper_material_path_evidence"), True),
        (("evidence_ownership_and_noninheritance", "taper_material_paths_integrated"), True),
        (("evidence_ownership_and_noninheritance", "st052m_parent_pde_receipt_belongs_to_taper_child"), True),
        (("evidence_ownership_and_noninheritance", "rank3_response_selects_production_tau"), True),
        (("evidence_ownership_and_noninheritance", "central_map_identity_proves_central_velocity_identity"), True),
        (("evidence_ownership_and_noninheritance", "sampled_order64_energy_scale_proves_CR001_energy_acceptance"), True),
        (("evidence_ownership_and_noninheritance", "custom_fd_divergence_below_1e5_proves_CR001_divergence_acceptance"), True),
        (("evidence_ownership_and_noninheritance", "target_free_morphology_pass_proves_OpenAI_correspondence"), True),
        (("composition_identity", "diagnostic_child_versioned_save_load_candidate"), True),
        (("composition_identity", "diagnostic_child_candidate_sha256_assigned"), True),
        (("composition_identity", "diagnostic_child_unified_velocity_api_registered"), True),
        (("composition_identity", "diagnostic_child_production_selected"), True),
        (("routing", "localized_taper_supersedes_governed_kappa_0p05_materialization_priority"), True),
        (("routing", "combined_taper_child_ready_for_pressure_or_restricted_force_reconstruction"), True),
        (("routing", "combined_taper_child_ready_for_production_selection"), True),
        (("independent_truth_states", "localized_taper_child_velocity_export_ready"), True),
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
        (("snapshot", "taper_pr_head"), "0" * 40),
        (("snapshot", "taper_pr_base_head"), "0" * 40),
        (("snapshot", "artifact_digest"), "sha256:" + "0" * 64),
        (("frozen_taper", "tau"), 0.075),
        (("frozen_taper", "window_abs_z_over_2"), [0.45, 0.82]),
        (("frozen_taper", "taper_energy_scale_reported"), 1.0),
        (("sampled_taper_evidence", "energy_quadrature_order"), 96),
        (("canonical_cr001", "validation_seed"), 9175291),
        (("canonical_cr001", "energy_quadrature_orders_per_axis"), [32, 64, 128]),
        (("canonical_cr001", "divergence_max"), 1e-4),
        (("canonical_cr001", "pde_residual_max"), 0.01),
        (("canonical_delivery", "candidate_family"), "ST052-M-taper"),
        (("canonical_delivery", "velocity_api"), "diagnostic:taper_velocity"),
    ],
)
def test_identity_threshold_and_protocol_drift_fail_closed(path, value):
    with pytest.raises(ValueError):
        _audit_mutated(path, value)


def test_future_selection_requires_new_identity_and_exact_child_replays():
    req = gov.load_contract()["future_selection_requirements"]
    assert req["assign_new_candidate_id_and_sha256_if_materialized"] is True
    assert req["create_reproducible_save_load_recipe_if_materialized"] is True
    assert req["register_explicit_velocity_x_y_z_t_api_or_adapter_if_materialized"] is True
    assert req["rerun_material_paths_on_exact_taper_child_before_material_path_claim"] is True
    assert req["run_public_reference_render_comparison_and_resolution_stability_before_visual_correspondence"] is True
    assert req["rebuild_or_refit_compatible_pressure_before_pde_claim"] is True
    assert req["rebuild_or_refit_only_existing_restricted_two_parameter_force_family"] is True
    assert req["do_not_use_free_or_residual_defined_force"] is True
    assert req["run_fresh_full_momentum_and_divergence_validation_on_exact_materialized_child"] is True
    assert req["run_registered_CR001_acceptance_contract_before_formal_pde_promotion"] is True
