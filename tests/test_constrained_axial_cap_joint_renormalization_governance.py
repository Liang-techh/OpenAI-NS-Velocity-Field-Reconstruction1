import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs" / "axial_cap_joint_renormalization_governance.json"
CONSTRAINTS_PATH = ROOT / "configs" / "constraints.json"
PROJECT_STATUS_PATH = ROOT / "project_status.json"


def _load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_axial_cap_joint_renorm_source_classes_stay_closed():
    contract = _load(CONTRACT_PATH)
    allowed = {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    assert set(contract["allowed_source_classes"]) == allowed
    assert set(contract["source_classification"].values()) <= allowed
    assert contract["source_classification"]["callable_velocity_delivery"] == "user_requirement"
    assert contract["source_classification"]["eq45_backbone"] == "public_source_fact"
    assert contract["source_classification"]["axial_cap_band_and_envelope"] == "autonomous_design"
    assert contract["source_classification"]["joint_global_energy_renormalization_rule"] == "autonomous_design"
    assert contract["source_classification"]["openai_hidden_numerical_profile"] == "pending_unknown"


def test_upstream_capacity_receipt_remains_diagnostic_and_unpromoted():
    contract = _load(CONTRACT_PATH)
    target = contract["audit_target"]
    evidence = contract["evidence_snapshot"]

    assert target["pull_request"] == 335
    assert target["upstream_head"] == "5e3a0c3d39dda3f2e40e0ef872c5dd25c289b307"
    assert target["mode"] == "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL"
    assert evidence["capacity_only"] is True
    assert evidence["trial_magnitude"] == 0.25
    assert evidence["negative_trial"]["common_velocity_scale"] == 0.9931700184515126
    assert evidence["positive_trial"]["common_velocity_scale"] == 0.9932358877752897
    assert evidence["negative_trial"]["q99_over_Zp"] == evidence["positive_trial"]["q99_over_Zp"] == 0.76171875
    assert evidence["baseline_q99_over_Zp"] == 0.38085937499999983
    assert evidence["held_out_pde_residual_evaluated"] is False
    assert evidence["production_coefficient_selected"] is False
    assert evidence["production_bound_selected"] is False
    assert evidence["candidate_sha_created"] is False


def test_joint_renormalization_is_new_candidate_transformation_not_energy_bookkeeping():
    semantics = _load(CONTRACT_PATH)["representation_semantics"]

    assert semantics["nontrivial_common_rescaling_changes_velocity_field"] is True
    assert semantics["nontrivial_common_rescaling_is_candidate_transformation"] is True
    assert semantics["renormalization_factor_is_part_of_candidate_definition"] is True
    assert semantics["requires_new_representation_family_identity_when_materialized"] is True
    assert semantics["requires_new_candidate_sha_when_materialized"] is True
    assert semantics["cannot_reuse_parent_pde_acceptance_results"] is True
    assert semantics["positive_common_scale_preserves_fixed_time_streamline_directions"] is True
    assert semantics["positive_common_scale_preserves_normalized_vorticity_location_fingerprints"] is True
    assert semantics["positive_common_scale_does_not_preserve_navier_stokes_balance"] is True
    assert semantics["positive_common_scale_does_not_certify_same_physical_time_pathlines"] is True
    assert semantics["positive_common_scale_does_not_transfer_pressure_or_forcing_consistency"] is True


def test_inherited_eq45_guard_and_sign_symmetric_tip_reach_cannot_select_production_child():
    contract = _load(CONTRACT_PATH)
    guard = contract["coefficient_and_guard_governance"]
    morphology = contract["morphology_semantics"]

    assert guard["plus_minus_0p25_is_diagnostic_only"] is True
    assert guard["trial_values_do_not_select_production_sign"] is True
    assert guard["trial_values_do_not_select_production_coefficient"] is True
    assert guard["trial_values_do_not_define_production_bound"] is True
    assert guard["inherited_eq45_coefficient_limit_is_representation_preflight_only"] is True
    assert guard["inherited_eq45_coefficient_limit_is_not_axial_cap_production_bound"] is True
    assert guard["passing_inherited_coefficient_preflight_is_not_candidate_acceptance"] is True
    assert guard["production_child_requires_explicit_autonomous_bound"] is True
    assert guard["production_child_requires_explicit_nonzero_value"] is True

    assert morphology["q90_q99_extension_is_target_free_capacity_evidence"] is True
    assert morphology["equal_q90_q99_for_both_trial_signs_does_not_select_sign"] is True
    assert morphology["enstrophy_based_tip_reach_can_be_even_in_coefficient_sign"] is True
    assert morphology["q90_q99_extension_does_not_establish_visual_correspondence"] is True
    assert morphology["q90_q99_extension_does_not_establish_openai_field_identity"] is True
    assert morphology["morphology_capacity_is_not_momentum_residual_jacobian_evidence"] is True
    assert morphology["morphology_capacity_is_not_pde_validation"] is True


def test_open_stacked_axial_cap_evidence_does_not_replace_live_integration_route():
    contract = _load(CONTRACT_PATH)
    project = _load(PROJECT_STATUS_PATH)
    integration = contract["integration_scope"]

    assert integration["axial_shoulder_dependency_integrated"] is True
    assert integration["axial_cap_mode_integrated_on_active_branch"] is False
    assert integration["joint_renormalization_evidence_integrated_on_active_branch"] is False
    assert integration["open_stacked_capacity_evidence_may_not_replace_live_route"] is True
    assert integration["current_live_route_remains_compact_poloidal_child_materialization"] is True
    assert integration["axial_cap_promotion_requires_explicit_integration_decision"] is True
    assert project["latest_integrated_poloidal_capacity"]["mode"] == "COMPACT_C4_ODD_Z_POLOIDAL"
    assert "COMPACT_C4_ODD_Z_POLOIDAL" in project["next_integration_task"]
    assert "additional_capacity_basis_growth_before_a_current_integrated_direction_is_materialized_and_screened" in project["integration_policy"]["avoid"]


def test_axial_cap_joint_renorm_contract_preserves_registered_cr001_protocol():
    contract = _load(CONTRACT_PATH)
    constraints = _load(CONSTRAINTS_PATH)
    snapshot = contract["canonical_invariants"]

    assert snapshot["nu"] == constraints["nu"] == 0.01
    assert snapshot["physical_domain"] == constraints["domain"]["physical"] == "R^3"
    assert snapshot["evaluation_box"] == constraints["domain"]["evaluation_box"]
    assert snapshot["support"] == constraints["domain"]["support"] == "r < 2 and abs(z) < 2"
    assert snapshot["time_interval"] == constraints["domain"]["time_interval"] == [0.25, 0.75]
    assert snapshot["reference_energy"] == constraints["nontriviality"]["reference_energy"] == 1.0
    assert snapshot["reference_energy_abs_tolerance"] == constraints["nontriviality"]["reference_energy_abs_tolerance"] == 0.001

    forcing = constraints["forcing"]
    assert snapshot["forcing_mode"] == forcing["mode"] == "restricted_two_parameter_family"
    assert snapshot["forcing_parameter_bounds"] == forcing["parameters"]
    assert "No residual-dependent basis or pointwise free force" in forcing["restriction"]

    validation = constraints["validation"]
    assert snapshot["validation_seed"] == validation["seed"] == 914027
    assert snapshot["held_out_points"] == validation["held_out_points"] == 4096
    assert snapshot["validation_times"] == validation["times"]
    assert snapshot["derivative_steps"] == validation["derivative_steps"] == [0.02, 0.01, 0.005]
    assert snapshot["divergence_max"] == validation["thresholds"]["divergence_max"] == 1e-5
    assert snapshot["divergence_L2"] == validation["thresholds"]["divergence_L2"] == 1e-5
    assert snapshot["pde_residual_max"] == validation["thresholds"]["pde_residual_max"] == 1e-3
    assert snapshot["pde_residual_L2"] == validation["thresholds"]["pde_residual_L2"] == 1e-3


def test_materialized_child_requires_fresh_full_replay_and_cannot_upgrade_truth_states():
    contract = _load(CONTRACT_PATH)
    replay = contract["acceptance_replay"]
    anti = contract["anti_shortcut_rules"]
    status = contract["current_truth_state"]
    project = _load(PROJECT_STATUS_PATH)
    required = set(replay["required_checks_after_materialization"])
    forbidden = "\n".join(contract["forbidden_inferences"]).lower()

    assert replay["restoring_reference_energy_is_not_acceptance"] is True
    assert replay["formal_validation_must_be_fresh_from_model_selection_data"] is True
    assert replay["coefficient_or_renormalization_selection_data_are_model_selection_evidence"] is True
    assert replay["repository_st006_superiority_requires_protocol_aligned_replay"] is True
    assert {
        "reference_energy",
        "validation_time_energy",
        "core_rotation_sign",
        "bipolar_radial_axial_direction_and_parity",
        "axis_regularity",
        "physical_support",
        "registered_numeric_divergence",
        "pressure_and_restricted_force_consistency",
        "full_per_component_momentum_residual",
        "actual_3d_streamline_and_vorticity_fingerprint",
    } <= required

    assert anti["free_or_residual_dependent_f_equals_R_forbidden"] is True
    assert anti["zero_or_collapsed_velocity_shortcut_forbidden"] is True
    assert anti["renormalization_may_not_change_forcing_family"] is True
    assert anti["renormalization_may_not_change_nu_domain_support_norms_derivative_ladder_or_thresholds"] is True
    assert anti["green_ci_may_not_promote_scientific_state"] is True
    assert anti["unintegrated_stacked_evidence_may_not_be_called_canonical_candidate"] is True

    assert project["states"]["velocity_export_ready"] is True
    assert status["canonical_velocity_changed_by_this_contract"] is False
    assert status["axial_cap_joint_renormalized_child_materialized"] is False
    assert status["axial_cap_production_coefficient_selected"] is False
    assert status["axial_cap_production_bound_selected"] is False
    assert status["visual_correspondence_verified"] is False
    assert status["pde_validated"] is False
    assert status["paper_exact"] is False
    assert status["openai_field_identified"] is False
    assert status["blowup_proved"] is False

    assert "pde validation" in forbidden
    assert "production bound" in forbidden
    assert "production coefficient sign" in forbidden
    assert "visual correspondence" in forbidden
    assert "openai field" in forbidden
    assert "paper-exact" in forbidden
    assert "canonical delivered velocity" in forbidden
    assert "green ci" in forbidden
