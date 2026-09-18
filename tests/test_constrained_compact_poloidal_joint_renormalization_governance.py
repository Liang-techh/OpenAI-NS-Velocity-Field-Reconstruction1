import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs" / "compact_poloidal_joint_renormalization_governance.json"
CONSTRAINTS_PATH = ROOT / "configs" / "constraints.json"


def _load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_joint_renormalization_source_classes_stay_closed():
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
    assert contract["source_classification"]["compact_poloidal_mode"] == "autonomous_design"
    assert contract["source_classification"]["joint_global_renormalization_rule"] == "autonomous_design"
    assert contract["source_classification"]["openai_hidden_numerical_profile"] == "pending_unknown"


def test_joint_renormalization_is_a_new_candidate_path_not_energy_bookkeeping():
    contract = _load(CONTRACT_PATH)
    semantics = contract["representation_semantics"]
    assert semantics["nontrivial_global_rescaling_changes_velocity_field"] is True
    assert semantics["nontrivial_global_rescaling_is_candidate_transformation"] is True
    assert semantics["requires_new_representation_family_identity"] is True
    assert semantics["requires_new_candidate_sha"] is True
    assert semantics["cannot_reuse_parent_acceptance_results"] is True
    assert semantics["fixed_parent_energy_envelope_remains_valid_for_direct_additive_materialization"] is True
    assert semantics["jointly_renormalized_child_is_a_distinct_materialization_path"] is True
    assert semantics["renormalization_factor_is_part_of_candidate_definition"] is True


def test_probe_coefficients_and_energy_roots_cannot_be_silently_promoted():
    contract = _load(CONTRACT_PATH)
    governance = contract["coefficient_governance"]
    assert governance["target_pr_probe_coefficients_are_diagnostic_only"] is True
    assert governance["probe_values_do_not_select_production_coefficient"] is True
    assert governance["probe_values_do_not_define_production_bound"] is True
    assert governance["energy_root_values_are_diagnostic_only_unless_separately_preregistered"] is True
    assert governance["production_child_requires_explicit_autonomous_bound"] is True
    assert governance["production_child_requires_explicit_nonzero_value"] is True


def test_energy_restoration_requires_full_fresh_acceptance_replay():
    contract = _load(CONTRACT_PATH)
    replay = contract["acceptance_replay"]
    required = set(replay["required_checks_after_materialization"])
    assert replay["restoring_reference_energy_is_not_acceptance"] is True
    assert replay["formal_validation_must_be_fresh_from_model_selection_data"] is True
    assert replay["selection_or_renormalization_tuning_data_are_model_selection_evidence"] is True
    assert {
        "reference_energy",
        "validation_time_energy",
        "core_rotation_sign",
        "bipolar_radial_axial_direction_and_parity",
        "axis_regularity",
        "physical_support",
        "registered_numeric_divergence",
        "full_per_component_momentum_residual",
    } <= required


def test_joint_renormalization_contract_preserves_registered_cr001_protocol():
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
    assert snapshot["derivative_steps"] == validation["derivative_steps"] == [0.02, 0.01, 0.005]
    assert snapshot["divergence_max"] == validation["thresholds"]["divergence_max"] == 1e-5
    assert snapshot["divergence_L2"] == validation["thresholds"]["divergence_L2"] == 1e-5
    assert snapshot["pde_residual_max"] == validation["thresholds"]["pde_residual_max"] == 1e-3
    assert snapshot["pde_residual_L2"] == validation["thresholds"]["pde_residual_L2"] == 1e-3


def test_joint_renormalization_cannot_upgrade_independent_truth_states():
    contract = _load(CONTRACT_PATH)
    anti = contract["anti_shortcut_rules"]
    status = contract["current_truth_state"]
    forbidden = "\n".join(contract["forbidden_inferences"]).lower()

    assert anti["free_or_residual_dependent_f_equals_R_forbidden"] is True
    assert anti["zero_or_collapsed_velocity_shortcut_forbidden"] is True
    assert anti["renormalization_may_not_change_forcing_family"] is True
    assert anti["renormalization_may_not_change_nu_domain_support_norms_derivative_ladder_or_thresholds"] is True

    assert status["joint_renormalized_nonzero_child_promoted"] is False
    assert status["pde_validated"] is False
    assert status["visual_correspondence_verified"] is False
    assert status["paper_exact"] is False
    assert status["openai_field_identified"] is False

    assert "pde validation" in forbidden
    assert "visual correspondence" in forbidden
    assert "paper-exact" in forbidden
    assert "openai field" in forbidden
    assert "green ci" in forbidden
