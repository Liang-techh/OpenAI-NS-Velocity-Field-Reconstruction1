import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GOVERNANCE = ROOT / "configs" / "axial_cap_temporal_ramp_governance.json"
CONSTRAINTS = ROOT / "configs" / "constraints.json"
PROJECT_STATUS = ROOT / "project_status.json"


def _load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_source_taxonomy_is_fail_closed_and_upstream_remains_capacity_only():
    governance = _load(GOVERNANCE)
    taxonomy = governance["source_taxonomy"]
    assert taxonomy["allowed_classes"] == [
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    ]
    assert set(taxonomy["classification"].values()) <= set(taxonomy["allowed_classes"])
    assert taxonomy["classification"]["linear_temporal_ramp_schedule"] == "autonomous_design"
    assert taxonomy["classification"]["diagnostic_terminal_magnitude_0p25"] == "autonomous_design"
    assert taxonomy["classification"]["openai_hidden_numeric_profile"] == "pending_unknown"
    assert governance["governed_upstream"]["evidence_status"] == "open_unintegrated_capacity_evidence"
    assert governance["governed_upstream"]["observations"]["held_out_pde_residual_evaluated"] is False
    assert governance["governed_upstream"]["observations"]["production_coefficient_selected"] is False
    assert governance["governed_upstream"]["observations"]["candidate_sha_created"] is False


def test_temporal_ramp_cannot_inherit_static_pde_evidence_or_select_sign_from_morphology():
    governance = _load(GOVERNANCE)
    upstream = governance["governed_upstream"]
    guards = governance["semantic_guards"]

    t0, t1 = upstream["registered_window"]
    amplitude = upstream["diagnostic_terminal_magnitude"]
    assert amplitude / (t1 - t0) == upstream["diagnostic_coefficient_time_derivative_magnitude"] == 0.5
    assert guards["time_derivative"]["diagnostic_abs_a_prime"] == 0.5
    assert guards["time_derivative"]["fresh_u_t_aware_validation_required"] is True
    assert guards["reference_replay"]["candidate_identity_must_change_if_materialized"] is True
    assert guards["reference_replay"]["candidate_sha_must_change_if_materialized"] is True

    morphology = guards["morphology"]
    assert upstream["observations"]["q99_far_tail_reach_increased"] is True
    assert upstream["observations"]["q90_core_reach_increased"] is False
    assert morphology["q90_no_gain_does_not_imply_no_pde_leverage"] is True
    assert morphology["q99_gain_does_not_establish_openai_correspondence"] is True
    assert morphology["near_sign_symmetry_does_not_select_production_sign"] is True
    assert morphology["production_sign_requires_signed_velocity_or_streamline_evidence_and_candidate_screening"] is True


def test_cr001_lock_matches_canonical_constraints_exactly():
    governance = _load(GOVERNANCE)
    constraints = _load(CONSTRAINTS)
    lock = governance["canonical_cr001_lock"]

    assert lock["nu"] == constraints["nu"] == 0.01
    assert lock["physical_domain"] == constraints["domain"]["physical"] == "R^3"
    assert lock["evaluation_box"] == constraints["domain"]["evaluation_box"]
    assert lock["support"] == constraints["domain"]["support"]
    assert lock["time_interval"] == constraints["domain"]["time_interval"]
    assert lock["forcing_mode"] == constraints["forcing"]["mode"]
    assert lock["forcing_parameter_bounds"] == constraints["forcing"]["parameters"]
    assert lock["reference_energy"] == constraints["nontriviality"]["reference_energy"]
    assert lock["reference_energy_abs_tolerance"] == constraints["nontriviality"]["reference_energy_abs_tolerance"]
    assert lock["minimum_energy_each_validation_time"] == constraints["nontriviality"]["minimum_energy_each_validation_time"]
    assert lock["maximum_energy_each_validation_time"] == constraints["nontriviality"]["maximum_energy_each_validation_time"]
    assert lock["validation_seed"] == constraints["validation"]["seed"]
    assert lock["held_out_points"] == constraints["validation"]["held_out_points"]
    assert lock["validation_times"] == constraints["validation"]["times"]
    assert lock["derivative_steps"] == constraints["validation"]["derivative_steps"]
    assert lock["divergence_max"] == constraints["validation"]["thresholds"]["divergence_max"]
    assert lock["divergence_L2"] == constraints["validation"]["thresholds"]["divergence_L2"]
    assert lock["pde_residual_max"] == constraints["validation"]["thresholds"]["pde_residual_max"]
    assert lock["pde_residual_L2"] == constraints["validation"]["thresholds"]["pde_residual_L2"]
    assert lock["free_or_residual_dependent_force_forbidden"] is True
    assert lock["collapsed_velocity_success_forbidden"] is True
    assert lock["threshold_relaxation_forbidden_without_new_experiment_version"] is True


def test_open_temporal_capacity_cannot_silently_replace_live_route_or_truth_states():
    governance = _load(GOVERNANCE)
    status = _load(PROJECT_STATUS)
    routing = governance["semantic_guards"]["routing"]
    states = governance["state_independence"]

    assert routing["current_live_next_task_must_remain_authoritative_until_explicitly_changed"] is True
    assert "COMPACT_C4_ODD_Z_POLOIDAL" in status["next_integration_task"]
    assert status["states"]["velocity_export_ready"] is True
    assert states["velocity_export_ready_can_be_true_while_pde_validated_is_false"] is True
    assert states["visual_correspondence_verified"] is False
    assert states["pde_validated"] is False
    assert states["paper_exact"] is False
    assert states["openai_field_identified"] is False
    assert states["blowup_proved"] is False
    assert governance["semantic_guards"]["benchmarking"]["st006_superiority_requires_protocol_aligned_replay"] is True
