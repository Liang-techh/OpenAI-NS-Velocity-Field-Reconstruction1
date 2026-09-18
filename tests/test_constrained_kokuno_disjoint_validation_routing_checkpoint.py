from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.kokuno_disjoint_validation_routing_checkpoint import (
    AGENT3_RECEIPT,
    SCHEMA,
    ST006_BASELINE,
    TASK_ID,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
)
from openai_ns_reconstruction.kokuno_multislice_covariance_preflight import HELD_OUT_CYCLE_TIMES
from openai_ns_reconstruction.kokuno_spacetime_covariance_preflight import SPATIAL_SCREEN_Z
from openai_ns_reconstruction.kokuno_spacetime_damping_holdout import VALIDATION_TIMES, VALIDATION_Z


def test_checkpoint_identity_and_hash_are_deterministic() -> None:
    first = build_checkpoint()
    second = build_checkpoint()
    assert first == second
    assert first["schema"] == SCHEMA
    assert first["task_id"] == TASK_ID
    assert first["checkpoint_sha256"] == checkpoint_sha256(first)
    validate_checkpoint(first)


def test_disjoint_validation_receipt_matches_executable_agent3_contract() -> None:
    assert AGENT3_RECEIPT["calibration_times"] == list(HELD_OUT_CYCLE_TIMES)
    assert AGENT3_RECEIPT["calibration_z"] == list(SPATIAL_SCREEN_Z)
    assert AGENT3_RECEIPT["validation_times"] == list(VALIDATION_TIMES)
    assert AGENT3_RECEIPT["validation_z"] == list(VALIDATION_Z)
    calibration = {
        (t, z) for z in AGENT3_RECEIPT["calibration_z"] for t in AGENT3_RECEIPT["calibration_times"]
    }
    validation = {
        (t, z) for z in AGENT3_RECEIPT["validation_z"] for t in AGENT3_RECEIPT["validation_times"]
    }
    assert calibration.isdisjoint(validation)
    assert AGENT3_RECEIPT["heldout_damping_reoptimized"] is False
    assert AGENT3_RECEIPT["validation_objective_used_for_selection"] is False
    assert AGENT3_RECEIPT["validation_transfer_preflight_passed"] is True


def test_local_damping_transfer_does_not_unlock_correction_or_pde() -> None:
    checkpoint = build_checkpoint()
    states = checkpoint["states"]
    a3 = checkpoint["upstream"]["agent3"]
    assert states["correction_calibration_validation_disjoint"] is True
    assert states["local_damping_transfer_preflight_passed"] is True
    assert a3["shared_selected_damping"] == pytest.approx(0.4926349922838171, rel=0, abs=0)
    assert a3["validation_fixed_damping_relative_stress_residual_rms"] == pytest.approx(
        0.6223757026222282, rel=0, abs=0
    )
    assert a3["validation_full_step_relative_stress_residual_rms"] == pytest.approx(
        0.9861371010802802, rel=0, abs=0
    )
    assert a3["rank_two_nodes"] == 0
    assert a3["required_rank_two_nodes"] == 25
    assert a3["upstream_bounded_inverse_passed"] is False
    assert a3["finite_correction_cycle_rerun_allowed"] is False
    assert states["correction_ready"] is False
    assert states["finite_correction_cycle_run"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False


def test_new_source_ingest_progress_is_recorded_without_ready_promotion() -> None:
    checkpoint = build_checkpoint()
    a1 = checkpoint["upstream"]["agent1"]
    a2 = checkpoint["upstream"]["agent2"]
    states = checkpoint["states"]
    assert a1["appendix_B_Xi_profile_executable"] is True
    assert a1["incoming_five_moment_discrepancy_bound"] is False
    assert a1["source_T_sh_lower_bound_verified"] is False
    assert a1["global_leading_ready"] is False
    assert a2["source_band_covering_schedule_executable"] is True
    assert a2["actual_positive_order_background_bound"] is False
    assert a2["actual_auxiliary_torus_mode_family_bound"] is False
    assert a2["public_q_scaled_by_beta_total_xyz_t_velocity_ready"] is False
    assert states["leading_ready"] is False
    assert states["oscillatory_ready"] is False


def test_st006_and_formal_gates_remain_unchanged() -> None:
    checkpoint = build_checkpoint()
    baseline = checkpoint["baseline_vs_kokuno"]["st006"]
    assert baseline == ST006_BASELINE
    assert baseline["momentum_sampled_max"] == pytest.approx(0.1082289305112118, rel=0, abs=0)
    assert baseline["momentum_volume_l2"] == pytest.approx(0.10758432876230622, rel=0, abs=0)
    assert checkpoint["formal_gates"]["held_out_normalized_momentum_max"] == 1.0e-3
    assert checkpoint["formal_gates"]["held_out_divergence_max"] == 1.0e-5
    assert checkpoint["formal_gates"]["changed_this_round"] is False
    assert checkpoint["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_truth_boundary_tampering_is_rejected() -> None:
    checkpoint = build_checkpoint()

    promoted = deepcopy(checkpoint)
    promoted["states"]["pde_validated"] = True
    promoted["checkpoint_sha256"] = checkpoint_sha256(promoted)
    with pytest.raises(ValueError, match="pde_validated"):
        validate_checkpoint(promoted)

    overlapped = deepcopy(checkpoint)
    overlapped["upstream"]["agent3"]["validation_times"] = [0.5]
    overlapped["upstream"]["agent3"]["validation_z"] = [0.08]
    overlapped["checkpoint_sha256"] = checkpoint_sha256(overlapped)
    with pytest.raises(ValueError, match="overlap"):
        validate_checkpoint(overlapped)

    retuned = deepcopy(checkpoint)
    retuned["upstream"]["agent3"]["heldout_damping_reoptimized"] = True
    retuned["checkpoint_sha256"] = checkpoint_sha256(retuned)
    with pytest.raises(ValueError, match="retune"):
        validate_checkpoint(retuned)

    relaxed = deepcopy(checkpoint)
    relaxed["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    relaxed["checkpoint_sha256"] = checkpoint_sha256(relaxed)
    with pytest.raises(ValueError, match="momentum gate"):
        validate_checkpoint(relaxed)

    instantiated = deepcopy(checkpoint)
    instantiated["candidate_artifact_contract"]["instantiated"] = True
    instantiated["checkpoint_sha256"] = checkpoint_sha256(instantiated)
    with pytest.raises(ValueError, match="candidate artifact"):
        validate_checkpoint(instantiated)
