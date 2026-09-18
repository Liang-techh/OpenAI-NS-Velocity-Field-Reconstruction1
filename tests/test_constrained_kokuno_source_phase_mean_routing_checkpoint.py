from __future__ import annotations

import json

import pytest

from openai_ns_reconstruction.kokuno_source_phase_mean_routing_checkpoint import (
    FORMAL_GATES,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
    write_checkpoint,
)


def test_checkpoint_promotes_only_audited_source_display_and_formal_identity() -> None:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    states = payload["states"]

    assert states["selected_source_rescaled_prefix_moments_executable"]
    assert states["displayed_source_phase_frame_bound"]
    assert states["displayed_source_phase_frame_independently_audited"]
    assert states["supplied_signed_complete_curl_physical_family_executable"]
    assert states["supplied_signed_complete_curl_physical_covariance_independently_audited"]
    assert states["formal_signed_mean_defect_identity_admitted"]

    for name in (
        "leading_ready",
        "actual_numeric_signed_mean_defect_materialized",
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_signed_auxiliary_rectangles_or_modes_bound",
        "actual_source_partition_labels_instantiated",
        "actual_source_physical_covariance_rank_two_assessed",
        "public_source_bound_xyz_t_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "signed_mean_inverse_input_ready",
        "real_candidate_defect_consumed",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        assert not states[name]


def test_checkpoint_preserves_independent_phase_frame_receipt_and_frozen_guards() -> None:
    payload = build_checkpoint()
    a4 = payload["upstream"]["agent4_ancestry"]

    assert a4["head"] == "31f7510545488818b3bcf979c5cac3cb5d826a49"
    assert a4["dedicated_run"] == 35398848549
    assert a4["standard_run"] == 35398848498
    assert a4["dedicated_status"] == "success"
    assert a4["standard_status"] == "success"
    assert a4["artifact_id"] == 10569677115
    assert a4["artifact_digest"] == "sha256:ff66fe2d26789dfa8fd5406659ad8d2ce0f9c09cab02cda15f63d742db5583af"
    assert a4["oracle_public_values_only"]
    assert not a4["oracle_reuses_agent2_fd_helper"]
    assert not a4["oracle_reuses_complete_curl_helper"]
    assert not a4["pressure_or_forcing_fit_used"]
    assert not a4["training_tensor_or_loss_used"]
    assert a4["local_structural_preflight_passed"]
    assert a4["value_relative_max"] == pytest.approx(6.661338147750939e-16)
    assert a4["frame_relative_rms"] == pytest.approx(2.7633777086542752e-17)
    assert a4["source_identity_relative_max"] == pytest.approx(6.685798127696402e-16)
    assert a4["angular_winding_relative_max"] == pytest.approx(3.3306690738754696e-16)
    assert a4["fd_relative_rms"][-1] == pytest.approx(4.3321017055747775e-12)
    assert min(a4["fd_refinement_ratios"]) >= 20.0
    assert a4["wrong_axial_normalization_relative_rms"] >= 0.20
    assert a4["wrong_tilt_mutation_relative_rms"] >= 0.30
    assert a4["carrier_k_values_observed"] == [2, 3, 4]


def test_checkpoint_records_agent3_numeric_mean_debt_blocker_without_fake_target() -> None:
    payload = build_checkpoint()
    a3 = payload["upstream"]["agent3_sibling"]

    assert a3["requested_cross_defect_identity_admitted"]
    assert a3["finite_head_band_admitted"]
    assert not a3["requested_cross_defect_theorem_machine_replayed"]
    assert not a3["actual_mean_cross_values_materialized"]
    assert not a3["missing_weight_values_materialized"]
    assert not a3["finite_head_mean_debt_materialized"]
    assert not a3["signed_mean_inverse_input_ready"]
    assert a3["physical_to_reference_rule"] == "H_ref * y = DeltaC / epsilon"
    assert not a3["real_candidate_defect_consumed"]
    assert not a3["finite_correction_cycle_rerun_allowed"]
    assert not a3["consumed_in_executable_ancestry"]


def test_checkpoint_records_agent1_progress_as_sibling_not_leading_ready() -> None:
    payload = build_checkpoint()
    a1 = payload["upstream"]["agent1_sibling"]

    assert a1["head"] == "f44182fe310cc87a010f8ae96a7c66d173a3229d"
    assert a1["dedicated_status"] == "success"
    assert a1["standard_status"] == "success"
    assert a1["selected_source_rescaled_prefix_moments_executable"]
    assert a1["shared_C_PA15_scaling_applied"]
    assert a1["PA16_input_tuple_formable"]
    assert not a1["actual_source_incoming_five_moment_discrepancy_bound"]
    assert not a1["source_T_sh_lower_bound_verified"]
    assert not a1["inner_to_outer_join_completed"]
    assert not a1["global_pressure_matched"]
    assert not a1["global_leading_profile_reconstructed"]
    assert not a1["consumed_in_executable_ancestry"]
    assert not payload["states"]["leading_ready"]


def test_checkpoint_keeps_source_hidden_values_and_formal_pde_gate_fail_closed() -> None:
    payload = build_checkpoint()
    a2 = payload["upstream"]["agent2_ancestry"]
    assert a2["displayed_source_phase_frame_bound"]
    assert a2["phase_and_n_Phi_derived_not_free_on_adapter_path"]
    assert not a2["actual_positive_order_background_bound"]
    assert not a2["actual_source_h_sigma_pulse_integrals_bound"]
    assert not a2["actual_signed_auxiliary_rectangles_or_modes_bound"]
    assert not a2["actual_source_partition_labels_instantiated"]
    assert not a2["public_source_bound_xyz_t_oscillatory_velocity_ready"]

    assert payload["formal_gates"] == FORMAL_GATES
    assert payload["formal_gates"]["held_out_normalized_momentum_max"] == 1.0e-3
    assert payload["formal_gates"]["held_out_divergence_max"] == 1.0e-5
    baseline = payload["baseline_vs_kokuno"]["st006"]
    assert baseline["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert baseline["momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_checkpoint_fails_closed_on_hidden_source_or_numeric_mean_promotion() -> None:
    payload = build_checkpoint()
    payload["states"]["oscillatory_ready"] = True
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="fail-closed state promoted"):
        validate_checkpoint(payload)

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent2_ancestry"]["actual_source_h_sigma_pulse_integrals_bound"] = True
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="displayed/source-hidden truth boundary changed"):
        validate_checkpoint(payload)

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent3_sibling"]["finite_head_mean_debt_materialized"] = True
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="opaque formal debt was promoted numerically"):
        validate_checkpoint(payload)

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent1_sibling"]["global_pressure_matched"] = True
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="leading sibling was over-promoted"):
        validate_checkpoint(payload)


def test_checkpoint_sha_roundtrip_and_receipt_copy_isolation(tmp_path) -> None:
    first = build_checkpoint()
    first["upstream"]["agent4_ancestry"]["value_relative_max"] = -1.0
    second = build_checkpoint()
    assert second["upstream"]["agent4_ancestry"]["value_relative_max"] > 0.0

    path = tmp_path / "checkpoint.json"
    payload = write_checkpoint(path)
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded == payload
    assert loaded["checkpoint_sha256"] == checkpoint_sha256(loaded)

    loaded["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    with pytest.raises(ValueError, match="checkpoint sha256 mismatch"):
        validate_checkpoint(loaded)
