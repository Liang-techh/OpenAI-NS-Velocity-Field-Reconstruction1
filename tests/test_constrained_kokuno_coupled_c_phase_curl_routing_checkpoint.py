from __future__ import annotations

import json

import pytest

from openai_ns_reconstruction.kokuno_coupled_c_phase_curl_routing_checkpoint import (
    AGENT1_COUPLED_C_RECEIPT,
    AGENT2_PHASE_CURL_RECEIPT,
    FORMAL_GATES,
    PREVIOUS_AGENT5_RECEIPT,
    TYPED_HANDOFF_CONTRACTS,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
    write_checkpoint,
)


def _rehash(payload):
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def test_v36_routes_new_sibling_seams_without_readiness_promotion() -> None:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    states = payload["states"]

    # Historical old selected path remains recorded as obstructed, while a
    # separately repropagated coupled-C path now clears the pointwise screen.
    assert states["selected_shared_C_pa10_path_obstructed"]
    assert states["fresh_coupled_C_pa10_pointwise_not_excluded"]
    assert not states["fresh_coupled_C_pa16_handoff_allowed"]

    assert states["candidate_signed_h_sigma_mass_executable"]
    assert states["candidate_signed_h_sigma_mass_independently_audited"]
    assert states["candidate_mass_bound_phase_curl_composition_ready"]
    assert states["candidate_h_sigma_spatial_derivatives_explicit"]

    for name in (
        "leading_ready",
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        assert not states[name]


def test_agent1_fresh_coupled_c_pass_is_only_a_negative_obstruction_removed() -> None:
    a1 = build_checkpoint()["upstream"]["agent1_coupled_c_sibling"]
    assert a1 == AGENT1_COUPLED_C_RECEIPT
    assert a1["head"] == "c1cef3ac2b3e4d58ab4b721b31a1d2de14abba65"
    assert a1["dedicated_run"] == 35409636108
    assert a1["dedicated_status"] == "success"
    assert a1["latest_standard_run_observed"] == 35409925427
    assert a1["latest_standard_status_observed"] == "in_progress"
    assert a1["log_C_factor"] == pytest.approx(32.0)
    assert a1["upstream_C_repropagated_before_core"]
    assert not a1["posthoc_C_multiplier_applied"]
    assert a1["pointwise_PA10_necessary_screen_passed"]
    assert a1["previous_pointwise_PA10_obstruction_removed"]
    assert a1["nontrivial_G_i_guard_passed"]
    assert not a1["passing_pointwise_screen_is_source_T_sh_certificate"]
    assert not a1["source_T_sh_lower_bound_verified"]
    assert not a1["selected_pa16_handoff_allowed"]
    assert not a1["global_pressure_matched"]
    assert not a1["global_leading_profile_reconstructed"]


def test_agent2_mass_bound_phase_curl_freezes_explicit_derivative_contract() -> None:
    a2 = build_checkpoint()["upstream"]["agent2_mass_bound_phase_curl_sibling"]
    assert a2 == AGENT2_PHASE_CURL_RECEIPT
    assert a2["head"] == "296ed78b79cade56547aa9b751735207685c80bc"
    assert a2["dedicated_run"] == 35409788333
    assert a2["standard_run"] == 35409788239
    assert a2["dedicated_status"] == "success"
    assert a2["standard_status"] == "success"
    assert a2["focused_passed"] == 20
    assert a2["candidate_signed_mass_consumed_directly"]
    assert a2["manual_h_sigma_reentry_removed"]
    assert a2["explicit_D_r_h_sigma_required"]
    assert a2["explicit_D_z_h_sigma_required"]
    assert a2["frozen_local_zero_derivative_mode_requires_exact_zero"]
    assert a2["q_scaled_by_sign_cartesian_arrays_executable"]
    assert a2["q_scaled_by_beta_cartesian_arrays_executable"]
    assert a2["q_scaled_total_cartesian_array_executable"]
    assert not a2["public_xyz_t_velocity_osc_provider_ready"]
    assert not a2["actual_positive_order_background_bound"]
    assert not a2["actual_auxiliary_torus_mode_family_bound"]
    assert not a2["pde_validated"]


def test_v36_freezes_typed_handoffs_and_fixed_gates() -> None:
    payload = build_checkpoint()
    contracts = payload["typed_handoff_contracts"]
    assert contracts == TYPED_HANDOFF_CONTRACTS
    assert not contracts["leading"]["ready"]
    assert contracts["leading"]["current_semantics"] == "pointwise_PA10_not_excluded_only"
    assert not contracts["oscillatory"]["ready"]
    assert contracts["oscillatory"]["missing_public_output"] == "velocity_osc(x,y,z,t)"
    assert "D_r h_sigma" in contracts["oscillatory"]["required_inputs"]
    assert "D_z h_sigma" in contracts["oscillatory"]["required_inputs"]
    assert not contracts["correction"]["ready"]
    assert not contracts["independent_validator"]["ready_for_final_gate"]

    assert payload["formal_gates"] == FORMAL_GATES
    assert payload["formal_gates"]["held_out_normalized_momentum_max"] == 1.0e-3
    assert payload["formal_gates"]["held_out_normalized_momentum_l2"] == 1.0e-3
    assert payload["formal_gates"]["held_out_divergence_max"] == 1.0e-5
    baseline = payload["baseline_vs_kokuno"]["st006"]
    assert baseline["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert baseline["momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_v36_preserves_previous_agent5_and_fails_closed_on_laundering() -> None:
    payload = build_checkpoint()
    previous = payload["upstream"]["previous_agent5_v35_ancestry"]
    assert previous == PREVIOUS_AGENT5_RECEIPT
    assert previous["head"] == "02d3c2b5663434101d9517a82d5cabe0fb667df6"
    assert previous["artifact_id"] == 10572867184
    assert previous["dedicated_status"] == "success"
    assert previous["standard_status"] == "success"

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent1_coupled_c_sibling"]["selected_pa16_handoff_allowed"] = True
    with pytest.raises(ValueError, match="receipt changed"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent2_mass_bound_phase_curl_sibling"]["explicit_D_r_h_sigma_required"] = False
    with pytest.raises(ValueError, match="receipt changed"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["states"]["public_provenance_labelled_xyz_t_oscillatory_velocity_ready"] = True
    with pytest.raises(ValueError, match="fail-closed state promoted"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="fail-closed state promoted"):
        validate_checkpoint(_rehash(payload))


def test_v36_roundtrip_is_deterministic(tmp_path) -> None:
    payload = build_checkpoint()
    path = tmp_path / "checkpoint.json"
    written = write_checkpoint(path)
    loaded = json.loads(path.read_text(encoding="utf-8"))

    assert loaded == written
    assert loaded["checkpoint_sha256"] == checkpoint_sha256(loaded)

    loaded["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    with pytest.raises(ValueError, match="checkpoint sha256 mismatch"):
        validate_checkpoint(loaded)
