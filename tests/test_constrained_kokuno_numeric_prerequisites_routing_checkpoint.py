from __future__ import annotations

import json

import pytest

from openai_ns_reconstruction.kokuno_numeric_prerequisites_routing_checkpoint import (
    AGENT1_PA10_SCREEN_RECEIPT,
    AGENT2_SIGNED_MASS_RECEIPT,
    AGENT3_AUTONOMOUS_MEAN_RECEIPT,
    FORMAL_GATES,
    PREVIOUS_AGENT5_RECEIPT,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
    write_checkpoint,
)


def _rehash(payload):
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def test_v35_records_three_numeric_prerequisites_without_readiness_promotion() -> None:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    states = payload["states"]

    assert states["autonomous_signed_rectangle_geometry_independently_audited"]
    assert states["selected_shared_C_pa10_path_obstructed"]
    assert states["leading_pa10_required_band_executable"]
    assert states["candidate_signed_h_sigma_mass_executable"]
    assert states["autonomous_finite_head_mean_factor_executable"]
    assert not states["same_cycle_requested_stress_materialized"]

    for name in (
        "leading_ready",
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        assert not states[name]


def test_agent1_pa10_screen_receipt_is_green_but_only_necessary() -> None:
    a1 = build_checkpoint()["upstream"]["agent1_pa10_screen_sibling"]
    assert a1 == AGENT1_PA10_SCREEN_RECEIPT
    assert a1["head"] == "864263ac2c9c1a030105dcf209a71b5602976d23"
    assert a1["dedicated_run"] == 35405746620
    assert a1["standard_run"] == 35405782809
    assert a1["dedicated_status"] == "success"
    assert a1["standard_status"] == "success"
    assert a1["pointwise_PA10_necessary_screen_executable"]
    assert a1["required_log_E_i_open_band_executable"]
    assert a1["selected_realization_excluded"]
    assert not a1["posthoc_C_multiplier_applied"]
    assert not a1["passing_pointwise_screen_is_source_T_sh_certificate"]
    assert not a1["source_T_sh_lower_bound_verified"]
    assert not a1["global_leading_profile_reconstructed"]


def test_agent2_candidate_signed_mass_receipt_does_not_claim_source_or_velocity() -> None:
    a2 = build_checkpoint()["upstream"]["agent2_signed_mass_sibling"]
    assert a2 == AGENT2_SIGNED_MASS_RECEIPT
    assert a2["head"] == "47f918e482981e7b124b48cef1781188302c4446"
    assert a2["dedicated_run"] == 35406299226
    assert a2["standard_run"] == 35406299071
    assert a2["focused_passed"] == 19
    assert a2["source_h_sigma_integral_formula_executable"]
    assert a2["candidate_h_sigma_numerically_materializable"]
    assert not a2["actual_source_h_sigma_pulse_integrals_bound"]
    assert not a2["actual_source_pulse_samples_recovered"]
    assert not a2["actual_positive_order_background_bound"]
    assert not a2["actual_auxiliary_torus_mode_family_bound"]
    assert not a2["public_provenance_labelled_xyz_t_oscillatory_velocity_ready"]
    assert not a2["pde_validated"]


def test_agent3_autonomous_mean_factor_is_not_formal_missing_weight_or_defect() -> None:
    a3 = build_checkpoint()["upstream"]["agent3_autonomous_mean_sibling"]
    assert a3 == AGENT3_AUTONOMOUS_MEAN_RECEIPT
    assert a3["head"] == "50d1b8ca965b35d5e0a7210c22e2302cd30c20e8"
    assert a3["dedicated_run"] == 35406792841
    assert a3["standard_run"] == 35406792929
    assert a3["physical_q"] == pytest.approx(0.040625)
    assert a3["active_mask_indices"] == [4, 5]
    assert a3["autonomous_missing_weight"] == pytest.approx(0.28410624176179694)
    assert a3["autonomous_finite_head_mean_factor_executable"]
    assert not a3["formal_theorem_machine_bump_identity_claimed"]
    assert not a3["formal_missing_weight_equality_claimed"]
    assert not a3["theorem_missing_weight_materialized"]
    assert not a3["requested_stress_actual_state_values_materialized"]
    assert not a3["finite_head_mean_debt_materialized"]
    assert not a3["real_candidate_defect_consumed"]
    assert not a3["finite_correction_cycle_rerun_allowed"]


def test_v35_preserves_previous_agent5_and_formal_gates() -> None:
    payload = build_checkpoint()
    previous = payload["upstream"]["previous_agent5_v34_ancestry"]
    assert previous == PREVIOUS_AGENT5_RECEIPT
    assert previous["head"] == "d244804a4b2322127f42239803d7023557be2325"
    assert previous["artifact_id"] == 10571992884
    assert previous["dedicated_status"] == "success"
    assert previous["standard_status"] == "success"

    assert payload["formal_gates"] == FORMAL_GATES
    assert payload["formal_gates"]["held_out_normalized_momentum_max"] == 1.0e-3
    assert payload["formal_gates"]["held_out_normalized_momentum_l2"] == 1.0e-3
    assert payload["formal_gates"]["held_out_divergence_max"] == 1.0e-5
    baseline = payload["baseline_vs_kokuno"]["st006"]
    assert baseline["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert baseline["momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_v35_fails_closed_on_truth_laundering_or_readiness_promotion() -> None:
    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent1_pa10_screen_sibling"]["passing_pointwise_screen_is_source_T_sh_certificate"] = True
    with pytest.raises(ValueError, match="source T_sh"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent2_signed_mass_sibling"]["actual_source_pulse_samples_recovered"] = True
    with pytest.raises(ValueError, match="recovered source"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent3_autonomous_mean_sibling"]["theorem_missing_weight_materialized"] = True
    with pytest.raises(ValueError, match="theorem missingWeight"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["states"]["same_cycle_requested_stress_materialized"] = True
    with pytest.raises(ValueError, match="requestedStress"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["states"]["candidate_artifact_instantiated"] = True
    with pytest.raises(ValueError, match="candidate_artifact_instantiated"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="pde_validated"):
        validate_checkpoint(_rehash(payload))


def test_v35_roundtrip_is_deterministic(tmp_path) -> None:
    payload = build_checkpoint()
    path = tmp_path / "checkpoint.json"
    written = write_checkpoint(path)
    loaded = json.loads(path.read_text(encoding="utf-8"))

    assert loaded == written
    assert loaded["checkpoint_sha256"] == checkpoint_sha256(loaded)

    loaded["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    with pytest.raises(ValueError, match="checkpoint sha256 mismatch"):
        validate_checkpoint(loaded)
