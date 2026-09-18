from __future__ import annotations

import json

import pytest

from openai_ns_reconstruction.kokuno_signed_covariance_routing_checkpoint import (
    FORMAL_GATES,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
    write_checkpoint,
)


def test_checkpoint_routes_signed_source_axis_without_physical_rank_promotion() -> None:
    payload = build_checkpoint()
    validate_checkpoint(payload)

    states = payload["states"]
    assert states["displayed_source_signed_covariance_reference_executable"]
    assert states["signed_covariance_reference_independently_audited"]
    assert states["signed_covariance_physical_to_reference_unit_bridge_ready"]
    assert not states["physical_complete_curl_signed_family_bound"]
    assert not states["genuinely_independent_second_covariance_column_ready"]
    assert not states["correction_ready"]
    assert not states["pde_validated"]

    a2 = payload["upstream"]["agent2"]
    assert a2["source_covariance_axis"] == "auxiliary_rectangle_sign_sigma_plus_minus"
    assert not a2["source_covariance_axis_is_dyadic_band_contrast"]
    assert not a2["source_covariance_axis_is_fourier_harmonic_sign"]

    a4 = payload["upstream"]["agent4"]
    assert a4["local_structural_preflight_passed"]
    assert a4["reference_value_max_relative_error"] < 5.0e-13
    assert a4["directional_derivative_relative_rms"][-1] < 1.0e-6
    assert min(a4["directional_derivative_refinement_ratios"]) > 8.0
    assert a4["semantic_sigma_swap_relative_max"] == 0.0
    assert not a4["physical_complete_curl_signed_family_audited"]


def test_checkpoint_preserves_epsilon_unit_bridge_truth_boundary() -> None:
    payload = build_checkpoint()
    a3 = payload["upstream"]["agent3"]
    assert a3["physical_covariance_relation"] == "DeltaC = epsilon * H_ref * y"
    assert a3["reference_inverse_relation"] == "H_ref * y = DeltaC / epsilon"
    assert a3["epsilon_division_required"]
    assert a3["structural_calibration_only"]
    assert not a3["real_candidate_defect_consumed"]
    assert a3["structural_cone_margin_fraction"] == pytest.approx(0.8)
    assert a3["structural_reference_condition_number"] == pytest.approx(1.5)


def test_checkpoint_keeps_formal_gates_and_st006_comparison_pending() -> None:
    payload = build_checkpoint()
    assert payload["formal_gates"] == FORMAL_GATES
    assert payload["formal_gates"]["held_out_normalized_momentum_max"] == 1.0e-3
    assert payload["formal_gates"]["held_out_divergence_max"] == 1.0e-5

    baseline = payload["baseline_vs_kokuno"]["st006"]
    assert baseline["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert baseline["momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert not baseline["pde_validated"]
    assert payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_checkpoint_fails_closed_on_rank_promotion_or_missing_epsilon_rule() -> None:
    payload = build_checkpoint()
    payload["states"]["genuinely_independent_second_covariance_column_ready"] = True
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="fail-closed state promoted"):
        validate_checkpoint(payload)

    # build_checkpoint intentionally exposes immutable receipt constants by value
    # convention; deep-copy here so this negative-control mutation cannot leak
    # into later tests through the module-level receipt object.
    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent3"]["epsilon_division_required"] = False
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="unit/truth boundary changed"):
        validate_checkpoint(payload)


def test_checkpoint_sha_and_roundtrip(tmp_path) -> None:
    path = tmp_path / "checkpoint.json"
    payload = write_checkpoint(path)
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded == payload
    assert loaded["checkpoint_sha256"] == checkpoint_sha256(loaded)

    loaded["upstream"]["agent4"]["physical_complete_curl_signed_family_audited"] = True
    with pytest.raises(ValueError, match="checkpoint sha256 mismatch"):
        validate_checkpoint(loaded)
