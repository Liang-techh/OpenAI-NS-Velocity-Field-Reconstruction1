import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_optimized_artifact import load_checked_candidate
from openai_ns_reconstruction.constrained_eq45_optimized_mixed_taper_capacity import (
    TASK_ID,
    _with_existing_mixed_delta,
    audit_eq45_optimized_mixed_taper_capacity,
)


def test_optimized_mixed_taper_capacity_calibration():
    report = audit_eq45_optimized_mixed_taper_capacity()

    assert report["task_id"] == TASK_ID
    assert report["candidate_sha256"] == (
        "b34902ee28b3965f245eceb5867df7409498063cb0ea80f1758931723b303c98"
    )
    assert report["baseline_parameter_values"] == pytest.approx([0.1, 0.05])
    assert report["mixed_numerical_rank"] == 2
    assert report["mixed_singular_values"] == pytest.approx(
        [0.007716335842830073, 0.002551698317418577], rel=2e-6, abs=1e-12
    )
    assert report["mixed_condition_number"] == pytest.approx(3.02400005132123, rel=2e-6)
    assert report["step_refinement_relative_change"] < 1e-5

    assert report["existing_radial_span_rank"] == 2
    assert report["combined_10_12_rank"] == 4
    assert report["combined_10_12_singular_values"] == pytest.approx(
        [
            0.15125077953195082,
            0.07982551057835283,
            0.006928954374582059,
            0.002389289889501695,
        ],
        rel=2e-6,
        abs=1e-12,
    )
    assert 60.0 < report["combined_10_12_condition_number"] < 66.0

    phi, swirl = report["channel_summary"]
    assert phi["parameter"] == "Phi(1,2)"
    assert swirl["parameter"] == "F(1,2)"
    assert 0.88 < phi["novelty_fraction_outside_existing_10_span"] < 0.91
    assert 0.97 < swirl["novelty_fraction_outside_existing_10_span"] < 1.0
    assert 0.74 < phi["axial_reach_response_fraction"] < 0.77
    assert 0.16 < swirl["axial_reach_response_fraction"] < 0.19
    assert phi["radial_thickness_response_fraction"] == pytest.approx(0.0, abs=1e-14)
    assert swirl["radial_thickness_response_fraction"] == pytest.approx(0.0, abs=1e-14)

    assert report["frozen_holdout"]["baseline_curl_rms_after_force"] == pytest.approx(
        9.335838456072842, rel=1e-10
    )
    residual_rows = report["residual_sensitivity_at_finest_coefficient_step"]
    assert max(abs(row["best_fractional_holdout_change"]) for row in residual_rows) < 2e-5
    assert max(abs(row["worst_fractional_holdout_change"]) for row in residual_rows) < 2e-5
    assert max(row["plus_velocity_relative_rms_change"] for row in residual_rows) < 2e-5
    assert max(row["minus_velocity_relative_rms_change"] for row in residual_rows) < 2e-5


def test_mixed_mode_perturbation_is_temporary_and_bounded():
    candidate = load_checked_candidate()
    basis = candidate.profile_basis
    index = basis.mode_indices.index((1, 2))

    plus = _with_existing_mixed_delta(candidate, channel="phi", delta=0.02)
    assert plus.profile_basis.phi_coefficients[index] == pytest.approx(
        basis.phi_coefficients[index] + 0.02
    )
    assert candidate.profile_basis.phi_coefficients[index] == pytest.approx(0.1)
    assert plus.sha256 != candidate.sha256

    with pytest.raises(ValueError, match="registered coefficient bound"):
        _with_existing_mixed_delta(candidate, channel="phi", delta=4.0)
    with pytest.raises(ValueError, match="channel"):
        _with_existing_mixed_delta(candidate, channel="bad", delta=0.01)
    with pytest.raises(ValueError, match="finite"):
        _with_existing_mixed_delta(candidate, channel="phi", delta=np.nan)


def test_mixed_taper_truth_boundary_stays_nonpromotional():
    report = audit_eq45_optimized_mixed_taper_capacity(coefficient_steps=(0.03, 0.015))
    assert report["task_id"] == TASK_ID
    assert report["candidate_sha256"] == load_checked_candidate().sha256
    assert report["baseline_parameter_values"] == pytest.approx([0.1, 0.05])
    assert report["mixed_numerical_rank"] == 2
    assert report["combined_10_12_rank"] == 4
    assert report["frozen_holdout"]["force_refit"] is False
    for channel in report["channel_summary"]:
        assert 0.0 <= channel["novelty_fraction_outside_existing_10_span"] <= 1.0
    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["production_coefficients_changed"] is False
    assert truth["mixed_mode_activated_in_production"] is False
    assert truth["existing_candidate_artifact_modified"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert truth["blowup_proved"] is False
