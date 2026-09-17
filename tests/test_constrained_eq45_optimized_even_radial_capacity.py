import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_optimized_artifact import (
    OPTIMIZED_CANDIDATE_SHA256,
    load_checked_candidate,
)
from openai_ns_reconstruction.constrained_eq45_optimized_even_radial_capacity import (
    PARAMETER_LABELS,
    TASK_ID,
    _with_existing_radial_delta,
    audit_eq45_optimized_even_radial_capacity,
)


def test_optimized_existing_radial_capacity_is_reproducible_and_truth_bounded(capsys):
    report = audit_eq45_optimized_even_radial_capacity(
        coefficient_steps=(0.04, 0.02),
        n_angles=12,
    )

    assert report["task_id"] == TASK_ID
    assert report["candidate_sha256"] == OPTIMIZED_CANDIDATE_SHA256
    assert report["parameter_labels"] == list(PARAMETER_LABELS)
    assert report["baseline_parameter_values"] == pytest.approx([-0.3, -0.2], abs=1e-15)

    response = np.asarray(report["finest_morphology_response_per_unit_coefficient"])
    assert response.shape == (12, 2)
    assert np.all(np.isfinite(response))
    assert report["numerical_rank"] == 2
    assert np.isfinite(report["condition_number"])
    assert report["condition_number"] < 1.0e6
    assert report["step_refinement_relative_change"] < 1.0e-2
    assert all(np.linalg.norm(response[:, column]) > 1.0e-5 for column in range(2))

    holdout = report["frozen_holdout"]
    assert holdout["seed"] == 914117
    assert holdout["sample_count"] == 16
    assert holdout["spatial_and_time_step"] == pytest.approx(0.005)
    assert holdout["force_refit"] is False
    assert holdout["baseline_curl_rms_after_force"] == pytest.approx(
        9.335838456067934, rel=1e-9, abs=1e-9
    )

    sensitivities = report["residual_sensitivity_at_finest_coefficient_step"]
    assert [row["parameter"] for row in sensitivities] == list(PARAMETER_LABELS)
    for row in sensitivities:
        for key in (
            "morphology_response_l2_per_unit",
            "plus_holdout_rms",
            "minus_holdout_rms",
            "best_fractional_holdout_change",
            "worst_fractional_holdout_change",
            "central_holdout_rms_derivative",
            "plus_velocity_relative_rms_change",
            "minus_velocity_relative_rms_change",
        ):
            assert np.isfinite(row[key])
        assert row["morphology_response_l2_per_unit"] > 0.0
        assert row["plus_velocity_relative_rms_change"] > 0.0
        assert row["minus_velocity_relative_rms_change"] > 0.0

    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["production_coefficients_changed"] is False
    assert truth["forcing_refit_on_validation"] is False
    assert truth["new_basis_added"] is False
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False

    with capsys.disabled():
        print("AGENT7_OPTIMIZED_EVEN_RADIAL_REPORT=" + json.dumps(report, sort_keys=True))


def test_optimized_even_radial_audit_fails_closed_on_identity_bounds_and_bad_ladders():
    candidate = load_checked_candidate()

    mutated = _with_existing_radial_delta(candidate, channel="phi", delta=0.01)
    with pytest.raises(ValueError, match="checked optimized Eq45 candidate identity"):
        audit_eq45_optimized_even_radial_capacity(mutated)

    with pytest.raises(ValueError, match="coefficient_steps"):
        audit_eq45_optimized_even_radial_capacity(coefficient_steps=(0.02, 0.04))

    with pytest.raises(ValueError, match="coefficient bound"):
        _with_existing_radial_delta(candidate, channel="phi", delta=4.31)

    with pytest.raises(ValueError, match="channel"):
        _with_existing_radial_delta(candidate, channel="bad", delta=0.01)
