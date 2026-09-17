import json

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
    pytest.fail(json.dumps(report, sort_keys=True))


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
    assert report["mixed_numerical_rank"] >= 1
    assert report["combined_10_12_rank"] >= report["existing_radial_span_rank"]
    assert report["frozen_holdout"]["force_refit"] is False
    for channel in report["channel_summary"]:
        assert 0.0 <= channel["novelty_fraction_outside_existing_10_span"] <= 1.0
    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["production_coefficients_changed"] is False
    assert truth["mixed_mode_activated_in_production"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert truth["blowup_proved"] is False
