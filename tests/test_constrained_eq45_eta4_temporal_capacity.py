import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_eta4_temporal_capacity import (
    audit_eq45_eta4_temporal_capacity,
)


CALIBRATION = (
    Path(__file__).resolve().parents[1]
    / "artifacts"
    / "constrained"
    / "eq45_eta4_temporal_capacity_seed.json"
)


def test_one_affine_phi_eta4_slope_adds_stable_temporal_aspect_capacity():
    candidate = Eq45VelocityCandidate.seed()
    original_sha = candidate.sha256
    checked = json.loads(CALIBRATION.read_text(encoding="utf-8"))

    report = audit_eq45_eta4_temporal_capacity(
        candidate,
        grid_size=17,
        coefficient_steps=(0.04, 0.02),
    )

    assert candidate.sha256 == original_sha
    assert report["candidate_sha256"] == original_sha == checked["candidate_sha256"]
    assert report["parameter_labels"] == checked["parameter_labels"]
    assert report["existing_static_parameter_count"] == 2
    assert report["proposed_new_parameter_count"] == 1
    assert np.allclose(report["normalized_time_tau"], [-1.0, 0.0, 1.0])

    response = np.asarray(report["finest_relative_aspect_response"], dtype=float)
    expected = np.asarray(checked["finest_relative_aspect_response"], dtype=float)
    assert response.shape == (3, 3)
    assert np.allclose(response, expected, rtol=5e-8, atol=1e-12)

    tau = np.asarray(report["normalized_time_tau"], dtype=float)
    assert np.allclose(response[:, 2], tau * response[:, 0], rtol=5e-8, atol=1e-12)
    assert report["numerical_rank"] == checked["numerical_rank"] == 3
    assert report["static_numerical_rank"] == 2
    assert np.allclose(
        report["singular_values"], checked["singular_values"], rtol=5e-8, atol=1e-12
    )
    assert np.allclose(
        report["static_singular_values"],
        checked["static_singular_values"],
        rtol=5e-8,
        atol=1e-12,
    )
    assert np.isclose(
        report["condition_number"], checked["condition_number"], rtol=5e-8
    )
    assert np.isclose(
        report["static_condition_number"],
        checked["static_condition_number"],
        rtol=5e-8,
    )
    assert report["condition_number"] < 100.0
    assert report["condition_number"] > report["static_condition_number"]

    novelty = report["temporal_slope_novelty_fraction"]
    assert np.isclose(novelty, checked["temporal_slope_novelty_fraction"], rtol=5e-8)
    assert 0.15 < novelty < 0.30
    assert report["coarse_to_finest_response_frobenius_change"] < 1.0e-3

    truth = report["truth_boundary"]
    assert truth == checked["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["time_dependent_lift_implemented"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert truth["blowup_proved"] is False


def test_eta4_temporal_capacity_fails_closed_on_invalid_contract_inputs():
    candidate = Eq45VelocityCandidate.seed()

    with pytest.raises(ValueError, match="odd"):
        audit_eq45_eta4_temporal_capacity(candidate, grid_size=18)
    with pytest.raises(ValueError, match="at least three"):
        audit_eq45_eta4_temporal_capacity(candidate, times=(0.25, 0.75))
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_eq45_eta4_temporal_capacity(candidate, times=(0.25, 0.5, 0.5))
    with pytest.raises(ValueError, match="delivery interval"):
        audit_eq45_eta4_temporal_capacity(candidate, times=(0.2, 0.5, 0.75))
    with pytest.raises(ValueError, match="strictly decreasing"):
        audit_eq45_eta4_temporal_capacity(candidate, coefficient_steps=(0.01, 0.02))
    with pytest.raises(ValueError, match="coefficient bound"):
        audit_eq45_eta4_temporal_capacity(
            candidate,
            coefficient_steps=(candidate.profile_basis.coefficient_limit + 0.1, 0.01),
        )
