import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_eta4_vorticity_capacity import (
    audit_eq45_eta4_vorticity_capacity,
    vorticity_core_rms,
)


def test_eq45_eta4_has_stable_two_channel_vorticity_core_capacity():
    candidate = Eq45VelocityCandidate.seed()
    original_sha = candidate.sha256

    report = audit_eq45_eta4_vorticity_capacity(
        candidate,
        grid_size=17,
        coefficient_steps=(0.04, 0.02),
    )

    assert candidate.sha256 == original_sha
    assert report["candidate_sha256"] == original_sha
    assert report["numerical_rank"] == 2
    assert np.isfinite(report["condition_number"])
    assert report["condition_number"] < 50.0
    assert report["coarse_to_finest_response_frobenius_change"] < 1.0e-3

    aspect_response = np.asarray(report["finest_relative_aspect_response"])
    assert aspect_response.shape == (3, 2)
    assert np.all(aspect_response[:, 0] < 0.0)
    assert np.all(aspect_response[:, 1] > 0.0)
    assert np.linalg.norm(aspect_response[:, 0]) > 2.0 * np.linalg.norm(
        aspect_response[:, 1]
    )

    best = report["max_mean_aspect_corner"]
    assert best["phi_eta4"] == -candidate.profile_basis.coefficient_limit
    assert best["swirl_eta4"] == candidate.profile_basis.coefficient_limit
    assert min(best["relative_aspect_change"]) > 0.04

    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert truth["blowup_proved"] is False


def test_vorticity_core_rms_is_finite_for_public_seed_velocity():
    metric = vorticity_core_rms(
        Eq45VelocityCandidate.seed(),
        0.5,
        grid_size=17,
        box_half_width=2.0,
    )
    assert metric.radial_rms > 0.0
    assert metric.axial_rms > 0.0
    assert metric.aspect_ratio_z_over_r > 0.0
    assert metric.omega_max > 0.0
    assert metric.enstrophy_integral > 0.0
    assert np.all(
        np.isfinite(
            [
                metric.radial_rms,
                metric.axial_rms,
                metric.aspect_ratio_z_over_r,
                metric.omega_max,
                metric.enstrophy_integral,
            ]
        )
    )


def test_eta4_vorticity_capacity_fails_closed_on_invalid_contract_inputs():
    candidate = Eq45VelocityCandidate.seed()

    with pytest.raises(ValueError, match="odd"):
        vorticity_core_rms(candidate, 0.5, grid_size=18)
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_eq45_eta4_vorticity_capacity(candidate, times=(0.25, 0.5, 0.5))
    with pytest.raises(ValueError, match="strictly decreasing"):
        audit_eq45_eta4_vorticity_capacity(
            candidate, coefficient_steps=(0.01, 0.02)
        )
    with pytest.raises(ValueError, match="coefficient bound"):
        audit_eq45_eta4_vorticity_capacity(
            candidate,
            coefficient_steps=(candidate.profile_basis.coefficient_limit + 0.1, 0.01),
        )
