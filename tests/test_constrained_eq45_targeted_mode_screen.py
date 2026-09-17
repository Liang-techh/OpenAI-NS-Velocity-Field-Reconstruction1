import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_targeted_mode_screen import (
    TRUTH_BOUNDARY,
    extend_eq45_basis,
    screen_eq45_targeted_modes,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "artifacts/constrained/eq45_targeted_mode_screen_seed.json"


def _extended_candidate(candidate):
    return Eq45VelocityCandidate(
        profile_basis=extend_eq45_basis(candidate.profile_basis),
        h=candidate.h,
        time_start=candidate.time_start,
        time_end=candidate.time_end,
    )


def test_zero_extension_preserves_frozen_public_velocity():
    candidate = Eq45VelocityCandidate.seed()
    extended = _extended_candidate(candidate)
    points = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [0.35, -0.15, 0.20],
            [-0.70, 0.45, -0.40],
            [1.10, 0.20, 0.65],
            [-0.50, -0.80, -0.70],
        ]
    )
    times = np.asarray([0.25, 0.33, 0.50, 0.62, 0.75])
    np.testing.assert_allclose(
        extended.at_points(points, times),
        candidate.at_points(points, times),
        rtol=0.0,
        atol=0.0,
    )


def test_seed_targeted_mode_screen_has_four_stable_independent_directions():
    report = screen_eq45_targeted_modes(Eq45VelocityCandidate.seed())
    assert report.mode_labels == (
        "phi_radial_i2_j0",
        "phi_axial_i0_j4",
        "swirl_radial_i2_j0",
        "swirl_axial_i0_j4",
    )
    assert report.feature_names == (
        "radial_thickness_ratio",
        "axial_reach_ratio",
        "shoulder_swirl_energy_fraction",
        "shoulder_axial_poloidal_fraction",
    )
    assert report.numerical_rank == 4
    assert report.condition_number < 2.5
    assert report.max_abs_column_cosine < 0.65
    assert report.derivative_refinement_relative_change < 1e-5
    assert report.response_matrix.shape == (12, 4)

    phi_axial = report.response_matrix[:, 1].reshape(3, 4)
    swirl_axial = report.response_matrix[:, 3].reshape(3, 4)

    # The even eta^4 additions vanish on the z=0 radial-thickness probe.
    np.testing.assert_allclose(phi_axial[:, 0], 0.0, atol=1e-12)
    np.testing.assert_allclose(swirl_axial[:, 0], 0.0, atol=1e-12)

    # Phi eta^4 independently moves axial reach and axial-poloidal content.
    assert np.all(phi_axial[:, 1] > 0.0)
    assert np.all(phi_axial[:, 3] > 0.0)

    # F eta^4 independently restores/changes shoulder swirl without changing
    # the poloidal split itself.
    assert np.all(swirl_axial[:, 2] > 0.0)
    np.testing.assert_allclose(swirl_axial[:, 3], 0.0, atol=1e-10)


def test_checked_seed_screen_artifact_replays_default_diagnostic():
    checked = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    report = screen_eq45_targeted_modes(Eq45VelocityCandidate.seed())

    assert checked["schema"] == "eq45_targeted_mode_screen_v1"
    assert checked["truth_boundary"] == TRUTH_BOUNDARY
    assert report.to_dict()["truth_boundary"] == TRUTH_BOUNDARY
    np.testing.assert_allclose(
        report.baseline_features,
        np.asarray(checked["baseline_features"]),
        rtol=2e-10,
        atol=2e-12,
    )
    np.testing.assert_allclose(
        report.response_matrix,
        np.asarray(checked["response_matrix"]),
        rtol=2e-8,
        atol=2e-10,
    )
    np.testing.assert_allclose(
        report.singular_values,
        np.asarray(checked["singular_values"]),
        rtol=2e-8,
        atol=2e-10,
    )
    assert report.numerical_rank == checked["numerical_rank"] == 4
    assert report.condition_number == pytest.approx(
        checked["condition_number"], rel=2e-8
    )
    assert report.derivative_refinement_relative_change == pytest.approx(
        checked["derivative_refinement_relative_change"], rel=5e-5
    )


def test_targeted_mode_screen_fails_closed_on_invalid_growth_controls():
    candidate = Eq45VelocityCandidate.seed()
    with pytest.raises(ValueError):
        extend_eq45_basis(
            candidate.profile_basis,
            perturbation=("phi", (0, 0), 0.1),
        )
    with pytest.raises(ValueError):
        screen_eq45_targeted_modes(candidate, times=(0.25, 0.50))
    with pytest.raises(ValueError):
        screen_eq45_targeted_modes(candidate, coefficient_delta=0.0)
    with pytest.raises(ValueError):
        screen_eq45_targeted_modes(candidate, coefficient_delta=1e-3, refinement_delta=1e-3)
