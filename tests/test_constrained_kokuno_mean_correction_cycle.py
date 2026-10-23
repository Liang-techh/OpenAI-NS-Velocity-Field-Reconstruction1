import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_mean_correction_cycle import (
    StressShapedMeanSwirlCorrection,
    _independent_divergence,
    _sample_cylindrical_points,
)
from openai_ns_reconstruction.kokuno_radial_stress import (
    CompactRadialStressInverse,
    RadialMeanDefectProfile,
)


def _manufactured_inverse():
    radii = np.linspace(0.12, 0.72, 65)
    s = 2.0 * (radii - 0.12) / 0.60 - 1.0
    gate = np.maximum(1.0 - s * s, 0.0) ** 5
    values = gate * (1.0 + 0.3 * s)
    profile = RadialMeanDefectProfile(
        radii=radii,
        raw_values=values,
        values=values,
        exponent=2,
        component="theta",
        time=0.5,
        z=0.0,
        angular_count=16,
        annulus=(0.12, 0.72),
    )
    return CompactRadialStressInverse(profile)


def test_stress_shaped_mean_swirl_is_compact_nonzero_and_divergence_free():
    correction = StressShapedMeanSwirlCorrection(_manufactured_inverse(), amplitude=0.02)
    points = np.array([
        [0.30, 0.00, 0.00],
        [0.00, 0.45, 0.10],
        [0.08, 0.00, 0.00],
        [0.80, 0.00, 0.00],
        [0.35, 0.00, 0.70],
    ])
    values = correction.at_points(points, 0.5)
    assert np.isfinite(values).all()
    assert np.linalg.norm(values[0]) > 0.0
    assert np.linalg.norm(values[1]) > 0.0
    assert np.all(values[2:] == 0.0)

    audit_points = _sample_cylindrical_points(
        seed=101, count=48, radial_bounds=(0.20, 0.64), z_half_width=0.35
    )
    divergence = _independent_divergence(correction, audit_points, step=1.0e-5)
    assert np.max(np.abs(divergence)) < 2.0e-8


def test_sampling_is_deterministic_bounded_and_seed_separated():
    a = _sample_cylindrical_points(seed=11, count=32)
    b = _sample_cylindrical_points(seed=11, count=32)
    c = _sample_cylindrical_points(seed=12, count=32)
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)
    radius = np.hypot(a[:, 0], a[:, 1])
    assert np.all((radius >= 0.16) & (radius <= 0.68))
    assert np.all(np.abs(a[:, 2]) <= 0.44)


def test_preregistered_amplitude_bound_fails_closed():
    inverse = _manufactured_inverse()
    with pytest.raises(ValueError, match="amplitude"):
        StressShapedMeanSwirlCorrection(inverse, amplitude=0.041)
