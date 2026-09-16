import numpy as np
import pytest

from openai_ns_reconstruction.constrained_offgrid_sampling import (
    assert_held_out_seed,
    sobol_validation_sample,
)


def test_sobol_sample_is_reproducible_power_of_two_and_bounded():
    bounds = np.array([
        [-2.0, 2.0],
        [-2.0, 2.0],
        [-2.0, 2.0],
        [0.25, 0.75],
    ])
    a = sobol_validation_sample(bounds, m=7, seed=914027)
    b = sobol_validation_sample(bounds, m=7, seed=914027)

    assert a.count == 128
    assert a.dimension == 4
    np.testing.assert_array_equal(a.points, b.points)
    assert np.all(a.points >= bounds[:, 0])
    assert np.all(a.points <= bounds[:, 1])
    assert not a.points.flags.writeable


def test_sampler_fail_closes_and_rejects_seed_reuse():
    with pytest.raises(ValueError):
        sobol_validation_sample(np.array([[1.0, 1.0]]), m=4, seed=1)
    with pytest.raises(ValueError):
        sobol_validation_sample(np.array([[0.0, np.inf]]), m=4, seed=1)
    with pytest.raises(ValueError):
        sobol_validation_sample(np.array([[0.0, 1.0]]), m=-1, seed=1)
    with pytest.raises(ValueError):
        assert_held_out_seed(training_seed=20260916, validation_seed=20260916)

    assert_held_out_seed(training_seed=20260916, validation_seed=914027)


def test_offgrid_sobol_probe_can_expose_narrow_hotspot_missed_by_coarse_grid():
    # Deliberately place a narrow diagnostic hotspot midway between a coarse
    # 9x9 structured grid. This does NOT prove universal detection; it is a
    # regression showing why a disjoint off-grid probe is useful in addition
    # to a structured validation grid.
    center = np.array([0.3125, 0.5625])
    radius = 0.025

    axis = np.linspace(0.0, 1.0, 9)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    coarse = np.column_stack([xx.ravel(), yy.ravel()])

    def defect(points):
        d = np.linalg.norm(points - center, axis=1)
        return (d < radius).astype(float)

    assert defect(coarse).max() == 0.0

    held_out = sobol_validation_sample(
        np.array([[0.0, 1.0], [0.0, 1.0]]),
        m=10,
        seed=914027,
    )
    assert defect(held_out.points).max() == 1.0
