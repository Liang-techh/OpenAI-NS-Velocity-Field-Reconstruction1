import json
from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_force import RestrictedForce
from openai_ns_reconstruction.constrained_localized_swirl import (
    COEFFICIENT_COUNT,
    FAMILY_ID,
    SPATIAL_BASIS,
    ZERO_COEFFICIENTS,
    LocalizedSwirlCandidate,
)
from openai_ns_reconstruction.constrained_momentum_budget import angular_moment
from openai_ns_reconstruction.constrained_outer_momentum import AngularMomentumCandidate
from openai_ns_reconstruction.constrained_tensor_candidate import TensorCandidate


def _parent(order=96):
    base = TensorCandidate(
        radial_shape=0.1,
        axial_shape=-0.1,
        swirl_radial_shape=0.2,
        swirl_axial_shape=-0.2,
    ).normalized(order=order)
    return AngularMomentumCandidate(
        base,
        RestrictedForce(a=0.8, c=1.1),
        order=order,
        outer_shape=(-0.6, 0.35),
        pressure_coefficients=(0.15,) * 18,
    )


def test_initial_and_core_values_are_preserved():
    parent = _parent()
    candidate = LocalizedSwirlCandidate(parent, np.linspace(-0.9, 0.9, 18))
    points = np.array(
        [
            [0.1, 0.0, 0.3],
            [0.2, 0.1, -0.4],
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 0.0, 2.0],
        ]
    )
    np.testing.assert_array_equal(
        candidate.velocity(points, 0.25), parent.velocity(points, 0.25)
    )
    np.testing.assert_array_equal(candidate.velocity(points[:3], 0.5), parent.velocity(points[:3], 0.5))
    np.testing.assert_array_equal(candidate.pressure(points, 0.5), parent.pressure(points, 0.5))
    assert candidate.force is parent.force
    assert candidate.correction_basis(points[0], 0.25).shape == (3, COEFFICIENT_COUNT)
    np.testing.assert_array_equal(candidate.correction_basis(points[0], 0.25), 0.0)


def test_each_gaussian_ring_has_zero_angular_moment_at_order_96():
    candidate = LocalizedSwirlCandidate(_parent(), np.zeros(COEFFICIENT_COUNT))
    # At w=1/2, the first nonzero cubic Bernstein factor is 3/8.  Exposing
    # one column at a time checks every spatial ring independently.
    for spatial_index in range(len(SPATIAL_BASIS)):
        column = 3 * spatial_index

        def mode(points, time, column=column):
            return candidate.correction_basis(points, time)[..., :, column] / (3.0 / 8.0)

        assert abs(angular_moment(mode, 0.5, 96)) < 2e-11


def test_localized_swirl_is_cartesian_divergence_free():
    candidate = LocalizedSwirlCandidate(_parent(), np.linspace(-0.85, 0.85, 18))
    points = np.random.default_rng(741).uniform(-1.15, 1.15, size=(80, 3))
    points = np.vstack((points, [[0.0, 0.0, 0.2], [0.0, 0.0, 0.0]]))
    h = 1e-5
    divergence = sum(
        (
            candidate.velocity(points + np.eye(3)[axis] * h, 0.47)[:, axis]
            - candidate.velocity(points - np.eye(3)[axis] * h, 0.47)[:, axis]
        )
        / (2.0 * h)
        for axis in range(3)
    )
    assert np.max(np.abs(divergence)) < 4e-7


def test_save_load_roundtrip_and_bounded_coefficients(tmp_path):
    candidate = LocalizedSwirlCandidate(_parent(order=48), np.linspace(-0.7, 0.7, 18), order=48)
    path = tmp_path / "localized_swirl.json"
    candidate.save(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["family"] == FAMILY_ID
    assert payload["parameters"]["parent"]["family"] == "angular_momentum_outer_v1"
    assert len(payload["parameters"]["coefficients"]) == COEFFICIENT_COUNT
    loaded = LocalizedSwirlCandidate.load(path)
    assert loaded == candidate
    assert isinstance(loaded.coefficients, tuple)
    points = np.array([[0.8, 0.1, 0.3], [0.2, -0.4, -0.5]])
    np.testing.assert_allclose(loaded.velocity(points, 0.58), candidate.velocity(points, 0.58))
    np.testing.assert_allclose(loaded.pressure(points, 0.58), candidate.pressure(points, 0.58))
    assert loaded.energy(0.5, order=32) == pytest.approx(candidate.energy(0.5, order=32))

    for bad in (np.nan, np.inf, -np.inf, -1.01, 1.01):
        values = np.zeros(COEFFICIENT_COUNT)
        values[4] = bad
        with pytest.raises(ValueError):
            LocalizedSwirlCandidate(_parent(order=24), values)
    with pytest.raises(ValueError):
        LocalizedSwirlCandidate(_parent(order=24), np.zeros(COEFFICIENT_COUNT - 1))
    with pytest.raises(FrozenInstanceError):
        LocalizedSwirlCandidate(_parent(order=24)).coefficients = ZERO_COEFFICIENTS
