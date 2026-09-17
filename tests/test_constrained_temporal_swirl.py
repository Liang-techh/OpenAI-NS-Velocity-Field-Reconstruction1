import json
from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_force import RestrictedForce
from openai_ns_reconstruction.constrained_momentum_budget import angular_moment
from openai_ns_reconstruction.constrained_outer_momentum import AngularMomentumCandidate
from openai_ns_reconstruction.constrained_temporal_swirl import (
    COEFFICIENT_COUNT,
    FAMILY_ID,
    SPATIAL_BASIS,
    TemporalSwirlCandidate,
)
from openai_ns_reconstruction.constrained_tensor_candidate import TensorCandidate


def _parent(order=48):
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
    candidate = TemporalSwirlCandidate(parent, np.linspace(-0.9, 0.9, 15))
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
    np.testing.assert_array_equal(
        candidate.velocity(points[:3], 0.5), parent.velocity(points[:3], 0.5)
    )
    np.testing.assert_array_equal(candidate.pressure(points, 0.5), parent.pressure(points, 0.5))
    assert candidate.force is parent.force

    times = np.array([0.25, 0.4, 0.6, 0.75])
    pointwise = np.column_stack(
        (0.1 * np.ones_like(times), 0.05 * np.ones_like(times), 0.2 * np.ones_like(times))
    )
    assert candidate.correction_basis(pointwise, times).shape == (4, 3, COEFFICIENT_COUNT)
    np.testing.assert_array_equal(candidate.correction_basis(points[0], 0.25), 0.0)


def test_each_spatial_mode_has_zero_angular_moment_at_selected_quadrature():
    candidate = TemporalSwirlCandidate(_parent(order=64), np.zeros(COEFFICIENT_COUNT), order=64)
    # At w=1/2, the k=1 Bernstein factor is 3/8.  Dividing one column by
    # that scalar exposes each spatial mode for the independent moment check.
    temporal_factor = 3.0 / 8.0
    for spatial_index in range(len(SPATIAL_BASIS)):
        column = 3 * spatial_index

        def mode(points, time, column=column):
            return candidate.correction_basis(points, time)[..., :, column] / temporal_factor

        assert abs(angular_moment(mode, 0.5, candidate.order)) < 2e-11


def test_temporal_swirl_is_cartesian_divergence_free():
    candidate = TemporalSwirlCandidate(_parent(order=48), np.linspace(-0.85, 0.85, 15))
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


def test_save_load_is_self_contained_and_uses_temporal_family(tmp_path):
    candidate = TemporalSwirlCandidate(_parent(order=32), np.linspace(-0.7, 0.7, 15))
    path = tmp_path / "temporal_swirl.json"
    candidate.save(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["family"] == FAMILY_ID
    assert payload["parameters"]["parent"]["family"] == "angular_momentum_outer_v1"
    assert len(payload["parameters"]["coefficients"]) == COEFFICIENT_COUNT

    loaded = TemporalSwirlCandidate.load(path)
    assert loaded == candidate
    assert isinstance(loaded.coefficients, tuple)
    points = np.array([[0.8, 0.1, 0.3], [0.2, -0.4, -0.5]])
    np.testing.assert_allclose(loaded.velocity(points, 0.58), candidate.velocity(points, 0.58))
    np.testing.assert_allclose(loaded.pressure(points, 0.58), candidate.pressure(points, 0.58))
    assert loaded.energy(0.5, order=32) == pytest.approx(candidate.energy(0.5, order=32))


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, -1.01, 1.01])
def test_coefficients_are_bounded_finite_and_immutable(bad):
    values = np.zeros(COEFFICIENT_COUNT)
    values[4] = bad
    with pytest.raises(ValueError):
        TemporalSwirlCandidate(_parent(order=24), values)

    with pytest.raises(ValueError):
        TemporalSwirlCandidate(_parent(order=24), np.zeros(COEFFICIENT_COUNT - 1))
    with pytest.raises(ValueError):
        TemporalSwirlCandidate(_parent(order=24), np.zeros((COEFFICIENT_COUNT, 1)))

    candidate = TemporalSwirlCandidate(_parent(order=24))
    with pytest.raises(FrozenInstanceError):
        candidate.coefficients = (0.0,) * COEFFICIENT_COUNT


def test_invalid_time_and_quadrature_inputs_are_rejected():
    candidate = TemporalSwirlCandidate(_parent(order=24))
    with pytest.raises(ValueError):
        candidate.velocity([0.0, 0.0, 0.0], 0.2)
    with pytest.raises(ValueError):
        candidate.energy(0.5, order=1)
