import json
from dataclasses import FrozenInstanceError, asdict

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_candidate import CompactCandidate
from openai_ns_reconstruction.constrained_tensor_candidate import (
    FAMILY_ID,
    TENSOR_BASIS,
    TensorCandidate,
)


def _base_candidate() -> CompactCandidate:
    return CompactCandidate(
        radial_shape=0.13,
        axial_shape=-0.17,
        swirl_radial_shape=0.4,
        swirl_axial_shape=-0.6,
        swirl_radial_time=1.2,
        swirl_axial_time=-0.7,
        poloidal_radial_time=1.1,
        poloidal_axial_time=-0.9,
    ).normalized(order=64)


def test_zero_tensor_reproduces_base_and_normalizes():
    base = _base_candidate()
    candidate = TensorCandidate(**asdict(base))
    points = np.random.default_rng(25).uniform(-1.5, 1.5, (80, 3))
    points = np.vstack((points, [[0.0, 0.0, 0.3], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]))

    np.testing.assert_array_equal(candidate.velocity(points, 0.4), base.velocity(points, 0.4))
    np.testing.assert_array_equal(
        candidate.vector_potential(points, 0.4), base.vector_potential(points, 0.4)
    )
    normalized = candidate.normalized(order=64)
    assert abs(normalized.energy(order=96) - 1.0) < 1e-8


def test_nonzero_tensor_agrees_with_legacy_numerical_curl():
    coefficients = np.zeros(27)
    coefficients[TENSOR_BASIS.index((1, 1, 1))] = 0.6
    coefficients[TENSOR_BASIS.index((2, 0, 2))] = -0.35
    swirl = np.zeros(27)
    swirl[TENSOR_BASIS.index((0, 1, 2))] = -0.4
    candidate = TensorCandidate(
        **asdict(_base_candidate()),
        poloidal_coefficients=coefficients,
        swirl_coefficients=swirl,
    )
    points = np.random.default_rng(45).uniform(-1.2, 1.2, (24, 3))
    points = np.vstack((points, [[0.0, 0.0, 0.3], [0.0, 0.0, 0.0]]))
    legacy = candidate.as_legacy_local_field()
    numerical = np.array(
        [legacy.velocity(*point, 0.5, eps=1e-5) for point in points]
    )

    np.testing.assert_allclose(
        numerical,
        candidate.velocity(points, 0.5),
        atol=2e-8,
        rtol=2e-6,
    )


def test_tensor_velocity_is_divergence_free_in_cartesian_coordinates():
    coefficients = np.linspace(-0.8, 0.8, 27)
    swirl = np.linspace(0.7, -0.7, 27)
    candidate = TensorCandidate(
        **asdict(_base_candidate()),
        poloidal_coefficients=coefficients,
        swirl_coefficients=swirl,
    )
    points = np.random.default_rng(89).uniform(-1.1, 1.1, (90, 3))
    points = np.vstack((points, [[0.0, 0.0, 0.2], [0.0, 0.0, 0.0]]))
    h = 1e-5
    divergence = sum(
        (
            candidate.velocity(points + np.eye(3)[axis] * h, 0.4)[:, axis]
            - candidate.velocity(points - np.eye(3)[axis] * h, 0.4)[:, axis]
        )
        / (2.0 * h)
        for axis in range(3)
    )
    assert np.max(np.abs(divergence)) < 3e-7


def test_tensor_artifact_round_trip_has_new_family_and_tuples(tmp_path):
    coefficients = np.zeros(27)
    coefficients[0] = 0.25
    coefficients[-1] = -0.5
    candidate = TensorCandidate(
        **asdict(_base_candidate()),
        poloidal_coefficients=coefficients,
        swirl_coefficients=tuple(-entry for entry in coefficients),
    )
    path = tmp_path / "tensor_candidate.json"
    candidate.save(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["family"] == FAMILY_ID
    assert len(payload["parameters"]["poloidal_coefficients"]) == 27

    loaded = TensorCandidate.load(path)
    assert loaded == candidate
    assert isinstance(loaded.poloidal_coefficients, tuple)
    assert isinstance(loaded.swirl_coefficients, tuple)


@pytest.mark.parametrize("field", ["poloidal_coefficients", "swirl_coefficients"])
def test_tensor_coefficients_are_bounded_finite_and_immutable(field):
    for bad in (np.nan, np.inf, -np.inf, -1.01, 1.01):
        values = np.zeros(27)
        values[4] = bad
        with pytest.raises(ValueError):
            TensorCandidate(**{field: values})
    with pytest.raises(ValueError):
        TensorCandidate(**{field: np.zeros(26)})
    with pytest.raises(ValueError):
        TensorCandidate(**{field: np.zeros((27, 1))})

    candidate = TensorCandidate()
    assert isinstance(getattr(candidate, field), tuple)
    with pytest.raises(FrozenInstanceError):
        setattr(candidate, field, (0.0,) * 27)

