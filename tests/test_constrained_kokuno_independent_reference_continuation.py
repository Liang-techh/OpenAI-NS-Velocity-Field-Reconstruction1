from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_independent_reference_continuation import (
    evaluate_fd2,
    sample_strata,
)
from openai_ns_reconstruction.kokuno_reference_continuation import (
    KokunoReferenceContinuationCandidate,
)


def test_fd2_recovers_manufactured_affine_operator() -> None:
    a, b, c, k = 0.02, -0.01, 0.03, 0.04

    def velocity(x, y, z, t):
        x, y, z = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
        )
        return np.stack((a * x, b * y, c * z), axis=-1)

    def pressure(x, y, z, t):
        return k * np.asarray(x, dtype=float)

    points = np.array(
        [
            [0.23, -0.17, 0.11],
            [-0.31, 0.19, -0.07],
            [0.14, 0.27, 0.16],
        ]
    )
    result = evaluate_fd2(velocity, pressure, points, 0.5, 0.003, nu=0.01)
    expected = np.column_stack(
        (
            a * a * points[:, 0] + k,
            b * b * points[:, 1],
            c * c * points[:, 2],
        )
    )
    np.testing.assert_allclose(result["residual"], expected, rtol=2e-11, atol=2e-11)
    np.testing.assert_allclose(
        result["divergence"], a + b + c, rtol=0.0, atol=2e-13
    )


def test_predeclared_strata_hit_transition_post_core_and_axis_near() -> None:
    candidate = KokunoReferenceContinuationCandidate(quadrature_points=16)
    samples = sample_strata(candidate, time=0.5, count=5, seed=9_172_991)

    inner = samples["inner_core"]
    transition = samples["reference_transition"]
    post = samples["post_core_reference"]
    axis = samples["axis_near"]

    assert np.max(inner["X"]) < candidate.X1
    assert np.min(transition["X"]) > candidate.X1
    assert np.max(transition["X"]) < candidate.X2
    assert np.min(post["X"]) > candidate.max_source_transition_X

    axis_radius = np.hypot(axis["points"][:, 0], axis["points"][:, 1])
    assert np.all(axis_radius >= 0.002)
    assert np.all(axis_radius <= 0.008)
    assert np.all(np.abs(axis["points"][:, 2]) >= 0.12)
    assert np.all(np.abs(axis["points"][:, 2]) <= 0.22)


def test_fd2_rejects_invalid_resolution() -> None:
    candidate = KokunoReferenceContinuationCandidate(quadrature_points=16)
    points = np.array([[0.2, 0.1, 0.0]])
    with pytest.raises(ValueError):
        evaluate_fd2(candidate.velocity, candidate.pressure, points, 0.5, 0.0)
