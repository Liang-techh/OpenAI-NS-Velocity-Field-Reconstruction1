import numpy as np

from openai_ns_reconstruction.kokuno_public_candidate_velocity import velocity_osc
from openai_ns_reconstruction.kokuno_public_oscillatory_independent_audit import (
    _divergence_metrics,
    _fd4_jacobian,
)


def test_independent_fd4_operator_on_analytic_divergence_free_rotation():
    points = np.asarray(
        [
            [0.41, -0.27, 0.13],
            [-0.62, 0.18, -0.31],
            [0.22, 0.53, 0.44],
        ],
        dtype=float,
    )
    times = np.asarray((0.31, 0.47, 0.69), dtype=float)

    def rotation(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((-y, x, np.zeros_like(z)), axis=-1)

    jac = _fd4_jacobian(points, times, rotation, 0.01)
    metrics = _divergence_metrics(jac)
    assert metrics["absolute_max"] < 1.0e-12
    assert metrics["relative_rms"] < 1.0e-12


def test_public_provider_black_box_smoke_and_radial_fail_safe():
    inside = np.asarray(
        velocity_osc(
            np.asarray((0.55, 0.71)),
            np.asarray((0.21, -0.28)),
            np.asarray((0.17, -0.33)),
            np.asarray((0.37, 0.61)),
        ),
        dtype=float,
    )
    assert inside.shape == (2, 3)
    assert np.all(np.isfinite(inside))
    assert np.linalg.norm(inside) > 0.0

    outside = np.asarray(
        velocity_osc(
            np.asarray((0.0, 0.05, 1.50)),
            np.asarray((0.0, 0.02, 0.0)),
            np.asarray((0.0, 0.2, -0.4)),
            np.asarray((0.31, 0.47, 0.63)),
        ),
        dtype=float,
    )
    assert np.array_equal(outside, np.zeros_like(outside))
