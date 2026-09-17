import numpy as np

from openai_ns_reconstruction.constrained_candidate import CompactCandidate
from openai_ns_reconstruction.constrained_derivatives import reference_velocity_derivatives


def test_derivative_oracle_consumes_active_compact_candidate_interface():
    candidate = CompactCandidate()
    points = np.array(
        [
            [0.18, -0.11, 0.07],
            [-0.21, 0.14, -0.09],
        ]
    )
    bundle = reference_velocity_derivatives(
        candidate,
        points,
        np.array([0.45, 0.55]),
        spatial_step=5.0e-5,
        time_step=5.0e-5,
    )

    assert bundle.value.shape == (2, 3)
    assert bundle.time.shape == (2, 3)
    assert bundle.gradient.shape == (2, 3, 3)
    assert bundle.hessian.shape == (2, 3, 3, 3)
    assert bundle.laplacian.shape == (2, 3)
    assert np.all(np.isfinite(bundle.value))
    assert np.all(np.isfinite(bundle.time))
    assert np.all(np.isfinite(bundle.gradient))
    assert np.all(np.isfinite(bundle.hessian))
    assert np.all(np.isfinite(bundle.laplacian))
    assert np.max(np.linalg.norm(bundle.value, axis=-1)) > 1.0e-8
    assert np.max(np.abs(bundle.divergence)) < 1.0e-5
