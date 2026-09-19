import numpy as np

from openai_ns_reconstruction.kokuno_agent4_radial_stress_spacetime_generalization import (
    ANGULAR_COUNT,
    CARTESIAN_FD8_STEP,
    PARENT_AGENT3_HEAD,
    PRIOR_ANCHOR_STATE,
    RADIAL_COUNTS,
    SPACETIME_STATES,
    _fd8_axis_samples,
)


def _polynomial_velocity(x, y, z, t):
    x = np.asarray(x, dtype=float)
    return np.stack((x**8, x**7, x**6), axis=-1)


def test_fd8_cartesian_operator_on_polynomial():
    x = np.array([0.17, 0.31, 0.53], dtype=float)
    zero = np.zeros_like(x)
    first, second = _fd8_axis_samples(
        _polynomial_velocity,
        x,
        zero,
        zero,
        zero,
        axis=0,
        step=0.01,
    )
    expected_first = np.stack((8.0 * x**7, 7.0 * x**6, 6.0 * x**5), axis=-1)
    expected_second = np.stack((56.0 * x**6, 42.0 * x**5, 30.0 * x**4), axis=-1)
    np.testing.assert_allclose(first, expected_first, rtol=2e-9, atol=2e-10)
    np.testing.assert_allclose(second, expected_second, rtol=2e-7, atol=2e-8)


def test_repaired_protocol_is_frozen_and_stencil_feasible():
    assert PARENT_AGENT3_HEAD == "a97f0882ab0b1fcc0d4b9575b5366fb44b35c233"
    assert RADIAL_COUNTS == (25, 49, 97)
    assert ANGULAR_COUNT == 16
    assert CARTESIAN_FD8_STEP == 0.0025
    assert SPACETIME_STATES == (
        {"t": 0.37, "z": -0.31},
        {"t": 0.63, "z": 0.31},
    )
    assert PRIOR_ANCHOR_STATE == {"t": 0.50, "z": 0.08}
    support_width = 1.2
    finest_spacing = support_width / (RADIAL_COUNTS[-1] - 1)
    assert 4.0 * CARTESIAN_FD8_STEP < finest_spacing
