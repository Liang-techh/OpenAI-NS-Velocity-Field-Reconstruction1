import numpy as np

from openai_ns_reconstruction.kokuno_agent4_radial_stress_spacetime_generalization import (
    CARTESIAN_FD8_STEP,
    PARENT_AGENT3_HEAD,
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


def test_protocol_is_frozen_and_distinct_from_prior_single_state_audit():
    assert PARENT_AGENT3_HEAD == "a97f0882ab0b1fcc0d4b9575b5366fb44b35c233"
    assert RADIAL_COUNTS == (33, 65, 129)
    assert CARTESIAN_FD8_STEP == 0.0025
    assert len(SPACETIME_STATES) == 3
    assert SPACETIME_STATES[0] == {"t": 0.37, "z": -0.31}
    assert SPACETIME_STATES[1] == {"t": 0.50, "z": 0.08}
    assert SPACETIME_STATES[2] == {"t": 0.63, "z": 0.31}
