import math

import numpy as np

from openai_ns_reconstruction.kokuno_public_oscillatory_time_derivative import (
    FD6_STEPS,
    evaluate_velocity_osc_dt,
    target_log_amplitude_derivative,
    verification_receipt,
    velocity_osc_dt,
)
from openai_ns_reconstruction.kokuno_public_z_pullback_velocity import default_field


def _autonomous_amplitudes(t: np.ndarray) -> np.ndarray:
    """Independent scalar reconstruction of the frozen sign-amplitude time factors."""
    field = default_field()
    tt = np.asarray(t, dtype=float)
    duration = field.time_max - field.time_min
    angle = 2.0 * math.pi * (tt - field.time_min) / duration
    normal = field.normal_target_ratio * (
        1.0 + field.time_modulation * np.sin(angle)
    )
    cross = field.cross_target_ratio * (
        1.0 + field.time_modulation * np.cos(angle)
    )
    # Constant positive h factors cancel from d(log a)/dt, so choose h_+=2,h_-=3.
    return np.stack(
        (
            np.sqrt((normal - cross) / 4.0),
            np.sqrt((normal + cross) / 6.0),
        ),
        axis=-1,
    )


def test_target_log_amplitude_derivative_matches_independent_scalar_fd():
    times = np.asarray((0.31, 0.43, 0.50, 0.62, 0.69))
    step = 1.0e-6
    finite_difference = (
        _autonomous_amplitudes(times + step) - _autonomous_amplitudes(times - step)
    ) / (2.0 * step)
    reference = _autonomous_amplitudes(times)
    expected = finite_difference / reference
    actual = target_log_amplitude_derivative(times)
    np.testing.assert_allclose(actual, expected, rtol=2.0e-8, atol=2.0e-9)


def test_public_time_derivative_shapes_aggregation_and_support():
    x = np.asarray((0.62, -0.71, 0.0, 1.8, 0.8))
    y = np.asarray((0.17, 0.24, 0.0, 0.0, 0.0))
    z = np.asarray((0.35, -0.82, 0.0, 0.0, 2.1))
    t = np.asarray((0.33, 0.57, 0.50, 0.50, 0.50))
    out = evaluate_velocity_osc_dt(x, y, z, t)

    total = np.asarray(out["velocity_dt_cartesian_total"])
    by_beta = np.asarray(out["velocity_dt_cartesian_by_beta"])
    by_sign = np.asarray(out["velocity_dt_cartesian_by_beta_sign"])
    assert total.shape == (5, 3)
    assert by_beta.shape[:1] == (5,)
    assert by_sign.shape[:1] == (5,)
    assert by_sign.shape[-2:] == (2, 3)
    np.testing.assert_allclose(np.sum(by_sign, axis=-2), by_beta, rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(np.sum(by_beta, axis=-2), total, rtol=0.0, atol=1.0e-12)
    assert np.linalg.norm(total[0]) > 0.0
    assert np.linalg.norm(total[1]) > 0.0
    np.testing.assert_array_equal(total[2:], np.zeros((3, 3)))
    np.testing.assert_array_equal(velocity_osc_dt(0.0, 0.0, 0.0, 0.5), np.zeros(3))


def test_exact_time_derivative_matches_independent_fd6_public_velocity():
    receipt = verification_receipt()
    assert receipt["fd6_steps"] == list(FD6_STEPS)
    assert receipt["sample_count"] == 36
    assert receipt["analytic_derivative_rms"] > 0.0
    assert receipt["analytic_derivative_sampled_max"] > 0.0
    assert receipt["support_exterior_absolute_max"] == 0.0
    rows = receipt["fd6_comparison"]
    assert len(rows) == 3
    assert rows[-1]["relative_rms_error"] < 2.0e-6
    assert rows[-1]["relative_max_error"] < 5.0e-6
    assert rows[-1]["relative_rms_error"] < rows[0]["relative_rms_error"]
    truth = receipt["truth_boundary"]
    assert truth["velocity_candidate_changed"] is False
    assert truth["source_formula_changed"] is False
    assert truth["autonomous_time_law_derivative_executable"] is True
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_time_derivative_rejects_out_of_interval_time():
    for bad in (0.249, 0.751):
        try:
            velocity_osc_dt(0.7, 0.0, 0.0, bad)
        except ValueError as exc:
            assert "outside" in str(exc)
        else:
            raise AssertionError("out-of-interval time must fail closed")
