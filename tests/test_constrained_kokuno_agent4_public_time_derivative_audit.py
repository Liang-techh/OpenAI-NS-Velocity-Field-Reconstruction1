import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_agent4_public_time_derivative_audit import (
    FD4_STEPS,
    GUARDS,
    _fd4_time_derivative,
    run_audit,
)


def test_fd4_time_derivative_calibrates_on_quartic_vector_field():
    x = np.asarray((0.2, -0.4, 0.7), dtype=float)
    y = np.asarray((0.1, 0.3, -0.2), dtype=float)
    z = np.asarray((-0.5, 0.6, 0.2), dtype=float)
    t = np.asarray((0.33, 0.47, 0.64), dtype=float)

    def manufactured(xv, yv, zv, tv):
        xv = np.asarray(xv, dtype=float)
        yv = np.asarray(yv, dtype=float)
        zv = np.asarray(zv, dtype=float)
        tv = np.asarray(tv, dtype=float)
        return np.stack(
            (
                tv**4 + 2.0 * xv,
                -3.0 * tv**3 + yv,
                5.0 * tv**2 - zv,
            ),
            axis=-1,
        )

    expected = np.stack((4.0 * t**3, -9.0 * t**2, 10.0 * t), axis=-1)
    actual = _fd4_time_derivative(manufactured, x, y, z, t, FD4_STEPS[0])
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2.0e-12)


def test_fd4_time_derivative_rejects_invalid_step():
    one = np.ones(2)
    with pytest.raises(ValueError):
        _fd4_time_derivative(lambda x, y, z, t: np.stack((t, t, t), axis=-1), one, one, one, one, 0.0)


def test_public_time_derivative_audit_is_fail_closed_and_truthful():
    receipt = run_audit()
    assert receipt["audited_agent2_pr"] == 579
    assert receipt["audited_agent2_head"] == "6c8e71c800a17c0e9df1802042f9feebf9dbce04"
    assert receipt["seed"] == 9173261
    assert receipt["fd4_steps"] == list(FD4_STEPS)
    assert receipt["guards_frozen_before_actions"] == GUARDS
    assert receipt["local_point_count"] == 36
    assert receipt["integral_point_count"] == 12
    assert len(receipt["fd4_comparison"]) == 3
    assert len(receipt["fundamental_theorem"]["quadrature_comparison"]) == 3
    assert receipt["support_exterior_absolute_max"] <= 1.0e-12
    assert receipt["public_derivative_rms"] > 1.0e-8
    assert receipt["fundamental_theorem"]["sign_flip_mutation_relative_rms"] >= 0.5
    assert receipt["public_time_derivative_preflight_passed"] == (len(receipt["failed_guards"]) == 0)
    assert receipt["independence"]["public_values_only"] is True
    assert receipt["independence"]["agent2_fd6_verifier_reused"] is False
    assert receipt["truth_boundary"]["heldout_ns_momentum_residual_assessed"] is False
    assert receipt["truth_boundary"]["formal_full_domain_pde_gate_assessed"] is False
    assert receipt["truth_boundary"]["pde_validated"] is False
