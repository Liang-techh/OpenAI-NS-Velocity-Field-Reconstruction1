from __future__ import annotations

import numpy as np

from openai_ns_reconstruction.kokuno_damped_quadratic_covariance_gain import (
    evaluate_damped_quadratic_covariance_gain,
    solve_damped_quadratic_covariance_gain,
)


def _measurement(linear: np.ndarray, quadratic: np.ndarray, *, l1: float = 0.2):
    linear = np.asarray(linear, dtype=float)
    quadratic = np.asarray(quadratic, dtype=float)
    return {
        "linear_covariance_change_theta_axial": linear.tolist(),
        "self_covariance_theta_axial": quadratic.tolist(),
        "exact_covariance_change_theta_axial": (linear + quadratic).tolist(),
        "quadratic_identity_passed": True,
        "aggregate_l1_update": float(l1),
    }


def test_exact_quartic_damping_recovers_known_half_step():
    # Along the proposed full update Delta C(lambda)=lambda+lambda^2 in the first
    # channel.  The target 0.75 is attained exactly at lambda=0.5, while the
    # undamped full step overshoots to 2.0.
    target = np.array([[0.75, 0.0], [0.75, 0.0]], dtype=float)
    linear = np.array([[1.0, 0.0], [1.0, 0.0]], dtype=float)
    quadratic = np.array([[1.0, 0.0], [1.0, 0.0]], dtype=float)
    solved = solve_damped_quadratic_covariance_gain(
        target,
        np.array([True, True]),
        _measurement(linear, quadratic, l1=0.2),
    )
    assert np.isclose(solved["selected_damping"], 0.5, rtol=0.0, atol=2e-12)
    assert solved["best_damped_relative_stress_residual_rms"] <= 5e-15
    assert solved["full_step_relative_stress_residual_rms"] > 1.0
    assert solved["damping_reduces_step"] is True
    assert np.isclose(solved["damped_aggregate_l1_update"], 0.1, rtol=0.0, atol=2e-12)


def test_damped_local_fit_can_pass_without_relaxing_existing_tolerance():
    target = np.array([[0.75, 0.0], [0.75, 0.0]], dtype=float)
    linear = np.array([[1.0, 0.0], [1.0, 0.0]], dtype=float)
    quadratic = np.array([[1.0, 0.0], [1.0, 0.0]], dtype=float)
    result = evaluate_damped_quadratic_covariance_gain(
        target,
        np.array([True, True]),
        _measurement(linear, quadratic),
        upstream_bounded_inverse_passed=True,
        spacetime_preflight_passed=True,
        radial_force_retained=True,
        public_velocity_correction_materialized=True,
    )
    assert result["best_damped_exact_stress_fit_within_existing_algebraic_tolerance"] is True
    assert result["damped_quadratic_covariance_preflight_passed"] is True
    assert result["finite_correction_cycle_rerun_allowed"] is True
    assert np.isclose(result["selected_damping"], 0.5, rtol=0.0, atol=2e-12)


def test_zero_response_selects_zero_step_and_remains_fail_closed():
    target = np.array([[1.0, -0.25], [0.5, 0.125]], dtype=float)
    zeros = np.zeros_like(target)
    result = evaluate_damped_quadratic_covariance_gain(
        target,
        np.array([True, True]),
        _measurement(zeros, zeros),
        upstream_bounded_inverse_passed=True,
        spacetime_preflight_passed=True,
        radial_force_retained=True,
        public_velocity_correction_materialized=True,
    )
    assert result["selected_damping"] == 0.0
    assert result["best_damped_relative_stress_residual_rms"] == 1.0
    assert result["positive_nontrivial_damped_step"] is False
    assert result["damped_quadratic_covariance_preflight_passed"] is False
    assert result["finite_correction_cycle_rerun_allowed"] is False


def test_exact_damped_fit_still_requires_radial_and_materialized_prerequisites():
    target = np.array([[0.75, 0.0], [0.75, 0.0]], dtype=float)
    linear = np.array([[1.0, 0.0], [1.0, 0.0]], dtype=float)
    quadratic = np.array([[1.0, 0.0], [1.0, 0.0]], dtype=float)
    blocked = evaluate_damped_quadratic_covariance_gain(
        target,
        np.array([True, True]),
        _measurement(linear, quadratic),
        upstream_bounded_inverse_passed=True,
        spacetime_preflight_passed=True,
        radial_force_retained=False,
        public_velocity_correction_materialized=True,
    )
    assert blocked["damped_quadratic_covariance_preflight_passed"] is True
    assert blocked["finite_correction_cycle_rerun_allowed"] is False

    blocked_materialization = evaluate_damped_quadratic_covariance_gain(
        target,
        np.array([True, True]),
        _measurement(linear, quadratic),
        upstream_bounded_inverse_passed=True,
        spacetime_preflight_passed=True,
        radial_force_retained=True,
        public_velocity_correction_materialized=False,
    )
    assert blocked_materialization["damped_quadratic_covariance_preflight_passed"] is True
    assert blocked_materialization["finite_correction_cycle_rerun_allowed"] is False
