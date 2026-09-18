from __future__ import annotations

import numpy as np

from openai_ns_reconstruction.kokuno_bounded_family_coefficients import (
    KokunoBoundedFamilyCoefficientCoordinates,
)
from openai_ns_reconstruction.kokuno_quadratic_covariance_gain_guard import (
    evaluate_quadratic_covariance_gain,
    measure_bounded_coordinate_quadratic_change,
)


def _synthetic_family(points: np.ndarray, time: float, phase: float):
    points = np.asarray(points, dtype=float)
    del time
    radius = np.hypot(points[:, 0], points[:, 1])
    theta = np.arctan2(points[:, 1], points[:, 0]) + float(phase)
    w0 = np.stack(
        (
            0.16 + 0.02 * radius + 0.01 * np.cos(theta),
            0.09 + 0.03 * np.cos(theta),
            0.05 + 0.02 * np.sin(theta),
        ),
        axis=-1,
    )
    w1 = np.stack(
        (
            -0.05 + 0.01 * radius + 0.02 * np.sin(theta),
            0.04 - 0.02 * np.sin(theta),
            -0.03 + 0.01 * np.cos(theta),
        ),
        axis=-1,
    )
    by_beta = np.stack((w0, w1), axis=-2)
    return {
        "beta_labels": ((5, (0, 0, 0)), (6, (1, 0, 0))),
        "velocity_physical_cylindrical_by_beta": by_beta,
        "velocity_physical_cylindrical_total": np.sum(by_beta, axis=-2),
    }


def test_exact_quadratic_decomposition_closes_for_common_coordinate():
    radii = np.linspace(0.10, 0.30, 9)
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=0.2)
    delta = 0.1
    measured = measure_bounded_coordinate_quadratic_change(
        _synthetic_family,
        radii,
        coordinate_contract=contract,
        delta_common=delta,
        delta_band=0.0,
        angular_count=8,
        phase_count=4,
    )
    base = np.asarray(measured["base_covariance_theta_axial"], dtype=float)
    linear = np.asarray(measured["linear_covariance_change_theta_axial"], dtype=float)
    self_term = np.asarray(measured["self_covariance_theta_axial"], dtype=float)
    exact = np.asarray(measured["exact_covariance_change_theta_axial"], dtype=float)

    np.testing.assert_allclose(linear, 2.0 * delta * base, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(self_term, delta * delta * base, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(exact, (2.0 * delta + delta * delta) * base, rtol=0.0, atol=3e-15)
    assert measured["quadratic_identity_passed"] is True
    assert measured["quadratic_identity_relative_rms"] <= 5e-12
    assert np.isclose(measured["self_to_linear_vector_rms_ratio"], abs(delta) / 2.0)


def test_band_coordinate_is_measured_from_physical_family_not_label_only():
    radii = np.linspace(0.10, 0.30, 9)
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=0.2)
    measured = measure_bounded_coordinate_quadratic_change(
        _synthetic_family,
        radii,
        coordinate_contract=contract,
        delta_common=0.0,
        delta_band=0.08,
        angular_count=8,
        phase_count=4,
    )
    assert measured["band_contrast_active"] is True
    assert measured["self_covariance_vector_rms"] > 0.0
    assert measured["linear_change_vector_rms"] > 0.0
    assert measured["quadratic_identity_passed"] is True


def test_quadratic_guard_rejects_linear_perfect_fit_spoiled_by_self_covariance():
    target = np.array([[1.0, 0.0], [0.5, -0.25]], dtype=float)
    self_term = 0.1 * target
    measurement = {
        "linear_covariance_change_theta_axial": target.tolist(),
        "self_covariance_theta_axial": self_term.tolist(),
        "exact_covariance_change_theta_axial": (target + self_term).tolist(),
        "quadratic_identity_passed": True,
    }
    result = evaluate_quadratic_covariance_gain(
        target,
        np.array([True, True]),
        measurement,
        upstream_bounded_inverse_passed=True,
        spacetime_preflight_passed=True,
        radial_force_retained=True,
        public_velocity_correction_materialized=True,
    )
    assert result["linear_predicted_relative_stress_residual_rms"] == 0.0
    assert np.isclose(result["exact_quadratic_relative_stress_residual_rms"], 0.1)
    assert result["exact_quadratic_stress_fit_within_existing_algebraic_tolerance"] is False
    assert result["quadratic_covariance_preflight_passed"] is False
    assert result["finite_correction_cycle_rerun_allowed"] is False


def test_exact_local_fit_still_requires_all_finite_cycle_prerequisites():
    target = np.array([[1.0, -0.2], [0.5, 0.3]], dtype=float)
    measurement = {
        "linear_covariance_change_theta_axial": target.tolist(),
        "self_covariance_theta_axial": np.zeros_like(target).tolist(),
        "exact_covariance_change_theta_axial": target.tolist(),
        "quadratic_identity_passed": True,
    }
    blocked = evaluate_quadratic_covariance_gain(
        target,
        np.array([True, True]),
        measurement,
        upstream_bounded_inverse_passed=True,
        spacetime_preflight_passed=True,
        radial_force_retained=False,
        public_velocity_correction_materialized=True,
    )
    assert blocked["quadratic_covariance_preflight_passed"] is True
    assert blocked["finite_correction_cycle_rerun_allowed"] is False

    allowed = evaluate_quadratic_covariance_gain(
        target,
        np.array([True, True]),
        measurement,
        upstream_bounded_inverse_passed=True,
        spacetime_preflight_passed=True,
        radial_force_retained=True,
        public_velocity_correction_materialized=True,
    )
    assert allowed["quadratic_covariance_preflight_passed"] is True
    assert allowed["finite_correction_cycle_rerun_allowed"] is True
