from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_family_covariance_cross_terms import (
    generate_actual_core_report,
    measure_family_covariance_jacobian,
)


def _constant_family(by_beta_values: np.ndarray, *, bad_total: bool = False):
    values = np.asarray(by_beta_values, dtype=float)
    labels = tuple(f"b{index}" for index in range(len(values)))

    def evaluate(points: np.ndarray, time: float, phase_offset: float):
        del time, phase_offset
        by_beta = np.broadcast_to(values, (len(points),) + values.shape).copy()
        total = np.sum(by_beta, axis=1)
        if bad_total:
            total = total.copy()
            total[:, 0] += 0.25
        return {
            "beta_labels": labels,
            "velocity_physical_cylindrical_by_beta": by_beta,
            "velocity_physical_cylindrical_total": total,
        }

    return evaluate


def test_two_label_cross_terms_and_euler_identity_are_retained():
    family = _constant_family(np.array([[1.0, 2.0, 3.0], [4.0, -1.0, 0.5]]))
    report = measure_family_covariance_jacobian(
        family,
        np.array([0.1, 0.2, 0.3]),
        angular_count=8,
        phase_count=4,
    )

    total = np.asarray(report["total_covariance_theta_axial"])
    diagonal = np.asarray(report["diagonal_self_covariance_theta_axial"])
    cross = np.asarray(report["cross_covariance_theta_axial"])
    jacobian = np.asarray(report["per_label_covariance_jacobian_theta_axial"])

    np.testing.assert_allclose(total, np.tile([5.0, 17.5], (3, 1)))
    np.testing.assert_allclose(diagonal, np.tile([-2.0, 5.0], (3, 1)))
    np.testing.assert_allclose(cross, np.tile([7.0, 12.5], (3, 1)))
    np.testing.assert_allclose(jacobian[:, 0, :], np.tile([11.0, 18.5], (3, 1)))
    np.testing.assert_allclose(jacobian[:, 1, :], np.tile([-1.0, 16.5], (3, 1)))
    np.testing.assert_allclose(np.sum(jacobian, axis=1), 2.0 * total)
    assert report["quadratic_euler_identity_checked"] is True
    assert report["quadratic_euler_relative_error"] < 1.0e-14
    assert report["cross_covariance_relative_vector_rms"] > 0.5
    assert report["bounded_inverse_ready"] is False


def test_custom_tangent_scales_preserve_coefficient_units():
    family = _constant_family(np.array([[1.0, 2.0, 3.0], [4.0, -1.0, 0.5]]))
    report = measure_family_covariance_jacobian(
        family,
        np.array([0.1, 0.2, 0.3]),
        coefficient_tangent_scales=[2.0, -0.5],
        angular_count=8,
        phase_count=4,
    )
    jacobian = np.asarray(report["per_label_covariance_jacobian_theta_axial"])
    np.testing.assert_allclose(jacobian[:, 0, :], np.tile([22.0, 37.0], (3, 1)))
    np.testing.assert_allclose(jacobian[:, 1, :], np.tile([0.5, -8.25], (3, 1)))
    assert report["coefficient_tangent_scales"] == [2.0, -0.5]
    assert report["quadratic_euler_identity_checked"] is False
    assert report["quadratic_euler_relative_error"] is None
    assert "coefficient units" in report["coefficient_semantics"]


def test_declared_total_must_equal_sum_of_by_beta_velocities():
    family = _constant_family(
        np.array([[1.0, 2.0, 3.0], [4.0, -1.0, 0.5]]),
        bad_total=True,
    )
    with pytest.raises(ValueError, match="does not equal the sum"):
        measure_family_covariance_jacobian(
            family,
            np.array([0.1, 0.2, 0.3]),
            angular_count=8,
            phase_count=4,
        )


def test_one_label_family_has_no_cross_covariance():
    family = _constant_family(np.array([[1.5, -2.0, 0.25]]))
    report = measure_family_covariance_jacobian(
        family,
        np.array([0.1, 0.2, 0.3]),
        angular_count=8,
        phase_count=4,
    )
    cross = np.asarray(report["cross_covariance_theta_axial"])
    np.testing.assert_allclose(cross, 0.0, atol=1.0e-15)
    total = np.asarray(report["total_covariance_theta_axial"])
    jacobian = np.asarray(report["per_label_covariance_jacobian_theta_axial"])
    np.testing.assert_allclose(jacobian[:, 0, :], 2.0 * total)


def test_real_one_label_calibration_reproduces_existing_agent3_response(tmp_path):
    report = generate_actual_core_report(
        output=tmp_path / "family.json",
        target_output=tmp_path / "target.json",
    )
    calibration = report["one_label_real_calibration"]
    assert calibration["relative_error_to_existing_current_response"] < 2.0e-12
    assert calibration["cross_covariance_vector_rms"] < 1.0e-15
    assert report["inputs"]["surrogate_defect_used"] is False
    truth = report["truth_boundary"]
    assert truth["family_cross_term_contract_executable"] is True
    assert truth["actual_source_multilabel_family_consumed"] is False
    assert truth["bounded_inverse_rerun_allowed"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["pde_validated"] is False
