import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_radial_aware_damping_holdout import (
    damped_axial_covariance_stress,
    evaluate_fixed_damping_radial_cell,
    evaluate_radial_aware_holdout_transfer,
)


def _measurement(linear_axial, quadratic_axial):
    linear = np.column_stack((np.zeros(len(linear_axial)), np.asarray(linear_axial, dtype=float)))
    quadratic = np.column_stack(
        (np.zeros(len(quadratic_axial)), np.asarray(quadratic_axial, dtype=float))
    )
    return {
        "linear_covariance_change_theta_axial": linear,
        "self_covariance_theta_axial": quadratic,
        "exact_covariance_change_theta_axial": linear + quadratic,
        "quadratic_identity_passed": True,
    }


def _radial_result(relative, *, improves=True, stable=True):
    return {
        "radial_relative_residual_rms_after_source_sign": relative,
        "radial_source_sign_improves": improves,
        "radial_source_sign_not_worse": relative <= 1.0,
        "radial_force_derivative_audit": {
            "radial_force_derivative_stability_preflight_passed": stable,
        },
    }


def _tangential_validation():
    return {
        "validation_transfer_preflight_passed": True,
        "heldout_damping_reoptimized": False,
        "validation_objective_used_for_selection": False,
        "cells": [
            {"time": 0.3125, "z": 0.07},
            {"time": 0.6875, "z": 0.09},
        ],
    }


def test_damped_axial_stress_uses_linear_plus_lambda_squared_self():
    measurement = _measurement(
        np.array([2.0, -1.0, 0.5]),
        np.array([4.0, 2.0, -2.0]),
    )
    actual = damped_axial_covariance_stress(measurement, 0.5)
    expected = 0.5 * np.array([2.0, -1.0, 0.5]) + 0.25 * np.array([4.0, 2.0, -2.0])
    assert np.allclose(actual, expected)


def test_radial_cell_uses_source_plus_divergence_sign_and_can_cancel_defect():
    n = 9
    defect = np.ones(n)
    steps = (0.02, 0.01, 0.005)
    derivatives = tuple(-np.ones(n) for _ in steps)
    result = evaluate_fixed_damping_radial_cell(defect, derivatives, z_steps=steps)

    assert result["source_radial_divergence_sign"] == "+d_z_sigma_1"
    assert result["gated_radial_defect_rms_before"] == pytest.approx(1.0)
    assert result["gated_radial_defect_rms_after_source_sign"] < 1e-14
    assert result["radial_source_sign_improves"] is True
    assert result["radial_source_sign_not_worse"] is True
    assert result["radial_force_derivative_audit"][
        "radial_force_derivative_stability_preflight_passed"
    ] is True


def test_radial_cell_rejects_source_sign_worsening_without_retuning():
    n = 9
    defect = np.ones(n)
    steps = (0.02, 0.01, 0.005)
    derivatives = tuple(0.25 * np.ones(n) for _ in steps)
    result = evaluate_fixed_damping_radial_cell(defect, derivatives, z_steps=steps)

    assert result["radial_relative_residual_rms_after_source_sign"] == pytest.approx(1.25)
    assert result["radial_source_sign_improves"] is False
    assert result["radial_source_sign_not_worse"] is False


def test_radial_aware_transfer_requires_same_disjoint_cells_and_nonworsening():
    radial = [
        {"time": 0.3125, "z": 0.07, "radial_result": _radial_result(0.98)},
        {"time": 0.6875, "z": 0.09, "radial_result": _radial_result(1.0, improves=False)},
    ]
    result = evaluate_radial_aware_holdout_transfer(_tangential_validation(), radial)
    assert result["tangential_validation_transfer_passed"] is True
    assert result["radial_source_sign_not_worse_every_validation_cell"] is True
    assert result["radial_aware_damping_transfer_preflight_passed"] is True
    assert result["heldout_damping_reoptimized"] is False


def test_radial_aware_transfer_fails_when_any_radial_cell_worsens():
    radial = [
        {"time": 0.3125, "z": 0.07, "radial_result": _radial_result(0.98)},
        {"time": 0.6875, "z": 0.09, "radial_result": _radial_result(1.001, improves=False)},
    ]
    result = evaluate_radial_aware_holdout_transfer(_tangential_validation(), radial)
    assert result["radial_source_sign_not_worse_every_validation_cell"] is False
    assert result["radial_aware_damping_transfer_preflight_passed"] is False


def test_radial_aware_transfer_rejects_cell_mismatch_and_unstable_derivative():
    mismatch = [
        {"time": 0.3125, "z": 0.07, "radial_result": _radial_result(0.98)},
        {"time": 0.6875, "z": 0.08, "radial_result": _radial_result(0.99)},
    ]
    with pytest.raises(ValueError, match="match exactly"):
        evaluate_radial_aware_holdout_transfer(_tangential_validation(), mismatch)

    unstable = [
        {
            "time": 0.3125,
            "z": 0.07,
            "radial_result": _radial_result(0.98, stable=False),
        },
        {"time": 0.6875, "z": 0.09, "radial_result": _radial_result(0.99)},
    ]
    with pytest.raises(ValueError, match="derivative stability"):
        evaluate_radial_aware_holdout_transfer(_tangential_validation(), unstable)
