from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from openai_ns_reconstruction.kokuno_heat_discrepancy_repair import (
    KokunoHeatDiscrepancyRepair,
)


def _direct_normalized_difference(
    repair: KokunoHeatDiscrepancyRepair,
    coefficients: np.ndarray,
    *,
    f_eta: float,
) -> np.ndarray:
    nodes, weights = leggauss(repair.quadrature_points)
    x_min, x_max = 0.80, 1.70
    half = 0.5 * (x_max - x_min)
    x = 0.5 * (x_max + x_min) + half * nodes
    weights = half * weights
    base = repair.background_over_e_star(x, f_eta)
    delta = repair.correction_over_e_star(x, coefficients)
    corrected = base + delta

    def moments(E: np.ndarray) -> np.ndarray:
        return np.asarray(
            [
                np.sum(weights * E * E / (2.0 * x)),
                -0.5 * np.sum(weights * E * E),
                np.sum(weights * np.sqrt(2.0 * x) * E),
            ]
        )

    return moments(corrected) - moments(base)


def test_exact_quadratic_rows_match_direct_full_moment_difference() -> None:
    repair = KokunoHeatDiscrepancyRepair(lambda_outer=0.07, quadrature_points=192)
    coefficients = np.asarray([0.003, -0.002, 0.0015])
    f_eta = float(repair.source_f(0.23))
    exact = repair.normalized_increment(coefficients, f_eta=f_eta)
    direct = _direct_normalized_difference(repair, coefficients, f_eta=f_eta)
    np.testing.assert_allclose(exact, direct, rtol=3e-11, atol=2e-13)


def test_source_linear_weights_and_analytic_coefficient_jacobian() -> None:
    repair = KokunoHeatDiscrepancyRepair(lambda_outer=0.08, quadrature_points=192)
    f_eta = 0.87
    x = np.asarray([0.93, 1.21, 1.49])
    weights = repair.linear_weights(x, f_eta=f_eta)
    expected = np.stack(
        [
            f_eta * x ** (-1.5 - repair.lambda_outer),
            -f_eta * x ** (-0.5 - repair.lambda_outer),
            np.sqrt(2.0) * np.sqrt(x),
        ]
    )
    np.testing.assert_array_equal(weights, expected)

    coefficients = np.asarray([0.002, -0.0015, 0.001])
    jac = repair.coefficient_jacobian(coefficients, f_eta=f_eta)
    eps = 2.0e-7
    fd = np.empty((3, 3))
    for column in range(3):
        plus = coefficients.copy()
        minus = coefficients.copy()
        plus[column] += eps
        minus[column] -= eps
        fd[:, column] = (
            repair.normalized_increment(plus, f_eta=f_eta)
            - repair.normalized_increment(minus, f_eta=f_eta)
        ) / (2.0 * eps)
    np.testing.assert_allclose(jac, fd, rtol=2e-7, atol=2e-10)

    linear = repair.coefficient_jacobian(np.zeros(3), f_eta=f_eta)
    assert np.linalg.matrix_rank(linear) == 3
    assert np.linalg.cond(linear) < 2.0e4


def test_bounded_nonlinear_inverse_closes_manufactured_small_discrepancy() -> None:
    repair = KokunoHeatDiscrepancyRepair(
        lambda_outer=0.05, quadrature_points=192, coefficient_limit=0.05
    )
    eta = 0.2
    f_eta = float(repair.source_f(eta))
    hidden_only_for_regression = np.asarray([0.003, -0.002, 0.0015])
    target = repair.normalized_increment(hidden_only_for_regression, f_eta=f_eta)
    result = repair.solve(target, f_eta=f_eta)
    assert result.success
    assert result.nfev <= 12
    assert result.max_abs_residual < 2.0e-11
    np.testing.assert_allclose(result.coefficients, hidden_only_for_regression, rtol=2e-8, atol=2e-9)
    np.testing.assert_allclose(result.achieved_normalized, target, rtol=0.0, atol=2e-11)


def test_compact_profile_derivative_and_physical_scaling() -> None:
    repair = KokunoHeatDiscrepancyRepair()
    coefficients = np.asarray([0.004, -0.003, 0.002])
    outside = np.asarray([0.70, 1.80])
    np.testing.assert_array_equal(repair.correction_over_e_star(outside, coefficients), 0.0)
    np.testing.assert_array_equal(
        repair.correction_derivative_over_e_star(outside, coefficients), 0.0
    )

    x = 1.24
    step = 1.0e-6
    fd = (
        repair.correction_over_e_star(x + step, coefficients)
        - repair.correction_over_e_star(x - step, coefficients)
    ) / (2.0 * step)
    analytic = repair.correction_derivative_over_e_star(x, coefficients)
    assert float(analytic) == pytest.approx(float(fd), rel=2e-6, abs=2e-8)

    normalized = np.asarray([1.0, -2.0, 3.0]) * 1.0e-4
    physical = repair.physical_increment(normalized, X_star=3.5, e_star=1.2)
    expected = normalized * np.asarray([1.2**2, 3.5 * 1.2**2, 3.5**1.5 * 1.2])
    np.testing.assert_array_equal(physical, expected)

    physical_profile = repair.physical_profile_correction(
        np.asarray([3.5 * 1.24]), X_star=3.5, e_star=1.2, coefficients=coefficients
    )
    assert abs(float(physical_profile["delta_E"][0])) > 0.0
    expected_dx = (
        1.2
        * float(repair.correction_derivative_over_e_star(1.24, coefficients))
        / 3.5
    )
    assert float(physical_profile["delta_E_X"][0]) == pytest.approx(expected_dx)


def test_serialization_truth_boundary_and_local_size_guard_fail_closed(tmp_path) -> None:
    repair = KokunoHeatDiscrepancyRepair(
        lambda_outer=0.04, quadrature_points=128, coefficient_limit=0.04
    )
    path = repair.save_json(tmp_path / "repair.json")
    replay = KokunoHeatDiscrepancyRepair.load_json(path)
    assert replay.sha256 == repair.sha256
    assert replay.to_payload() == repair.to_payload()

    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["truth_boundary"]["core_to_heat_matching_completed"] = True
    unsigned = {key: value for key, value in payload.items() if key != "sha256"}
    payload["sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    with pytest.raises(ValueError):
        KokunoHeatDiscrepancyRepair.from_payload(payload)

    with pytest.raises(ValueError):
        repair.normalized_increment(np.asarray([0.05, 0.0, 0.0]), f_eta=1.0)
    with pytest.raises(ValueError):
        repair.solve(np.asarray([0.5, 0.5, 0.5]), f_eta=1.0)
    with pytest.raises(ValueError):
        repair.source_f(1.1)
