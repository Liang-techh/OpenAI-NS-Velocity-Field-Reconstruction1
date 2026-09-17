from __future__ import annotations

import json
import math

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from openai_ns_reconstruction.kokuno_five_moment_repair import (
    X_SCALED_MAX,
    X_SCALED_MIN,
    KokunoFiveMomentRepair,
)


def _direct_moments(
    repair: KokunoFiveMomentRepair,
    coefficients: np.ndarray,
    *,
    eta: float,
    P_star: float,
    f_eta: float,
    order: int = 256,
) -> np.ndarray:
    nodes, weights = leggauss(order)
    half = 0.5 * (X_SCALED_MAX - X_SCALED_MIN)
    center = 0.5 * (X_SCALED_MAX + X_SCALED_MIN)
    x = center + half * nodes
    weights = half * weights
    corrected = repair.corrected_scaled_profiles(
        x, eta, P_star, f_eta, coefficients
    )
    base = repair.corrected_scaled_profiles(
        x, eta, P_star, f_eta, np.zeros(5)
    )

    def moments(values: dict[str, np.ndarray]) -> np.ndarray:
        U = values["U"]
        E = values["E"]
        H = np.sqrt(2.0 * x) * E
        return np.asarray(
            [
                np.sum(weights * U),
                np.sum(weights * H),
                np.sum(weights * U * H),
                np.sum(weights * (U * U - 0.5 * E * E)),
                np.sum(weights * E * E / (2.0 * x)),
            ]
        )

    return moments(corrected) - moments(base)


def test_pa14_exact_increment_matches_direct_five_moment_difference() -> None:
    repair = KokunoFiveMomentRepair(quadrature_points=192)
    coefficients = np.array([0.003, -0.002, 0.002, -0.0015, 0.001])
    kwargs = dict(eta=0.21, P_star=1.2, f_eta=0.9)
    exact = repair.moment_increment(coefficients, **kwargs)
    direct = _direct_moments(repair, coefficients, **kwargs)
    np.testing.assert_allclose(exact, direct, rtol=2e-10, atol=2e-12)


def test_analytic_moment_jacobian_and_pa16_block_structure() -> None:
    repair = KokunoFiveMomentRepair(quadrature_points=192)
    eta = 0.23
    kwargs = dict(eta=eta, P_star=1.15, f_eta=0.85)
    coefficients = np.array([0.002, -0.001, 0.0015, -0.0012, 0.0008])
    jac = repair.moment_jacobian(coefficients, **kwargs)

    eps = 2.0e-7
    fd = np.empty((5, 5))
    for column in range(5):
        plus = coefficients.copy()
        minus = coefficients.copy()
        plus[column] += eps
        minus[column] -= eps
        fd[:, column] = (
            repair.moment_increment(plus, **kwargs)
            - repair.moment_increment(minus, **kwargs)
        ) / (2.0 * eps)
    np.testing.assert_allclose(jac, fd, rtol=2e-7, atol=2e-10)

    linear = repair.transformed_jacobian(np.zeros(5), **kwargs)
    assert np.linalg.matrix_rank(linear) == 5
    assert np.linalg.cond(linear) < 2.0e5
    # PA.16 is exactly block diagonal after J->J-4eta I and S->S-8eta M.
    np.testing.assert_allclose(linear[:2, 2:], 0.0, atol=3e-18)
    np.testing.assert_allclose(linear[2:, :2], 0.0, atol=3e-18)


def test_bounded_nonlinear_solve_closes_a_manufactured_small_discrepancy() -> None:
    repair = KokunoFiveMomentRepair(quadrature_points=192, coefficient_limit=0.05)
    kwargs = dict(eta=0.2, P_star=1.2, f_eta=0.9)
    hidden_only_for_regression = np.array([0.003, -0.002, 0.002, -0.0015, 0.001])
    target = repair.moment_increment(hidden_only_for_regression, **kwargs)
    result = repair.solve(target, **kwargs)

    assert result.success
    assert result.nfev <= 20
    assert result.max_abs_residual < 2.0e-11
    np.testing.assert_allclose(result.achieved_scaled_moments, target, rtol=0.0, atol=2e-11)
    np.testing.assert_allclose(result.coefficients, hidden_only_for_regression, rtol=2e-8, atol=2e-9)
    assert max(abs(value) for value in result.coefficients) < repair.coefficient_limit


def test_compact_profile_correction_and_radial_derivative() -> None:
    repair = KokunoFiveMomentRepair()
    coefficients = np.array([0.004, -0.003, 0.002, -0.001, 0.0015])

    outside = np.array([math.exp(-6.2), math.exp(-4.8)])
    u, e = repair.perturbation(outside, coefficients)
    ux, ex = repair.perturbation_dx(outside, coefficients)
    np.testing.assert_array_equal(u, 0.0)
    np.testing.assert_array_equal(e, 0.0)
    np.testing.assert_array_equal(ux, 0.0)
    np.testing.assert_array_equal(ex, 0.0)

    # This point lies where the first U bump and first E bump both have support.
    x = math.exp(-5.775)
    step = 2.0e-8
    plus_u, plus_e = repair.perturbation(x + step, coefficients)
    minus_u, minus_e = repair.perturbation(x - step, coefficients)
    analytic_u, analytic_e = repair.perturbation_dx(x, coefficients)
    assert float(analytic_u) == pytest.approx(
        float((plus_u - minus_u) / (2.0 * step)), rel=3e-5, abs=1e-7
    )
    assert float(analytic_e) == pytest.approx(
        float((plus_e - minus_e) / (2.0 * step)), rel=3e-5, abs=1e-7
    )

    physical = repair.physical_profile_correction(
        np.array([0.1, 0.2]), 100.0, coefficients
    )
    # x=X/X_R lies entirely below exp(-6), so the physical corrections vanish.
    np.testing.assert_array_equal(physical["delta_U"], 0.0)
    np.testing.assert_array_equal(physical["delta_E"], 0.0)


def test_pa15_physical_moment_scaling_is_exact() -> None:
    scaled = np.array([1.0, -2.0, 3.0, -4.0, 5.0]) * 1e-4
    X_R = 3.5
    physical = KokunoFiveMomentRepair.physical_moment_increment(scaled, X_R)
    expected = scaled * np.array(
        [X_R, X_R ** 1.5, X_R ** 1.5, X_R, 1.0]
    )
    np.testing.assert_array_equal(physical, expected)


def test_serialization_truth_boundary_and_bounds_fail_closed(tmp_path) -> None:
    repair = KokunoFiveMomentRepair(quadrature_points=128, coefficient_limit=0.04)
    path = repair.save_json(tmp_path / "repair.json")
    replay = KokunoFiveMomentRepair.load_json(path)
    assert replay.sha256 == repair.sha256
    assert replay.to_payload() == repair.to_payload()

    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["truth_boundary"]["global_leading_profile_reconstructed"] = True
    unsigned = {key: value for key, value in payload.items() if key != "sha256"}
    import hashlib

    payload["sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    with pytest.raises(ValueError):
        KokunoFiveMomentRepair.from_payload(payload)

    with pytest.raises(ValueError):
        repair.moment_increment(
            np.array([0.05, 0.0, 0.0, 0.0, 0.0]), eta=0.0, P_star=1.0, f_eta=1.0
        )
    with pytest.raises(ValueError):
        repair.corrected_scaled_profiles(
            np.array([math.exp(-5.5)]), 0.0, 0.01, 0.01,
            np.array([0.0, 0.0, -0.04, -0.04, -0.04]),
        )
