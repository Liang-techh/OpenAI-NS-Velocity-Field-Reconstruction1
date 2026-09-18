import json
import math

import numpy as np
from numpy.polynomial.legendre import leggauss

from openai_ns_reconstruction.kokuno_shear_moment_repair import (
    KokunoShearMomentRepair,
)


def _f(eta: float) -> float:
    return 1.0 / (1.0 + eta * eta)


def test_pa17_zero_jacobian_matches_source_power_blocks():
    order = 192
    repair = KokunoShearMomentRepair(quadrature_points=order)
    eta = 0.23
    f = _f(eta)
    jac = repair.coefficient_jacobian(np.zeros(5), f_eta=f)

    nodes, weights = leggauss(order)
    lo, hi = 0.80, 1.70
    xi = 0.5 * (hi - lo) * nodes + 0.5 * (hi + lo)
    w = 0.5 * (hi - lo) * weights
    basis = repair.basis_values(xi)
    lam = repair.lambda_outer
    u_weights = np.stack(
        [np.ones_like(xi), np.sqrt(2.0) * f * xi ** (-lam)], axis=0
    )
    e_weights = np.stack(
        [
            np.sqrt(2.0 * xi),
            -f * xi ** (-0.5 - lam),
            f * xi ** (-1.5 - lam),
        ],
        axis=0,
    )
    expected = np.zeros((5, 5))
    for column in range(2):
        expected[0, column] = np.sum(w * u_weights[0] * basis[column])
        expected[1, column] = np.sum(w * u_weights[1] * basis[column])
    for column in range(2, 5):
        expected[2, column] = np.sum(w * e_weights[0] * basis[column])
        expected[3, column] = np.sum(w * e_weights[1] * basis[column])
        expected[4, column] = np.sum(w * e_weights[2] * basis[column])

    assert np.allclose(jac, expected, rtol=3.0e-13, atol=3.0e-14)
    singular = np.linalg.svd(jac, compute_uv=False)
    assert singular[-1] > 0.0
    assert singular[0] / singular[-1] < 1.0e8


def test_exact_nonlinear_map_jacobian_and_manufactured_inverse():
    repair = KokunoShearMomentRepair(quadrature_points=192)
    eta = 0.2
    f = _f(eta)
    coefficients = np.array([0.003, -0.002, 0.002, -0.0015, 0.001])
    jac = repair.coefficient_jacobian(coefficients, f_eta=f)

    step = 2.0e-7
    fd = np.empty_like(jac)
    for column in range(5):
        plus = coefficients.copy()
        minus = coefficients.copy()
        plus[column] += step
        minus[column] -= step
        fd[:, column] = (
            repair.normalized_increment(plus, f_eta=f)
            - repair.normalized_increment(minus, f_eta=f)
        ) / (2.0 * step)
    assert np.allclose(jac, fd, rtol=2.0e-8, atol=2.0e-10)

    target = repair.normalized_increment(coefficients, f_eta=f)
    solved = repair.solve(target, f_eta=f)
    assert solved.success
    assert solved.max_abs_residual < 2.0e-12
    assert np.allclose(solved.coefficients, coefficients, rtol=0.0, atol=2.0e-9)


def test_i1_binding_profile_and_velocity_correction_are_executable():
    repair = KokunoShearMomentRepair(quadrature_points=128)
    report = repair.log_scale_report()
    low, high = report["I1_log_interval"]
    support_low = report["log_X_1"] + math.log(report["xi_support"][0])
    support_high = report["log_X_1"] + math.log(report["xi_support"][1])
    assert low < support_low < support_high < high
    assert math.isfinite(report["X_1"]) and report["X_1"] > 0.0
    assert math.isfinite(report["e_1"]) and report["e_1"] > 0.0
    for value in report["physical_moment_scale_logs"].values():
        assert math.isfinite(value)

    coefficients = np.array([0.003, -0.002, 0.002, -0.0015, 0.001])
    coefficient_eta = np.array([0.001, -0.0005, -0.0007, 0.0003, 0.0002])
    xi = 1.18
    eta = 0.2
    values = repair.profile_correction_logX(
        report["log_X_1"] + math.log(xi),
        eta,
        coefficients=coefficients,
        coefficient_eta=coefficient_eta,
    )
    for key in (
        "delta_E",
        "delta_E_X",
        "delta_E_eta",
        "delta_F",
        "delta_U",
        "delta_U_X",
        "delta_U_eta",
        "delta_M",
        "delta_M_eta",
        "delta_v0",
    ):
        assert np.all(np.isfinite(values[key]))
    assert abs(float(values["delta_E"])) > 0.0

    X = repair.X_1 * xi
    r = math.sqrt(2.0 * X)
    velocity = repair.velocity_correction(
        r,
        0.0,
        0.0,
        0.0,
        coefficients=coefficients,
        coefficient_eta=np.zeros(5),
    )
    assert velocity.shape == (3,)
    assert np.all(np.isfinite(velocity))
    assert np.linalg.norm(velocity) > 0.0


def test_physical_scaling_and_fail_closed_serialization(tmp_path):
    repair = KokunoShearMomentRepair(quadrature_points=96)
    normalized = np.array([1.0e-5, -2.0e-5, 3.0e-5, -4.0e-5, 5.0e-5])
    physical = repair.physical_increment(normalized)
    assert physical.shape == (5,)
    assert np.all(np.isfinite(physical))
    assert np.all(np.sign(physical) == np.sign(normalized))

    path = repair.save_json(tmp_path / "pa17.json")
    loaded = KokunoShearMomentRepair.load_json(path)
    assert loaded.sha256 == repair.sha256
    assert loaded.to_payload() == repair.to_payload()

    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["truth_boundary"]["I1_cone_repair_applied_to_actual_modulation"] = True
    path.write_text(json.dumps(payload), encoding="utf-8")
    try:
        KokunoShearMomentRepair.load_json(path)
    except ValueError as exc:
        assert "hash or content mismatch" in str(exc)
    else:
        raise AssertionError("tampered truth metadata must fail closed")


def test_outside_i1_refuses_and_eta_derivative_is_explicit():
    repair = KokunoShearMomentRepair()
    low, _ = repair.outer_schedule.reserved_log_intervals()["I1"]
    coeff = np.zeros(5)
    try:
        repair.profile_correction_logX(
            low - 0.1,
            0.0,
            coefficients=coeff,
            coefficient_eta=coeff,
        )
    except ValueError as exc:
        assert "open I1" in str(exc)
    else:
        raise AssertionError("evaluation outside source I1 must fail closed")
