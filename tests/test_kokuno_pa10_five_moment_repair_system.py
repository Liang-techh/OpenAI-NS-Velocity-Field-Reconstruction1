from __future__ import annotations

import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_five_moment_repair_system import (
    AUTONOMOUS_E_BUMP_INTERVALS,
    AUTONOMOUS_U_BUMP_INTERVALS,
    COEFFICIENT_TRUST_RADIUS,
    LOG_X_SUPPORT,
    KokunoPA10FiveMomentRepairSystem,
)


def _independent_rule(order: int = 512) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    y0, y1 = LOG_X_SUPPORT
    y = 0.5 * (y1 - y0) * nodes + 0.5 * (y0 + y1)
    dy_weight = 0.5 * (y1 - y0) * weights
    x = np.exp(y)
    return x, dy_weight * x


def _independent_pa14(
    repair: KokunoPA10FiveMomentRepairSystem,
    coefficients: np.ndarray,
    eta: float,
    *,
    order: int = 512,
) -> np.ndarray:
    x, w = _independent_rule(order)
    basis = repair.basis_values(x)
    u = basis["U_bumps"] @ coefficients[:2]
    e = basis["E_bumps"] @ coefficients[2:]
    f = 1.0 / (1.0 + eta * eta)
    U0 = np.full_like(x, 4.0 * eta)
    E0 = repair.P_star * f * x**0.1
    root = np.sqrt(2.0 * x)
    H0 = root * E0
    return np.asarray(
        [
            np.sum(w * u),
            np.sum(w * root * e),
            np.sum(w * (H0 * u + U0 * root * e + root * u * e)),
            np.sum(w * (2.0 * U0 * u + u * u - E0 * e - 0.5 * e * e)),
            np.sum(w * (E0 * e / x + 0.5 * e * e / x)),
        ]
    )


def test_autonomous_basis_has_source_count_order_support_and_unit_dx_integral() -> None:
    repair = KokunoPA10FiveMomentRepairSystem()
    assert len(AUTONOMOUS_U_BUMP_INTERVALS) == 2
    assert len(AUTONOMOUS_E_BUMP_INTERVALS) == 3
    for family in (AUTONOMOUS_U_BUMP_INTERVALS, AUTONOMOUS_E_BUMP_INTERVALS):
        assert all(LOG_X_SUPPORT[0] < a < b < LOG_X_SUPPORT[1] for a, b in family)
        assert all(family[j][1] < family[j + 1][0] for j in range(len(family) - 1))

    x, w = _independent_rule(768)
    basis = repair.basis_values(x)
    assert np.all(basis["U_bumps"] >= 0.0)
    assert np.all(basis["E_bumps"] >= 0.0)
    np.testing.assert_allclose(w @ basis["U_bumps"], np.ones(2), rtol=0.0, atol=1e-5)
    np.testing.assert_allclose(w @ basis["E_bumps"], np.ones(3), rtol=0.0, atol=1e-5)

    # The autonomous supports stop strictly before the public join log x=-5.
    join_basis = repair.basis_values(math.exp(-5.01))
    np.testing.assert_array_equal(join_basis["U_bumps"], np.zeros(2))
    np.testing.assert_array_equal(join_basis["E_bumps"], np.zeros(3))


def test_exact_pa14_normalized_increment_map_matches_independent_quadrature() -> None:
    repair = KokunoPA10FiveMomentRepairSystem()
    coefficients = np.array([8.0e-7, -6.0e-7, 5.0e-7, -4.0e-7, 3.0e-7])
    eta = 0.3
    actual = repair.moment_increment(coefficients, eta, order=512)
    expected = _independent_pa14(repair, coefficients, eta, order=512)
    np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=2e-15)
    assert float(np.max(np.abs(actual))) > 0.0


def test_analytic_coefficient_jacobian_against_centered_difference() -> None:
    repair = KokunoPA10FiveMomentRepairSystem()
    coefficients = np.array([5.0e-7, -4.0e-7, 3.0e-7, -2.0e-7, 1.0e-7])
    eta = -0.25
    jac = repair.coefficient_jacobian(coefficients, eta)
    h = 1.0e-9
    fd = np.zeros((5, 5))
    for j in range(5):
        direction = np.zeros(5)
        direction[j] = h
        fd[:, j] = (
            repair.moment_increment(coefficients + direction, eta)
            - repair.moment_increment(coefficients - direction, eta)
        ) / (2.0 * h)
    scale = np.maximum(np.abs(fd), 1.0)
    assert float(np.max(np.abs(jac - fd) / scale)) < 3e-8


def test_linearized_five_bump_system_is_full_rank_across_eta_probes() -> None:
    repair = KokunoPA10FiveMomentRepairSystem()
    for eta in (-0.9, -0.5, 0.0, 0.5, 0.9):
        linear = repair.linearized_system(eta)
        assert linear["rank"] == 5
        assert abs(linear["determinant"]) > 1e-5
        assert linear["condition_number_2"] < 3.0e5
        assert float(np.min(linear["singular_values"])) > 5e-4


def test_local_newton_solver_recovers_manufactured_exact_moment_delta() -> None:
    repair = KokunoPA10FiveMomentRepairSystem()
    eta = 0.3
    known = np.array([8.0e-7, -6.0e-7, 5.0e-7, -4.0e-7, 3.0e-7])
    # Generate the target with an implementation-distinct 512-node quadrature.
    target = _independent_pa14(repair, known, eta, order=512)
    solved = repair.solve_moment_delta(target, eta)
    assert solved["converged"] is True
    assert solved["iterations"] <= 5
    assert solved["residual_inf"] <= 2e-13
    np.testing.assert_allclose(solved["achieved_delta"], target, rtol=0.0, atol=2e-13)
    assert float(np.max(np.abs(solved["coefficients"]))) < COEFFICIENT_TRUST_RADIUS
    # Independent quadrature checks the solved coefficients reproduce the same target.
    replay = _independent_pa14(repair, solved["coefficients"], eta, order=768)
    scale = np.maximum(np.abs(target), 1e-10)
    assert float(np.max(np.abs(replay - target) / scale)) < 3e-5


def test_repaired_pair_is_nontrivial_inside_and_exactly_ideal_near_join() -> None:
    repair = KokunoPA10FiveMomentRepairSystem()
    coeff = np.array([8.0e-7, -6.0e-7, 5.0e-7, -4.0e-7, 3.0e-7])
    x = np.exp(np.array([-5.82, -5.52, -5.20, -5.01]))
    out = repair.repaired_values(x, 0.2, coeff)
    assert float(np.max(np.abs(out["delta_U"]))) > 0.0
    assert float(np.max(np.abs(out["delta_E"]))) > 0.0
    assert out["delta_U"][-1] == 0.0
    assert out["delta_E"][-1] == 0.0
    assert out["U_repaired"][-1] == out["U_ideal"][-1]
    assert out["E_repaired"][-1] == out["E_ideal"][-1]


def test_serialization_hash_and_truth_boundary(tmp_path) -> None:
    repair = KokunoPA10FiveMomentRepairSystem()
    path = tmp_path / "five_moment_repair.json"
    repair.save_configuration(path)
    loaded = KokunoPA10FiveMomentRepairSystem.load_configuration(path)
    assert loaded.configuration() == repair.configuration()
    assert loaded.semantic_sha256 == repair.semantic_sha256

    truth = repair.truth_boundary
    assert truth["public_pa14_exact_moment_increment_map_executable"] is True
    assert truth["public_two_u_three_e_support_class_respected"] is True
    assert truth["repository_autonomous_bump_shapes_materialized"] is True
    assert truth["local_small_discrepancy_solver_executable"] is True
    assert truth["source_bump_shapes_numerically_identified"] is False
    assert truth["source_repair_coefficients_identified"] is False
    assert truth["actual_upstream_five_moment_discrepancy_materialized"] is False
    assert truth["eta_dependent_repair_coefficient_functions_materialized"] is False
    assert truth["actual_inner_to_outer_bridge_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_fail_closed_outside_local_scope() -> None:
    repair = KokunoPA10FiveMomentRepairSystem()
    zeros = np.zeros(5)
    with pytest.raises(ValueError):
        repair.repaired_values(math.exp(-6.1), 0.0, zeros)
    with pytest.raises(ValueError):
        repair.repaired_values(math.exp(-4.9), 0.0, zeros)
    with pytest.raises(ValueError):
        repair.moment_increment(zeros, 1.01)
    with pytest.raises(ValueError):
        repair.solve_moment_delta(np.ones(5), 0.0)
    with pytest.raises(ValueError):
        repair.moment_increment(np.full(5, 2.0 * COEFFICIENT_TRUST_RADIUS), 0.0)
