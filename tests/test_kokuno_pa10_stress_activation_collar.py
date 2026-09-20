from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_reference_continuation_collar import (
    SELECTED_T1,
    smooth_step,
)
from openai_ns_reconstruction.kokuno_pa10_stress_activation_collar import (
    QUADRATURE_ORDER,
    SELECTED_KAPPA0,
    KokunoPA10StressActivationCollar,
)


def _independent_activation_integrals(
    collar: KokunoPA10StressActivationCollar,
    X: np.ndarray,
    eta: np.ndarray,
    *,
    order: int = 160,
) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    y = np.log(X / collar.X_0)
    half = 0.5 * y
    y_nodes = half[:, None] * (1.0 + nodes[None, :])
    y_weights = half[:, None] * weights[None, :]
    X_nodes = collar.X_0 * np.exp(y_nodes)
    eta_nodes = np.broadcast_to(eta[:, None], X_nodes.shape)
    stress = collar.reference_stress_data(X_nodes, eta_nodes)
    e_a = (1.0 - SELECTED_KAPPA0) * smooth_step(y_nodes / SELECTED_T1)
    return (
        np.sum(y_weights * 0.5 * e_a * stress["p1_reference"], axis=1),
        np.sum(
            y_weights * 0.5 * e_a * X_nodes * stress["ns_reference"],
            axis=1,
        ),
    )


def test_source_activation_schedule_and_autonomous_kappa0_boundary() -> None:
    collar = KokunoPA10StressActivationCollar()
    assert SELECTED_KAPPA0 == 0.1
    assert QUADRATURE_ORDER == 64
    assert 0.0 < SELECTED_KAPPA0 < 0.5
    assert collar.X_0 < collar.X_1

    eta = np.array([-0.7, 0.0, 0.7])
    start = collar.activation_schedule(np.full_like(eta, collar.X_0), eta)
    end = collar.activation_schedule(np.full_like(eta, collar.X_1), eta)
    np.testing.assert_array_equal(start["e_a"], np.zeros_like(eta))
    np.testing.assert_array_equal(start["kappa"], np.ones_like(eta))
    np.testing.assert_array_equal(
        end["e_a"], np.full_like(eta, 1.0 - SELECTED_KAPPA0)
    )
    np.testing.assert_array_equal(
        end["kappa"], np.full_like(eta, SELECTED_KAPPA0)
    )


def test_reference_stress_data_replays_stress_free_natural_identities() -> None:
    collar = KokunoPA10StressActivationCollar()
    fractions = np.array([0.0, 0.2, 0.55, 0.9, 1.0])
    X = collar.X_0 * np.exp(SELECTED_T1 * fractions)
    eta = np.array([-0.7, -0.3, 0.0, 0.35, 0.72])
    stress = collar.reference_stress_data(X, eta)
    natural = collar.reference.source_normalized.values(X, eta)
    natural_d = collar.reference.source_normalized.derivatives(X, eta)

    expected_p1 = -2.0 * X * natural_d["Phi_0_X"] / natural["Phi_0"]
    expected_ns = -2.0 * natural_d["U_0_X"]
    expected_p2 = X * expected_ns / natural["E_0"]
    expected_v = expected_p1 + expected_p2**2 / expected_p1
    np.testing.assert_allclose(stress["p1_reference"], expected_p1, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(stress["ns_reference"], expected_ns, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(stress["p2_reference"], expected_p2, rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(stress["v_reference"], expected_v, rtol=3e-15, atol=0.0)
    assert np.all(stress["p1_reference"] > 0.0)
    assert np.all(np.isfinite(stress["v_reference"]))


def test_activation_values_match_independent_high_order_quadrature() -> None:
    collar = KokunoPA10StressActivationCollar()
    fractions = np.array([0.18, 0.42, 0.68, 0.94])
    X = collar.X_0 * np.exp(SELECTED_T1 * fractions)
    eta = np.array([-0.65, -0.2, 0.25, 0.7])
    actual = collar.values(X, eta)
    reference = collar.reference.values(X, eta)
    delta_log_F, delta_U = _independent_activation_integrals(collar, X, eta)
    expected_F = reference["F_reference"] * np.exp(delta_log_F)
    expected_U = reference["U_reference"] + delta_U
    np.testing.assert_allclose(actual["F_activation"], expected_F, rtol=4e-12, atol=1e-14)
    np.testing.assert_allclose(actual["U_activation"], expected_U, rtol=4e-12, atol=1e-14)
    np.testing.assert_allclose(
        actual["E_activation"], np.sqrt(2.0 * X) * actual["F_activation"],
        rtol=2e-15, atol=0.0
    )
    assert float(np.max(np.abs(actual["delta_log_F_from_reference"]))) > 0.0


def test_activation_starts_exactly_from_reference_and_satisfies_difference_odes() -> None:
    collar = KokunoPA10StressActivationCollar()
    eta0 = np.linspace(-0.8, 0.8, 9)
    X0 = np.full_like(eta0, collar.X_0)
    actual0 = collar.values(X0, eta0)
    reference0 = collar.reference.values(X0, eta0)
    np.testing.assert_allclose(actual0["F_activation"], reference0["F_reference"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(actual0["U_activation"], reference0["U_reference"], rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(actual0["delta_log_F_from_reference"], np.zeros_like(eta0))
    np.testing.assert_array_equal(actual0["delta_U_from_reference"], np.zeros_like(eta0))

    fractions = np.array([0.2, 0.45, 0.7, 0.9])
    X = collar.X_0 * np.exp(SELECTED_T1 * fractions)
    eta = np.array([-0.6, -0.15, 0.3, 0.68])
    actual = collar.values(X, eta)
    actual_d = collar.radial_derivatives(X, eta)
    reference = collar.reference.values(X, eta)
    reference_d = collar.reference.radial_derivatives(X, eta)
    stress = collar.reference_stress_data(X, eta)
    schedule = collar.activation_schedule(X, eta)

    lhs_log_F = X * (
        actual_d["F_activation_X"] / actual["F_activation"]
        - reference_d["F_reference_X"] / reference["F_reference"]
    )
    rhs_log_F = 0.5 * schedule["e_a"] * stress["p1_reference"]
    lhs_U = X * (actual_d["U_activation_X"] - reference_d["U_reference_X"])
    rhs_U = 0.5 * schedule["e_a"] * X * stress["ns_reference"]
    np.testing.assert_allclose(lhs_log_F, rhs_log_F, rtol=2e-12, atol=2e-13)
    np.testing.assert_allclose(lhs_U, rhs_U, rtol=2e-12, atol=2e-13)


def test_analytic_radial_derivatives_replay_centered_difference() -> None:
    collar = KokunoPA10StressActivationCollar()
    fractions = np.array([0.28, 0.52, 0.78])
    X = collar.X_0 * np.exp(SELECTED_T1 * fractions)
    eta = np.array([-0.55, 0.08, 0.62])
    analytic = collar.radial_derivatives(X, eta)
    rel_h = 1.5e-6
    Xp = X * (1.0 + rel_h)
    Xm = X * (1.0 - rel_h)
    plus = collar.values(Xp, eta)
    minus = collar.values(Xm, eta)
    denom = Xp - Xm
    fd_F = (plus["F_activation"] - minus["F_activation"]) / denom
    fd_U = (plus["U_activation"] - minus["U_activation"]) / denom
    scale_F = np.maximum(np.abs(fd_F), 1.0)
    scale_U = np.maximum(np.abs(fd_U), 1.0)
    assert float(np.max(np.abs(analytic["F_activation_X"] - fd_F) / scale_F)) < 4e-6
    assert float(np.max(np.abs(analytic["U_activation_X"] - fd_U) / scale_U)) < 4e-6


def test_actual_shear_schedule_matches_public_formula() -> None:
    collar = KokunoPA10StressActivationCollar()
    X = collar.X_0 * np.exp(SELECTED_T1 * np.array([0.2, 0.6, 1.0]))
    eta = np.array([-0.5, 0.1, 0.65])
    shear = collar.shear_schedule(X, eta)
    stress = collar.reference_stress_data(X, eta)
    actual = collar.values(X, eta)
    schedule = collar.activation_schedule(X, eta)
    expected_first = schedule["kappa"] * stress["p1_reference"]
    expected_second = (
        schedule["kappa"]
        * stress["p2_reference"]
        * stress["E_reference"]
        / actual["E_activation"]
    )
    np.testing.assert_allclose(shear["shear_first"], expected_first, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(shear["shear_second"], expected_second, rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(
        shear["t_s"], expected_second / expected_first, rtol=3e-15, atol=0.0
    )


def test_serialization_truth_boundary_and_fail_closed_domain(tmp_path) -> None:
    collar = KokunoPA10StressActivationCollar()
    path = tmp_path / "activation_collar.json"
    collar.save_configuration(path)
    loaded = KokunoPA10StressActivationCollar.load_configuration(path)
    assert loaded.configuration() == collar.configuration()
    assert loaded.semantic_sha256 == collar.semantic_sha256

    truth = collar.truth_boundary
    assert truth["public_stress_activation_formula_executable"] is True
    assert truth["stress_activation_materialized"] is True
    assert truth["selected_kappa0_is_repository_autonomous"] is True
    assert truth["selected_kappa0_chosen_from_ns_residual"] is False
    assert truth["selected_kappa0_global_source_smallness_admitted"] is False
    assert truth["stress_activation_cone_admissibility_independently_certified"] is False
    assert truth["kappa0_continuation_to_X100_materialized"] is False
    assert truth["actual_inner_to_outer_bridge_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    with pytest.raises(ValueError):
        collar.values(collar.X_0 * (1.0 - 1e-6), 0.0)
    with pytest.raises(ValueError):
        collar.values(collar.X_1 * (1.0 + 1e-6), 0.0)
    payload = collar.configuration()
    payload["selected_kappa0"] = 0.11
    with pytest.raises(ValueError):
        KokunoPA10StressActivationCollar.from_configuration(payload)
