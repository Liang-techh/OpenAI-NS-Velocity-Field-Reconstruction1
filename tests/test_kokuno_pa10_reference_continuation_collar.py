from __future__ import annotations

import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_reference_continuation_collar import (
    QUADRATURE_ORDER,
    SELECTED_T1,
    SOURCE_ANALYTIC_Y_UPPER_TIMES_LAMBDA,
    SOURCE_X0_TIMES_LAMBDA,
    SOURCE_X_B_REF,
    SOURCE_X_I,
    KokunoPA10ReferenceContinuationCollar,
    smooth_step,
)


def _independent_transition_integrals(
    collar: KokunoPA10ReferenceContinuationCollar,
    X: np.ndarray,
    eta: np.ndarray,
    *,
    order: int = 160,
) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    y = np.log(X / collar.X_0)
    half = 0.5 * (y - SELECTED_T1)
    mid = 0.5 * (y + SELECTED_T1)
    y_nodes = mid[:, None] + half[:, None] * nodes[None, :]
    y_weights = half[:, None] * weights[None, :]
    X_nodes = collar.X_0 * np.exp(y_nodes)
    eta_nodes = np.broadcast_to(eta[:, None], X_nodes.shape)
    natural = collar.source_normalized.values(X_nodes, eta_nodes)
    natural_d = collar.source_normalized.derivatives(X_nodes, eta_nodes)
    gate = 1.0 - smooth_step((y_nodes - SELECTED_T1) / SELECTED_T1)
    return (
        np.sum(
            y_weights
            * gate
            * X_nodes
            * natural_d["F_0_X"]
            / natural["F_0"],
            axis=1,
        ),
        np.sum(y_weights * gate * X_nodes * natural_d["U_0_X"], axis=1),
    )


def test_public_constants_exact_smooth_step_and_geometric_bound() -> None:
    collar = KokunoPA10ReferenceContinuationCollar()
    assert SOURCE_X0_TIMES_LAMBDA == 4.0
    assert SOURCE_X_B_REF == 100.0
    assert SOURCE_X_I == 110.0
    assert SELECTED_T1 == 0.005
    assert QUADRATURE_ORDER == 64
    assert collar.geometric_source_bound_satisfied is True
    assert 0.0 < SELECTED_T1 < collar.geometric_t1_upper_bound
    assert 4.0 * math.exp(2.0 * SELECTED_T1) < SOURCE_ANALYTIC_Y_UPPER_TIMES_LAMBDA
    assert collar.X_0 == pytest.approx(4.0 / collar.Lambda, rel=0.0, abs=0.0)

    s = np.linspace(-0.2, 1.2, 1401)
    sigma = smooth_step(s)
    assert np.all(np.diff(sigma) >= 0.0)
    np.testing.assert_array_equal(sigma[s <= 0.0], np.zeros(np.sum(s <= 0.0)))
    np.testing.assert_array_equal(sigma[s >= 1.0], np.ones(np.sum(s >= 1.0)))
    interior = np.linspace(0.01, 0.99, 201)
    np.testing.assert_allclose(
        smooth_step(interior) + smooth_step(1.0 - interior),
        np.ones_like(interior),
        rtol=0.0,
        atol=3e-15,
    )


def test_reference_matches_natural_profile_through_y_t1() -> None:
    collar = KokunoPA10ReferenceContinuationCollar()
    X = collar.X_0 * np.exp(np.array([-0.2, 0.0, 0.5 * SELECTED_T1, SELECTED_T1]))
    eta = np.array([-0.6, -0.2, 0.25, 0.7])
    actual = collar.values(X, eta)
    natural = collar.source_normalized.values(X, eta)
    np.testing.assert_allclose(actual["F_reference"], natural["F_0"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(actual["U_reference"], natural["U_0"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        actual["E_reference"], natural["E_0"], rtol=2e-15, atol=0.0
    )

    actual_d = collar.radial_derivatives(X, eta)
    natural_d = collar.source_normalized.derivatives(X, eta)
    np.testing.assert_allclose(
        actual_d["F_reference_X"], natural_d["F_0_X"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        actual_d["U_reference_X"], natural_d["U_0_X"], rtol=0.0, atol=0.0
    )


def test_transition_ode_matches_independent_high_order_quadrature() -> None:
    collar = KokunoPA10ReferenceContinuationCollar()
    fractions = np.array([0.2, 0.45, 0.7, 0.95])
    y = SELECTED_T1 * (1.0 + fractions)
    X = collar.X_0 * np.exp(y)
    eta = np.array([-0.7, -0.2, 0.3, 0.75])

    actual = collar.values(X, eta)
    start = collar.source_normalized.values(np.full_like(X, collar.X_1), eta)
    log_F_increment, U_increment = _independent_transition_integrals(
        collar, X, eta
    )
    expected_F = start["F_0"] * np.exp(log_F_increment)
    expected_U = start["U_0"] + U_increment
    np.testing.assert_allclose(
        actual["F_reference"], expected_F, rtol=3e-12, atol=1e-14
    )
    np.testing.assert_allclose(
        actual["U_reference"], expected_U, rtol=3e-12, atol=1e-14
    )
    assert float(np.max(np.abs(actual["F_reference"] - start["F_0"]))) > 0.0
    assert float(np.max(np.abs(actual["U_reference"] - start["U_0"]))) > 0.0


def test_radial_derivative_replays_centered_difference_inside_transition() -> None:
    collar = KokunoPA10ReferenceContinuationCollar()
    y = SELECTED_T1 * np.array([1.25, 1.55, 1.85])
    X = collar.X_0 * np.exp(y)
    eta = np.array([-0.55, 0.1, 0.65])
    analytic = collar.radial_derivatives(X, eta)

    relative_h = 2.0e-6
    X_plus = X * (1.0 + relative_h)
    X_minus = X * (1.0 - relative_h)
    plus = collar.values(X_plus, eta)
    minus = collar.values(X_minus, eta)
    denominator = X_plus - X_minus
    fd_F = (plus["F_reference"] - minus["F_reference"]) / denominator
    fd_U = (plus["U_reference"] - minus["U_reference"]) / denominator

    scale_F = np.maximum(np.abs(fd_F), 1.0)
    scale_U = np.maximum(np.abs(fd_U), 1.0)
    assert (
        float(np.max(np.abs(analytic["F_reference_X"] - fd_F) / scale_F))
        < 2e-6
    )
    assert (
        float(np.max(np.abs(analytic["U_reference_X"] - fd_U) / scale_U))
        < 2e-6
    )


def test_endpoint_handoff_is_value_continuous_and_zero_radial_jet() -> None:
    collar = KokunoPA10ReferenceContinuationCollar()
    eta = np.linspace(-0.9, 0.9, 13)

    left = collar.values(np.full_like(eta, collar.X_1), eta)
    natural_left = collar.source_normalized.values(
        np.full_like(eta, collar.X_1), eta
    )
    np.testing.assert_allclose(left["F_reference"], natural_left["F_0"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(left["U_reference"], natural_left["U_0"], rtol=0.0, atol=0.0)

    handoff = collar.handoff_at_collar_end(eta)
    assert np.all(np.isfinite(handoff["F_reference"]))
    assert np.all(np.isfinite(handoff["U_reference"]))
    np.testing.assert_array_equal(
        handoff["F_reference_X"], np.zeros_like(handoff["F_reference_X"])
    )
    np.testing.assert_array_equal(
        handoff["U_reference_X"], np.zeros_like(handoff["U_reference_X"])
    )


def test_serialization_hash_and_truth_boundary(tmp_path) -> None:
    collar = KokunoPA10ReferenceContinuationCollar()
    path = tmp_path / "reference_collar.json"
    collar.save_configuration(path)
    loaded = KokunoPA10ReferenceContinuationCollar.load_configuration(path)
    assert loaded.configuration() == collar.configuration()
    assert loaded.semantic_sha256 == collar.semantic_sha256

    truth = collar.truth_boundary
    assert truth["public_first_reference_collar_formula_executable"] is True
    assert truth["exact_public_smooth_step_executable"] is True
    assert truth["selected_t1_is_repository_autonomous"] is True
    assert truth["selected_t1_chosen_from_ns_residual"] is False
    assert truth["selected_t1_full_later_smallness_admitted"] is False
    assert truth["stress_activation_materialized"] is False
    assert truth["kappa0_continuation_to_X100_materialized"] is False
    assert truth["final_interpolation_to_Xi_materialized"] is False
    assert truth["actual_G_i_at_Xi_materialized"] is False
    assert truth["actual_upstream_five_moment_discrepancy_materialized"] is False
    assert truth["actual_inner_to_outer_bridge_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_fail_closed_beyond_first_collar_and_on_mutated_configuration() -> None:
    collar = KokunoPA10ReferenceContinuationCollar()
    with pytest.raises(ValueError):
        collar.values(collar.X_2 * 1.000001, 0.0)
    with pytest.raises(ValueError):
        collar.values(0.5 * collar.X_0, collar.eta_interval[1] + 1.0e-3)
    payload = collar.configuration()
    payload["selected_t1"] = 0.006
    with pytest.raises(ValueError):
        KokunoPA10ReferenceContinuationCollar.from_configuration(payload)
