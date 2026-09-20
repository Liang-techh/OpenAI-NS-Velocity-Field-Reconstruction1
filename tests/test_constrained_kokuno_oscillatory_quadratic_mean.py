from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_oscillatory_quadratic_mean as meanmod


def test_analytic_cylindrical_mean_projection_recovers_m0_mode() -> None:
    receipt = meanmod.analytic_regression_receipt()
    assert receipt["schema"] == meanmod.SCHEMA
    assert receipt["task"] == meanmod.TASK
    assert receipt["real_agent2_numeric_mean_observed_in_this_receipt"] is False
    assert receipt["analytic_projection_absolute_max_error"] < 2.0e-14

    witness = receipt["analytic_witness"]
    assert witness["angular_orders"] == [32, 64, 128]
    assert max(witness["successive_mean_relative_differences"]) < 2.0e-14
    mean = np.asarray(witness["mean_cylindrical"], dtype=float)
    expected = np.broadcast_to(np.asarray([2.0, -0.5, 0.25]), mean.shape)
    np.testing.assert_allclose(mean, expected, rtol=0.0, atol=2.0e-14)


def test_projection_is_rotating_cylindrical_mean_not_naive_cartesian_mean() -> None:
    radius = np.asarray([0.8])
    z = np.asarray([0.1])
    t = np.asarray([0.5])
    mean, _ = meanmod._mean_for_order(
        meanmod._analytic_regression_evaluator,
        radius,
        z,
        t,
        128,
    )
    np.testing.assert_allclose(mean[0], [2.0, -0.5, 0.25], atol=2.0e-14, rtol=0.0)

    theta = 2.0 * np.pi * np.arange(128) / 128.0
    x = radius[..., None] * np.cos(theta)
    y = radius[..., None] * np.sin(theta)
    diagnostic = meanmod._analytic_regression_evaluator(
        x,
        y,
        z[..., None] + np.zeros_like(theta),
        t[..., None] + np.zeros_like(theta),
        spatial_step=meanmod.AGENT2_SPATIAL_STEP,
    )
    naive = np.mean(np.asarray(diagnostic["self_advection"]), axis=-2)
    # A radial/tangential axisymmetric vector rotates in Cartesian coordinates;
    # averaging Cartesian x/y would therefore erase the physically relevant m=0 profile.
    np.testing.assert_allclose(naive[0, :2], [0.0, 0.0], atol=2.0e-14, rtol=0.0)
    assert abs(float(naive[0, 2]) - 0.25) < 2.0e-14


def test_exact_backend_binding_rejects_same_name_from_wrong_source_blob() -> None:
    def fake(x, y, z, t, *, spatial_step=meanmod.AGENT2_SPATIAL_STEP):
        del x, y, z, t, spatial_step
        return {"self_advection": np.zeros((3,), dtype=float)}

    fake.__module__ = meanmod.AGENT2_MODULE
    fake.__name__ = meanmod.AGENT2_FUNCTION
    with pytest.raises(ValueError, match="source blob"):
        meanmod.ExactAgent2SelfAdvectionBackend.bind(fake)


def test_ring_domain_and_payload_fail_closed() -> None:
    with pytest.raises(ValueError, match="0 < r < 2"):
        meanmod._materialize_with_evaluator(
            meanmod._analytic_regression_evaluator,
            radius=0.0,
            z=0.0,
            t=0.5,
            backend=None,
            backend_kind="test",
        )
    with pytest.raises(ValueError, match="time"):
        meanmod._materialize_with_evaluator(
            meanmod._analytic_regression_evaluator,
            radius=0.5,
            z=0.0,
            t=0.8,
            backend=None,
            backend_kind="test",
        )

    def missing_payload(x, y, z, t, *, spatial_step):
        del x, y, z, t, spatial_step
        return {}

    with pytest.raises(TypeError, match="self_advection"):
        meanmod._materialize_with_evaluator(
            missing_payload,
            radius=0.5,
            z=0.0,
            t=0.5,
            backend=None,
            backend_kind="test",
        )


def test_public_api_has_no_scientific_escape_hatch_inputs() -> None:
    forbidden = {
        "residual",
        "defect",
        "stress",
        "debt",
        "pressure",
        "forcing",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "quadrature_order",
        "angular_order",
        "spatial_step",
        "delta_y",
        "delta_a",
        "alpha",
    }
    for function in (
        meanmod.materialize_agent2_oscillatory_quadratic_mean,
        meanmod.ExactAgent2SelfAdvectionBackend.bind,
    ):
        parameters = set(inspect.signature(function).parameters)
        assert not (parameters & forbidden)


def test_truth_boundary_keeps_quadratic_subterm_separate_from_ns_defect() -> None:
    boundary = meanmod.truth_boundary()
    assert boundary["agent2_self_advection_reimplemented_by_agent3"] is False
    assert boundary["cylindrical_m0_projection_executable"] is True
    assert boundary["projection_uses_naive_cartesian_vector_mean"] is False
    assert boundary["quadratic_mean_is_complete_ns_defect"] is False
    assert boundary["leading_cross_terms_included"] is False
    assert boundary["pressure_gradient_included"] is False
    assert boundary["viscous_term_included"] is False
    assert boundary["restricted_forcing_included"] is False
    assert boundary["caller_supplied_mean_values_allowed"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5
