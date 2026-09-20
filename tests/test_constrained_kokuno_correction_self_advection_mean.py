from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_correction_self_advection_mean as meanmod


def test_analytic_projection_recovers_cylindrical_mean() -> None:
    receipt = meanmod.analytic_regression_receipt()
    assert receipt["schema"] == meanmod.SCHEMA
    assert receipt["task"] == meanmod.TASK
    assert receipt["real_agent2_numeric_correction_self_mean_observed_in_this_receipt"] is False
    assert receipt["analytic_projection_absolute_max_error"] < 2.0e-14

    witness = receipt["analytic_witness"]
    assert witness["angular_orders"] == [32, 64, 128]
    assert max(witness["successive_mean_relative_differences"]) < 2.0e-14
    expected = np.broadcast_to(np.asarray([0.12, -0.04, 0.18]), (3, 3))
    np.testing.assert_allclose(
        np.asarray(witness["mean_correction_self_advection_cylindrical"], dtype=float),
        expected,
        rtol=0.0,
        atol=2.0e-14,
    )


def test_exact_backend_binding_rejects_same_name_from_wrong_source_blob() -> None:
    def fake(correction, x, y, z, t, **kwargs):
        del correction, x, y, z, t, kwargs
        return {}

    fake.__module__ = meanmod.AGENT2_MODULE
    fake.__name__ = meanmod.AGENT2_FUNCTION
    with pytest.raises(ValueError, match="source blob"):
        meanmod.ExactAgent2CorrectionSelfAdvectionBackend.bind(fake)


def test_payload_shape_and_finiteness_fail_closed() -> None:
    def missing(correction, x, y, z, t, **kwargs):
        del correction, x, y, z, t, kwargs
        return {}

    with pytest.raises(TypeError, match="correction_self_advection"):
        meanmod._materialize_with_evaluator(
            missing,
            object(),
            radius=0.5,
            z=0.0,
            t=0.5,
            backend=None,
            backend_kind="test",
        )

    def wrong_shape(correction, x, y, z, t, **kwargs):
        del correction, y, z, t, kwargs
        x = np.asarray(x, dtype=float)
        return {"correction_self_advection": np.zeros(x.shape + (2,), dtype=float)}

    with pytest.raises(ValueError, match="unexpected ring shape"):
        meanmod._materialize_with_evaluator(
            wrong_shape,
            object(),
            radius=0.5,
            z=0.0,
            t=0.5,
            backend=None,
            backend_kind="test",
        )

    def nonfinite(correction, x, y, z, t, **kwargs):
        del correction, y, z, t, kwargs
        x = np.asarray(x, dtype=float)
        value = np.zeros(x.shape + (3,), dtype=float)
        value[..., 0] = np.nan
        return {"correction_self_advection": value}

    with pytest.raises(ValueError, match="non-finite"):
        meanmod._materialize_with_evaluator(
            nonfinite,
            object(),
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
        "inverse",
        "pressure",
        "forcing",
        "mean",
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
        "damping_grid",
    }
    for function in (
        meanmod.materialize_agent2_correction_self_advection_mean,
        meanmod.ExactAgent2CorrectionSelfAdvectionBackend.bind,
    ):
        parameters = set(inspect.signature(function).parameters)
        assert not (parameters & forbidden)


def test_ring_domain_is_inherited_and_fail_closed() -> None:
    with pytest.raises(ValueError, match="0 < r < 2"):
        meanmod._materialize_with_evaluator(
            meanmod._analytic_regression_evaluator,
            object(),
            radius=0.0,
            z=0.0,
            t=0.5,
            backend=None,
            backend_kind="test",
        )


def test_truth_boundary_keeps_self_mean_separate_from_full_defect() -> None:
    boundary = meanmod.truth_boundary()
    assert boundary["agent2_correction_self_advection_reimplemented_by_agent3"] is False
    assert boundary["cylindrical_m0_correction_self_projection_executable"] is True
    assert boundary["projection_uses_naive_cartesian_vector_mean"] is False
    assert boundary["correction_self_mean_is_complete_ns_defect"] is False
    assert boundary["oscillatory_self_advection_mean_included_here"] is False
    assert boundary["oscillation_correction_cross_mean_included_here"] is False
    assert boundary["correction_nonlinear_mean_aggregate_complete"] is False
    assert boundary["leading_cross_terms_included"] is False
    assert boundary["pressure_gradient_included"] is False
    assert boundary["viscous_term_included"] is False
    assert boundary["restricted_forcing_included"] is False
    assert boundary["real_agent3_delta_a_bound"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5
