from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_oscillatory_correction_cross_mean as meanmod


def test_analytic_cross_projection_recovers_termwise_cylindrical_means() -> None:
    receipt = meanmod.analytic_regression_receipt()
    assert receipt["schema"] == meanmod.SCHEMA
    assert receipt["task"] == meanmod.TASK
    assert receipt["real_agent2_numeric_cross_mean_observed_in_this_receipt"] is False
    assert receipt["analytic_projection_absolute_max_error"] < 2.0e-14

    witness = receipt["analytic_witness"]
    assert witness["angular_orders"] == [32, 64, 128]
    assert max(witness["successive_cross_mean_relative_differences"]) < 2.0e-14
    assert witness["termwise_mean_closure_absolute_max"] < 2.0e-14
    expected = np.broadcast_to(np.asarray([0.25, -0.05, 0.30]), (3, 3))
    np.testing.assert_allclose(
        np.asarray(witness["mean_cross_cylindrical"], dtype=float),
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
        meanmod.ExactAgent2CrossAdvectionBackend.bind(fake)


def test_payload_and_termwise_closure_fail_closed() -> None:
    def missing(correction, x, y, z, t, **kwargs):
        del correction, x, y, z, t, kwargs
        return {}

    with pytest.raises(TypeError, match="both mixed terms"):
        meanmod._materialize_with_evaluator(
            missing,
            object(),
            radius=0.5,
            z=0.0,
            t=0.5,
            backend=None,
            backend_kind="test",
        )

    def inconsistent(correction, x, y, z, t, **kwargs):
        del correction, z, t, kwargs
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        zeros = np.zeros(x.shape + (3,), dtype=float)
        bad = np.ones(y.shape + (3,), dtype=float)
        return {
            "oscillation_advects_correction": zeros,
            "correction_advects_oscillation": zeros,
            "cross_advection": bad,
        }

    with pytest.raises(ValueError, match="termwise closure"):
        meanmod._materialize_with_evaluator(
            inconsistent,
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
        meanmod.materialize_agent2_oscillatory_correction_cross_mean,
        meanmod.ExactAgent2CrossAdvectionBackend.bind,
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


def test_truth_boundary_keeps_cross_mean_separate_from_full_defect() -> None:
    boundary = meanmod.truth_boundary()
    assert boundary["agent2_cross_advection_reimplemented_by_agent3"] is False
    assert boundary["cylindrical_m0_cross_projection_executable"] is True
    assert boundary["projection_uses_naive_cartesian_vector_mean"] is False
    assert boundary["cross_mean_is_complete_ns_defect"] is False
    assert boundary["oscillatory_self_advection_mean_included_here"] is False
    assert boundary["correction_self_advection_included"] is False
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
