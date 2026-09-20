from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_inner_leading_oscillatory_cross_mean import (
    AGENT1_HEAD,
    AGENT1_PR,
    AGENT2_HEAD,
    AGENT2_INNER_LEADING_FD4_STEP,
    AGENT2_OSCILLATORY_FD6_STEP,
    AGENT2_PR,
    ANGULAR_ORDERS,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    SCHEMA,
    TASK,
    ExactAgent2InnerLeadingOscillatoryCrossBackend,
    _analytic_regression_evaluator,
    _materialize_with_evaluator,
    analytic_regression_receipt,
    materialize_agent2_inner_leading_oscillatory_cross_mean,
    truth_boundary,
)


class DummyInner:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.stack((1.0 + 0.0 * x, 0.0 * y, 0.0 * z + 0.0 * t), axis=-1)


def test_analytic_projection_recovers_two_term_and_total_means() -> None:
    receipt = analytic_regression_receipt()
    assert receipt["schema"] == SCHEMA
    assert receipt["task"] == TASK
    assert receipt["projection_ladder"] == list(ANGULAR_ORDERS)
    assert receipt["agent2_pr"] == AGENT2_PR
    assert receipt["agent2_expected_head"] == AGENT2_HEAD
    assert receipt["agent1_inner_velocity_pr"] == AGENT1_PR
    assert receipt["agent1_inner_velocity_expected_head"] == AGENT1_HEAD
    assert receipt["analytic_projection_absolute_max_error"] < 2.0e-14
    witness = receipt["analytic_witness"]
    assert witness["termwise_mean_closure_absolute_max"] < 2.0e-14
    assert receipt["real_agent2_numeric_inner_leading_cross_mean_observed_in_this_receipt"] is False


def test_projected_terms_close_to_projected_cross() -> None:
    witness = _materialize_with_evaluator(
        _analytic_regression_evaluator,
        DummyInner(),
        radius=np.asarray([0.51, 0.74]),
        z=np.asarray([-0.2, 0.3]),
        t=np.asarray([0.47, 0.53]),
        backend=None,
        backend_kind="test",
    )
    term_sum = (
        witness.mean_inner_advects_oscillation_cylindrical
        + witness.mean_oscillation_advects_inner_cylindrical
    )
    assert np.max(np.abs(term_sum - witness.mean_cross_cylindrical)) < 2.0e-14
    assert witness.cross_mean_rms <= witness.full_ring_cross_rms + 1.0e-14


def test_payload_fails_closed_on_step_drift() -> None:
    def bad_lead_step(inner, x, y, z, t, **kwargs):
        out = _analytic_regression_evaluator(inner, x, y, z, t, **kwargs)
        return SimpleNamespace(**{**out.__dict__, "inner_leading_spatial_step": 2.0e-3})

    with pytest.raises(ValueError, match="inner-leading FD4 step"):
        _materialize_with_evaluator(
            bad_lead_step,
            DummyInner(),
            [0.6],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )

    def bad_osc_step(inner, x, y, z, t, **kwargs):
        out = _analytic_regression_evaluator(inner, x, y, z, t, **kwargs)
        return SimpleNamespace(**{**out.__dict__, "oscillatory_spatial_step": 2.0e-3})

    with pytest.raises(ValueError, match="oscillatory FD6 step"):
        _materialize_with_evaluator(
            bad_osc_step,
            DummyInner(),
            [0.6],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )


def test_payload_fails_closed_on_shape_nonfinite_and_termwise_drift() -> None:
    def bad_shape(inner, x, y, z, t, **kwargs):
        out = _analytic_regression_evaluator(inner, x, y, z, t, **kwargs)
        return SimpleNamespace(
            **{**out.__dict__, "cross_advection": out.cross_advection[..., :2]}
        )

    with pytest.raises(ValueError, match="unexpected ring shape"):
        _materialize_with_evaluator(
            bad_shape,
            DummyInner(),
            [0.6],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )

    def nonfinite(inner, x, y, z, t, **kwargs):
        out = _analytic_regression_evaluator(inner, x, y, z, t, **kwargs)
        cross = np.array(out.cross_advection, copy=True)
        cross.reshape(-1)[0] = np.nan
        return SimpleNamespace(**{**out.__dict__, "cross_advection": cross})

    with pytest.raises(ValueError, match="non-finite"):
        _materialize_with_evaluator(
            nonfinite,
            DummyInner(),
            [0.6],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )

    def bad_closure(inner, x, y, z, t, **kwargs):
        out = _analytic_regression_evaluator(inner, x, y, z, t, **kwargs)
        cross = np.array(out.cross_advection, copy=True)
        cross[..., 0] += 1.0e-4
        return SimpleNamespace(**{**out.__dict__, "cross_advection": cross})

    with pytest.raises(ValueError, match="termwise closure"):
        _materialize_with_evaluator(
            bad_closure,
            DummyInner(),
            [0.6],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )


def test_exact_backend_rejects_non_agent2_callable() -> None:
    with pytest.raises(ValueError, match="module identity"):
        ExactAgent2InnerLeadingOscillatoryCrossBackend.bind(_analytic_regression_evaluator)


def test_public_api_exposes_no_scientific_escape_hatches() -> None:
    parameters = inspect.signature(
        materialize_agent2_inner_leading_oscillatory_cross_mean
    ).parameters
    assert set(parameters) == {
        "backend",
        "inner_leading_backend",
        "radius",
        "z",
        "t",
    }
    forbidden = {
        "residual",
        "defect",
        "mean",
        "stress",
        "debt",
        "inverse",
        "pressure",
        "forcing",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "quadrature_order",
        "angular_order",
        "viscosity",
        "nu",
        "inner_leading_spatial_step",
        "oscillatory_spatial_step",
        "delta_y",
        "delta_a",
        "alpha",
        "damping_grid",
    }
    assert forbidden.isdisjoint(parameters)


def test_ring_validation_and_inner_backend_protocol_fail_closed() -> None:
    with pytest.raises(ValueError):
        _materialize_with_evaluator(
            _analytic_regression_evaluator,
            DummyInner(),
            [0.5, 0.7],
            [0.0, 0.1, 0.2],
            [0.5, 0.5],
            backend=None,
            backend_kind="test",
        )
    with pytest.raises(ValueError):
        _materialize_with_evaluator(
            _analytic_regression_evaluator,
            DummyInner(),
            [0.0],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )

    class NoVelocity:
        pass

    # Public binder type-check happens before any evaluator call.
    fake = object.__new__(ExactAgent2InnerLeadingOscillatoryCrossBackend)
    object.__setattr__(fake, "evaluator", _analytic_regression_evaluator)
    object.__setattr__(fake, "source_file", "test")
    object.__setattr__(fake, "source_blob_sha", "wrong")
    with pytest.raises(ValueError, match="provenance drifted"):
        materialize_agent2_inner_leading_oscillatory_cross_mean(
            fake, NoVelocity(), [0.6], [0.0], [0.5]
        )


def test_truth_boundary_keeps_inner_cross_mean_below_complete_defect_claim() -> None:
    boundary = truth_boundary()
    assert boundary["agent2_inner_leading_oscillatory_cross_reimplemented_by_agent3"] is False
    assert boundary["cylindrical_m0_inner_leading_oscillatory_cross_projection_executable"] is True
    assert boundary["uses_agent1_pa10_inner_contraction_center_only"] is True
    assert boundary["inner_leading_is_final_corrected_fixed_point"] is False
    assert boundary["outer_global_leading_join_materialized"] is False
    assert boundary["global_leading_velocity_materialized"] is False
    assert boundary["inner_leading_oscillatory_cross_mean_is_complete_ns_defect"] is False
    assert boundary["inner_leading_self_advection_included_here"] is False
    assert boundary["inner_leading_velocity_dt_included_here"] is False
    assert boundary["oscillatory_self_advection_mean_included_here"] is False
    assert boundary["correction_transport_mean_included_here"] is False
    assert boundary["pressure_gradient_included"] is False
    assert boundary["restricted_forcing_included"] is False
    assert boundary["real_agent3_delta_a_bound"] is False
    assert boundary["real_full_candidate_bound"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["caller_supplied_residual_allowed"] is False
    assert boundary["caller_supplied_mean_values_allowed"] is False
    assert boundary["caller_supplied_derivative_resolution_allowed"] is False
    assert boundary["caller_supplied_scientific_threshold_allowed"] is False
    assert boundary["pde_validated"] is False
    assert boundary["paper_exact"] is False
    assert boundary["blowup_proved"] is False
    assert boundary["final_normalized_momentum_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE
    assert boundary["final_normalized_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE


def test_frozen_derivative_steps_remain_exact() -> None:
    assert AGENT2_INNER_LEADING_FD4_STEP == 1.0e-3
    assert AGENT2_OSCILLATORY_FD6_STEP == 1.0e-3
