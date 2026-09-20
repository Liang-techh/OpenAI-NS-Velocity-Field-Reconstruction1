from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_correction_transport_mean import (
    AGENT2_VISCOSITY,
    ANGULAR_ORDERS,
    COMPONENT_KEYS,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    SCHEMA,
    TASK,
    ExactAgent2CorrectionTransportBackend,
    _analytic_regression_evaluator,
    _materialize_with_evaluator,
    analytic_regression_receipt,
    materialize_agent2_correction_transport_mean,
    truth_boundary,
)


def test_analytic_projection_recovers_component_and_total_means() -> None:
    receipt = analytic_regression_receipt()
    assert receipt["schema"] == SCHEMA
    assert receipt["task"] == TASK
    assert receipt["projection_ladder"] == list(ANGULAR_ORDERS)
    errors = receipt["analytic_projection_absolute_max_errors"]
    assert max(errors.values()) < 2.0e-14

    witness = receipt["analytic_witness"]
    assert witness["pointwise_component_closure_absolute_max"] < 2.0e-14
    assert witness["projected_component_closure_absolute_max"] < 2.0e-14
    assert witness["projected_component_closure_relative_rms"] < 2.0e-14
    assert witness["viscosity"] == AGENT2_VISCOSITY
    assert receipt["real_agent2_numeric_correction_transport_mean_observed_in_this_receipt"] is False


def test_component_projection_closes_to_projected_total() -> None:
    witness = _materialize_with_evaluator(
        _analytic_regression_evaluator,
        correction=object(),
        radius=np.asarray([0.51, 0.74]),
        z=np.asarray([-0.2, 0.3]),
        t=np.asarray([0.47, 0.53]),
        backend=None,
        backend_kind="test",
    )
    component_sum = (
        witness.mean_correction_velocity_dt_cylindrical
        + witness.mean_cross_advection_cylindrical
        + witness.mean_correction_self_advection_cylindrical
        + witness.mean_correction_viscous_term_cylindrical
    )
    assert np.max(
        np.abs(component_sum - witness.mean_correction_transport_cylindrical)
    ) < 2.0e-14
    assert witness.correction_transport_mean_rms <= (
        witness.full_ring_correction_transport_rms + 1.0e-14
    )


def test_payload_contract_fails_closed_on_missing_or_drifted_metadata() -> None:
    def missing_key(correction, x, y, z, t):
        payload = _analytic_regression_evaluator(correction, x, y, z, t)
        del payload["correction_self_advection"]
        return payload

    with pytest.raises(TypeError, match="correction_self_advection"):
        _materialize_with_evaluator(
            missing_key,
            object(),
            [0.6],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )

    def bad_viscosity(correction, x, y, z, t):
        payload = _analytic_regression_evaluator(correction, x, y, z, t)
        payload["viscosity"] = 0.02
        return payload

    with pytest.raises(ValueError, match="viscosity"):
        _materialize_with_evaluator(
            bad_viscosity,
            object(),
            [0.6],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )

    def bad_handoff(correction, x, y, z, t):
        payload = _analytic_regression_evaluator(correction, x, y, z, t)
        payload["handoff_contract"]["complete_ns_defect"] = True
        return payload

    with pytest.raises(ValueError, match="complete_ns_defect"):
        _materialize_with_evaluator(
            bad_handoff,
            object(),
            [0.6],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )


def test_payload_contract_fails_closed_on_shape_and_nonfinite_values() -> None:
    def bad_shape(correction, x, y, z, t):
        payload = _analytic_regression_evaluator(correction, x, y, z, t)
        payload["cross_advection"] = payload["cross_advection"][..., :2]
        return payload

    with pytest.raises(ValueError, match="unexpected ring shape"):
        _materialize_with_evaluator(
            bad_shape,
            object(),
            [0.6],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )

    def nonfinite(correction, x, y, z, t):
        payload = _analytic_regression_evaluator(correction, x, y, z, t)
        value = np.array(payload["correction_transport_increment"], copy=True)
        value.reshape(-1)[0] = np.nan
        payload["correction_transport_increment"] = value
        return payload

    with pytest.raises(ValueError, match="non-finite"):
        _materialize_with_evaluator(
            nonfinite,
            object(),
            [0.6],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )


def test_exact_backend_rejects_non_agent2_callable() -> None:
    with pytest.raises(ValueError, match="module identity"):
        ExactAgent2CorrectionTransportBackend.bind(_analytic_regression_evaluator)


def test_public_api_exposes_no_scientific_escape_hatches() -> None:
    parameters = inspect.signature(materialize_agent2_correction_transport_mean).parameters
    assert set(parameters) == {"backend", "correction", "radius", "z", "t"}
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
        "oscillatory_spatial_step",
        "correction_spatial_step",
        "correction_laplacian_step",
        "delta_y",
        "delta_a",
        "alpha",
        "damping_grid",
    }
    assert forbidden.isdisjoint(parameters)


def test_ring_validation_is_inherited_and_fail_closed() -> None:
    with pytest.raises(ValueError):
        _materialize_with_evaluator(
            _analytic_regression_evaluator,
            object(),
            [0.5, 0.7],
            [0.0, 0.1, 0.2],
            [0.5, 0.5],
            backend=None,
            backend_kind="test",
        )
    with pytest.raises(ValueError):
        _materialize_with_evaluator(
            _analytic_regression_evaluator,
            object(),
            [0.0],
            [0.0],
            [0.5],
            backend=None,
            backend_kind="test",
        )


def test_truth_boundary_keeps_transport_mean_below_complete_defect_claim() -> None:
    boundary = truth_boundary()
    assert boundary["agent2_correction_transport_reimplemented_by_agent3"] is False
    assert boundary["cylindrical_m0_correction_transport_projection_executable"] is True
    assert boundary["componentwise_projected_closure_checked"] is True
    assert boundary["correction_transport_mean_is_complete_ns_defect"] is False
    assert boundary["includes_correction_time_derivative"] is True
    assert boundary["includes_oscillation_correction_cross_advection"] is True
    assert boundary["includes_correction_self_advection"] is True
    assert boundary["includes_correction_viscosity"] is True
    assert boundary["leading_cross_terms_included"] is False
    assert boundary["pressure_gradient_included"] is False
    assert boundary["restricted_forcing_included"] is False
    assert boundary["real_agent3_delta_a_bound"] is False
    assert boundary["real_full_candidate_bound"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["caller_supplied_residual_allowed"] is False
    assert boundary["caller_supplied_mean_values_allowed"] is False
    assert boundary["caller_supplied_viscosity_allowed"] is False
    assert boundary["caller_supplied_derivative_resolution_allowed"] is False
    assert boundary["caller_supplied_scientific_threshold_allowed"] is False
    assert boundary["pde_validated"] is False
    assert boundary["paper_exact"] is False
    assert boundary["blowup_proved"] is False
    assert boundary["final_normalized_momentum_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE
    assert boundary["final_normalized_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE


def test_component_keys_are_frozen() -> None:
    assert COMPONENT_KEYS == (
        "correction_velocity_dt",
        "cross_advection",
        "correction_self_advection",
        "correction_viscous_term",
    )
