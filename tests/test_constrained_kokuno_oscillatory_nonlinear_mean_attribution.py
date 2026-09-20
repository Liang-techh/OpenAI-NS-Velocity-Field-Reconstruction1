from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_nonlinear_mean_attribution import (
    ANGULAR_ORDERS,
    _analytic_provider,
    _materialize_from_provider,
    analytic_regression_receipt,
    materialize_oscillatory_nonlinear_mean_attribution,
    truth_boundary,
)


def _witness():
    return _materialize_from_provider(
        _analytic_provider,
        radius=np.asarray([0.35, 0.65, 0.95]),
        z=np.asarray([-0.2, 0.0, 0.25]),
        t=np.asarray([0.44, 0.50, 0.56]),
        backend=None,
        backend_kind="test-mechanics",
    )


def test_registered_rotating_cylindrical_means_are_recovered():
    witness = _witness()
    expected_a = np.broadcast_to(np.asarray([0.21, -0.06, 0.14]), (3, 3))
    expected_b = np.broadcast_to(np.asarray([-0.08, 0.02, 0.11]), (3, 3))
    expected_q = np.broadcast_to(np.asarray([0.17, -0.03, 0.09]), (3, 3))
    expected_mixed = expected_a + expected_b
    expected_total = expected_mixed + expected_q

    np.testing.assert_allclose(
        witness.mean_inner_advects_oscillation_cylindrical, expected_a, atol=2e-14, rtol=0
    )
    np.testing.assert_allclose(
        witness.mean_oscillation_advects_inner_cylindrical, expected_b, atol=2e-14, rtol=0
    )
    np.testing.assert_allclose(
        witness.mean_mixed_cross_cylindrical, expected_mixed, atol=2e-14, rtol=0
    )
    np.testing.assert_allclose(
        witness.mean_oscillatory_self_advection_cylindrical, expected_q, atol=2e-14, rtol=0
    )
    np.testing.assert_allclose(
        witness.mean_aggregate_nonlinear_cylindrical, expected_total, atol=2e-14, rtol=0
    )

    assert witness.angular_orders == ANGULAR_ORDERS
    assert witness.pointwise_three_piece_closure_absolute_max == 0.0
    assert witness.pointwise_mixed_closure_absolute_max == 0.0
    assert witness.projected_three_piece_closure_absolute_max < 2e-14
    assert witness.projected_mixed_closure_absolute_max < 2e-14
    assert witness.quadratic_mean_rms > 0.0
    assert witness.aggregate_mean_rms > 0.0
    assert witness.quadratic_mean_to_full_rms_ratio <= 1.0 + 1e-12
    assert witness.mixed_mean_to_full_rms_ratio <= 1.0 + 1e-12
    assert witness.aggregate_mean_to_full_rms_ratio <= 1.0 + 1e-12


def test_corrupted_aggregate_fails_closed():
    def bad_provider(x, y, z, t):
        a, b, q, mixed, aggregate = _analytic_provider(x, y, z, t)
        corrupted = np.array(aggregate, copy=True)
        corrupted[..., 0] += 1e-5
        return a, b, q, mixed, corrupted

    with pytest.raises(ValueError, match="three-piece nonlinear closure"):
        _materialize_from_provider(
            bad_provider,
            radius=np.asarray([0.4, 0.7]),
            z=np.asarray([0.0, 0.1]),
            t=np.asarray([0.48, 0.52]),
            backend=None,
            backend_kind="negative-control",
        )


def test_corrupted_mixed_term_fails_closed():
    def bad_provider(x, y, z, t):
        a, b, q, mixed, aggregate = _analytic_provider(x, y, z, t)
        corrupted = np.array(mixed, copy=True)
        corrupted[..., 2] -= 1e-5
        return a, b, q, corrupted, aggregate

    with pytest.raises(ValueError, match="mixed nonlinear closure"):
        _materialize_from_provider(
            bad_provider,
            radius=np.asarray([0.4, 0.7]),
            z=np.asarray([0.0, 0.1]),
            t=np.asarray([0.48, 0.52]),
            backend=None,
            backend_kind="negative-control",
        )


def test_malformed_or_nonfinite_payload_fails_closed():
    def bad_shape(x, y, z, t):
        a, b, q, mixed, aggregate = _analytic_provider(x, y, z, t)
        return a[..., :2], b, q, mixed, aggregate

    with pytest.raises(ValueError, match="returned shape"):
        _materialize_from_provider(
            bad_shape,
            radius=np.asarray([0.4]),
            z=np.asarray([0.0]),
            t=np.asarray([0.5]),
            backend=None,
            backend_kind="bad-shape",
        )

    def bad_finite(x, y, z, t):
        a, b, q, mixed, aggregate = _analytic_provider(x, y, z, t)
        a = np.array(a, copy=True)
        a[..., 0] = np.nan
        return a, b, q, mixed, aggregate

    with pytest.raises(ValueError, match="non-finite"):
        _materialize_from_provider(
            bad_finite,
            radius=np.asarray([0.4]),
            z=np.asarray([0.0]),
            t=np.asarray([0.5]),
            backend=None,
            backend_kind="bad-finite",
        )


def test_public_api_has_no_surrogate_or_threshold_knobs():
    signature = inspect.signature(materialize_oscillatory_nonlinear_mean_attribution)
    assert set(signature.parameters) == {"handoff_backend", "radius", "z", "t"}
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "inverse", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "angular_order",
        "spatial_step", "time_step", "viscosity", "nu", "delta_y", "delta_a",
        "normalized_score", "scientific_threshold",
    }
    assert forbidden.isdisjoint(signature.parameters)


def test_truth_boundary_stays_fail_closed():
    truth = truth_boundary()
    assert truth["current_typed_nonlinear_decomposition_consumed"] is True
    assert truth["oscillatory_quadratic_self_mean_materialized"] is True
    assert truth["agent2_curl_or_advection_reimplemented_by_agent3"] is False
    assert truth["radial_inverse_performed_in_this_increment"] is False
    assert truth["pressure_gradient_included"] is False
    assert truth["restricted_forcing_included"] is False
    assert truth["global_corrected_leading_join_materialized"] is False
    assert truth["complete_ns_defect"] is False
    assert truth["scoped_nonlinear_mean_authorized_as_correction_target"] is False
    assert truth["mean_correction_velocity_materialized"] is False
    assert truth["real_candidate_finite_correction_cycle_run"] is False
    assert truth["heldout_normalized_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["same_protocol_comparable_to_st006"] is False
    assert truth["pde_validated"] is False
    assert truth["final_normalized_momentum_gate"] == 1e-3
    assert truth["final_normalized_divergence_gate"] == 1e-5


def test_analytic_receipt_is_mechanics_only():
    receipt = analytic_regression_receipt()
    assert receipt["mechanics_only"] is True
    assert receipt["candidate_residual_evidence"] is False
    assert max(receipt["projection_errors"].values()) < 2e-14


def test_wrong_backend_type_is_rejected():
    class Fake:
        pass

    with pytest.raises(TypeError, match="ExactAgent2OscillatoryNonlinearHandoff"):
        materialize_oscillatory_nonlinear_mean_attribution(
            Fake(), np.asarray([0.4]), np.asarray([0.0]), np.asarray([0.5])  # type: ignore[arg-type]
        )
