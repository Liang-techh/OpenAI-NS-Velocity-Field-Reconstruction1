from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_transport_mean import (
    ANGULAR_ORDERS,
    REPOSITORY_VISCOSITY,
    ExactAgent2OscillatoryTransportHandoff,
    _materialize_from_provider,
    materialize_oscillatory_transport_mean,
    truth_boundary,
)


def _cylindrical_field(mean: tuple[float, float, float], mode: float):
    def field(x, y, z, t):
        xb, yb, zb, tb = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        del zb, tb
        theta = np.arctan2(yb, xb)
        c = np.cos(theta)
        s = np.sin(theta)
        m2 = mode * np.cos(2.0 * theta)
        radial = mean[0] + m2
        tangential = mean[1] - 0.35 * m2
        axial = mean[2] + 0.6 * mode * np.sin(2.0 * theta)
        return np.stack(
            (radial * c - tangential * s, radial * s + tangential * c, axial),
            axis=-1,
        )
    return field


def _mechanics_provider(*, break_before_after: bool = False):
    time_field = _cylindrical_field((0.10, -0.02, 0.03), 0.21)
    nonlinear_field = _cylindrical_field((0.20, -0.04, 0.05), -0.17)
    viscous_field = _cylindrical_field((-0.03, 0.01, -0.02), 0.09)
    inner_field = _cylindrical_field((0.40, 0.10, -0.20), 0.13)

    def provider(x, y, z, t):
        time_part = time_field(x, y, z, t)
        nonlinear = nonlinear_field(x, y, z, t)
        viscous = viscous_field(x, y, z, t)
        delta = time_part + nonlinear + viscous
        inner = inner_field(x, y, z, t)
        combined = inner + delta
        if break_before_after:
            theta = np.arctan2(np.asarray(y, dtype=float), np.asarray(x, dtype=float))
            combined = combined + np.stack(
                (1.0e-3 * np.cos(theta), 1.0e-3 * np.sin(theta), np.zeros_like(theta)),
                axis=-1,
            )
        return (
            time_part,
            nonlinear,
            viscous,
            delta,
            inner,
            combined,
            REPOSITORY_VISCOSITY,
        )
    return provider


def test_mechanics_projection_recovers_registered_means_and_two_closures():
    witness = _materialize_from_provider(
        _mechanics_provider(),
        np.array([0.6, 0.9]),
        np.array([-0.2, 0.25]),
        np.array([0.48, 0.52]),
        backend=None,
        backend_kind="manufactured-mechanics-only",
    )
    assert witness.angular_orders == ANGULAR_ORDERS
    np.testing.assert_allclose(
        witness.mean_time_increment_cylindrical,
        np.array([[0.10, -0.02, 0.03], [0.10, -0.02, 0.03]]),
        atol=2.0e-14,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        witness.mean_nonlinear_increment_cylindrical,
        np.array([[0.20, -0.04, 0.05], [0.20, -0.04, 0.05]]),
        atol=2.0e-14,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        witness.mean_viscous_increment_cylindrical,
        np.array([[-0.03, 0.01, -0.02], [-0.03, 0.01, -0.02]]),
        atol=2.0e-14,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        witness.mean_transport_increment_cylindrical,
        np.array([[0.27, -0.05, 0.06], [0.27, -0.05, 0.06]]),
        atol=2.0e-14,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        witness.mean_inner_plus_oscillatory_transport_cylindrical
        - witness.mean_inner_transport_cylindrical,
        witness.mean_transport_increment_cylindrical,
        atol=2.0e-14,
        rtol=0.0,
    )
    assert witness.pointwise_decomposition_closure_absolute_max < 2.0e-15
    assert witness.pointwise_before_after_closure_absolute_max < 2.0e-15
    assert witness.projected_decomposition_closure_absolute_max < 2.0e-14
    assert witness.projected_before_after_closure_absolute_max < 2.0e-14
    assert max(witness.successive_transport_mean_relative_differences) < 2.0e-14
    assert witness.transport_increment_mean_rms > 1.0e-12
    assert witness.transport_increment_mean_to_full_rms_ratio <= 1.0 + 1.0e-12


def test_before_after_transport_difference_is_not_silently_ignored():
    witness = _materialize_from_provider(
        _mechanics_provider(break_before_after=True),
        np.array([0.75]),
        np.array([0.1]),
        np.array([0.5]),
        backend=None,
        backend_kind="negative-control",
    )
    assert witness.pointwise_before_after_closure_absolute_max > 9.0e-4
    assert witness.projected_before_after_closure_absolute_max > 9.0e-4


def test_public_materializer_rejects_untyped_surrogate_backend():
    class Fake:
        pass

    with pytest.raises(TypeError, match="ExactAgent2OscillatoryTransportHandoff"):
        materialize_oscillatory_transport_mean(
            Fake(), np.array([0.7]), np.array([0.0]), np.array([0.5])
        )


def test_exact_binding_rejects_wrong_agent2_class_identity():
    class FakeHandoff:
        handoff_sha256 = "a" * 64
        source_contract_sha256 = "b" * 64

        def evaluate(self, x, y, z, t):
            raise AssertionError("not reached")

        def semantic_payload(self):
            return {}

    with pytest.raises(ValueError, match="module identity"):
        ExactAgent2OscillatoryTransportHandoff.bind(FakeHandoff())


def test_public_api_has_no_scientific_tuning_or_surrogate_inputs():
    forbidden = {
        "residual",
        "defect",
        "mean",
        "stress",
        "inverse",
        "pressure",
        "forcing",
        "target",
        "gain",
        "alpha",
        "damping",
        "nu",
        "viscosity",
        "delta_y",
        "delta_a",
        "angular_order",
        "spatial_step",
        "time_step",
        "scientific_threshold",
    }
    params = set(inspect.signature(materialize_oscillatory_transport_mean).parameters)
    assert not (params & forbidden)
    assert params == {"handoff_backend", "radius", "z", "t"}


def test_truth_boundary_stays_below_complete_defect_and_cycle():
    truth = truth_boundary()
    assert truth["oscillatory_transport_increment_materialized"] is True
    assert truth["cylindrical_m0_projection_materialized"] is True
    assert truth["pressure_gradient_included"] is False
    assert truth["restricted_forcing_included"] is False
    assert truth["global_corrected_leading_join_materialized"] is False
    assert truth["complete_ns_defect"] is False
    assert truth["radial_inverse_performed"] is False
    assert truth["scoped_transport_mean_authorized_as_correction_target"] is False
    assert truth["real_candidate_finite_correction_cycle_run"] is False
    assert truth["heldout_normalized_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["same_protocol_comparable_to_st006"] is False
    assert truth["pde_validated"] is False
