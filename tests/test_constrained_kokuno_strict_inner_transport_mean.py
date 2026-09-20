from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_strict_inner_transport_mean import (
    AGENT1_HEAD,
    AGENT1_SOURCE_BLOB_SHA,
    AGENT2_HEAD,
    AGENT2_SOURCE_BLOB_SHA,
    ANGULAR_ORDERS,
    REPOSITORY_VISCOSITY,
    _materialize_from_transport_provider,
    _mechanics_provider,
    build_mechanics_report,
    materialize_strict_inner_transport_mean,
    truth_boundary,
)


def _mechanics_witness():
    return _materialize_from_transport_provider(
        _mechanics_provider(),
        np.array([0.7, 1.1]),
        np.array([-0.2, 0.3]),
        np.array([0.48, 0.52]),
        agent1_backend=None,
        agent2_backend=None,
        backend_kind="manufactured-mechanics-only",
    )


def test_rotating_frame_recovers_registered_transport_components_and_closure():
    witness = _mechanics_witness()
    expected_time = np.array([0.07, 0.02, 0.11])
    expected_advection = np.array([0.36, -0.03, 0.35])
    expected_viscous = np.array([-0.015, 0.01, -0.02])
    expected_transport = np.array([0.415, 0.0, 0.44])

    np.testing.assert_allclose(
        witness.mean_time_derivative_cylindrical,
        np.broadcast_to(expected_time, witness.mean_time_derivative_cylindrical.shape),
        atol=2.0e-14,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        witness.mean_advection_cylindrical,
        np.broadcast_to(expected_advection, witness.mean_advection_cylindrical.shape),
        atol=2.0e-14,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        witness.mean_viscous_cylindrical,
        np.broadcast_to(expected_viscous, witness.mean_viscous_cylindrical.shape),
        atol=2.0e-14,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        witness.mean_transport_cylindrical,
        np.broadcast_to(expected_transport, witness.mean_transport_cylindrical.shape),
        atol=2.0e-14,
        rtol=0.0,
    )
    assert witness.angular_orders == ANGULAR_ORDERS == (32, 64, 128)
    assert max(witness.successive_transport_mean_relative_differences) < 2.0e-14
    assert witness.pointwise_component_closure_absolute_max < 2.0e-14
    assert witness.projected_component_closure_absolute_max < 2.0e-14
    assert witness.transport_mean_rms > 0.0
    assert witness.full_ring_transport_rms >= witness.transport_mean_rms
    assert witness.transport_mean_to_full_rms_ratio <= 1.0 + 1.0e-12
    assert witness.viscosity == REPOSITORY_VISCOSITY == 0.01


def test_malformed_nonfinite_or_wrong_viscosity_payloads_fail_closed():
    base = _mechanics_provider()

    def wrong_shape(x, y, z, t):
        time_part, advection, viscous, transport, viscosity = base(x, y, z, t)
        return time_part[..., :2], advection, viscous, transport, viscosity

    with pytest.raises(ValueError, match="expected"):
        _materialize_from_transport_provider(
            wrong_shape,
            [0.8],
            [0.0],
            [0.5],
            agent1_backend=None,
            agent2_backend=None,
            backend_kind="negative-control",
        )

    def nonfinite(x, y, z, t):
        time_part, advection, viscous, transport, viscosity = base(x, y, z, t)
        transport = np.array(transport, copy=True)
        transport[..., 1] = np.nan
        return time_part, advection, viscous, transport, viscosity

    with pytest.raises(ValueError, match="non-finite"):
        _materialize_from_transport_provider(
            nonfinite,
            [0.8],
            [0.0],
            [0.5],
            agent1_backend=None,
            agent2_backend=None,
            backend_kind="negative-control",
        )

    def wrong_viscosity(x, y, z, t):
        time_part, advection, viscous, transport, _ = base(x, y, z, t)
        return time_part, advection, viscous, transport, 0.02

    with pytest.raises(ValueError, match="0.01"):
        _materialize_from_transport_provider(
            wrong_viscosity,
            [0.8],
            [0.0],
            [0.5],
            agent1_backend=None,
            agent2_backend=None,
            backend_kind="negative-control",
        )


def test_inconsistent_total_is_exposed_by_both_pointwise_and_projected_closure():
    base = _mechanics_provider()

    def inconsistent(x, y, z, t):
        time_part, advection, viscous, transport, viscosity = base(x, y, z, t)
        offset = np.zeros_like(transport)
        offset[..., 2] = 1.0e-4
        return time_part, advection, viscous, transport + offset, viscosity

    witness = _materialize_from_transport_provider(
        inconsistent,
        [0.8],
        [0.0],
        [0.5],
        agent1_backend=None,
        agent2_backend=None,
        backend_kind="negative-control",
    )
    assert witness.pointwise_component_closure_absolute_max >= 9.9e-5
    assert witness.projected_component_closure_absolute_max >= 9.9e-5


def test_public_materializer_exposes_no_retuning_or_residual_escape_hatches():
    parameters = inspect.signature(materialize_strict_inner_transport_mean).parameters
    assert list(parameters) == [
        "transport_backend",
        "inner_backend",
        "radius",
        "z",
        "t",
    ]
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
        "normalized_score",
        "alpha",
        "damping",
        "angular_order",
        "spatial_step",
        "time_step",
        "viscosity",
        "nu",
        "delta_y",
        "delta_a",
        "scientific_threshold",
    }
    assert forbidden.isdisjoint(parameters)


def test_truth_boundary_and_exact_upstream_pins_remain_fail_closed():
    boundary = truth_boundary()
    assert AGENT1_HEAD == "e96ee90144976a992b62d76d4361a94eb16bd91e"
    assert AGENT1_SOURCE_BLOB_SHA == "74b0e185e8e0bbb08695d493e0a091c305a03006"
    assert AGENT2_HEAD == "016e3c152d7d11e4fedcd79b819df337f8942983"
    assert AGENT2_SOURCE_BLOB_SHA == "d6082ebd346ba8150321bf3d62f104c107a70300"
    assert boundary["strict_inner_pressure_forcing_free_transport_mean_materialized"] is True
    assert boundary["time_derivative_included"] is True
    assert boundary["convective_term_included"] is True
    assert boundary["base_viscous_term_included"] is True
    assert boundary["complete_ns_defect"] is False
    assert boundary["pressure_gradient_included"] is False
    assert boundary["restricted_forcing_included"] is False
    assert boundary["correction_transport_included"] is False
    assert boundary["global_corrected_leading_join_materialized"] is False
    assert boundary["mean_correction_velocity_materialized_from_this_witness"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5


def test_mechanics_report_is_explicitly_not_real_candidate_evidence():
    report = build_mechanics_report()
    assert report["mechanics_only"] is True
    assert report["backend_kind"] == "manufactured-mechanics-only"
    np.testing.assert_allclose(report["expected_time_mean"], [0.07, 0.02, 0.11])
    np.testing.assert_allclose(report["expected_advection_mean"], [0.36, -0.03, 0.35])
    np.testing.assert_allclose(report["expected_viscous_mean"], [-0.015, 0.01, -0.02])
    np.testing.assert_allclose(report["expected_transport_mean"], [0.415, 0.0, 0.44])
    assert report["truth_boundary"]["complete_ns_defect"] is False
