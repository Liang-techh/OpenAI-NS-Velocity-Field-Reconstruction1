from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_cartesian_center_transport_viscous import (
    VISCOSITY,
    KokunoPA10CartesianCenterTransportViscous,
    _assemble_transport,
    build_report,
)


def test_transport_assembly_signs_and_fixed_viscosity() -> None:
    dt = np.array([[1.0, -2.0, 3.0], [0.5, 0.25, -0.75]])
    adv = np.array([[4.0, 5.0, -6.0], [-1.0, 2.0, 3.0]])
    lap = np.array([[10.0, -20.0, 30.0], [40.0, 50.0, -60.0]])
    expected = dt + adv - 0.01 * lap
    actual = _assemble_transport(dt, adv, lap)
    assert VISCOSITY == 0.01
    assert np.array_equal(actual, expected)


def test_transport_report_frozen_scoped_gates() -> None:
    report = build_report()
    checks = report["machine_checks"]

    assert checks["viscosity"] == VISCOSITY
    assert checks["assembly_closure_max_abs"] < 1.0e-12
    assert checks["transport_nontrivial_rms"] > 1.0e-10
    assert checks["time_term_rms"] > 0.0
    assert checks["self_advection_rms"] > 0.0
    assert checks["viscous_term_rms"] > 0.0
    assert checks["inner_boundary_stencil_fails_closed"] is True
    assert checks["all_probe_values_finite"] is True

    truth = report["truth_boundary"]
    assert truth["inner_cartesian_center_transport_viscous_executable"] is True
    assert truth["transport_viscosity_is_repository_choice"] is True
    assert truth["transport_viscosity_caller_tunable"] is False
    assert truth["matched_pressure_gradient_included"] is False
    assert truth["restricted_forcing_included"] is False
    assert truth["complete_momentum_residual"] is False
    assert truth["source_center_is_final_corrected_fixed_point"] is False
    assert truth["global_cartesian_spacetime_leading_velocity_materialized"] is False
    assert truth["outer_join_localization_materialized"] is False
    assert truth["complete_candidate_api_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    gates = report["scientific_gates"]
    assert gates["momentum_max_l2"] == 1.0e-3
    assert gates["divergence_max_l2"] == 1.0e-5
    assert gates["free_residual_defined_forcing_forbidden"] is True


def test_configuration_round_trip_preserves_transport_semantics(tmp_path) -> None:
    candidate = KokunoPA10CartesianCenterTransportViscous()
    path = tmp_path / "transport.json"
    payload = candidate.save_configuration(path)
    reloaded = KokunoPA10CartesianCenterTransportViscous.load_configuration(path)

    assert reloaded.configuration() == payload
    assert reloaded.field_sha256 == candidate.field_sha256
    assert reloaded.temporal_derivative_sha256 == candidate.temporal_derivative_sha256
    assert reloaded.spatial_derivative_sha256 == candidate.spatial_derivative_sha256
    assert reloaded.self_advection_sha256 == candidate.self_advection_sha256
    assert reloaded.laplacian_sha256 == candidate.laplacian_sha256
    assert reloaded.transport_viscous_sha256 == candidate.transport_viscous_sha256

    x0, x1 = candidate.spatial.source_X_interval
    p = candidate.field.cartesian_from_similarity(
        x0 + 0.23 * (x1 - x0), 0.08, 0.5, 0.7
    )
    args = (p["x"], p["y"], p["z"], p["t"])
    before = candidate.transport_viscous(*args)
    after = reloaded.transport_viscous(*args)
    assert np.array_equal(before, after)


def test_inner_boundary_transport_fails_closed_with_laplacian_stencil() -> None:
    candidate = KokunoPA10CartesianCenterTransportViscous()
    _, x1 = candidate.spatial.source_X_interval
    p = candidate.field.cartesian_from_similarity(x1, 0.0, 0.5, 0.0)
    with pytest.raises(ValueError, match="source inner X interval"):
        candidate.transport_viscous(p["x"], p["y"], p["z"], p["t"])


def test_configuration_fails_closed_on_viscosity_or_operator_drift() -> None:
    candidate = KokunoPA10CartesianCenterTransportViscous()

    payload = copy.deepcopy(candidate.configuration())
    payload["viscosity"] = 0.02
    with pytest.raises(ValueError, match="viscosity drifted"):
        KokunoPA10CartesianCenterTransportViscous.from_configuration(payload)

    payload = copy.deepcopy(candidate.configuration())
    payload["operator"] = "residual-defined-forcing"
    with pytest.raises(ValueError, match="unsupported transport-viscous operator"):
        KokunoPA10CartesianCenterTransportViscous.from_configuration(payload)


def test_transport_assembly_rejects_bad_component_axis_and_nonfinite_values() -> None:
    with pytest.raises(ValueError, match="component axis"):
        _assemble_transport(np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2)))

    bad = np.zeros((1, 3))
    bad[0, 1] = np.nan
    with pytest.raises(ValueError, match="finite"):
        _assemble_transport(bad, np.zeros((1, 3)), np.zeros((1, 3)))
