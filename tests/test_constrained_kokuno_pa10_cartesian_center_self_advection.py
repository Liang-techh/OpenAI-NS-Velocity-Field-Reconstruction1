from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_cartesian_center_self_advection import (
    KokunoPA10CartesianCenterSelfAdvection,
    build_report,
)


def test_self_advection_is_exact_jacobian_velocity_contraction() -> None:
    candidate = KokunoPA10CartesianCenterSelfAdvection()
    x0, x1 = candidate.spatial.source_X_interval
    X = x0 + (x1 - x0) * np.array([0.19, 0.27, 0.35])
    eta = np.array([-0.28, 0.03, 0.31])
    t = np.array([0.46, 0.50, 0.54])
    theta = np.array([0.37, 1.29, 2.41])
    p = candidate.field.cartesian_from_similarity(X, eta, t, theta)

    velocity = candidate.velocity(p["x"], p["y"], p["z"], p["t"])
    jacobian = candidate.velocity_jacobian(p["x"], p["y"], p["z"], p["t"])
    expected = np.einsum("...ij,...j->...i", jacobian, velocity)
    actual = candidate.self_advection(p["x"], p["y"], p["z"], p["t"])

    assert actual.shape == (3, 3)
    assert np.array_equal(actual, expected)
    assert np.all(np.isfinite(actual))
    assert float(np.sqrt(np.mean(actual * actual))) > 1.0e-10


def test_self_advection_report_frozen_engineering_gates() -> None:
    report = build_report()
    checks = report["machine_checks"]

    assert checks["production_matches_explicit_jacobian_contraction_abs_max"] == 0.0
    assert checks["self_advection_vs_directional_fd4_fine_relative_max"] < 5.0e-5
    assert checks["directional_fd4_fine_vs_coarse_relative_max"] < 5.0e-5
    assert checks["rotation_covariance_relative_max"] < 5.0e-10
    assert checks["axis_transverse_self_advection_abs_max"] < 1.0e-12
    assert checks["self_advection_nontrivial_rms"] > 1.0e-10
    assert checks["all_probe_values_finite"] is True

    truth = report["truth_boundary"]
    assert truth["inner_cartesian_center_self_advection_executable"] is True
    assert truth["self_advection_uses_finite_difference_in_production"] is False
    assert truth["source_center_is_final_corrected_fixed_point"] is False
    assert truth["global_cartesian_spacetime_leading_velocity_materialized"] is False
    assert truth["outer_join_localization_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["viscous_laplacian_materialized_here"] is False
    assert truth["complete_candidate_api_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    gates = report["scientific_gates"]
    assert gates["momentum_max_l2"] == 1.0e-3
    assert gates["divergence_max_l2"] == 1.0e-5
    assert gates["free_residual_defined_forcing_forbidden"] is True


def test_configuration_round_trip_preserves_semantics(tmp_path) -> None:
    candidate = KokunoPA10CartesianCenterSelfAdvection()
    path = tmp_path / "self_advection.json"
    payload = candidate.save_configuration(path)
    reloaded = KokunoPA10CartesianCenterSelfAdvection.load_configuration(path)

    assert reloaded.configuration() == payload
    assert reloaded.field_sha256 == candidate.field_sha256
    assert reloaded.temporal_derivative_sha256 == candidate.temporal_derivative_sha256
    assert reloaded.spatial_derivative_sha256 == candidate.spatial_derivative_sha256
    assert reloaded.self_advection_sha256 == candidate.self_advection_sha256

    p = candidate.field.cartesian_from_similarity(0.2 * candidate.spatial.source_X_interval[1], 0.1, 0.5, 0.7)
    before = candidate.self_advection(p["x"], p["y"], p["z"], p["t"])
    after = reloaded.self_advection(p["x"], p["y"], p["z"], p["t"])
    assert np.array_equal(before, after)


def test_configuration_fails_closed_on_realization_drift() -> None:
    candidate = KokunoPA10CartesianCenterSelfAdvection()
    payload = copy.deepcopy(candidate.configuration())
    payload["self_advection_realization"] = "finite-difference"
    with pytest.raises(ValueError, match="unsupported self-advection realization"):
        KokunoPA10CartesianCenterSelfAdvection.from_configuration(payload)
