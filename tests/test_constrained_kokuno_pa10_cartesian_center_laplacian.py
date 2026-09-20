from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_cartesian_center_laplacian import (
    PRODUCTION_SPATIAL_STEP,
    KokunoPA10CartesianCenterLaplacian,
    _fd6_laplacian,
    build_report,
)


def test_fd6_laplacian_manufactured_polynomial() -> None:
    def velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack(
            (
                x**6 + 2.0 * y**4 - 3.0 * z**2 + t,
                x**4 * y**2 + y**6 + z**4 - 0.5 * t,
                2.0 * x**2 + 3.0 * y**2 + 5.0 * z**2 + x * y,
            ),
            axis=-1,
        )

    x = np.array([0.17, -0.23, 0.31])
    y = np.array([-0.19, 0.27, 0.13])
    z = np.array([0.22, -0.14, 0.29])
    t = np.array([0.4, 0.5, 0.6])
    actual = _fd6_laplacian(velocity, x, y, z, t, step=5.0e-3)
    expected = np.stack(
        (
            30.0 * x**4 + 24.0 * y**2 - 6.0,
            12.0 * x**2 * y**2 + 2.0 * x**4 + 30.0 * y**4 + 12.0 * z**2,
            np.full_like(x, 20.0),
        ),
        axis=-1,
    )
    scale = np.maximum(np.abs(expected), 1.0)
    assert np.max(np.abs(actual - expected) / scale) < 2.0e-7


def test_laplacian_report_frozen_engineering_gates() -> None:
    report = build_report()
    checks = report["machine_checks"]

    assert checks["production_step"] == PRODUCTION_SPATIAL_STEP
    assert checks["laplacian_vs_independent_fd4_fine_relative_max"] < 5.0e-3
    assert checks["independent_fd4_fine_vs_coarse_relative_max"] < 5.0e-3
    assert checks["rotation_covariance_relative_max"] < 5.0e-4
    assert checks["laplacian_nontrivial_rms"] > 1.0e-10
    assert checks["inner_boundary_stencil_fails_closed"] is True
    assert checks["all_probe_values_finite"] is True

    truth = report["truth_boundary"]
    assert truth["inner_cartesian_center_velocity_laplacian_executable"] is True
    assert truth["laplacian_uses_finite_difference_in_production"] is True
    assert truth["laplacian_is_public_source_formula"] is False
    assert truth["production_laplacian_step_caller_tunable"] is False
    assert truth["outside_inner_domain_zero_extension_used"] is False
    assert truth["source_center_is_final_corrected_fixed_point"] is False
    assert truth["global_cartesian_spacetime_leading_velocity_materialized"] is False
    assert truth["outer_join_localization_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["complete_candidate_api_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    gates = report["scientific_gates"]
    assert gates["momentum_max_l2"] == 1.0e-3
    assert gates["divergence_max_l2"] == 1.0e-5
    assert gates["free_residual_defined_forcing_forbidden"] is True


def test_configuration_round_trip_preserves_laplacian_semantics(tmp_path) -> None:
    candidate = KokunoPA10CartesianCenterLaplacian()
    path = tmp_path / "laplacian.json"
    payload = candidate.save_configuration(path)
    reloaded = KokunoPA10CartesianCenterLaplacian.load_configuration(path)

    assert reloaded.configuration() == payload
    assert reloaded.field_sha256 == candidate.field_sha256
    assert reloaded.temporal_derivative_sha256 == candidate.temporal_derivative_sha256
    assert reloaded.spatial_derivative_sha256 == candidate.spatial_derivative_sha256
    assert reloaded.self_advection_sha256 == candidate.self_advection_sha256
    assert reloaded.laplacian_sha256 == candidate.laplacian_sha256

    x0, x1 = candidate.spatial.source_X_interval
    p = candidate.field.cartesian_from_similarity(
        x0 + 0.22 * (x1 - x0), 0.1, 0.5, 0.7
    )
    before = candidate.velocity_laplacian(p["x"], p["y"], p["z"], p["t"])
    after = reloaded.velocity_laplacian(p["x"], p["y"], p["z"], p["t"])
    assert np.array_equal(before, after)


def test_inner_boundary_laplacian_fails_closed_instead_of_zero_extension() -> None:
    candidate = KokunoPA10CartesianCenterLaplacian()
    _, x1 = candidate.spatial.source_X_interval
    p = candidate.field.cartesian_from_similarity(x1, 0.0, 0.5, 0.0)
    with pytest.raises(ValueError, match="source inner X interval"):
        candidate.velocity_laplacian(p["x"], p["y"], p["z"], p["t"])


def test_configuration_fails_closed_on_step_or_realization_drift() -> None:
    candidate = KokunoPA10CartesianCenterLaplacian()

    payload = copy.deepcopy(candidate.configuration())
    payload["production_spatial_step"] = 2.0 * PRODUCTION_SPATIAL_STEP
    with pytest.raises(ValueError, match="step drifted"):
        KokunoPA10CartesianCenterLaplacian.from_configuration(payload)

    payload = copy.deepcopy(candidate.configuration())
    payload["laplacian_realization"] = "residual-defined"
    with pytest.raises(ValueError, match="unsupported Laplacian realization"):
        KokunoPA10CartesianCenterLaplacian.from_configuration(payload)
