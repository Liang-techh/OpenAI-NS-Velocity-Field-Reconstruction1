from __future__ import annotations

from dataclasses import replace
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_nonlinear_mean_attribution import (
    _analytic_provider,
    _materialize_from_provider,
)
from openai_ns_reconstruction.kokuno_oscillatory_quadratic_nonlinear_radial_stress import (
    MEAN_PIECE_CLOSURE_ABSOLUTE_GATE,
    PARENT_AGENT3_HEAD,
    PARENT_AGENT3_SOURCE_BLOB,
    RADIAL_OPERATOR_AGENT3_HEAD,
    RADIAL_OPERATOR_SOURCE_BLOB,
    STRESS_PIECE_CLOSURE_RELATIVE_GATE,
    _materialize_from_mean_witness,
    build_mechanics_report,
    materialize_oscillatory_quadratic_nonlinear_radial_stress,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)


def _fixture():
    radii = np.linspace(0.02, 0.80, 79)
    geometry = StrictInnerTransportRadialGeometry(
        time=0.50,
        axial_z=0.13,
        radii=tuple(float(v) for v in radii),
        bump_center=0.50,
        bump_halfwidth=0.28,
    )
    mean = _materialize_from_provider(
        _analytic_provider,
        radii,
        np.full_like(radii, geometry.axial_z),
        np.full_like(radii, geometry.time),
        backend=None,
        backend_kind="analytic-projection-regression-only",
    )
    result = _materialize_from_mean_witness(
        mean,
        geometry,
        backend_kind="manufactured-mechanics-only",
    )
    return radii, geometry, mean, result


def test_quadratic_mean_radial_stress_preserves_piece_attribution() -> None:
    _, _, mean, result = _fixture()

    assert result.mean_piece_closure_absolute_max <= MEAN_PIECE_CLOSURE_ABSOLUTE_GATE
    assert result.theta_stress_piece_closure_relative_max <= STRESS_PIECE_CLOSURE_RELATIVE_GATE
    assert result.axial_stress_piece_closure_relative_max <= STRESS_PIECE_CLOSURE_RELATIVE_GATE

    np.testing.assert_allclose(
        mean.mean_aggregate_nonlinear_cylindrical,
        mean.mean_mixed_cross_cylindrical
        + mean.mean_oscillatory_self_advection_cylindrical,
        rtol=0.0,
        atol=MEAN_PIECE_CLOSURE_ABSOLUTE_GATE,
    )

    assert np.max(np.abs(result.radial_quadratic_mean_profile)) > 1.0e-8
    assert np.max(np.abs(result.radial_mixed_mean_profile)) > 1.0e-8
    assert np.max(np.abs(result.radial_aggregate_mean_profile)) > 1.0e-8

    for stress in (
        result.theta_quadratic_stress,
        result.theta_mixed_stress,
        result.theta_aggregate_stress,
        result.axial_quadratic_stress,
        result.axial_mixed_stress,
        result.axial_aggregate_stress,
    ):
        assert stress["first_cell_polynomial_degree"] == 2
        assert stress["quadrature_lower_limit"] == 0.0
        assert abs(float(stress["bump_weighted_integral"]) - 1.0) <= 1.0e-11
        assert abs(float(stress["moment_complement_weighted_moment"])) <= (
            1.0e-9 * max(1.0, abs(float(stress["weighted_moment"])))
        )
        assert abs(float(stress["stress_outer_edge"])) <= (
            1.0e-9 * max(1.0, float(stress["stress_max_abs"]))
        )


def test_receipt_and_truth_boundary_stay_scoped() -> None:
    report = build_mechanics_report()
    assert report["mechanics_only"] is True
    assert report["candidate_residual_evidence"] is False
    assert report["registered_quadratic_mean"] == [0.17, -0.03, 0.09]
    assert report["registered_mixed_mean"] == [0.13, -0.04, 0.25]
    assert report["registered_aggregate_mean"] == [0.30, -0.07, 0.34]

    boundary = truth_boundary()
    assert boundary["parent_agent3_head"] == PARENT_AGENT3_HEAD
    assert boundary["parent_agent3_source_blob"] == PARENT_AGENT3_SOURCE_BLOB
    assert boundary["radial_operator_agent3_head"] == RADIAL_OPERATOR_AGENT3_HEAD
    assert boundary["radial_operator_source_blob"] == RADIAL_OPERATOR_SOURCE_BLOB
    assert boundary["scoped_quadratic_nonlinear_radial_inverse_performed"] is True
    assert boundary["radial_component_recorded_but_not_inverted"] is True
    assert boundary["stress_level_quadratic_mixed_attribution_checked"] is True
    assert boundary["agent2_curl_or_advection_reimplemented_by_agent3"] is False

    for key in (
        "pressure_gradient_included",
        "restricted_forcing_included",
        "global_corrected_leading_join_materialized",
        "complete_ns_defect",
        "scoped_quadratic_radial_stress_authorized_as_correction_target",
        "mean_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert boundary[key] is False


def test_public_materializer_has_no_scientific_tuning_surface() -> None:
    parameters = set(
        inspect.signature(
            materialize_oscillatory_quadratic_nonlinear_radial_stress
        ).parameters
    )
    assert parameters == {"handoff_backend", "geometry"}
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "viscosity",
        "nu", "angular_order", "spatial_step", "time_step", "delta_a",
        "normalized_score", "scientific_threshold",
    }
    assert parameters.isdisjoint(forbidden)
    assert truth_boundary()["forbidden_public_parameters_absent"] is True


def test_mean_piece_attribution_drift_fails_closed() -> None:
    _, geometry, mean, _ = _fixture()
    mutated = replace(
        mean,
        mean_aggregate_nonlinear_cylindrical=(
            np.asarray(mean.mean_aggregate_nonlinear_cylindrical, dtype=float)
            + np.array([0.0, 2.0e-6, 0.0])
        ),
    )
    with pytest.raises(ValueError, match="nonlinear mean attribution closure"):
        _materialize_from_mean_witness(
            mutated,
            geometry,
            backend_kind="deliberately-mutated-negative-control",
        )


def test_geometry_identity_drift_fails_closed() -> None:
    radii, geometry, mean, _ = _fixture()
    shifted = StrictInnerTransportRadialGeometry(
        time=geometry.time,
        axial_z=geometry.axial_z,
        radii=tuple(float(v) for v in (radii + 1.0e-5)),
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    with pytest.raises(ValueError, match="radial grid"):
        _materialize_from_mean_witness(
            mean,
            shifted,
            backend_kind="deliberately-shifted-negative-control",
        )
