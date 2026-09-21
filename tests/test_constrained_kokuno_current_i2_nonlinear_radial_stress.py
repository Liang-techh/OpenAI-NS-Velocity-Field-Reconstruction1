from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_i2_nonlinear_radial_stress import (
    PARENT_AGENT3_SOURCE_BLOB,
    _materialize_current_mean_witness,
    materialize_current_i2_nonlinear_radial_stress,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_oscillatory_nonlinear_mean_attribution import (
    _analytic_provider,
    _materialize_from_provider,
)
from openai_ns_reconstruction.kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)


def _mechanics_fixture():
    radii = np.linspace(0.02, 0.78, 39)
    geometry = StrictInnerTransportRadialGeometry(
        time=0.50,
        axial_z=0.13,
        radii=tuple(float(v) for v in radii),
        bump_center=0.50,
        bump_halfwidth=0.22,
    )
    witness = _materialize_from_provider(
        _analytic_provider,
        radii,
        np.full_like(radii, geometry.axial_z),
        np.full_like(radii, geometry.time),
        backend=None,
        backend_kind="analytic-current-i2-radial-stress-routing-mechanics-only",
    )
    routed = _materialize_current_mean_witness(
        witness,
        geometry,
        parent_source_blob=PARENT_AGENT3_SOURCE_BLOB,
    )
    return routed


def test_current_i2_mean_routes_through_inherited_compact_radial_operator():
    routed = _mechanics_fixture()
    inherited = routed.inherited

    assert routed.mean_witness is inherited.mean_witness
    assert routed.mean_piece_closure_absolute_max <= 1.0e-12
    assert routed.theta_stress_piece_closure_relative_max <= 2.0e-11
    assert routed.axial_stress_piece_closure_relative_max <= 2.0e-11

    reports = (
        inherited.theta_quadratic_stress,
        inherited.theta_mixed_stress,
        inherited.theta_aggregate_stress,
        inherited.axial_quadratic_stress,
        inherited.axial_mixed_stress,
        inherited.axial_aggregate_stress,
    )
    for report in reports:
        scale = max(1.0, abs(float(report["weighted_moment"])))
        assert abs(float(report["moment_complement_weighted_moment"])) <= 2.0e-12 * scale
        assert abs(float(report["stress_outer_edge"])) <= 2.0e-12 * scale


def test_receipt_records_exact_pr1073_current_i2_lineage():
    receipt = _mechanics_fixture().to_receipt()
    boundary = receipt["truth_boundary"]

    assert receipt["backend_kind"] == (
        "exact-current-a3-pr1073-mean-to-inherited-compact-radial-stress"
    )
    assert receipt["parent_agent3_1073_source_blob"] == PARENT_AGENT3_SOURCE_BLOB
    assert boundary["parent_agent3_pr"] == 1073
    assert boundary["parent_agent3_head"] == "a798056d7ddabe0bf4020082861807bc5ac1efbd"
    assert boundary["agent2_composite_pr"] == 1071
    assert boundary["agent2_composite_head"] == "48d69e37e78e7f7f0e4e9936f28ff7974719288d"
    assert boundary["agent1_leading_pr"] == 1061
    assert boundary["agent1_leading_head"] == "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3"
    assert boundary["current_I2_nonlinear_radial_inverse_performed"] is True
    assert boundary["current_I2_quadratic_radial_stress_materialized"] is True
    assert boundary["public_float64_I2_delta_verified"] is False
    assert boundary["radial_force_materialized_in_this_increment"] is False
    assert boundary["complete_ns_defect"] is False
    assert boundary["residual_reduction_claimed"] is False


def test_parent_source_blob_is_fail_closed():
    routed = _mechanics_fixture()
    with pytest.raises(ValueError, match="exact A3 #1073"):
        _materialize_current_mean_witness(
            routed.mean_witness,
            routed.inherited.geometry,
            parent_source_blob="0" * 40,
        )


def test_public_contract_has_no_scientific_tuning_knobs():
    params = inspect.signature(materialize_current_i2_nonlinear_radial_stress).parameters
    assert list(params) == ["backend", "geometry"]
    for forbidden in (
        "residual", "defect", "mean", "source", "stress", "pressure", "forcing",
        "gain", "damping", "spatial_step", "derivative_step", "angular_order",
        "viscosity", "correction", "stage_budget", "scientific_threshold",
    ):
        assert forbidden not in params


def test_materializer_rejects_untyped_backend_before_any_residual_work():
    geometry = _mechanics_fixture().inherited.geometry
    with pytest.raises(TypeError, match="ExactCurrentI2NonlinearBackend"):
        materialize_current_i2_nonlinear_radial_stress(object(), geometry)


def test_truth_boundary_keeps_force_full_ns_and_cycle_claims_closed():
    boundary = truth_boundary()
    assert boundary["current_I2_leading_plus_oscillatory_identity_consumed"] is True
    assert boundary["current_I2_mixed_nonlinear_mean_consumed"] is True
    assert boundary["current_I2_quadratic_oscillatory_mean_consumed"] is True
    assert boundary["current_I2_aggregate_nonlinear_mean_consumed"] is True
    assert boundary["moment_complement_preserved_by_inherited_operator"] is True
    assert boundary["polynomial_first_cell_axis_extrapolation_inherited"] is True
    assert boundary["outer_edge_stress_closure_inherited"] is True
    assert boundary["independent_third_radial_moment_inverse_required"] is False
    assert boundary["radial_force_requires_later_partial_z_sigma_1_step"] is True
    assert boundary["agent2_curl_or_jacobian_reimplemented_by_agent3"] is False
    assert boundary["forbidden_public_parameters_absent"] is True
    assert boundary["public_parameters"] == ("backend", "geometry")

    for key in (
        "current_I3_overlay_consumed",
        "current_I4_overlay_consumed",
        "outer_global_velocity_consumed",
        "radial_force_materialized_in_this_increment",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect",
        "scoped_current_I2_radial_stress_authorized_as_complete_ns_defect",
        "scoped_current_I2_radial_stress_authorized_as_correction_target",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
        "blowup_proof_claimed",
    ):
        assert boundary[key] is False, key
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5
