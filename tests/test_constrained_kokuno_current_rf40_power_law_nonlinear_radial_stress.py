from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_rf40_power_law_nonlinear_radial_stress import (
    PARENT_AGENT3_SOURCE_BLOB,
    _materialize_current_mean_witness,
    materialize_current_rf40_power_law_nonlinear_radial_stress,
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
        backend_kind="analytic-rf40-power-law-routing-mechanics-only",
    )
    routed = _materialize_current_mean_witness(
        witness,
        geometry,
        parent_source_blob=PARENT_AGENT3_SOURCE_BLOB,
    )
    return routed


def test_rf40_power_law_mean_routes_through_inherited_compact_radial_operator():
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


def test_receipt_records_exact_pr1028_lineage_and_power_law_scope():
    receipt = _mechanics_fixture().to_receipt()
    boundary = receipt["truth_boundary"]

    assert receipt["backend_kind"] == (
        "exact-current-a3-pr1028-mean-to-inherited-compact-radial-stress"
    )
    assert receipt["parent_agent3_1028_source_blob"] == PARENT_AGENT3_SOURCE_BLOB
    assert boundary["current_RF40_power_law_nonlinear_radial_inverse_performed"] is True
    assert boundary["current_RF40_power_law_quadratic_radial_stress_materialized"] is True
    assert boundary["partial_domain_through_RF40_power_law_only"] is True
    assert boundary["velocity_after_RF40_power_law_materialized"] is False
    assert boundary["complete_ns_defect"] is False
    assert (
        boundary[
            "scoped_current_RF40_power_law_radial_stress_authorized_as_correction_target"
        ]
        is False
    )
    assert boundary["residual_reduction_claimed"] is False


def test_parent_source_blob_is_fail_closed():
    routed = _mechanics_fixture()
    witness = routed.mean_witness
    geometry = routed.inherited.geometry
    with pytest.raises(ValueError, match="exact A3 #1028"):
        _materialize_current_mean_witness(
            witness,
            geometry,
            parent_source_blob="0" * 40,
        )


def test_public_contract_has_no_scientific_tuning_knobs():
    params = inspect.signature(
        materialize_current_rf40_power_law_nonlinear_radial_stress
    ).parameters
    assert list(params) == ["backend", "geometry"]
    for forbidden in (
        "residual", "defect", "mean", "source", "stress", "pressure", "forcing",
        "gain", "damping", "spatial_step", "derivative_step", "angular_order",
        "viscosity", "correction", "stage_budget", "scientific_threshold",
    ):
        assert forbidden not in params


def test_truth_boundary_keeps_full_ns_and_correction_claims_closed():
    boundary = truth_boundary()
    assert (
        boundary["current_leading_plus_oscillation_through_RF40_power_law_consumed"]
        is True
    )
    assert boundary["moment_complement_preserved_by_inherited_operator"] is True
    assert boundary["polynomial_first_cell_axis_extrapolation_inherited"] is True
    assert boundary["independent_third_radial_moment_inverse_required"] is False
    assert boundary["radial_force_requires_later_partial_z_sigma_1_step"] is True
    assert boundary["agent2_curl_or_jacobian_reimplemented_by_agent3"] is False
    assert boundary["forbidden_public_parameters_absent"] is True
    for key in (
        "velocity_after_RF40_power_law_materialized",
        "full_post_XR_RF40_current_lineage_materialized",
        "current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed",
        "outer_global_leading_velocity_materialized",
        "pressure_gradient_included",
        "matched_global_pressure_materialized",
        "restricted_forcing_included",
        "complete_ns_defect",
        "scoped_current_RF40_power_law_radial_stress_authorized_as_correction_target",
        "mean_correction_velocity_materialized",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
        "blowup_proved",
    ):
        assert boundary[key] is False, key
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5
