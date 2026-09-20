from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_strict_inner_transport_mean import (
    _materialize_from_transport_provider,
)
from openai_ns_reconstruction.kokuno_strict_inner_transport_radial_stress import (
    PARENT_AGENT3_HEAD,
    RADIAL_OPERATOR_AGENT3_HEAD,
    StrictInnerTransportRadialGeometry,
    _materialize_from_mean_witness,
    _mechanics_transport_provider,
    build_mechanics_report,
    materialize_strict_inner_transport_radial_stress,
    truth_boundary,
)


def _mechanics_case():
    radii = np.linspace(0.20, 0.80, 49)
    geometry = StrictInnerTransportRadialGeometry(
        time=0.50,
        axial_z=0.13,
        radii=tuple(float(value) for value in radii),
        bump_center=0.50,
        bump_halfwidth=0.30,
    )
    mean_witness = _materialize_from_transport_provider(
        _mechanics_transport_provider(),
        radii,
        np.full_like(radii, geometry.axial_z),
        np.full_like(radii, geometry.time),
        agent1_backend=None,
        agent2_backend=None,
        backend_kind="manufactured-mechanics-only",
    )
    stress_witness = _materialize_from_mean_witness(
        mean_witness,
        geometry,
        backend_kind="manufactured-mechanics-only",
    )
    return radii, geometry, mean_witness, stress_witness


def test_scoped_transport_mean_routes_theta_and_axial_profiles_into_compact_stress():
    radii, _, mean_witness, stress_witness = _mechanics_case()
    expected_radial = 0.05 * radii
    expected_theta = 0.20 * (1.0 - ((radii - 0.50) / 0.30) ** 2) ** 2
    expected_axial = 0.12 + 0.03 * radii

    np.testing.assert_allclose(
        stress_witness.radial_mean_profile, expected_radial, atol=5.0e-14, rtol=0.0
    )
    np.testing.assert_allclose(
        stress_witness.theta_mean_profile, expected_theta, atol=5.0e-14, rtol=0.0
    )
    np.testing.assert_allclose(
        stress_witness.axial_mean_profile, expected_axial, atol=5.0e-14, rtol=0.0
    )
    np.testing.assert_allclose(
        stress_witness.theta_mean_profile,
        mean_witness.mean_transport_cylindrical[:, 1],
        atol=0.0,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        stress_witness.axial_mean_profile,
        mean_witness.mean_transport_cylindrical[:, 2],
        atol=0.0,
        rtol=0.0,
    )

    theta = stress_witness.theta_stress
    axial = stress_witness.axial_stress
    assert theta["exponent"] == 2
    assert axial["exponent"] == 1
    assert abs(float(theta["bump_weighted_integral"]) - 1.0) < 2.0e-14
    assert abs(float(axial["bump_weighted_integral"]) - 1.0) < 2.0e-14
    assert abs(float(theta["moment_complement_weighted_moment"])) < 2.0e-14
    assert abs(float(axial["moment_complement_weighted_moment"])) < 2.0e-14
    assert abs(float(theta["stress_inner_edge"])) < 2.0e-14
    assert abs(float(theta["stress_outer_edge"])) < 2.0e-14
    assert abs(float(axial["stress_inner_edge"])) < 2.0e-14
    assert abs(float(axial["stress_outer_edge"])) < 2.0e-14
    assert float(theta["stress_rms"]) > 1.0e-6
    assert float(axial["stress_rms"]) > 1.0e-6


def test_geometry_and_mean_identity_mismatches_fail_closed():
    radii, geometry, mean_witness, _ = _mechanics_case()

    with pytest.raises(ValueError, match="at least nine"):
        StrictInnerTransportRadialGeometry(
            time=0.5,
            axial_z=0.0,
            radii=tuple(np.linspace(0.2, 0.8, 8)),
            bump_center=0.5,
            bump_halfwidth=0.3,
        )

    with pytest.raises(ValueError, match="strictly increasing"):
        StrictInnerTransportRadialGeometry(
            time=0.5,
            axial_z=0.0,
            radii=(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.5, 0.8),
            bump_center=0.5,
            bump_halfwidth=0.3,
        )

    with pytest.raises(ValueError, match="support"):
        StrictInnerTransportRadialGeometry(
            time=0.5,
            axial_z=0.0,
            radii=tuple(np.linspace(0.3, 0.7, 17)),
            bump_center=0.5,
            bump_halfwidth=0.3,
        )

    shifted = StrictInnerTransportRadialGeometry(
        time=geometry.time,
        axial_z=geometry.axial_z,
        radii=tuple(float(value) for value in (radii + 1.0e-4)),
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    with pytest.raises(ValueError, match="radial grid"):
        _materialize_from_mean_witness(
            mean_witness,
            shifted,
            backend_kind="negative-control",
        )


def test_public_materializer_has_no_surrogate_mean_stress_or_retuning_inputs():
    parameters = inspect.signature(
        materialize_strict_inner_transport_radial_stress
    ).parameters
    assert list(parameters) == ["transport_backend", "inner_backend", "geometry"]
    forbidden = {
        "residual",
        "defect",
        "mean",
        "source",
        "stress",
        "inverse",
        "pressure",
        "forcing",
        "target",
        "gain",
        "alpha",
        "damping",
        "angular_order",
        "spatial_step",
        "time_step",
        "viscosity",
        "nu",
        "delta_y",
        "delta_a",
        "normalized_score",
        "scientific_threshold",
    }
    assert forbidden.isdisjoint(parameters)


def test_truth_boundary_keeps_scoped_stress_out_of_correction_cycle():
    boundary = truth_boundary()
    assert PARENT_AGENT3_HEAD == "7c3e5c23df9125722f055df4696f593ea17dd0e8"
    assert RADIAL_OPERATOR_AGENT3_HEAD == "56024553b981833a283f98648c8599011d6f38ab"
    assert boundary["strict_inner_transport_radial_stress_materialized"] is True
    assert boundary["mean_recomputed_from_exact_upstream_in_public_api"] is True
    assert boundary["caller_supplied_mean_allowed"] is False
    assert boundary["caller_supplied_stress_allowed"] is False
    assert boundary["compact_support_geometry_checked"] is True
    assert boundary["moment_complement_recorded"] is True
    assert boundary["theta_exponent"] == 2
    assert boundary["axial_exponent"] == 1
    assert boundary["complete_ns_defect"] is False
    assert boundary["pressure_gradient_included"] is False
    assert boundary["restricted_forcing_included"] is False
    assert boundary["global_corrected_leading_join_materialized"] is False
    assert boundary["correction_transport_included"] is False
    assert boundary["scoped_transport_stress_authorized_as_correction_target"] is False
    assert boundary["finite_head_debt_authorized_from_scoped_transport"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["same_protocol_comparable_to_st006"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5


def test_mechanics_report_is_explicitly_not_candidate_residual_evidence():
    report = build_mechanics_report()
    assert report["mechanics_only"] is True
    evaluation = report["evaluation"]
    assert evaluation["backend_kind"] == "manufactured-mechanics-only"
    assert evaluation["theta_e2"]["stress_rms"] > 0.0
    assert evaluation["axial_e1"]["stress_rms"] > 0.0
    boundary = report["truth_boundary"]
    assert boundary["complete_ns_defect"] is False
    assert boundary["scoped_transport_stress_authorized_as_correction_target"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
