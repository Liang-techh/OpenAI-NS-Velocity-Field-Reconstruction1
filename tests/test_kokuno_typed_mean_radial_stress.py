from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_compact_radial_stress_adapter import (
    _compact_radial_stress,
)
from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)
from openai_ns_reconstruction.kokuno_typed_mean_radial_stress import (
    RadialMeanStressGeometry,
    deterministic_receipt,
    materialize_actual_mean_radial_stress,
    truth_boundary,
)


def _fixture():
    identity = CycleIdentity("typed-radial-test-v1", 0, "t0")
    center = 0.50
    halfwidth = 0.30

    def psi(x: float, y: float) -> float:
        radius = math.hypot(x, y)
        s = (radius - center) / halfwidth
        if abs(s) >= 1.0:
            return 0.0
        return math.cos(0.5 * math.pi * s) ** 8

    velocity = VectorFieldProvider(
        identity,
        "test:u=t*psi*(-y,x,1)",
        lambda x, y, z, t: (-t * y * psi(x, y), t * x * psi(x, y), t * psi(x, y)),
    )
    velocity_dt = VectorFieldProvider(
        identity,
        "test:u_t=psi*(-y,x,1)",
        lambda x, y, z, t: (-y * psi(x, y), x * psi(x, y), psi(x, y)),
    )
    pressure = ScalarFieldProvider(identity, "test:p=0", lambda x, y, z, t: 0.0)
    forcing = RestrictedForcingProvider(
        identity,
        "test:f=0",
        "zero forcing fixed independently of residual",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    radii = tuple(float(value) for value in np.linspace(0.20, 0.80, 49))
    geometry = RadialMeanStressGeometry(
        time=0.0,
        axial_z=0.11,
        radii=radii,
        bump_center=center,
        bump_halfwidth=halfwidth,
        angular_count=32,
        phase=0.023,
    )
    expected_psi = np.asarray([psi(radius, 0.0) for radius in radii], dtype=float)
    return velocity, velocity_dt, pressure, forcing, geometry, expected_psi


def test_actual_mean_profiles_feed_existing_radial_operator_without_surrogate():
    velocity, velocity_dt, pressure, forcing, geometry, expected_psi = _fixture()
    result = materialize_actual_mean_radial_stress(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        geometry,
        viscosity=0.01,
        spatial_step=0.005,
    )
    radii = np.asarray(geometry.radii, dtype=float)
    expected_theta = radii * expected_psi

    assert np.max(np.abs(result.radial_mean_profile)) <= 2.0e-12
    assert np.max(np.abs(result.theta_mean_profile - expected_theta)) <= 2.0e-12
    assert np.max(np.abs(result.axial_mean_profile - expected_psi)) <= 2.0e-12

    direct_theta = _compact_radial_stress(
        radii,
        expected_theta,
        exponent=2,
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    direct_axial = _compact_radial_stress(
        radii,
        expected_psi,
        exponent=1,
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    np.testing.assert_allclose(result.theta_stress["stress"], direct_theta["stress"], rtol=0.0, atol=2.0e-12)
    np.testing.assert_allclose(result.axial_stress["stress"], direct_axial["stress"], rtol=0.0, atol=2.0e-12)
    assert abs(float(result.theta_stress["moment_complement_weighted_moment"])) <= 2.0e-12
    assert abs(float(result.axial_stress["moment_complement_weighted_moment"])) <= 2.0e-12
    assert abs(float(result.theta_stress["stress_inner_edge"])) <= 2.0e-12
    assert abs(float(result.theta_stress["stress_outer_edge"])) <= 2.0e-12
    assert abs(float(result.axial_stress["stress_inner_edge"])) <= 2.0e-12
    assert abs(float(result.axial_stress["stress_outer_edge"])) <= 2.0e-12


def test_public_api_cannot_accept_precomputed_success_inputs():
    parameters = set(inspect.signature(materialize_actual_mean_radial_stress).parameters)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "source",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
    }
    assert parameters.isdisjoint(forbidden)
    boundary = truth_boundary()
    assert boundary["raw_same_cycle_providers_required"] is True
    assert boundary["caller_supplied_precomputed_defect_allowed"] is False
    assert boundary["caller_supplied_mean_source_allowed"] is False
    assert boundary["caller_supplied_stress_allowed"] is False
    assert boundary["actual_mean_channels_fed_directly_to_existing_radial_operator"] is True
    assert boundary["radial_operator_revalidated_here"] is False


def test_same_cycle_identity_mismatch_fails_closed():
    velocity, _, pressure, forcing, geometry, _ = _fixture()
    wrong = CycleIdentity("other-cycle", 0, "t0")
    velocity_dt = VectorFieldProvider(
        wrong,
        "wrong:u_t",
        lambda x, y, z, t: (-y, x, 1.0),
    )
    with pytest.raises(ValueError, match="same-cycle identity"):
        materialize_actual_mean_radial_stress(
            velocity,
            velocity_dt,
            pressure,
            forcing,
            geometry,
            viscosity=0.01,
            spatial_step=0.005,
        )


def test_geometry_rejects_noncompact_or_invalid_radial_sampling():
    radii = tuple(float(value) for value in np.linspace(0.20, 0.80, 49))
    with pytest.raises(ValueError, match="support"):
        RadialMeanStressGeometry(
            time=0.0,
            axial_z=0.0,
            radii=radii,
            bump_center=0.50,
            bump_halfwidth=0.40,
        )
    with pytest.raises(ValueError, match="strictly increasing"):
        RadialMeanStressGeometry(
            time=0.0,
            axial_z=0.0,
            radii=(0.20, 0.25, 0.30, 0.35, 0.40, 0.40, 0.50, 0.60, 0.70),
            bump_center=0.45,
            bump_halfwidth=0.20,
        )


def test_receipt_keeps_full_candidate_and_validation_claims_false():
    receipt = deterministic_receipt()
    analytic = receipt["analytic_regression"]
    assert analytic["maximum_theta_mean_error"] <= 2.0e-12
    assert analytic["maximum_axial_mean_error"] <= 2.0e-12
    assert analytic["maximum_radial_mean_error"] <= 2.0e-12
    boundary = receipt["truth_boundary"]
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["full_same_cycle_composite_requested_stress_materialized"] is False
    assert boundary["candidate_finite_head_mean_debt_materialized"] is False
    assert boundary["signed_mean_inverse_input_ready"] is False
    assert boundary["public_velocity_correction_materialized"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
