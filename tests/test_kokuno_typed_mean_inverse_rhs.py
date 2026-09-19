from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)
from openai_ns_reconstruction.kokuno_typed_mean_inverse_rhs import (
    MeanInverseReferenceScale,
    deterministic_receipt,
    materialize_typed_mean_inverse_rhs,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_typed_mean_radial_stress import RadialMeanStressGeometry


def _fixture():
    identity = CycleIdentity("typed-mean-inverse-rhs-test-v1", 0, "t0")
    center = 0.50
    halfwidth = 0.30
    radii = tuple(float(value) for value in np.linspace(0.20, 0.80, 49))

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
    geometry = RadialMeanStressGeometry(
        time=0.0,
        axial_z=0.11,
        radii=radii,
        bump_center=center,
        bump_halfwidth=halfwidth,
        angular_count=32,
        phase=0.023,
    )
    reference = MeanInverseReferenceScale(
        identity=identity,
        radii=radii,
        epsilon=tuple(0.23 + 0.04 * radius for radius in radii),
        producer_kind="analytic-regression",
        provenance="fixed before defect evaluation",
        source_reference_scale_certified=False,
    )
    return velocity, velocity_dt, pressure, forcing, geometry, reference


def test_typed_actual_debt_maps_profilewise_to_delta_c_over_epsilon():
    velocity, velocity_dt, pressure, forcing, geometry, reference = _fixture()
    result = materialize_typed_mean_inverse_rhs(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        geometry,
        reference,
        viscosity=0.01,
        spatial_step=0.005,
    )
    epsilon = np.asarray(reference.epsilon, dtype=float)
    np.testing.assert_allclose(result.theta_rhs * epsilon, result.debt.theta_debt, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(result.axial_rhs * epsilon, result.debt.axial_debt, rtol=0.0, atol=1.0e-15)
    assert result.rhs_rms > 0.0
    assert result.rhs_max_abs > 0.0
    assert result.reconstruction_closure_max_abs <= 1.0e-14


def test_public_api_has_no_naked_epsilon_or_precomputed_scientific_target_inputs():
    parameters = set(inspect.signature(materialize_typed_mean_inverse_rhs).parameters)
    forbidden = {
        "epsilon",
        "residual",
        "defect",
        "mean_source",
        "source",
        "stress",
        "requested_stress",
        "debt",
        "delta_c",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "h_ref",
        "inverse_matrix",
        "correction",
    }
    assert parameters.isdisjoint(forbidden)
    boundary = truth_boundary()
    assert boundary["naked_epsilon_public_parameter_allowed"] is False
    assert boundary["provenance_bearing_reference_scale_required"] is True
    assert boundary["caller_supplied_precomputed_defect_allowed"] is False
    assert boundary["caller_supplied_precomputed_debt_allowed"] is False
    assert boundary["caller_supplied_inverse_rhs_allowed"] is False
    assert boundary["profile_level_delta_c_over_epsilon_materialized"] is True


def test_reference_scale_identity_and_grid_are_fail_closed():
    velocity, velocity_dt, pressure, forcing, geometry, reference = _fixture()
    wrong_identity = MeanInverseReferenceScale(
        identity=CycleIdentity("other-cycle", 0, "t0"),
        radii=reference.radii,
        epsilon=reference.epsilon,
        producer_kind="analytic-regression",
        provenance="wrong cycle mutation",
    )
    with pytest.raises(ValueError, match="same correction-cycle identity"):
        materialize_typed_mean_inverse_rhs(
            velocity,
            velocity_dt,
            pressure,
            forcing,
            geometry,
            wrong_identity,
            viscosity=0.01,
            spatial_step=0.005,
        )

    shifted = MeanInverseReferenceScale(
        identity=reference.identity,
        radii=tuple(radius + 1.0e-4 for radius in reference.radii),
        epsilon=reference.epsilon,
        producer_kind="analytic-regression",
        provenance="wrong radial grid mutation",
    )
    with pytest.raises(ValueError, match="radial grid"):
        materialize_typed_mean_inverse_rhs(
            velocity,
            velocity_dt,
            pressure,
            forcing,
            geometry,
            shifted,
            viscosity=0.01,
            spatial_step=0.005,
        )


def test_reference_scale_rejects_nonpositive_epsilon():
    _, _, _, _, _, reference = _fixture()
    bad = list(reference.epsilon)
    bad[len(bad) // 2] = 0.0
    with pytest.raises(ValueError, match="strictly positive"):
        MeanInverseReferenceScale(
            identity=reference.identity,
            radii=reference.radii,
            epsilon=tuple(bad),
            producer_kind="analytic-regression",
            provenance="nonpositive epsilon mutation",
        )


def test_receipt_keeps_inverse_rhs_separate_from_inverse_and_scientific_claims():
    receipt = deterministic_receipt()
    evaluation = receipt["analytic_regression"]["evaluation"]
    rhs = evaluation["mean_inverse_rhs"]
    assert rhs["rhs_rms"] > 0.0
    assert rhs["rhs_max_abs"] > 0.0
    assert rhs["reconstruction_closure_max_abs"] <= 1.0e-14
    assert evaluation["reference_scale"]["source_reference_scale_certified"] is False

    boundary = receipt["truth_boundary"]
    assert boundary["formal_missing_weight_materialized"] is False
    assert boundary["h_ref_materialized"] is False
    assert boundary["signed_mean_inverse_solved"] is False
    assert boundary["signed_mean_inverse_input_ready"] is False
    assert boundary["public_velocity_correction_materialized"] is False
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["candidate_finite_head_mean_debt_materialized"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
