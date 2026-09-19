from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_autonomous_finite_head_factor_adapter import (
    EXPECTED_AUTONOMOUS_MISSING_WEIGHT,
    replay_autonomous_finite_head_factor,
)
from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)
from openai_ns_reconstruction.kokuno_typed_mean_finite_head_debt import (
    deterministic_receipt,
    materialize_typed_mean_finite_head_debt,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_typed_mean_radial_stress import (
    RadialMeanStressGeometry,
)


def _fixture():
    identity = CycleIdentity("a5-typed-mean-debt-test-v1", 0, "t0")
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
    geometry = RadialMeanStressGeometry(
        time=0.0,
        axial_z=0.11,
        radii=tuple(float(value) for value in np.linspace(0.20, 0.80, 49)),
        bump_center=center,
        bump_halfwidth=halfwidth,
        angular_count=32,
        phase=0.023,
    )
    return velocity, velocity_dt, pressure, forcing, geometry


def test_pure_factor_adapter_replays_frozen_586_weight_without_a2_graph():
    factor = replay_autonomous_finite_head_factor()
    assert factor["autonomous_missing_weight"] == pytest.approx(
        EXPECTED_AUTONOMOUS_MISSING_WEIGHT, abs=5.0e-15
    )
    assert factor["origin_pr"] == 514
    assert factor["replay_pr"] == 586
    assert factor["integration_adapter_only"] is True
    assert factor["formal_theorem_machine_bump_identity_claimed"] is False
    assert factor["formal_missing_weight_equality_claimed"] is False


def test_typed_actual_stress_maps_to_frozen_finite_head_debt_without_surrogate():
    velocity, velocity_dt, pressure, forcing, geometry = _fixture()
    result = materialize_typed_mean_finite_head_debt(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        geometry,
        viscosity=0.01,
        spatial_step=0.005,
    )
    weight = float(result.autonomous_factor["autonomous_missing_weight"])
    assert weight == pytest.approx(EXPECTED_AUTONOMOUS_MISSING_WEIGHT, abs=5.0e-15)

    theta_stress = np.asarray(result.radial_stress.theta_stress["stress"], dtype=float)
    axial_stress = np.asarray(result.radial_stress.axial_stress["stress"], dtype=float)
    np.testing.assert_allclose(result.theta_debt, -weight * theta_stress, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(result.axial_debt, -weight * axial_stress, rtol=0.0, atol=1.0e-15)
    assert result.requested_stress_rms > 0.0
    assert result.debt_rms > 0.0
    assert result.debt_rms / result.requested_stress_rms == pytest.approx(weight, rel=2.0e-15)
    assert result.identity_closure_max_abs <= 1.0e-14
    assert result.autonomous_factor["formal_theorem_machine_bump_identity_claimed"] is False
    assert result.autonomous_factor["formal_missing_weight_equality_claimed"] is False


def test_public_api_rejects_precomputed_defect_stress_factor_and_success_inputs():
    parameters = set(inspect.signature(materialize_typed_mean_finite_head_debt).parameters)
    forbidden = {
        "residual", "defect", "mean_source", "source", "stress", "requested_stress",
        "target", "factor", "missing_weight", "gain", "normalized_score",
        "scientific_threshold", "prepared_n", "band", "coordinate_q",
    }
    assert parameters.isdisjoint(forbidden)
    boundary = truth_boundary()
    assert boundary["raw_same_cycle_providers_required"] is True
    assert boundary["caller_supplied_precomputed_defect_allowed"] is False
    assert boundary["caller_supplied_mean_source_allowed"] is False
    assert boundary["caller_supplied_stress_allowed"] is False
    assert boundary["caller_supplied_finite_head_factor_allowed"] is False
    assert boundary["forbidden_public_parameters_absent"] is True
    assert boundary["typed_actual_stress_fed_directly_to_frozen_finite_head_factor"] is True
    assert boundary["typed_actual_state_finite_head_debt_executable"] is True
    assert boundary["integration_adapter_only"] is True


def test_same_cycle_identity_mismatch_still_fails_closed():
    velocity, _, pressure, forcing, geometry = _fixture()
    wrong = CycleIdentity("wrong-cycle", 0, "t0")
    velocity_dt = VectorFieldProvider(wrong, "wrong:u_t", lambda x, y, z, t: (-y, x, 1.0))
    with pytest.raises(ValueError, match="same-cycle identity"):
        materialize_typed_mean_finite_head_debt(
            velocity,
            velocity_dt,
            pressure,
            forcing,
            geometry,
            viscosity=0.01,
            spatial_step=0.005,
        )


def test_receipt_keeps_autonomous_factor_and_full_candidate_truth_separate():
    receipt = deterministic_receipt()
    analytic = receipt["analytic_regression"]
    assert analytic["autonomous_missing_weight"] == pytest.approx(
        EXPECTED_AUTONOMOUS_MISSING_WEIGHT, abs=5.0e-15
    )
    assert analytic["debt_to_requested_stress_rms_ratio"] == pytest.approx(
        EXPECTED_AUTONOMOUS_MISSING_WEIGHT, rel=2.0e-15
    )
    debt = analytic["evaluation"]["finite_head_mean_debt"]
    assert debt["requested_stress_rms"] > 0.0
    assert debt["debt_rms"] > 0.0
    assert debt["identity_closure_max_abs"] <= 1.0e-14

    boundary = receipt["truth_boundary"]
    assert boundary["autonomous_factor_is_repository_engineering"] is True
    assert boundary["integration_adapter_only"] is True
    assert boundary["formal_theorem_machine_bump_identity_claimed"] is False
    assert boundary["formal_missing_weight_materialized"] is False
    assert boundary["formal_missing_weight_equality_claimed"] is False
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["full_same_cycle_composite_requested_stress_materialized"] is False
    assert boundary["candidate_finite_head_mean_debt_materialized"] is False
    assert boundary["signed_mean_inverse_input_ready"] is False
    assert boundary["public_velocity_correction_materialized"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
