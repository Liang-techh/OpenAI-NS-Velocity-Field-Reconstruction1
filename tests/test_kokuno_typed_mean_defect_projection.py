from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)
from openai_ns_reconstruction.kokuno_typed_mean_defect_projection import (
    CylindricalRing,
    deterministic_receipt,
    evaluate_actual_mean_defect_ring,
    evaluate_disjoint_actual_mean_defect,
    truth_boundary,
)


def _providers(identity: CycleIdentity):
    velocity = VectorFieldProvider(
        identity,
        "analytic:u=t*(-y,x,0)",
        lambda x, y, z, t: (-t * y, t * x, 0.0),
    )
    velocity_dt = VectorFieldProvider(
        identity,
        "analytic:u_t=(-y,x,0)",
        lambda x, y, z, t: (-y, x, 0.0),
    )
    pressure = ScalarFieldProvider(
        identity, "analytic:p=0", lambda x, y, z, t: 0.0
    )
    forcing = RestrictedForcingProvider(
        identity,
        "analytic:f=0",
        "zero forcing fixed independently of residual",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    return velocity, velocity_dt, pressure, forcing


def test_actual_defect_projects_to_exact_rotating_theta_mean() -> None:
    identity = CycleIdentity("projection-test", 0, "fixed")
    providers = _providers(identity)
    ring = CylindricalRing(time=0.0, radius=0.37, axial_z=-0.21, angular_count=32)
    result = evaluate_actual_mean_defect_ring(
        *providers, ring, viscosity=0.01, spatial_step=0.005
    )

    assert result.identity == identity
    assert result.theta_mean == pytest.approx(0.37, rel=0.0, abs=2.0e-12)
    assert result.radial_mean == pytest.approx(0.0, rel=0.0, abs=2.0e-12)
    assert result.axial_mean == pytest.approx(0.0, rel=0.0, abs=2.0e-12)
    assert result.theta_rms == pytest.approx(0.37, rel=0.0, abs=2.0e-12)
    assert result.raw_vector_rms == pytest.approx(0.37, rel=0.0, abs=2.0e-12)


def test_heldin_heldout_ring_points_are_disjoint_and_not_retuned() -> None:
    identity = CycleIdentity("projection-split", 3, "same-cycle")
    providers = _providers(identity)
    result = evaluate_disjoint_actual_mean_defect(
        *providers,
        held_in_rings=(
            CylindricalRing(0.0, 0.20, -0.17, 32),
            CylindricalRing(0.0, 0.34, 0.09, 32),
        ),
        held_out_rings=(
            CylindricalRing(0.0, 0.27, 0.23, 32, 0.031),
            CylindricalRing(0.0, 0.41, -0.29, 32, 0.047),
        ),
        viscosity=0.01,
        spatial_step=0.005,
    )

    assert result.identity == identity
    assert [item.theta_mean for item in result.held_in] == pytest.approx([0.20, 0.34])
    assert [item.theta_mean for item in result.held_out] == pytest.approx([0.27, 0.41])
    assert all(abs(item.radial_mean) < 2.0e-12 for item in result.held_in + result.held_out)
    assert all(abs(item.axial_mean) < 2.0e-12 for item in result.held_in + result.held_out)


def test_reusing_validation_ring_is_rejected() -> None:
    identity = CycleIdentity("projection-overlap", 0, "fixed")
    providers = _providers(identity)
    ring = CylindricalRing(0.0, 0.31, 0.12, 16)
    with pytest.raises(ValueError, match="held-in and held-out"):
        evaluate_disjoint_actual_mean_defect(
            *providers,
            held_in_rings=(ring,),
            held_out_rings=(ring,),
            viscosity=0.01,
            spatial_step=0.005,
        )


def test_cycle_identity_mismatch_is_rejected_by_source_defect_contract() -> None:
    identity = CycleIdentity("projection-identity", 0, "fixed")
    velocity, velocity_dt, pressure, forcing = _providers(identity)
    other = CycleIdentity("projection-identity", 1, "mutated")
    mismatched_pressure = ScalarFieldProvider(
        other, "analytic:p=0/mutated", lambda x, y, z, t: 0.0
    )
    ring = CylindricalRing(0.0, 0.29, -0.08, 16)
    with pytest.raises(ValueError, match="same-cycle identity"):
        evaluate_actual_mean_defect_ring(
            velocity,
            velocity_dt,
            mismatched_pressure,
            forcing,
            ring,
            viscosity=0.01,
            spatial_step=0.005,
        )


def test_public_projection_api_cannot_accept_precomputed_success_inputs() -> None:
    parameters = inspect.signature(evaluate_disjoint_actual_mean_defect).parameters
    forbidden = {
        "residual",
        "defect",
        "stress",
        "target",
        "gain",
        "normalized_score",
    }
    assert forbidden.isdisjoint(parameters)


def test_deterministic_receipt_records_real_nonzero_defect_and_truth_boundary() -> None:
    receipt = deterministic_receipt()
    regression = receipt["analytic_regression"]
    boundary = receipt["truth_boundary"]

    assert regression["maximum_theta_mean_error"] < 2.0e-12
    assert regression["maximum_radial_mean_error"] < 2.0e-12
    assert regression["maximum_axial_mean_error"] < 2.0e-12
    split = regression["split"]
    observed = [
        item["theta_mean"] for item in split["held_in"] + split["held_out"]
    ]
    assert np.asarray(observed) == pytest.approx([0.20, 0.34, 0.27, 0.41])
    assert boundary["raw_same_cycle_providers_required"] is True
    assert boundary["caller_supplied_precomputed_defect_allowed"] is False
    assert boundary["mean_channels_materialized_from_actual_defect"] is True
    assert boundary["ready_as_typed_input_for_radial_stress_reconstruction"] is True
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["pde_validated"] is False
