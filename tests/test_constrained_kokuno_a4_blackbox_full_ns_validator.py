from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_a4_blackbox_full_ns_validator as validator


class SolidRotationCandidate:
    """Analytic incompressible field with exactly balanced pressure.

    u=(-y,x,0), p=(x^2+y^2)/2, f=0 gives
    (u.grad)u=(-x,-y,0) and grad p=(x,y,0), so the full NS residual is zero.
    """

    def __init__(self, contract_sha: str, *, pressure_scale: float = 1.0, residual_defined_forcing: bool = False):
        self.contract_sha = contract_sha
        self.pressure_scale = float(pressure_scale)
        self.residual_defined_forcing = bool(residual_defined_forcing)

    def velocity(self, x, y, z, t):
        xb, yb, zb, tb = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.stack((-yb, xb, np.zeros_like(zb)), axis=-1)

    def pressure(self, x, y, z, t):
        xb, yb, _, _ = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return 0.5 * self.pressure_scale * (xb * xb + yb * yb)

    def forcing(self, x, y, z, t):
        xb, yb, zb, tb = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.zeros(xb.shape + (3,), dtype=float)

    def validation_metadata(self):
        return {
            "candidate_id": "solid-rotation-mechanics",
            "candidate_sha256": "solid-rotation-sha",
            "physical_contract_sha256": self.contract_sha,
            "stage": "leading_only",
            "global_leading_velocity_materialized": True,
            "outer_join_materialized": True,
            "matched_pressure_included": True,
            "restricted_forcing_included": True,
            "restricted_forcing_preregistered": True,
            "residual_defined_forcing": self.residual_defined_forcing,
            "whole_domain_support_materialized": True,
            "leading_included": True,
            "oscillatory_included": False,
            "correction_included": False,
        }


class TimePolynomialCandidate(SolidRotationCandidate):
    def velocity(self, x, y, z, t):
        xb, yb, zb, tb = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.stack((tb**4, np.zeros_like(tb), np.zeros_like(tb)), axis=-1)


class MissingPressureCandidate:
    def velocity(self, x, y, z, t):
        xb, yb, zb, tb = np.broadcast_arrays(x, y, z, t)
        return np.zeros(xb.shape + (3,), dtype=float)

    def forcing(self, x, y, z, t):
        xb, yb, zb, tb = np.broadcast_arrays(x, y, z, t)
        return np.zeros(xb.shape + (3,), dtype=float)

    def validation_metadata(self):
        return {}


def _contract(forcing_parameter_sha: str = "force-params-v1"):
    return validator.build_fixed_physical_contract(
        contract_id="test-fixed-contract",
        matched_pressure_identity_sha256="pressure-identity-v1",
        restricted_forcing_parameters_sha256=forcing_parameter_sha,
        restricted_forcing_preregistration_sha256="forcing-preregistered-v1",
    )


def _sample_points():
    return (
        np.asarray((0.17, -0.31, 0.44, -0.63, 0.21), dtype=float),
        np.asarray((-0.29, 0.52, 0.11, -0.18, -0.47), dtype=float),
        np.asarray((0.09, -0.24, 0.37, 0.28, -0.16), dtype=float),
        np.asarray((0.25, 0.3125, 0.4375, 0.5625, 0.75), dtype=float),
    )


def test_fd4_full_ns_operator_replays_analytic_solid_rotation_at_all_three_levels():
    contract = _contract()
    candidate = SolidRotationCandidate(contract.physical_contract_sha256)
    x, y, z, t = _sample_points()
    levels = validator.three_resolution_levels(candidate, x, y, z, t)
    assert tuple(level.derivative_step for level in levels) == validator.DERIVATIVE_STEPS
    for level in levels:
        np.testing.assert_allclose(level.residual, 0.0, atol=3e-11, rtol=0.0)
        np.testing.assert_allclose(level.divergence, 0.0, atol=3e-12, rtol=0.0)
        assert level.normalized_sampled_max < 5e-11
        assert level.normalized_sampled_rms < 5e-11


def test_pressure_mutation_is_detected_by_complete_residual():
    contract = _contract()
    candidate = SolidRotationCandidate(contract.physical_contract_sha256, pressure_scale=0.99)
    x, y, z, t = _sample_points()
    fine = validator.full_ns_level(
        candidate, x, y, z, t, derivative_step=validator.DERIVATIVE_STEPS[-1]
    )
    # The 1% pressure error leaves R=-0.01*(x,y,0), clearly above 1e-3 here.
    assert fine.normalized_sampled_max > validator.MOMENTUM_GATE
    assert fine.normalized_sampled_rms > validator.MOMENTUM_GATE


def test_time_endpoint_stencil_replays_quartic_derivative_without_out_of_window_sampling():
    contract = _contract()
    candidate = TimePolynomialCandidate(contract.physical_contract_sha256)
    x = np.asarray((0.1, 0.1, 0.1), dtype=float)
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    t = np.asarray((0.25, 0.50, 0.75), dtype=float)
    coords = validator._broadcast_xyzt(x, y, z, t)
    out = validator._fd4_time_first(candidate.velocity, coords, 0.005)
    expected = np.stack((4.0 * t**3, np.zeros_like(t), np.zeros_like(t)), axis=-1)
    np.testing.assert_allclose(out, expected, atol=2e-10, rtol=0.0)


def test_admission_fails_closed_when_candidate_surfaces_or_scientific_flags_are_missing():
    contract = _contract()
    missing = validator.candidate_admission_failures(
        MissingPressureCandidate(), contract, expected_stage="leading_only"
    )
    assert "missing_public_surface:pressure" in missing

    candidate = SolidRotationCandidate(
        contract.physical_contract_sha256, residual_defined_forcing=True
    )
    failures = validator.candidate_admission_failures(
        candidate, contract, expected_stage="leading_only"
    )
    assert "residual_defined_forcing_forbidden" in failures


def test_physical_contract_hash_changes_if_restricted_forcing_parameters_change():
    left = _contract("force-params-left")
    right = _contract("force-params-right")
    assert left.physical_contract_sha256 != right.physical_contract_sha256

    candidate = SolidRotationCandidate(left.physical_contract_sha256)
    failures = validator.candidate_admission_failures(
        candidate, right, expected_stage="leading_only"
    )
    assert "physical_contract_identity_mismatch" in failures


def test_scientific_viscosity_cannot_be_retuned_through_full_ns_level():
    contract = _contract()
    candidate = SolidRotationCandidate(contract.physical_contract_sha256)
    x, y, z, t = _sample_points()
    with pytest.raises(ValueError, match="viscosity is frozen"):
        validator.full_ns_level(
            candidate,
            x,
            y,
            z,
            t,
            derivative_step=validator.DERIVATIVE_STEPS[-1],
            nu=validator.NU * (1.0 + validator.NU_PERTURBATION_FRACTION),
        )


def test_frozen_cloud_and_axis_near_probes_are_deterministic_and_include_axis():
    first = validator.frozen_heldout_cloud()
    second = validator.frozen_heldout_cloud()
    for a, b in zip(first, second):
        np.testing.assert_array_equal(a, b)
    x, y, z, t, weights = first
    assert x.size == validator.HELDOUT_POINTS * len(validator.VALIDATION_TIMES)
    assert np.all(np.isfinite(weights))
    # Each time slice integrates a constant over the [-2,2]^3 box to volume 64.
    for time_value in validator.VALIDATION_TIMES:
        assert np.sum(weights[np.isclose(t, time_value)]) == pytest.approx(64.0)

    ax, ay, az, at = validator.frozen_axis_near_probes()
    radii = np.sqrt(ax * ax + ay * ay)
    assert np.any(radii == 0.0)
    assert np.any((radii > 0.0) & (radii <= 1e-9))
    assert set(np.unique(at)) == set(validator.VALIDATION_TIMES)


def test_current_strict_inner_artifact_is_machine_locked_as_ineligible_for_complete_validation():
    state = validator.current_strict_inner_ineligibility()
    assert state["latest_agent1_pr"] == 893
    assert state["latest_agent2_pr"] == 894
    assert state["latest_agent3_pr"] == 895
    assert state["complete_blackbox_candidate_available"] is False
    assert state["matched_pressure_included"] is False
    assert state["restricted_forcing_included"] is False
    assert state["agent3_correction_velocity_materialized"] is False
    assert state["pde_validated"] is False


def test_public_contract_locks_final_gates_and_quadrature_requirement():
    contract = validator.public_contract()
    assert contract["derivative_steps"] == [0.02, 0.01, 0.005]
    assert contract["frozen_gates"]["momentum_max"] == 1e-3
    assert contract["frozen_gates"]["momentum_volume_l2"] == 1e-3
    assert contract["frozen_gates"]["divergence_max"] == 1e-5
    assert contract["frozen_gates"]["divergence_volume_l2"] == 1e-5
    assert contract["truth_boundary"]["residual_defined_free_forcing_allowed"] is False
    assert contract["truth_boundary"]["quadrature_ladder_required_for_pde_promotion"] is True
    assert contract["truth_boundary"]["pde_validated"] is False
