from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_current_joined_profile_independent_audit import (
    DERIVATIVE_MAX_GATE,
    DERIVATIVE_RMS_GATE,
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_GATE,
    LOG_STEPS,
    SEED,
    XI_JET_GATE,
    XH_E_SLOPE_GATE,
    XH_F_SLOPE_GATE,
    XH_U_SLOPE_GATE,
    audit_rebound_public_profile,
    enforce_current_joined_profile_independent_audit,
    protocol,
)
from openai_ns_reconstruction.kokuno_pa16_current_joined_profile import (
    KokunoPA16CurrentJoinedProfile,
)


class _ManufacturedJoinedProfile:
    """Smooth public surface with the same endpoint dimensionless slopes."""

    X_i = 1.0
    log_X_R = 12.0
    X_h = math.exp(7.0)  # X_R * exp(-5)
    route_ready = True
    certificate = SimpleNamespace(selected_T_sh=2.0)

    def profile_values(self, X, eta):
        X = np.asarray(X, dtype=float)
        eta = np.asarray(eta, dtype=float)
        s = np.log(X)
        F = X ** (-0.4)
        E = math.sqrt(2.0) * X**0.1
        # Nonzero interior derivative but exact ideal U and zero slope at X_h.
        U = 4.0 * eta + 0.02 * (s - math.log(self.X_h)) ** 2
        H = np.sqrt(2.0 * X) * E
        return {
            "F_current_joined": F,
            "U_current_joined": U,
            "E_current_joined": E,
            "H_current_joined": H,
        }

    def radial_derivatives(self, X, eta):
        X = np.asarray(X, dtype=float)
        s = np.log(X)
        F = X ** (-0.4)
        E = math.sqrt(2.0) * X**0.1
        return {
            "F_current_joined_X": -0.4 * F / X,
            "U_current_joined_X": 0.04 * (s - math.log(self.X_h)) / X,
            "E_current_joined_X": 0.1 * E / X,
        }


class _Wrapper:
    def __init__(self, base):
        self.base = base

    def __getattr__(self, name):
        return getattr(self.base, name)


class _ProductionFDerivativeDrift(_Wrapper):
    def radial_derivatives(self, X, eta):
        out = dict(self.base.radial_derivatives(X, eta))
        out["F_current_joined_X"] = 0.99 * np.asarray(out["F_current_joined_X"])
        return out


class _EndpointValuePreservingETangent(_Wrapper):
    def profile_values(self, X, eta):
        out = dict(self.base.profile_values(X, eta))
        X_arr = np.asarray(X, dtype=float)
        factor = 1.0 + 1.0e-3 * np.log(X_arr / float(self.X_h))
        E = np.asarray(out["E_current_joined"], dtype=float) * factor
        out["E_current_joined"] = E
        # Preserve the public H=sqrt(2X)E identity so the derivative audit,
        # rather than a trivial algebraic inconsistency, must catch mutation.
        out["H_current_joined"] = np.sqrt(2.0 * X_arr) * E
        return out


class _BlockedRoute(_ManufacturedJoinedProfile):
    route_ready = False


def test_frozen_protocol_and_final_truth_boundary():
    p = protocol()
    assert p["seed"] == SEED == 9173651
    assert tuple(p["log_step_ladder"]) == LOG_STEPS == (4e-3, 2e-3, 1e-3)
    assert p["derivative_relative_rms_gate"] == DERIVATIVE_RMS_GATE == 5e-3
    assert p["derivative_relative_max_gate"] == DERIVATIVE_MAX_GATE == 2e-2
    assert p["Xi_left_right_jet_gate"] == XI_JET_GATE == 2e-2
    assert p["Xh_E_slope_gate"] == XH_E_SLOPE_GATE == 2e-4
    assert p["Xh_F_slope_gate"] == XH_F_SLOPE_GATE == 2e-4
    assert p["Xh_U_slope_gate"] == XH_U_SLOPE_GATE == 5e-3
    assert p["final_momentum_gate_unassessed"] == FINAL_MOMENTUM_GATE == 1e-3
    assert p["final_divergence_gate_unassessed"] == FINAL_DIVERGENCE_GATE == 1e-5
    assert p["residual_defined_forcing_forbidden"] is True


def test_manufactured_public_surface_passes_independent_logX_audit():
    receipt = audit_rebound_public_profile(_ManufacturedJoinedProfile())
    assert receipt["audit_pass"] is True, receipt
    assert receipt["status"] == "pass"
    assert receipt["gates"]["finest_derivative_pass"] is True
    assert receipt["gates"]["endpoint_pass"] is True
    assert receipt["truth_boundary"]["pde_validated"] is False
    assert receipt["truth_boundary"]["global_cartesian_leading_velocity_materialized"] is False


def test_production_derivative_drift_is_rejected():
    receipt = audit_rebound_public_profile(
        _ProductionFDerivativeDrift(_ManufacturedJoinedProfile())
    )
    assert receipt["audit_pass"] is False
    fine = receipt["derivative_levels"][-1]["metrics"]["F_X"]
    assert fine["relative_rms"] > DERIVATIVE_RMS_GATE


def test_endpoint_value_preserving_public_tangent_is_rejected():
    base = _ManufacturedJoinedProfile()
    mutated = _EndpointValuePreservingETangent(base)
    # Mutation is exactly value-preserving at X_h.
    eta = np.asarray([-0.4, 0.0, 0.4])
    base_exit = base.profile_values(np.full(3, base.X_h), eta)["E_current_joined"]
    mut_exit = mutated.profile_values(np.full(3, base.X_h), eta)["E_current_joined"]
    np.testing.assert_array_equal(mut_exit, base_exit)
    receipt = audit_rebound_public_profile(mutated)
    assert receipt["audit_pass"] is False
    assert (
        receipt["endpoints"]["X_h"]["E_dimensionless_slope_max_abs_error"]
        > XH_E_SLOPE_GATE
    )


def test_infeasible_current_tsh_route_fails_closed_without_outer_values():
    receipt = audit_rebound_public_profile(_BlockedRoute())
    assert receipt["status"] == "blocked_by_current_tsh_geometry"
    assert receipt["audit_pass"] is False
    assert receipt["route_ready"] is False
    assert "derivative_levels" not in receipt


def test_enforcer_rejects_failed_or_pde_promoted_receipt():
    good = audit_rebound_public_profile(_ManufacturedJoinedProfile())
    enforce_current_joined_profile_independent_audit(good)

    failed = dict(good)
    failed["audit_pass"] = False
    failed["status"] = "fail"
    with pytest.raises(RuntimeError, match="did not pass"):
        enforce_current_joined_profile_independent_audit(failed)

    promoted = dict(good)
    promoted["truth_boundary"] = dict(good["truth_boundary"])
    promoted["truth_boundary"]["pde_validated"] = True
    with pytest.raises(RuntimeError, match="cannot promote"):
        enforce_current_joined_profile_independent_audit(promoted)


def test_real_a1_configuration_roundtrip_preserves_semantic_identity(tmp_path):
    field = KokunoPA16CurrentJoinedProfile()
    path = tmp_path / "joined.json"
    field.save_configuration(path)
    rebound = KokunoPA16CurrentJoinedProfile.load_configuration(path)
    assert rebound.semantic_sha256 == field.semantic_sha256
    assert rebound.truth_boundary["cartesian_velocity_materialized"] is False
    assert rebound.truth_boundary["matched_global_pressure_materialized"] is False
    assert rebound.truth_boundary["restricted_forcing_materialized"] is False
    assert rebound.truth_boundary["pde_validated"] is False
