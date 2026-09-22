from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_current_lminus_h_bridge_independent_audit import (
    FD4_Y_STEPS,
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_GATE,
    FINAL_QUADRATURE,
    INHOMOGENEOUS_P_MIN,
    ODE_DEFECT_NORM_MAX_GATE,
    ODE_DEFECT_NORM_RMS_GATE,
    RK4_STEPS_PER_UNIT,
    SEED,
    TRANSPORT_NORM_MAX_GATE,
    TRANSPORT_NORM_RMS_GATE,
    _rk4_at_targets,
    audit_loaded_current_bridge,
    default_candidate,
    enforce_preregistered_gates,
    heldout_eta,
    heldout_y,
    public_api_has_no_scientific_tuning_knobs,
)
from openai_ns_reconstruction.kokuno_current_exterior_lminus_h_bridge import (
    KokunoCurrentExteriorLMinusHBridge,
)


def _audit_roundtrip(tmp_path):
    candidate = default_candidate()
    path = candidate.save_configuration(tmp_path / "candidate.json")
    loaded = KokunoCurrentExteriorLMinusHBridge.load_configuration(path)
    return audit_loaded_current_bridge(loaded, candidate)


def test_protocol_is_frozen_before_results() -> None:
    assert SEED == 9174001
    assert RK4_STEPS_PER_UNIT == (64, 128, 256)
    assert FD4_Y_STEPS == (2.0**-7, 2.0**-8, 2.0**-9)
    assert TRANSPORT_NORM_MAX_GATE == 2.0e-8
    assert TRANSPORT_NORM_RMS_GATE == 2.0e-8
    assert ODE_DEFECT_NORM_MAX_GATE == 2.0e-8
    assert ODE_DEFECT_NORM_RMS_GATE == 2.0e-8
    assert INHOMOGENEOUS_P_MIN == 1.0e-12
    assert FINAL_MOMENTUM_GATE == 1.0e-3
    assert FINAL_DIVERGENCE_GATE == 1.0e-5
    assert FINAL_QUADRATURE == (24, 48, 96)


def test_heldout_protocol_is_deterministic_offgrid_and_stencil_safe() -> None:
    candidate = default_candidate()
    eta_a = heldout_eta(candidate)
    eta_b = heldout_eta(candidate)
    y_a = heldout_y()
    y_b = heldout_y()
    assert np.array_equal(eta_a, eta_b)
    assert np.array_equal(y_a, y_b)
    assert eta_a.shape == (33,)
    assert y_a.shape == (24,)
    assert len(np.unique(eta_a)) == 33
    assert len(np.unique(y_a)) == 24
    assert np.min(y_a) > 2.0 * max(FD4_Y_STEPS)
    lo, hi = candidate.eta_interval
    assert np.all((eta_a > lo) & (eta_a < hi))
    scaled = (eta_a[9:] - float(lo)) / (float(hi) - float(lo))
    assert np.any(np.abs(scaled * 2048.0 - np.rint(scaled * 2048.0)) > 1.0e-8)


def test_independent_rk4_replays_manufactured_closed_form() -> None:
    q0 = np.asarray([1.3, 0.7, 2.1], dtype=float)
    p0 = np.asarray([0.2, -0.1, 0.035], dtype=float)
    y = np.asarray([0.2, 1.1, 3.0, 6.5], dtype=float)
    h = 0.005
    got = _rk4_at_targets(q0, p0, y, 256, h)
    yy = y[None, :]
    expected = np.exp(-(1.0 - h) * yy) * (
        q0[:, None] - p0[:, None] * (1.0 - np.exp(-h * yy))
    )
    assert np.allclose(got, expected, rtol=4.0e-11, atol=2.0e-13)


def test_roundtrip_audit_preserves_truth_boundary_and_three_resolutions(tmp_path) -> None:
    report = _audit_roundtrip(tmp_path)
    assert report["task"] == "K4-VAL-129"
    assert report["rk4_steps_per_unit"] == list(RK4_STEPS_PER_UNIT)
    assert report["fd4_y_steps"] == list(FD4_Y_STEPS)
    assert len(report["resolutions"]) == 3
    assert report["save_load"] == {
        "semantic_exact": True,
        "configuration_exact": True,
        "q_state_exact": True,
    }
    truth = report["truth_boundary"]
    assert truth["current_l_minus_h_bridge_independent_audit_executed"] is True
    assert truth["full_eta_interval_common_scalar_bridge_length_established"] is False
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["eta_dependent_cartesian_matching_boundary_materialized"] is False
    assert truth["cartesian_terminal_multiplier_composed"] is False
    assert truth["complete_ns_admission_ready"] is False
    assert truth["pde_validated"] is False
    assert report["final_project_admission_ready"] is False
    target = report["target_time_diagnostic"]
    assert target["diagnostic_only_not_a_scalar_bridge_construction"] is True
    assert target["independent_rk4_target_time_spread"] >= 0.0


def test_preregistered_scoped_gates_are_enforced_on_real_roundtrip(tmp_path) -> None:
    report = _audit_roundtrip(tmp_path)
    enforce_preregistered_gates(report)


def test_production_analytic_derivative_corruption_is_diagnostic_only(
    monkeypatch, tmp_path
) -> None:
    original = KokunoCurrentExteriorLMinusHBridge.state

    def corrupted(self, y, eta):
        out = dict(original(self, y, eta))
        out["D_y_Q_s"] = np.asarray(out["D_y_Q_s"], dtype=float) + 1.0e-3
        return out

    monkeypatch.setattr(KokunoCurrentExteriorLMinusHBridge, "state", corrupted)
    report = _audit_roundtrip(tmp_path)
    fine = report["resolutions"][-1]
    assert fine["analytic_derivative_diagnostic_normalized_max"] > 1.0e-5
    enforce_preregistered_gates(report)


def test_public_q_corruption_is_detected(monkeypatch, tmp_path) -> None:
    original = KokunoCurrentExteriorLMinusHBridge.state

    def corrupted(self, y, eta):
        out = dict(original(self, y, eta))
        out["Q_s"] = np.asarray(out["Q_s"], dtype=float) + 1.0e-3
        return out

    monkeypatch.setattr(KokunoCurrentExteriorLMinusHBridge, "state", corrupted)
    report = _audit_roundtrip(tmp_path)
    fine = report["resolutions"][-1]
    assert (
        fine["transport_normalized_max"] > TRANSPORT_NORM_MAX_GATE
        or fine["ode_defect_normalized_max"] > ODE_DEFECT_NORM_MAX_GATE
    )
    with pytest.raises(AssertionError):
        enforce_preregistered_gates(report)


def test_truth_promotion_fails_closed(tmp_path) -> None:
    report = _audit_roundtrip(tmp_path)
    report["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError, match="PDE truth promotion"):
        enforce_preregistered_gates(report)


def test_public_audit_api_exposes_no_scientific_tuning_knobs() -> None:
    sig = inspect.signature(audit_loaded_current_bridge)
    assert tuple(sig.parameters) == ("loaded", "pre_serialization_reference")
    assert public_api_has_no_scientific_tuning_knobs() is True
