from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_current_qs_release2_independent_audit import (
    CURRENT_MINUS_IDEAL_MIN,
    DERIVATIVE_CHANNEL_NORM_MAX_GATE,
    ETA_FD_STEPS,
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_GATE,
    FINAL_QUADRATURE,
    Q_RECON_NORM_MAX_GATE,
    SEED,
    _fd4,
    audit_loaded_current_qs_state,
    default_candidate,
    enforce_preregistered_gates,
    heldout_eta,
    public_api_has_no_scientific_tuning_knobs,
)
from openai_ns_reconstruction.kokuno_current_exterior_qs_release2_state import (
    KokunoCurrentExteriorQsRelease2State,
)


def _audit_roundtrip(tmp_path):
    candidate = default_candidate()
    path = candidate.save_configuration(tmp_path / "candidate.json")
    loaded = KokunoCurrentExteriorQsRelease2State.load_configuration(path)
    return audit_loaded_current_qs_state(loaded, candidate)


def test_protocol_is_frozen_before_results() -> None:
    assert SEED == 9173991
    assert ETA_FD_STEPS == (2.0**-8, 2.0**-9, 2.0**-10)
    assert Q_RECON_NORM_MAX_GATE == 1.0e-7
    assert DERIVATIVE_CHANNEL_NORM_MAX_GATE == 2.0e-6
    assert CURRENT_MINUS_IDEAL_MIN == 1.0e-12
    assert FINAL_MOMENTUM_GATE == 1.0e-3
    assert FINAL_DIVERGENCE_GATE == 1.0e-5
    assert FINAL_QUADRATURE == (24, 48, 96)


def test_heldout_eta_is_deterministic_offgrid_and_stencil_safe() -> None:
    candidate = default_candidate()
    a = heldout_eta(candidate)
    b = heldout_eta(candidate)
    assert np.array_equal(a, b)
    assert a.shape == (33,)
    assert len(np.unique(a)) == len(a)
    lo, hi = candidate.eta_interval
    h = max(ETA_FD_STEPS)
    assert np.all(a - 2.0 * h > lo)
    assert np.all(a + 2.0 * h < hi)
    # The seeded random subset should not collapse onto a regular grid.
    scaled = (a[9:] - float(lo)) / (float(hi) - float(lo))
    assert np.any(np.abs(scaled * 1024.0 - np.rint(scaled * 1024.0)) > 1.0e-8)


def test_a4_fd4_is_exact_on_quartic_manufactured_channel() -> None:
    x = np.asarray([-0.71, -0.31, 0.17, 0.63], dtype=float)
    h = 2.0**-9
    got = _fd4(lambda y: y**4 - 2.0 * y**3 + 0.7 * y, x, h)
    expected = 4.0 * x**3 - 6.0 * x**2 + 0.7
    assert np.allclose(got, expected, rtol=2.0e-10, atol=2.0e-11)


def test_roundtrip_audit_preserves_truth_boundary_and_three_resolutions(tmp_path) -> None:
    report = _audit_roundtrip(tmp_path)
    assert report["task"] == "K4-VAL-128"
    assert report["eta_fd_steps"] == list(ETA_FD_STEPS)
    assert len(report["resolutions"]) == 3
    assert report["save_load"] == {
        "semantic_exact": True,
        "configuration_exact": True,
        "q_state_exact": True,
    }
    truth = report["truth_boundary"]
    assert truth["current_q_s_release2_independent_audit_executed"] is True
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["cartesian_terminal_multiplier_composed"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["complete_ns_admission_ready"] is False
    assert truth["pde_validated"] is False
    assert report["final_project_admission_ready"] is False
    bridge = report["bridge_diagnostic"]
    assert bridge["diagnostic_only_not_a_bridge_construction"] is True
    assert bridge["eligible_probe_count"] + bridge["ineligible_probe_count"] == 33


def test_preregistered_scoped_gates_are_enforced_on_real_roundtrip(tmp_path) -> None:
    report = _audit_roundtrip(tmp_path)
    enforce_preregistered_gates(report)


def test_public_derivative_channel_corruption_is_detected(monkeypatch, tmp_path) -> None:
    original = KokunoCurrentExteriorQsRelease2State.state

    def corrupted(self, eta):
        out = dict(original(self, eta))
        out["I_eta_over_XH"] = np.asarray(out["I_eta_over_XH"], dtype=float) + 1.0e-3
        return out

    monkeypatch.setattr(KokunoCurrentExteriorQsRelease2State, "state", corrupted)
    report = _audit_roundtrip(tmp_path)
    fine = report["resolutions"][-1]
    assert fine["derivative_channel_normalized_max"] > DERIVATIVE_CHANNEL_NORM_MAX_GATE
    with pytest.raises(AssertionError, match="derivative-channel"):
        enforce_preregistered_gates(report)


def test_public_q_corruption_is_detected_by_independent_reconstruction(monkeypatch, tmp_path) -> None:
    original = KokunoCurrentExteriorQsRelease2State.state

    def corrupted(self, eta):
        out = dict(original(self, eta))
        out["Q_s_current_release2_end"] = (
            np.asarray(out["Q_s_current_release2_end"], dtype=float) + 1.0e-3
        )
        return out

    monkeypatch.setattr(KokunoCurrentExteriorQsRelease2State, "state", corrupted)
    report = _audit_roundtrip(tmp_path)
    fine = report["resolutions"][-1]
    assert fine["q_reconstruction_normalized_max"] > Q_RECON_NORM_MAX_GATE
    with pytest.raises(AssertionError, match="Q_s normalized max"):
        enforce_preregistered_gates(report)


def test_truth_promotion_fails_closed(tmp_path) -> None:
    report = _audit_roundtrip(tmp_path)
    report["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError, match="PDE truth promotion"):
        enforce_preregistered_gates(report)


def test_bridge_spread_is_diagnostic_and_not_a_public_tuning_knob() -> None:
    sig = inspect.signature(audit_loaded_current_qs_state)
    assert tuple(sig.parameters) == ("loaded", "pre_serialization_reference")
    assert public_api_has_no_scientific_tuning_knobs() is True
