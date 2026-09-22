from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_exterior_qs_release2_state import (
    KokunoCurrentExteriorQsRelease2State,
)


def _sample_eta() -> np.ndarray:
    return np.asarray([-0.81, -0.53, -0.21, 0.0, 0.27, 0.66], dtype=float)


def test_truth_boundary_keeps_current_state_separate_from_bridge_and_pde() -> None:
    current = KokunoCurrentExteriorQsRelease2State()
    truth = current.truth_boundary
    assert truth["current_q_s_release2_endpoint_materialized"] is True
    assert truth["current_q_s_uses_general_radial_identity"] is True
    assert truth["current_absolute_I_state_consumed"] is True
    assert truth["current_actual_M_M_eta_state_consumed"] is True
    assert truth["current_actual_J_J_eta_state_consumed"] is True
    assert truth["source_ideal_q_s_relabelled_as_current"] is False
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_relative_swirl_angular_state_retains_only_actual_row_residual() -> None:
    current = KokunoCurrentExteriorQsRelease2State()
    eta = _sample_eta()
    target, target_eta = current.split.target.target_with_eta(eta)
    sol, c1_eta, c2_eta = current.split.compensator.solve_target_jet(
        target, target_eta
    )
    row = current.split.compensator.angular_row_scaled
    residual = row[0] * np.asarray(sol.c1) + row[1] * np.asarray(sol.c2) - target
    residual_eta = row[0] * np.asarray(c1_eta) + row[1] * np.asarray(c2_eta) - target_eta
    attenuation = math.exp(
        -(1.0 - current.lambda_value)
        * (current.split.hold_length - current.split.compensator.y1)
    )
    r, r_eta = current._angular_relative_endpoint(eta)
    r_star = 1.0 / (1.0 - current.lambda_value)
    assert np.allclose(r - r_star, attenuation * residual, rtol=2e-12, atol=1e-300)
    assert np.allclose(r_eta, attenuation * residual_eta, rtol=2e-12, atol=1e-300)


def test_current_pulse_J_reconstructs_exact_1133_public_scaled_row() -> None:
    current = KokunoCurrentExteriorQsRelease2State()
    eta = _sample_eta()
    _, _, j_scaled = current._pulse_j_normalized_with_eta(eta)
    public = current.pulse_end.similarity_profile_values_logX(
        np.full_like(eta, current.pulse_end.log_X_pulse_end), eta
    )["public_J_row_scaled"]
    public = np.asarray(public, dtype=float)
    scale = np.maximum.reduce(
        [np.abs(j_scaled), np.abs(public), np.full_like(j_scaled, 1e-300)]
    )
    assert np.max(np.abs(j_scaled - public) / scale) <= 2e-11


def test_release2_current_Qs_is_exact_general_radial_identity_decomposition() -> None:
    current = KokunoCurrentExteriorQsRelease2State()
    state = current.state(_sample_eta())
    rebuilt = (
        state["term_minus_W"]
        + state["term_I"]
        + state["term_I_eta"]
        + state["term_J_eta"]
        + state["term_J"]
    )
    q = state["Q_s_current_release2_end"]
    assert np.array_equal(np.asarray(q), np.asarray(rebuilt))
    assert np.all(np.isfinite(q))
    assert np.allclose(
        q - state["Q_s_source_ideal_release2_end"],
        state["Q_s_current_minus_source_ideal"],
        rtol=0.0,
        atol=4e-15 * max(1.0, float(np.max(np.abs(q)))),
    )


def test_W_is_built_from_actual_1204_carried_M_state() -> None:
    current = KokunoCurrentExteriorQsRelease2State()
    eta = _sample_eta()
    state = current.state(eta)
    d = 1.0 - eta * eta
    expected = (
        1.0
        - 2.0 * current.D * eta * state["M_over_X"]
        - d * state["M_eta_over_X"]
    )
    assert np.array_equal(np.asarray(state["W"]), np.asarray(expected))


def test_angular_release2_transport_uses_actual_XH_ratio_not_ideal_reset() -> None:
    current = KokunoCurrentExteriorQsRelease2State()
    eta = _sample_eta()
    r_entry, r_entry_eta = current._angular_relative_endpoint(eta)
    r_end, r_end_eta = current._angular_release2_endpoint(eta)

    entry = current.split.split_profile_logX(
        np.full_like(eta, current.split.log_X_hold_end), eta
    )
    end = current.release2.profile_logX(
        np.full_like(eta, current.release2.log_X_release2_end), eta
    )
    log_xh_entry = (
        math.log(2.0)
        + 2.0 * current.split.log_X_hold_end
        + np.asarray(entry["log_F_base"], dtype=float)
    )
    log_xh_end = (
        math.log(2.0)
        + 2.0 * current.release2.log_X_release2_end
        + np.asarray(end["log_F"], dtype=float)
    )
    factor = np.exp(log_xh_entry - log_xh_end)
    r_star = 1.0 / (1.0 - current.lambda_value)
    r_ideal_end = (current.parent.q_in_source_ideal + 1.0) / (1.0 - current.h_value)
    assert np.allclose(
        r_end - r_ideal_end,
        (r_entry - r_star) * factor,
        rtol=2e-12,
        atol=1e-300,
    )
    assert np.allclose(r_end_eta, r_entry_eta * factor, rtol=2e-12, atol=1e-300)


def test_report_is_finite_but_does_not_materialize_matching_bridge() -> None:
    current = KokunoCurrentExteriorQsRelease2State()
    report = current.report()
    q = np.asarray(report["q_s_current_release2_end"], dtype=float)
    delta = np.asarray(report["q_s_current_minus_source_ideal"], dtype=float)
    assert np.all(np.isfinite(q))
    assert np.all(np.isfinite(delta))
    assert report["current_l_minus_h_matching_bridge_materialized"] is False
    assert math.isfinite(report["max_abs_current_minus_source_ideal"])
    assert report["q_p"] > 0.0


def test_save_load_replays_semantic_identity_and_state(tmp_path) -> None:
    current = KokunoCurrentExteriorQsRelease2State()
    path = current.save_configuration(tmp_path / "current_qs.json")
    loaded = KokunoCurrentExteriorQsRelease2State.load_configuration(path)
    assert loaded.semantic_sha256 == current.semantic_sha256
    eta = np.asarray([-0.47, 0.13, 0.59])
    before = current.state(eta)["Q_s_current_release2_end"]
    after = loaded.state(eta)["Q_s_current_release2_end"]
    assert np.array_equal(np.asarray(before), np.asarray(after))


def test_configuration_truth_mutation_fails_closed(tmp_path) -> None:
    current = KokunoCurrentExteriorQsRelease2State()
    path = current.save_configuration(tmp_path / "current_qs.json")
    import json

    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["truth_boundary"]["pde_validated"] = True
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="configuration/provenance mismatch"):
        KokunoCurrentExteriorQsRelease2State.load_configuration(path)


def test_public_state_api_exposes_no_scientific_tuning_knobs() -> None:
    sig = inspect.signature(KokunoCurrentExteriorQsRelease2State.state)
    assert tuple(sig.parameters) == ("self", "eta")
    forbidden = {
        "residual",
        "forcing",
        "pressure",
        "optimizer",
        "tolerance",
        "threshold",
        "gain",
        "h",
        "lambda_value",
        "precision",
    }
    assert forbidden.isdisjoint(sig.parameters)
    assert not hasattr(KokunoCurrentExteriorQsRelease2State, "velocity")
