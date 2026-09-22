from __future__ import annotations

import inspect
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_exterior_lminus_h_bridge import (
    KokunoCurrentExteriorLMinusHBridge,
)


def _eta() -> np.ndarray:
    return np.asarray([-0.77, -0.42, -0.13, 0.0, 0.31, 0.69], dtype=float)


def test_truth_boundary_materializes_transport_not_global_matching_bridge() -> None:
    bridge = KokunoCurrentExteriorLMinusHBridge()
    truth = bridge.truth_boundary
    assert truth["current_q_s_release2_endpoint_materialized"] is True
    assert truth["current_l_minus_h_q_s_transport_materialized"] is True
    assert truth["current_eta_resolved_pointwise_target_time_diagnostic_materialized"] is True
    assert truth["full_eta_interval_common_scalar_bridge_length_established"] is False
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["eta_dependent_cartesian_matching_boundary_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_bridge_entry_exactly_replays_1225_current_state() -> None:
    bridge = KokunoCurrentExteriorLMinusHBridge()
    eta = _eta()
    parent = bridge.parent.state(eta)
    state = bridge.state(0.0, eta)
    assert np.array_equal(state["Q_s"], parent["Q_s_current_release2_end"])
    assert np.array_equal(state["M_over_X"], parent["M_over_X"])
    assert np.array_equal(state["M_eta_over_X"], parent["M_eta_over_X"])
    assert np.allclose(state["W"], parent["W"], rtol=0.0, atol=2e-15)


def test_closed_form_replays_general_inhomogeneous_bridge_formula() -> None:
    bridge = KokunoCurrentExteriorLMinusHBridge()
    eta = _eta()
    y = np.asarray([0.1, 0.7, 2.0, 5.0, 8.0, 11.0])
    entry = bridge.entry_state(eta)
    state = bridge.state(y, eta)
    q_expected = np.exp(-(1.0 - bridge.h_value) * y) * (
        entry["Q_s_entry"] - entry["P_entry"] * (1.0 - np.exp(-bridge.h_value * y))
    )
    assert np.allclose(state["Q_s"], q_expected, rtol=3e-14, atol=1e-300)
    assert np.allclose(state["M_over_X"], entry["M_over_X_entry"] * np.exp(-y), rtol=2e-14)
    assert np.allclose(state["M_eta_over_X"], entry["M_eta_over_X_entry"] * np.exp(-y), rtol=2e-14)


def test_analytic_derivative_satisfies_public_specialized_Qs_ode() -> None:
    bridge = KokunoCurrentExteriorLMinusHBridge()
    eta = _eta()
    y = np.asarray([0.23, 0.91, 2.7, 4.4, 7.2, 9.8])
    state = bridge.state(y, eta)
    lhs = state["D_y_Q_s"] + (1.0 - bridge.h_value) * state["Q_s"]
    rhs = bridge.h_value * (state["W"] - 1.0)
    assert np.allclose(lhs, rhs, rtol=2e-13, atol=2e-14 * max(1.0, float(np.max(np.abs(rhs)))))

    step = 2.0e-5
    qp = bridge.state(y + step, eta)["Q_s"]
    qm = bridge.state(y - step, eta)["Q_s"]
    fd = (qp - qm) / (2.0 * step)
    scale = np.maximum(1.0, np.abs(state["D_y_Q_s"]))
    assert np.max(np.abs(fd - state["D_y_Q_s"]) / scale) < 2.0e-8


def test_pointwise_target_time_hits_existing_terminal_Qp_without_promoting_scalar_bridge() -> None:
    bridge = KokunoCurrentExteriorLMinusHBridge()
    eta = _eta()
    times = bridge.pointwise_target_time(eta)
    assert np.all(np.isfinite(times))
    assert np.all(times > 0.0)
    hit = bridge.state(times, eta)["Q_s"]
    scale = max(1.0, abs(bridge.q_p))
    assert np.max(np.abs(hit - bridge.q_p)) <= 2.0e-13 * scale
    report = bridge.pointwise_target_report(eta)
    assert report["finite_probe_only"] is True
    assert report["full_eta_interval_common_scalar_bridge_length_established"] is False
    assert report["current_l_minus_h_matching_bridge_materialized"] is False
    assert report["eta_dependent_cartesian_matching_boundary_materialized"] is False
    assert math.isfinite(report["pointwise_target_time_spread"])


def test_source_ideal_reference_time_is_reference_only() -> None:
    bridge = KokunoCurrentExteriorLMinusHBridge()
    report = bridge.pointwise_target_report(np.asarray([-0.5, 0.0, 0.5]))
    expected = math.log(
        bridge.parent.parent.q_in_source_ideal / bridge.q_p
    ) / (1.0 - bridge.h_value)
    assert report["source_ideal_reference_time"] == pytest.approx(expected, rel=0.0, abs=2e-15)
    assert report["truth_boundary"]["source_ideal_q_s_relabelled_as_current"] is False


def test_save_load_replays_semantic_identity_and_transport(tmp_path) -> None:
    bridge = KokunoCurrentExteriorLMinusHBridge()
    path = bridge.save_configuration(tmp_path / "bridge.json")
    loaded = KokunoCurrentExteriorLMinusHBridge.load_configuration(path)
    assert loaded.semantic_sha256 == bridge.semantic_sha256
    eta = np.asarray([-0.47, 0.13, 0.59])
    y = np.asarray([0.3, 3.0, 9.0])
    assert np.array_equal(loaded.state(y, eta)["Q_s"], bridge.state(y, eta)["Q_s"])


def test_configuration_truth_mutation_fails_closed(tmp_path) -> None:
    bridge = KokunoCurrentExteriorLMinusHBridge()
    path = bridge.save_configuration(tmp_path / "bridge.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["truth_boundary"]["current_l_minus_h_matching_bridge_materialized"] = True
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="configuration/provenance mismatch"):
        KokunoCurrentExteriorLMinusHBridge.load_configuration(path)


def test_public_apis_expose_no_scientific_or_bridge_selection_tuning_knobs() -> None:
    state_sig = inspect.signature(KokunoCurrentExteriorLMinusHBridge.state)
    root_sig = inspect.signature(KokunoCurrentExteriorLMinusHBridge.pointwise_target_time)
    assert tuple(state_sig.parameters) == ("self", "y", "eta")
    assert tuple(root_sig.parameters) == ("self", "eta")
    forbidden = {
        "residual",
        "forcing",
        "pressure",
        "optimizer",
        "tolerance",
        "threshold",
        "gain",
        "h",
        "precision",
        "iterations",
        "bridge_length",
        "eta_representative",
    }
    assert forbidden.isdisjoint(state_sig.parameters)
    assert forbidden.isdisjoint(root_sig.parameters)
    assert not hasattr(KokunoCurrentExteriorLMinusHBridge, "velocity")
