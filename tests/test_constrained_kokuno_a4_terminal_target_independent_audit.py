from __future__ import annotations

import copy
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_terminal_target_independent_audit import (
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_GATE,
    SIMPSON_SUBINTERVALS,
    _f_prime_ref,
    _q_p_ibp_reference,
    audit_loaded_terminal_target,
    default_target,
    enforce_preregistered_gates,
    public_api_has_no_scientific_tuning_knobs,
)
from openai_ns_reconstruction.kokuno_public_terminal_multiplier_target import (
    KokunoPublicTerminalMultiplierTarget,
)


def _roundtrip(tmp_path):
    target = default_target()
    path = tmp_path / "candidate.json"
    target.save_configuration(path)
    loaded = KokunoPublicTerminalMultiplierTarget.load_configuration(path)
    return target, loaded


def test_independent_terminal_target_passes_frozen_scoped_gates(tmp_path):
    target, loaded = _roundtrip(tmp_path)
    report = audit_loaded_terminal_target(loaded, target)
    enforce_preregistered_gates(report)

    q = report["q_p"]
    assert list(q["independent_by_resolution"]) == [str(n) for n in SIMPSON_SUBINTERVALS]
    assert q["fine_independent_vs_public_relative_difference"] <= 2.0e-10
    assert q["medium_to_fine_relative_change"] <= 5.0e-11
    assert report["slope"]["maximum"] < report["slope"]["h_over_4"]
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["final_project_admission_ready"] is False
    assert report["gates"]["final_project_momentum_gate_unchanged"] == FINAL_MOMENTUM_GATE
    assert report["gates"]["final_project_divergence_gate_unchanged"] == FINAL_DIVERGENCE_GATE


def test_independent_reference_is_three_resolution_and_nontrivial():
    target = default_target()
    h = target.h_value
    rho = target.rho_o
    values = [_q_p_ibp_reference(h, rho, n) for n in SIMPSON_SUBINTERVALS]
    assert all(np.isfinite(v) and 0.0 < v < h for v in values)
    assert abs(values[-1] - values[-2]) / values[-1] <= 5.0e-11


def test_wrong_terminal_slope_sign_is_detectable():
    target = default_target()
    y = np.asarray([2.0], dtype=float)
    correct = float(_f_prime_ref(y, target.rho_o)[0])
    wrong = float(_f_prime_ref(y, -target.rho_o)[0])
    assert correct > 0.0
    assert wrong < 0.0


@pytest.mark.parametrize("key,factor", [("h_current", 1.001), ("c_o_autonomous", 1.001), ("rho_o", 1.001)])
def test_serialized_parameter_drift_fails_closed(key, factor):
    target = default_target()
    payload = copy.deepcopy(target.configuration())
    payload[key] = float(payload[key]) * factor
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoPublicTerminalMultiplierTarget.from_configuration(payload)


def test_truth_promotion_fails_closed():
    target = default_target()
    payload = copy.deepcopy(target.configuration())
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoPublicTerminalMultiplierTarget.from_configuration(payload)


def test_scientific_report_tamper_is_rejected(tmp_path):
    target, loaded = _roundtrip(tmp_path)
    report = audit_loaded_terminal_target(loaded, target)
    bad = copy.deepcopy(report)
    bad["q_p"]["fine_independent_vs_public_relative_difference"] = 1.0e-3
    with pytest.raises(AssertionError):
        enforce_preregistered_gates(bad)

    bad = copy.deepcopy(report)
    bad["slope"]["maximum"] = bad["slope"]["h_over_4"]
    with pytest.raises(AssertionError):
        enforce_preregistered_gates(bad)

    bad = copy.deepcopy(report)
    bad["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError):
        enforce_preregistered_gates(bad)


def test_public_audit_surface_has_no_scientific_tuning_knobs():
    assert public_api_has_no_scientific_tuning_knobs() is True
    assert tuple(inspect.signature(audit_loaded_terminal_target).parameters) == (
        "loaded",
        "pre_serialization_reference",
    )
