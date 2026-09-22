from __future__ import annotations

import inspect
import json
from pathlib import Path

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_a4_public_ideal_qs_independent_audit as module
from openai_ns_reconstruction.kokuno_public_ideal_exterior_qs_schedule import (
    KokunoPublicIdealExteriorQsSchedule,
)


def _roundtrip(tmp_path: Path):
    original = module.default_schedule()
    path = tmp_path / "candidate.json"
    original.save_configuration(path)
    loaded = KokunoPublicIdealExteriorQsSchedule.load_configuration(path)
    return original, loaded


def test_independent_three_resolution_schedule_passes_frozen_scoped_gates(tmp_path) -> None:
    original, loaded = _roundtrip(tmp_path)
    report = module.audit_loaded_public_ideal_schedule(loaded, original)
    module.enforce_preregistered_gates(report)

    assert report["task"] == "K4-VAL-127"
    assert report["upstream_head"] == module.UPSTREAM_HEAD
    assert report["independent_operator"]["steps_per_unit"] == [2048, 4096, 8192]
    assert report["independent_operator"]["seed"] == 9173981
    assert report["save_load_semantic_exact"] is True
    assert report["fine_endpoint_relative_max_to_public"] <= module.ENDPOINT_PUBLIC_REL_GATE
    assert report["fine_offgrid_release_relative_max_to_public"] <= module.OFFGRID_PUBLIC_REL_GATE
    assert report["medium_to_fine_endpoint_relative_max"] <= module.MEDIUM_FINE_ENDPOINT_REL_GATE
    assert report["fine_hold_linear_relative_max_to_public"] <= module.HOLD_LINEAR_REL_GATE
    assert report["fine_bridge_endpoint_relative_error_to_q_p"] <= module.BRIDGE_ENDPOINT_REL_GATE
    assert report["truth_boundary"]["current_candidate_q_s_audited"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["final_project_admission_ready"] is False


def test_reference_rk4_is_distinct_from_production_gl192_and_converges() -> None:
    # This regression exercises only the A4-owned ODE path; no A1 q_release*
    # method or production quadrature is used here.
    lam = 0.05
    h = 0.005
    q0 = (lam - h) / (1.0 - lam)
    values = [
        module._rk4_to(module._release1_rhs_ref, q0, 1.0, n, lam, h)
        for n in module.RK4_STEPS_PER_UNIT
    ]
    assert all(np.isfinite(values))
    assert all(v > 0.0 for v in values)
    assert module._relative_difference(values[-2], values[-1]) <= module.MEDIUM_FINE_ENDPOINT_REL_GATE

    params = inspect.signature(module._rk4_to).parameters
    assert "nodes" not in params
    assert "weights" not in params
    assert "quadrature" not in params


def test_wrong_ode_sign_is_detectably_incompatible_with_public_schedule() -> None:
    schedule = module.default_schedule()
    lam = schedule.lambda_value
    h = schedule.h_value
    q0 = (lam - h) / (1.0 - lam)

    def wrong_rhs(s, q, lam, h):
        # Deliberately flip the source term sign while preserving the damping term.
        sig = float(module._sigma_ref(np.asarray(s)))
        ell = -lam - (1.0 - lam) * sig
        return +ell + h - (1.0 + ell) * q

    wrong = module._rk4_to(wrong_rhs, q0, 1.0, 2048, lam, h)
    assert module._relative_difference(wrong, schedule.q_release1_end_source_ideal) > 1.0e-2


def test_public_value_corruption_is_caught_after_reference_is_formed(tmp_path, monkeypatch) -> None:
    original, loaded = _roundtrip(tmp_path)
    actual = loaded.q_release2

    def corrupted(s):
        return np.asarray(actual(s), dtype=float) * 1.001

    monkeypatch.setattr(loaded, "q_release2", corrupted)
    report = module.audit_loaded_public_ideal_schedule(loaded, original)
    assert report["fine_offgrid_release_relative_max_to_public"] > module.OFFGRID_PUBLIC_REL_GATE
    with pytest.raises(AssertionError):
        module.enforce_preregistered_gates(report)


def test_configuration_truth_or_provenance_drift_fails_closed() -> None:
    schedule = module.default_schedule()
    payload = schedule.configuration()

    bad_truth = json.loads(json.dumps(payload))
    bad_truth["truth_boundary"]["current_q_s_release2_endpoint_materialized"] = True
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoPublicIdealExteriorQsSchedule.from_configuration(bad_truth)

    bad_source = json.loads(json.dumps(payload))
    bad_source["source_blob"] = "0" * 40
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoPublicIdealExteriorQsSchedule.from_configuration(bad_source)


def test_tampered_report_cannot_promote_current_or_pde_truth(tmp_path) -> None:
    original, loaded = _roundtrip(tmp_path)
    report = module.audit_loaded_public_ideal_schedule(loaded, original)

    tampered = json.loads(json.dumps(report))
    tampered["truth_boundary"]["current_candidate_q_s_audited"] = True
    with pytest.raises(AssertionError):
        module.enforce_preregistered_gates(tampered)

    tampered = json.loads(json.dumps(report))
    tampered["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError):
        module.enforce_preregistered_gates(tampered)


def test_public_audit_api_has_no_scientific_tuning_knobs() -> None:
    assert module.public_api_has_no_scientific_tuning_knobs() is True
    assert tuple(inspect.signature(module.audit_loaded_public_ideal_schedule).parameters) == (
        "loaded",
        "pre_serialization_reference",
    )
    assert module.FINAL_MOMENTUM_GATE == 1.0e-3
    assert module.FINAL_DIVERGENCE_GATE == 1.0e-5
    assert module.FINAL_QUADRATURE == (24, 48, 96)
