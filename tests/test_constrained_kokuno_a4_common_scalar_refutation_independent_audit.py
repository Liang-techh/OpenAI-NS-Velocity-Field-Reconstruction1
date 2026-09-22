from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_common_scalar_refutation_independent_audit import (
    FD4_Y_STEPS,
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_GATE,
    FINAL_QUADRATURE,
    FINE_NORMALIZED_SLOPE_MIN,
    PRODUCTION_ROOT_DISAGREEMENT_GATE,
    REFINEMENT_FLOOR,
    REFINEMENT_WORSEN_FACTOR,
    ROOT_SNAPSHOT_ITERATIONS,
    ROOT_TARGET_NORM_GATE,
    SEED,
    _independent_reference,
    _regula_falsi_snapshots,
    audit_loaded_common_scalar_refutation,
    default_candidate,
    enforce_preregistered_gates,
    heldout_eta,
    public_api_has_no_scientific_tuning_knobs,
)
from openai_ns_reconstruction.kokuno_current_bridge_common_scalar_refutation import (
    KokunoCurrentBridgeCommonScalarRefutation,
)


def _audit_roundtrip(tmp_path):
    candidate = default_candidate()
    path = candidate.save_configuration(tmp_path / "candidate.json")
    loaded = KokunoCurrentBridgeCommonScalarRefutation.load_configuration(path)
    return audit_loaded_common_scalar_refutation(loaded, candidate)


def test_protocol_is_frozen_before_results() -> None:
    assert SEED == 9174011
    assert FD4_Y_STEPS == (2.0**-6, 2.0**-7, 2.0**-8)
    assert ROOT_SNAPSHOT_ITERATIONS == (32, 48, 64)
    assert ROOT_TARGET_NORM_GATE == 2.0e-10
    assert PRODUCTION_ROOT_DISAGREEMENT_GATE == 2.0e-8
    assert FINE_NORMALIZED_SLOPE_MIN == 1.0e-8
    assert REFINEMENT_WORSEN_FACTOR == 1.25
    assert REFINEMENT_FLOOR == 5.0e-10
    assert FINAL_MOMENTUM_GATE == 1.0e-3
    assert FINAL_DIVERGENCE_GATE == 1.0e-5
    assert FINAL_QUADRATURE == (24, 48, 96)


def test_heldout_eta_is_deterministic_strict_interior_and_not_parent_grid() -> None:
    candidate = default_candidate()
    a = heldout_eta(candidate)
    b = heldout_eta(candidate)
    assert np.array_equal(a, b)
    assert a.shape == (41,)
    assert len(np.unique(a)) == 41
    lo, hi = candidate.eta_interval
    assert np.all((a > lo) & (a < hi))
    production = candidate.witness_eta()
    distances = np.min(np.abs(a[:, None] - production[None, :]), axis=1)
    assert np.count_nonzero(distances > 1.0e-10) >= 30


def test_regula_falsi_is_implementation_distinct_and_converges() -> None:
    # Positive-to-negative monotone manufactured residual with root log(3).
    residual = lambda y: 3.0 * np.exp(-float(y)) - 1.0
    snapshots = _regula_falsi_snapshots(
        residual,
        0.0,
        2.0,
        residual(0.0),
        residual(2.0),
    )
    assert set(snapshots) == set(ROOT_SNAPSHOT_ITERATIONS)
    mids = [snapshots[it]["midpoint"] for it in ROOT_SNAPSHOT_ITERATIONS]
    assert abs(mids[-1] - np.log(3.0)) < 1.0e-11
    assert snapshots[64]["lower_residual"] >= 0.0
    assert snapshots[64]["upper_residual"] <= 0.0


def test_independent_reference_does_not_call_a1_certificate_or_bisection(
    monkeypatch,
) -> None:
    candidate = default_candidate()

    def forbidden(*args, **kwargs):
        raise AssertionError("production witness path entered independent reference")

    monkeypatch.setattr(
        KokunoCurrentBridgeCommonScalarRefutation,
        "unique_root_certificate",
        forbidden,
    )
    monkeypatch.setattr(
        KokunoCurrentBridgeCommonScalarRefutation,
        "unique_root_bracket",
        forbidden,
    )
    report = _independent_reference(candidate)
    assert report["eta_count"] == 41
    assert report["root_found_count"] >= 1
    assert report["finite_sample_spread_is_diagnostic_only"] is True


def test_roundtrip_preserves_truth_boundary_and_three_slope_resolutions(tmp_path) -> None:
    report = _audit_roundtrip(tmp_path)
    assert report["task"] == "K4-VAL-130"
    assert report["fd4_y_steps"] == list(FD4_Y_STEPS)
    assert report["root_snapshot_iterations"] == list(ROOT_SNAPSHOT_ITERATIONS)
    assert report["production_report_read_only_after_reference"] is True
    assert report["save_load"] == {
        "semantic_exact": True,
        "configuration_exact": True,
        "entry_q_exact": True,
    }
    independent = report["independent_reference"]
    assert independent["eta_count"] == 41
    assert independent["root_found_count"] >= 1
    for record in independent["records"]:
        assert len(record["fd4_steps"]) == 3
        assert len(record["fd4_slopes"]) == 3
        assert len(record["root_by_resolution"]) == 3
    truth = report["truth_boundary"]
    assert truth["current_common_scalar_refutation_independent_audit_executed"] is True
    assert truth["full_eta_interval_target_totality_established"] is False
    assert truth["full_eta_interval_root_uniqueness_established"] is False
    assert truth["smooth_eta_target_time_map_established"] is False
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["eta_dependent_cartesian_matching_boundary_materialized"] is False
    assert truth["complete_ns_admission_ready"] is False
    assert truth["pde_validated"] is False
    assert report["final_project_admission_ready"] is False


def test_preregistered_gates_are_enforced_on_real_roundtrip(tmp_path) -> None:
    report = _audit_roundtrip(tmp_path)
    enforce_preregistered_gates(report)


def test_production_pair_corruption_is_detected(monkeypatch, tmp_path) -> None:
    original = KokunoCurrentBridgeCommonScalarRefutation.refutation_report

    def corrupted(self):
        out = dict(original(self))
        if out["common_scalar_bridge_refuted_by_disjoint_unique_root_brackets"]:
            pair = dict(out["disjoint_root_pair"])
            pair["earlier_root_bracket"] = [
                float(pair["earlier_root_bracket"][0]) + 0.1,
                float(pair["earlier_root_bracket"][1]) + 0.1,
            ]
            out["disjoint_root_pair"] = pair
        return out

    monkeypatch.setattr(
        KokunoCurrentBridgeCommonScalarRefutation,
        "refutation_report",
        corrupted,
    )
    report = _audit_roundtrip(tmp_path)
    pair = report["production_pair_independent_audit"]
    if pair["production_refutation_claim_present"]:
        assert max(pair["production_vs_independent_root_abs_disagreement"]) > PRODUCTION_ROOT_DISAGREEMENT_GATE
        with pytest.raises(AssertionError):
            enforce_preregistered_gates(report)
    else:
        assert pair["absence_of_production_pair_is_inconclusive"] is True


def test_absence_of_production_pair_is_not_promoted_to_common_bridge(
    monkeypatch, tmp_path
) -> None:
    original = KokunoCurrentBridgeCommonScalarRefutation.refutation_report

    def inconclusive(self):
        out = dict(original(self))
        out["common_scalar_bridge_refuted_by_disjoint_unique_root_brackets"] = False
        out["disjoint_root_pair"] = None
        out["disjoint_gap_lower_bound"] = 0.0
        return out

    monkeypatch.setattr(
        KokunoCurrentBridgeCommonScalarRefutation,
        "refutation_report",
        inconclusive,
    )
    report = _audit_roundtrip(tmp_path)
    pair = report["production_pair_independent_audit"]
    assert pair["production_refutation_claim_present"] is False
    assert pair["absence_of_production_pair_is_inconclusive"] is True
    assert report["truth_boundary"]["full_eta_interval_common_scalar_bridge_length_established"] is False
    enforce_preregistered_gates(report)


def test_truth_promotion_fails_closed(tmp_path) -> None:
    report = _audit_roundtrip(tmp_path)
    report["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError, match="PDE truth promotion"):
        enforce_preregistered_gates(report)


def test_public_audit_api_exposes_no_scientific_tuning_knobs() -> None:
    sig = inspect.signature(audit_loaded_common_scalar_refutation)
    assert tuple(sig.parameters) == ("loaded", "pre_serialization_reference")
    assert public_api_has_no_scientific_tuning_knobs() is True
