from __future__ import annotations

import copy
import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_a4_relative_swirl_split_independent_audit as m
from openai_ns_reconstruction.kokuno_pa16_current_cartesian_relative_swirl_split import (
    KokunoPA16CurrentCartesianRelativeSwirlSplit,
)


class _Delegate:
    def __init__(self, inner):
        self.inner = inner

    def __getattr__(self, name):
        return getattr(self.inner, name)


class _ZeroDelta(_Delegate):
    def velocity_split(self, x, y, z, t):
        base, delta = self.inner.velocity_split(x, y, z, t)
        return base, np.zeros_like(np.asarray(delta, dtype=float))


class _AngularDivergentDelta(_Delegate):
    def velocity_split(self, x, y, z, t):
        base, delta = self.inner.velocity_split(x, y, z, t)
        theta = np.arctan2(np.asarray(y, dtype=float), np.asarray(x, dtype=float))
        factor = 1.0 + 0.05 * np.sin(theta)
        return base, np.asarray(delta, dtype=float) * np.expand_dims(factor, axis=-1) if np.ndim(factor) else np.asarray(delta, dtype=float) * float(factor)


def test_real_save_load_black_box_audit_emits_three_resolutions_and_truth_boundary(tmp_path):
    original = m.default_candidate()
    path = tmp_path / "candidate.json"
    original.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianRelativeSwirlSplit.load_configuration(path)
    report = m.audit_loaded_split(loaded, original)

    assert report["schema"] == m.SCHEMA
    assert report["task"] == "K4-VAL-121"
    assert report["upstream_head"] == m.UPSTREAM_HEAD
    assert report["seed"] == 9173921
    assert len(report["resolutions"]) == 3
    assert [row["relative_step"] for row in report["resolutions"]] == list(m.RELATIVE_FD_LADDER)
    assert len(report["active_points"]) == 16
    assert len(report["control_points"]) == 12
    assert report["save_load_replay_relative_max"] <= m.SAVE_LOAD_REPLAY_GATE
    assert report["off_bump_delta_abs_max"] <= m.ZERO_ABSOLUTE_GATE
    assert report["exact_axis_delta_abs_max"] <= m.ZERO_ABSOLUTE_GATE
    assert report["resolutions"][-1]["active_ratio_max"] > 0.0
    assert report["resolutions"][-1]["active_ratio_max"] < np.finfo(float).eps
    assert report["truth_boundary"]["binary64_total_relative_swirl_sum_is_resolved"] is False
    assert report["truth_boundary"]["current_cartesian_relative_swirl_composed"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["final_project_admission_ready"] is False


def test_scoped_fd_is_implementation_distinct_and_canonical_absolute_ladder_is_only_diagnostic():
    candidate = m.default_candidate()
    active, _ = m._heldout_points(candidate)
    assert len(active) == 16
    for row in m._canonical_representability(active)["ladder"]:
        assert row["absolute_step"] in m.CANONICAL_ABSOLUTE_FD_LADDER
        assert 0.0 <= row["representable_point_fraction"] <= 1.0
    sig = inspect.signature(m.audit_loaded_split)
    assert list(sig.parameters) == ["loaded", "pre_serialization_reference"]
    forbidden = {
        "threshold", "residual", "forcing", "pressure", "viscosity", "nu",
        "step", "gain", "optimizer", "tolerance", "seed",
    }
    assert forbidden.isdisjoint(sig.parameters)


def test_zero_delta_mutation_is_detected_by_nontriviality_floor():
    candidate = m.default_candidate()
    active, _ = m._heldout_points(candidate)
    obs = m._point_observables(_ZeroDelta(candidate), active[0])
    assert obs["delta_norm"] == 0.0
    assert obs["abs_delta_over_abs_base_tangential"] == 0.0

    report = {
        "resolutions": [
            {"normalized_divergence_sampled_max": 0.0},
            {"normalized_divergence_sampled_max": 0.0},
            {
                "normalized_divergence_sampled_max": 0.0,
                "normalized_divergence_sample_rms": 0.0,
                "relative_leakage_max": 0.0,
                "active_ratio_max": 0.0,
            },
        ],
        "off_bump_delta_abs_max": 0.0,
        "exact_axis_delta_abs_max": 0.0,
        "save_load_replay_relative_max": 0.0,
    }
    assert "delta_channel_collapsed" in m._assess(report)


def test_angular_modulation_mutation_creates_detectable_normalized_divergence():
    candidate = m.default_candidate()
    active, _ = m._heldout_points(candidate)
    mutated = _AngularDivergentDelta(candidate)
    row = m._resolution_metrics(mutated, active[:4], m.RELATIVE_FD_LADDER[-1])
    assert row["normalized_divergence_sampled_max"] > m.NORMALIZED_DIVERGENCE_GATE


def test_truth_and_provenance_promotion_fail_closed(tmp_path):
    original = m.default_candidate()
    path = tmp_path / "candidate.json"
    payload = original.save_configuration(path)

    for key in (
        "binary64_total_relative_swirl_sum_is_resolved",
        "current_cartesian_relative_swirl_composed",
        "pde_validated",
    ):
        mutated = copy.deepcopy(payload)
        mutated["truth_boundary"][key] = True
        with pytest.raises(ValueError, match="configuration/provenance drift"):
            KokunoPA16CurrentCartesianRelativeSwirlSplit.from_configuration(mutated)


def test_gate_enforcement_rejects_posthoc_metric_promotion():
    base = {
        "resolutions": [
            {"normalized_divergence_sampled_max": 1e-8},
            {"normalized_divergence_sampled_max": 1e-8},
            {
                "normalized_divergence_sampled_max": 1e-8,
                "normalized_divergence_sample_rms": 1e-8,
                "relative_leakage_max": 1e-12,
                "active_ratio_max": 1e-37,
            },
        ],
        "off_bump_delta_abs_max": 0.0,
        "exact_axis_delta_abs_max": 0.0,
        "save_load_replay_relative_max": 0.0,
        "truth_boundary": {"pde_validated": False},
        "final_project_admission_ready": False,
    }
    m.enforce_preregistered_gates(base)

    bad = copy.deepcopy(base)
    bad["resolutions"][-1]["normalized_divergence_sampled_max"] = 1e-3
    with pytest.raises(AssertionError, match="fine_normalized_divergence_max"):
        m.enforce_preregistered_gates(bad)

    bad = copy.deepcopy(base)
    bad["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError, match="pde_validated"):
        m.enforce_preregistered_gates(bad)

    bad = copy.deepcopy(base)
    bad["final_project_admission_ready"] = True
    with pytest.raises(AssertionError, match="final project admission"):
        m.enforce_preregistered_gates(bad)
