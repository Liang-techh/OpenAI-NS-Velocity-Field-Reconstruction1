from __future__ import annotations

import inspect
from decimal import Decimal, localcontext

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_relative_swirl_decimal_total_independent_audit import (
    AUDIT_DECIMAL_DIGITS,
    CANONICAL_ABSOLUTE_FD_LADDER,
    FINAL_PROJECT_DIVERGENCE_GATE,
    FINAL_PROJECT_MOMENTUM_GATE,
    NORMALIZED_DIVERGENCE_GATE,
    RELATIVE_FD_LADDER,
    SEED,
    TASK,
    _canonical_representability,
    _configure_decimal,
    _fd4_jacobian_decimal,
    _heldout_points,
    audit_loaded_decimal_total,
    default_candidate,
    enforce_preregistered_gates,
    public_api_has_no_scientific_tuning_knobs,
)
from openai_ns_reconstruction.kokuno_pa16_current_cartesian_relative_swirl_decimal_composed import (
    KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed,
)


def _load_roundtrip(tmp_path):
    original = default_candidate()
    path = tmp_path / "candidate.json"
    original.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed.load_configuration(path)
    return original, loaded


def test_protocol_is_frozen_and_public_api_has_no_scientific_tuning_knobs():
    assert TASK == "K4-VAL-122"
    assert SEED == 9173931
    assert RELATIVE_FD_LADDER == (2.0e-6, 1.0e-6, 5.0e-7)
    assert CANONICAL_ABSOLUTE_FD_LADDER == (0.02, 0.01, 0.005)
    assert NORMALIZED_DIVERGENCE_GATE == 1.0e-5
    assert FINAL_PROJECT_MOMENTUM_GATE == 1.0e-3
    assert FINAL_PROJECT_DIVERGENCE_GATE == 1.0e-5
    assert public_api_has_no_scientific_tuning_knobs()
    assert set(inspect.signature(audit_loaded_decimal_total).parameters) == {
        "loaded",
        "pre_serialization_reference",
    }


def test_heldout_points_are_seeded_offgrid_and_cover_both_bumps(tmp_path):
    _, loaded = _load_roundtrip(tmp_path)
    a0, c0 = _heldout_points(loaded)
    a1, c1 = _heldout_points(loaded)
    assert a0 == a1
    assert c0 == c1
    assert len(a0) == 16
    assert len(c0) == 12
    roles = {str(p["role"]) for p in a0}
    assert any(role.startswith("bump1_interior") for role in roles)
    assert any(role.startswith("bump2_interior") for role in roles)
    assert all(float(p["x"]) != 0.0 and float(p["y"]) != 0.0 for p in a0)


class _ManufacturedSolenoidal:
    def velocity(self, x, y, z, t):
        with localcontext() as ctx:
            _configure_decimal(ctx)
            xd = Decimal.from_float(float(x))
            yd = Decimal.from_float(float(y))
            td = Decimal.from_float(float(t))
            return np.asarray((-td * yd, td * xd, Decimal(0)), dtype=object)


def test_independent_fd4_operator_calibrates_on_time_dependent_solenoidal_field():
    point = {
        "x": 0.37,
        "y": -0.23,
        "z": 0.19,
        "t": 0.47,
        "r": 1.0,
        "q": 1.0,
        "role": "manufactured",
    }
    jac = _fd4_jacobian_decimal(_ManufacturedSolenoidal(), point, 1.0e-5)
    with localcontext() as ctx:
        _configure_decimal(ctx)
        div = jac[0][0] + jac[1][1] + jac[2][2]
    assert abs(div) < Decimal("1e-80")
    assert AUDIT_DECIMAL_DIGITS >= 100


def test_canonical_absolute_stencil_collapse_is_detected_without_retuning():
    point = {"x": 1.0e40, "y": -1.0e40, "z": 0.0}
    result = _canonical_representability([point])
    assert result["all_representable"] is False
    assert any(not row["all_points_all_axes_representable"] for row in result["ladder"])


def test_real_save_load_public_total_audit_emits_three_resolution_receipt(tmp_path):
    original, loaded = _load_roundtrip(tmp_path)
    report = audit_loaded_decimal_total(loaded, original)
    assert report["task"] == TASK
    assert report["candidate_semantic_sha256"] == loaded.semantic_sha256
    assert report["relative_fd_ladder"] == list(RELATIVE_FD_LADDER)
    assert len(report["resolutions"]) == 3
    assert report["save_load_replay_exact"] is True
    assert report["decimal_string_export_replay_exact"] is True
    composition = report["composition_checks"]
    assert composition["active_total_minus_base_replays_delta_exact"] is True
    assert composition["active_binary64_base_plus_delta_erased"] is True
    assert composition["active_correction_ratio_max"] > 0.0
    assert composition["off_bump_total_equals_independent_base_exact"] is True
    assert composition["exact_axis_transverse_velocity_zero"] is True
    assert report["final_project_admission_ready"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    # Scientific PASS is deliberately not hard-coded here: exact-head workflow
    # uploads the raw report before enforcing the preregistered gate.
    assert isinstance(report["failures"], list)


class _DropDeltaMutation(KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed):
    def velocity(self, x, y, z, t):
        base, _ = self.split.velocity_split(x, y, z, t)
        base = np.asarray(base, dtype=float)
        out = np.empty(base.shape, dtype=object)
        with localcontext() as ctx:
            ctx.prec = 96
            for idx in np.ndindex(base.shape):
                out[idx] = +Decimal.from_float(float(base[idx]))
        return out


class _DivergentMutation(KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed):
    def velocity(self, x, y, z, t):
        out = np.asarray(super().velocity(x, y, z, t), dtype=object).copy()
        shape = out.shape[:-1]
        xb = np.broadcast_to(np.asarray(x, dtype=float), shape)
        with localcontext() as ctx:
            _configure_decimal(ctx)
            for idx in np.ndindex(shape):
                out[idx + (0,)] = out[idx + (0,)] + Decimal("1e-3") * Decimal.from_float(
                    float(xb[idx])
                )
        return out


def test_dropping_sub_epsilon_delta_fails_independent_composition_gate(tmp_path):
    reference = default_candidate()
    mutated = _DropDeltaMutation()
    report = audit_loaded_decimal_total(mutated, reference)
    assert "decimal_total_does_not_replay_parent_delta" in report["failures"]
    assert report["scoped_decimal_total_gate_passed"] is False
    with pytest.raises(AssertionError, match="K4-VAL-122"):
        enforce_preregistered_gates(report)


def test_explicit_divergent_total_mutation_is_detected(tmp_path):
    reference = default_candidate()
    mutated = _DivergentMutation()
    report = audit_loaded_decimal_total(mutated, reference)
    assert (
        "fine_total_normalized_divergence_max" in report["failures"]
        or "fine_total_normalized_divergence_rms" in report["failures"]
    )
    with pytest.raises(AssertionError, match="K4-VAL-122"):
        enforce_preregistered_gates(report)


def test_truth_promotion_fails_closed_before_scientific_execution(tmp_path):
    original, loaded = _load_roundtrip(tmp_path)

    class _TruthMutation(KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed):
        @property
        def truth_boundary(self):
            truth = dict(super().truth_boundary)
            truth["pde_validated"] = True
            return truth

    mutated = _TruthMutation()
    with pytest.raises(ValueError, match="truth promotion"):
        audit_loaded_decimal_total(mutated, original)
