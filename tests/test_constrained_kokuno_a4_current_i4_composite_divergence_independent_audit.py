from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_a4_current_i4_composite_divergence_independent_audit as audit


class _LinearSolenoidal:
    D = 1.0

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, float), np.asarray(y, float), np.asarray(z, float), np.asarray(t, float)
        )
        return np.stack((2.0 * x + y, -3.0 * y + z, x + y + z), axis=-1)


class _ZeroField:
    D = 1.0

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, float), np.asarray(y, float), np.asarray(z, float), np.asarray(t, float)
        )
        return np.stack((np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)), axis=-1)


def _metric(value: float = 1.0e-6, speed: float = 1.0) -> dict[str, float]:
    return {
        "sampled_max": value,
        "pooled_weighted_rms": value,
        "i4_entry_seam_sampled_max": value,
        "axis_axis_near_sampled_max": value,
        "speed_rms": speed,
    }


def _passing_receipt() -> dict:
    checks = {
        "save_reload_semantic_identity": True,
        "save_reload_velocity_replay_max_abs": 0.0,
        "velocity_mutation_detected": True,
        "offgrid_order_invariant_all_components": True,
        "rehash_consistent_leading_mutation_rejected": True,
        "rehash_consistent_oscillatory_mutation_rejected": True,
        "post_i4_fails_closed": True,
        "manufactured_solenoidal_calibration_pass": True,
    }
    resolutions = []
    for value in (4.0e-6, 2.0e-6, 1.0e-6):
        resolutions.append({"leading_plus_oscillatory": _metric(value)})
    return {
        "inner_oscillatory_signal_rms": 1.0e-6,
        "resolutions": resolutions,
        "checks": checks,
    }


def test_protocol_is_frozen_to_current_i4_heads_and_gates():
    assert audit.UPSTREAM_PR == 1080
    assert audit.UPSTREAM_HEAD == "c40d8ddecd2971544a6e07dab093436b423cf326"
    assert audit.AGENT1_PR == 1079
    assert audit.AGENT1_HEAD == "b06742ca6e189499192ede3cce40f62cdc1e35ca"
    assert audit.SEED == 9173831
    assert audit.SPATIAL_STEPS == (0.02, 0.01, 0.005)
    assert audit.DIVERGENCE_GATE == 1.0e-5


def test_truth_boundary_cannot_promote_full_pde_or_missing_source_corrections():
    truth = audit._TRUTH_BOUNDARY
    assert truth["current_i4_leading_plus_oscillatory_velocity_consumed"] is True
    assert truth["source_positive_order_i3_correction_materialized"] is False
    assert truth["source_i4_mean_correction_materialized"] is False
    assert truth["terminal_global_leading_completion_materialized"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_fd2_recovers_linear_solenoidal_divergence():
    points = np.asarray(((0.1, 0.2, -0.3), (-0.4, 0.5, 0.6)), dtype=float)
    times = np.asarray((0.31, 0.63), dtype=float)
    div, jac = audit._divergence(_LinearSolenoidal(), points, times, 0.005)
    assert jac.shape == (2, 3, 3)
    assert np.max(np.abs(div)) <= 2.0e-13


def test_manufactured_time_dependent_solenoidal_calibration():
    assert audit._manufactured_calibration() is True


def test_xlinear_mutation_adds_expected_divergence_signal():
    points = np.asarray(((0.2, -0.3, 0.1), (-0.4, 0.1, -0.2)), dtype=float)
    times = np.asarray((0.31, 0.63), dtype=float)
    base_div, _ = audit._divergence(_ZeroField(), points, times, 0.005)
    mutated_div, _ = audit._divergence(audit._XLinearMutation(_ZeroField()), points, times, 0.005)
    signal = np.max(np.abs(mutated_div - base_div))
    assert signal == pytest.approx(audit.MUTATION_EPSILON, rel=1.0e-12, abs=1.0e-14)
    assert signal >= audit.MUTATION_DETECTION_FLOOR


def test_weighted_metrics_match_direct_formula():
    div = np.asarray((1.0, -2.0, 3.0))
    weights = np.asarray((2.0, 1.0, 3.0))
    volume, l2, rms = audit._weighted_metrics(div, weights)
    expected_sq = 2.0 * 1.0**2 + 1.0 * 2.0**2 + 3.0 * 3.0**2
    assert volume == 6.0
    assert l2 == pytest.approx(np.sqrt(expected_sq))
    assert rms == pytest.approx(np.sqrt(expected_sq / volume))


def test_weighted_metrics_reject_nonpositive_volume():
    with pytest.raises(RuntimeError):
        audit._weighted_metrics(np.asarray((1.0,)), np.asarray((0.0,)))


def test_stability_policy_has_fixed_floor_and_factor():
    assert audit._stability_ok(1.0e-6, 1.0e-6)
    assert audit._stability_ok(0.0, audit.STABILITY_FLOOR)
    assert not audit._stability_ok(1.0e-6, 2.0e-6)


def test_synthetic_passing_receipt_passes_frozen_gate_logic():
    assert audit._derive_audit_pass(_passing_receipt()) is True


def test_each_scientific_gate_can_fail_without_threshold_rewrite():
    base = _passing_receipt()
    bad = copy.deepcopy(base)
    bad["resolutions"][-1]["leading_plus_oscillatory"]["sampled_max"] = 1.01e-5
    assert audit._derive_audit_pass(bad) is False

    bad = copy.deepcopy(base)
    bad["inner_oscillatory_signal_rms"] = 0.0
    assert audit._derive_audit_pass(bad) is False

    bad = copy.deepcopy(base)
    bad["checks"]["rehash_consistent_oscillatory_mutation_rejected"] = False
    assert audit._derive_audit_pass(bad) is False


def test_component_order_invariance_checks_all_aggregate_metrics():
    value = {
        "sampled_max": 1.0,
        "pooled_weighted_rms": 2.0,
        "pooled_volume_l2_estimate": 3.0,
    }
    assert audit._component_order_invariant(value, copy.deepcopy(value))
    bad = copy.deepcopy(value)
    bad["pooled_weighted_rms"] += 2.0 * audit.ORDER_INVARIANCE_ATOL
    assert not audit._component_order_invariant(value, bad)


def test_physical_volume_jacobian_positive_on_preregistered_chart():
    class _Field:
        D = 2.0

    eta = np.asarray((-0.45, 0.0, 0.45))
    t = np.asarray((0.31, 0.47, 0.71))
    jac = audit._physical_volume_jacobian(_Field(), eta, t)
    assert jac.shape == eta.shape
    assert np.all(np.isfinite(jac))
    assert np.all(jac > 0.0)
