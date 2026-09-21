from __future__ import annotations

import copy
import inspect
import math

import numpy as np

from openai_ns_reconstruction.kokuno_a4_current_i1_composite_divergence_independent_audit import (
    DIVERGENCE_GATE,
    MUTATION_DETECTION_FLOOR,
    OSCILLATORY_SIGNAL_RMS_MIN,
    SEED,
    SPATIAL_STEPS,
    UPSTREAM_HEAD,
    _TRUTH_BOUNDARY,
    _XLinearMutation,
    _derive_audit_pass,
    _sha256,
    independent_fd2_jacobian,
    make_heldout_probes,
    materialize_receipt,
)


class _ManufacturedSolenoidal:
    D = 1.0

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.stack((x + 0.2 * t, -y + 0.1 * z, 0.3 * x - 0.1 * t), axis=-1)


class _ProbeGeometryOnly:
    D = 1.0
    X_I1_start = 2.0
    X_I1_end = 4.0


def _div(field, points, times, h):
    jac = independent_fd2_jacobian(field, points, times, h)
    return np.trace(jac, axis1=1, axis2=2)


def test_manufactured_solenoidal_field_locks_fd_sign_and_trace():
    field = _ManufacturedSolenoidal()
    points = np.asarray(
        [[0.17, -0.31, 0.11], [-0.44, 0.29, -0.08], [0.03, 0.02, 0.51]],
        dtype=float,
    )
    times = np.asarray([0.31, 0.47, 0.71], dtype=float)
    for h in SPATIAL_STEPS:
        assert np.max(np.abs(_div(field, points, times, h))) < 2.0e-13


def test_x_linear_mutation_is_detected_independently():
    base = _ManufacturedSolenoidal()
    mutated = _XLinearMutation(base)
    points = np.asarray([[0.17, -0.31, 0.11], [-0.44, 0.29, -0.08]], dtype=float)
    times = np.asarray([0.31, 0.63], dtype=float)
    delta = _div(mutated, points, times, SPATIAL_STEPS[-1]) - _div(
        base, points, times, SPATIAL_STEPS[-1]
    )
    assert np.max(np.abs(delta)) >= MUTATION_DETECTION_FLOOR
    np.testing.assert_allclose(delta, 1.0e-3, rtol=0.0, atol=2.0e-13)


def test_fresh_probe_protocol_is_deterministic_and_covers_three_i1_zones():
    field = _ProbeGeometryOnly()
    a = make_heldout_probes(field)
    b = make_heldout_probes(field)
    np.testing.assert_array_equal(a.points, b.points)
    np.testing.assert_array_equal(a.times, b.times)
    np.testing.assert_array_equal(a.weights, b.weights)
    np.testing.assert_array_equal(a.region, b.region)
    assert SEED == 9173811
    labels = set(a.region.tolist())
    assert labels == {"pre_i1_control", "i1_zone_0", "i1_zone_1", "i1_zone_2"}
    assert a.points.shape == (64, 3)
    assert a.seam_points.shape == (16, 3)
    assert a.axis_points.shape == (5, 3)
    assert np.all(a.weights > 0.0)


def test_materializer_exposes_no_scientific_tuning_knob():
    assert list(inspect.signature(materialize_receipt).parameters) == []


def test_identity_hash_is_canonical_under_key_reordering():
    left = {"b": [1.0, 2.0], "a": {"x": 3.0}}
    right = {"a": {"x": 3.0}, "b": [1.0, 2.0]}
    assert _sha256(left) == _sha256(right)


def _metric(maximum=2.0e-6, rms=1.0e-6):
    return {
        "sampled_max": maximum,
        "pooled_weighted_rms": rms,
        "pooled_volume_l2_estimate": rms,
        "i1_entry_seam_sampled_max": 2.0e-6,
        "axis_axis_near_sampled_max": 2.0e-6,
        "speed_rms": 1.0e-2,
    }


def _passing_synthetic_receipt():
    resolution = {
        "step": SPATIAL_STEPS[0],
        "leading": _metric(),
        "leading_plus_oscillatory": _metric(),
        "oscillatory_increment": _metric(),
    }
    resolutions = [copy.deepcopy(resolution) for _ in SPATIAL_STEPS]
    for entry, step in zip(resolutions, SPATIAL_STEPS):
        entry["step"] = step
    return {
        "upstream_head": UPSTREAM_HEAD,
        "inner_oscillatory_signal_rms": max(1.0e-6, OSCILLATORY_SIGNAL_RMS_MIN),
        "resolutions": resolutions,
        "checks": {
            "save_reload_semantic_identity": True,
            "save_reload_velocity_replay_max_abs": 0.0,
            "velocity_mutation_detected": True,
            "offgrid_order_invariant_all_components": True,
            "nested_agent1_gate_mutation_rejected": True,
            "rehash_consistent_oscillatory_mutation_rejected": True,
            "post_i1_fails_closed": True,
            "runtime_closure_gate_is_frozen": True,
        },
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
    }


def test_preregistered_divergence_gate_is_not_posthoc_relaxable():
    receipt = _passing_synthetic_receipt()
    assert _derive_audit_pass(receipt)
    receipt["resolutions"][-1]["leading_plus_oscillatory"]["sampled_max"] = DIVERGENCE_GATE * 1.001
    assert not _derive_audit_pass(receipt)


def test_resolution_degradation_is_rejected():
    receipt = _passing_synthetic_receipt()
    receipt["resolutions"][1]["leading_plus_oscillatory"]["pooled_weighted_rms"] = 1.0e-6
    receipt["resolutions"][2]["leading_plus_oscillatory"]["pooled_weighted_rms"] = 4.0e-6
    assert not _derive_audit_pass(receipt)


def test_zero_oscillatory_signal_cannot_vacuously_pass():
    receipt = _passing_synthetic_receipt()
    receipt["inner_oscillatory_signal_rms"] = 0.0
    assert not _derive_audit_pass(receipt)


def test_mutation_firewall_is_required_for_admission():
    for key in (
        "save_reload_semantic_identity",
        "velocity_mutation_detected",
        "offgrid_order_invariant_all_components",
        "nested_agent1_gate_mutation_rejected",
        "rehash_consistent_oscillatory_mutation_rejected",
        "post_i1_fails_closed",
        "runtime_closure_gate_is_frozen",
    ):
        receipt = _passing_synthetic_receipt()
        receipt["checks"][key] = False
        assert not _derive_audit_pass(receipt), key


def test_upstream_and_final_truth_boundaries_are_frozen():
    assert UPSTREAM_HEAD == "109527f520abb29bbe10372b0517eda44bcad0b6"
    assert _TRUTH_BOUNDARY["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert _TRUTH_BOUNDARY["after_correction_ns_residual_assessed"] is False
    assert _TRUTH_BOUNDARY["pde_validated"] is False
    assert math.isclose(DIVERGENCE_GATE, 1.0e-5, rel_tol=0.0, abs_tol=0.0)
