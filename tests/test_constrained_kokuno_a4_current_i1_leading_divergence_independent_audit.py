from __future__ import annotations

import copy
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_current_i1_leading_divergence_independent_audit import (
    DIVERGENCE_GATE,
    MUTATION_DETECTION_FLOOR,
    SEED,
    SPATIAL_STEPS,
    UPSTREAM_HEAD,
    _XLinearMutation,
    _derive_audit_pass,
    _gate_mutation_rejected,
    _nested_gate_mutation_rejected,
    _post_i1_fails_closed,
    independent_fd2_jacobian,
    make_heldout_probes,
    materialize_receipt,
)
from openai_ns_reconstruction.kokuno_pa16_current_cartesian_i1_frozen_gate import (
    I1_CLOSURE_TOLERANCE,
    KokunoPA16CurrentCartesianI1FrozenGate,
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
        # trace(grad u) = 1 - 1 + 0 = 0, with explicit time dependence.
        return np.stack((x + 0.2 * t, -y + 0.1 * z, 0.3 * x - 0.1 * t), axis=-1)


def _div(field, points, times, h):
    jac = independent_fd2_jacobian(field, points, times, h)
    return np.trace(jac, axis1=1, axis2=2)


def test_manufactured_solenoidal_field_locks_fd_sign_and_trace():
    field = _ManufacturedSolenoidal()
    points = np.asarray(
        [
            [0.17, -0.31, 0.11],
            [-0.44, 0.29, -0.08],
            [0.03, 0.02, 0.51],
        ],
        dtype=float,
    )
    times = np.asarray([0.31, 0.47, 0.71], dtype=float)
    for h in SPATIAL_STEPS:
        assert np.max(np.abs(_div(field, points, times, h))) < 2.0e-13


def test_x_linear_mutation_is_detected_independently():
    base = _ManufacturedSolenoidal()
    mutated = _XLinearMutation(base)
    points = np.asarray(
        [[0.17, -0.31, 0.11], [-0.44, 0.29, -0.08]], dtype=float
    )
    times = np.asarray([0.31, 0.63], dtype=float)
    delta = _div(mutated, points, times, SPATIAL_STEPS[-1]) - _div(
        base, points, times, SPATIAL_STEPS[-1]
    )
    assert np.max(np.abs(delta)) >= MUTATION_DETECTION_FLOOR
    np.testing.assert_allclose(delta, 1.0e-3, rtol=0.0, atol=2.0e-13)


def test_fresh_probe_protocol_is_deterministic_and_covers_i1():
    field = KokunoPA16CurrentCartesianI1FrozenGate()
    a = make_heldout_probes(field)
    b = make_heldout_probes(field)
    np.testing.assert_array_equal(a.points, b.points)
    np.testing.assert_array_equal(a.times, b.times)
    np.testing.assert_array_equal(a.weights, b.weights)
    np.testing.assert_array_equal(a.region, b.region)
    assert SEED == 9173801
    labels = set(a.region.tolist())
    assert "pre_i1_control" in labels
    assert {"i1_zone_0", "i1_zone_1", "i1_zone_2"}.issubset(labels)
    assert a.seam_points.shape == (16, 3)
    assert a.axis_points.shape == (5, 3)
    assert np.all(a.weights > 0.0)


def test_frozen_gate_mutations_fail_closed():
    field = KokunoPA16CurrentCartesianI1FrozenGate()
    assert field.closure_tolerance == I1_CLOSURE_TOLERANCE == 5.0e-7
    assert _gate_mutation_rejected(field)
    assert _nested_gate_mutation_rejected(field)


def test_post_i1_probe_fails_closed():
    field = KokunoPA16CurrentCartesianI1FrozenGate()
    assert _post_i1_fails_closed(field)


def test_materializer_exposes_no_scientific_tuning_knob():
    assert list(inspect.signature(materialize_receipt).parameters) == []


def _passing_synthetic_receipt():
    resolution = {
        "sampled_max": 2.0e-6,
        "pooled_weighted_rms": 1.0e-6,
        "i1_entry_seam_sampled_max": 2.0e-6,
        "axis_axis_near_sampled_max": 2.0e-6,
        "speed_rms": 1.0e-2,
    }
    return {
        "upstream_head": UPSTREAM_HEAD,
        "resolutions": [copy.deepcopy(resolution) for _ in SPATIAL_STEPS],
        "checks": {
            "save_reload_semantic_identity": True,
            "save_reload_velocity_replay_max_abs": 0.0,
            "velocity_mutation_detected": True,
            "offgrid_order_invariant": True,
            "frozen_gate_mutation_rejected": True,
            "nested_parent_gate_mutation_rejected": True,
            "post_i1_fails_closed": True,
        },
    }


def test_preregistered_gate_is_not_posthoc_relaxable():
    receipt = _passing_synthetic_receipt()
    assert _derive_audit_pass(receipt)
    receipt["resolutions"][-1]["sampled_max"] = DIVERGENCE_GATE * 1.001
    assert not _derive_audit_pass(receipt)


def test_resolution_degradation_is_rejected():
    receipt = _passing_synthetic_receipt()
    receipt["resolutions"][1]["pooled_weighted_rms"] = 1.0e-6
    receipt["resolutions"][2]["pooled_weighted_rms"] = 4.0e-6
    assert not _derive_audit_pass(receipt)


def test_wrapper_rejects_serialized_gate_drift_directly():
    field = KokunoPA16CurrentCartesianI1FrozenGate()
    payload = copy.deepcopy(field.configuration())
    payload["frozen_i1_closure_gate"]["closure_tolerance"] = 1.0e-5
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianI1FrozenGate.from_configuration(payload)


def test_upstream_identity_is_frozen():
    assert UPSTREAM_HEAD == "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5"
