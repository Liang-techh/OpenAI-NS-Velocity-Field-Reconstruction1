from __future__ import annotations

import copy
import hashlib
import json

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_a4_current_partial_leading_oscillatory_divergence_independent_audit as audit


class _LinearField:
    def __init__(self, matrix: np.ndarray):
        self.matrix = np.asarray(matrix, dtype=float)

    def velocity(self, x, y, z, t):
        xx, yy, zz, _tt = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        pts = np.stack((xx, yy, zz), axis=-1)
        return np.einsum("ij,...j->...i", self.matrix, pts)


class _SumField:
    def __init__(self, leading, oscillatory):
        self.leading = leading
        self.oscillatory = oscillatory

    def velocity(self, x, y, z, t):
        return self.leading.velocity(x, y, z, t) + self.oscillatory.velocity(x, y, z, t)


class _SemanticField:
    def __init__(self):
        self._payload = {
            "schema": "manufactured-semantic-field",
            "upstream": {"head": "abc", "version": 1},
            "truth_boundary": {"pde_validated": False},
        }
        self.semantic_sha256 = audit._sha256(self._payload)

    def semantic_payload(self):
        return copy.deepcopy(self._payload)


LEADING_MATRIX = np.asarray(
    [
        [1.0, 2.0, 0.0],
        [0.0, -0.5, 1.0],
        [0.3, 0.0, -0.5],
    ],
    dtype=float,
)
OSC_MATRIX = np.asarray(
    [
        [0.0, 0.2, 0.0],
        [0.0, 0.0, 0.1],
        [-0.2, 0.0, 0.0],
    ],
    dtype=float,
)


def test_independent_fd2_recovers_linear_jacobian_and_zero_divergence():
    field = _LinearField(LEADING_MATRIX)
    points = np.asarray(((0.3, -0.2, 0.4), (-0.6, 0.5, -0.1)), dtype=float)
    times = np.asarray((0.31, 0.63), dtype=float)
    jac = audit.independent_fd2_jacobian(field, points, times, 0.005)
    expected = np.broadcast_to(LEADING_MATRIX, jac.shape)
    assert np.max(np.abs(jac - expected)) <= 1.0e-12
    assert np.max(np.abs(np.trace(jac, axis1=1, axis2=2))) <= 1.0e-12


def test_heldout_probes_are_deterministic_physical_and_axis_explicit():
    first = audit.make_heldout_probes()
    second = audit.make_heldout_probes()
    assert np.array_equal(first.points, second.points)
    assert np.array_equal(first.times, second.times)
    assert np.array_equal(first.weights, second.weights)
    assert first.points.shape == (len(audit.TIMES) * audit.POINTS_PER_TIME, 3)
    assert first.axis_points.shape == (5, 3)
    radius = np.hypot(first.points[:, 0], first.points[:, 1])
    assert np.all(radius >= audit.RADIAL_INTERVAL[0])
    assert np.all(radius <= audit.RADIAL_INTERVAL[1])
    assert np.all(first.points[:, 2] >= audit.AXIAL_INTERVAL[0])
    assert np.all(first.points[:, 2] <= audit.AXIAL_INTERVAL[1])
    assert np.all(first.weights > 0.0)
    axis_radius = np.hypot(first.axis_points[:, 0], first.axis_points[:, 1])
    assert axis_radius[0] == 0.0
    assert axis_radius[1] == 0.0
    assert np.all(axis_radius[2:] > 0.0)
    assert np.max(axis_radius[2:]) <= 1.0e-6 * 1.000001


def test_staged_metrics_recover_solenoidal_leading_total_and_increment():
    leading = _LinearField(LEADING_MATRIX)
    oscillatory = _LinearField(OSC_MATRIX)
    total = _SumField(leading, oscillatory)
    probes = audit.make_heldout_probes()
    metrics = audit._metrics_for_step(total, leading, probes, audit.SPATIAL_STEPS[-1])
    assert metrics["leading"]["sampled_max"] <= 1.0e-11
    assert metrics["total"]["sampled_max"] <= 1.0e-11
    assert metrics["oscillatory_increment"]["sampled_max"] <= 1.0e-11
    assert metrics["leading"]["pooled_weighted_rms"] <= 1.0e-11
    assert metrics["total"]["pooled_weighted_rms"] <= 1.0e-11
    assert metrics["oscillatory_increment"]["pooled_weighted_rms"] <= 1.0e-11
    assert metrics["oscillatory_speed_rms"] > audit.NONTRIVIAL_OSC_SPEED_RMS


def test_cartesian_divergence_mutation_is_detected():
    leading = _LinearField(LEADING_MATRIX)
    total = _SumField(leading, _LinearField(OSC_MATRIX))
    mutated = audit._XLinearMutation(total, epsilon=audit.MUTATION_EPSILON)
    probes = audit.make_heldout_probes()
    metrics = audit._metrics_for_step(mutated, leading, probes, audit.SPATIAL_STEPS[-1])
    assert metrics["total"]["sampled_max"] >= audit.MUTATION_DETECTION_FLOOR


def test_semantic_snapshot_replay_is_deterministic():
    field = _SemanticField()
    ok, observed = audit._semantic_snapshot_replay(field)
    assert ok is True
    assert observed == field.semantic_sha256


def _metric_group(value: float) -> dict[str, object]:
    return {
        "sampled_max": value,
        "pooled_weighted_rms": value,
        "pooled_volume_l2_estimate": value,
        "normalized_sampled_max": value,
        "normalized_weighted_rms": value,
        "per_time": [],
        "worst_divergence": value,
        "worst_jacobian_frobenius": 1.0,
        "worst_normalized_divergence": value,
    }


def _resolution(step: float, value: float) -> dict[str, object]:
    return {
        "step": step,
        "total": _metric_group(value),
        "leading": _metric_group(value),
        "oscillatory_increment": _metric_group(value),
        "axis_axis_near": {
            "total_sampled_max": value,
            "leading_sampled_max": value,
            "oscillatory_increment_sampled_max": value,
        },
        "total_speed_rms": 1.0,
        "oscillatory_speed_rms": 0.1,
        "worst_total_witness": {
            "point": [0.6, 0.0, 0.0],
            "time": 0.31,
            "divergence": value,
            "normalized_divergence": value,
        },
        "worst_oscillatory_witness": {
            "point": [0.6, 0.0, 0.0],
            "time": 0.31,
            "divergence_increment": value,
            "normalized_divergence": value,
        },
    }


def _synthetic_receipt() -> dict[str, object]:
    receipt: dict[str, object] = {
        "schema": audit.SCHEMA,
        "upstream": {
            "pr": audit.UPSTREAM_PR,
            "head": audit.UPSTREAM_HEAD,
            "source_blob": audit.UPSTREAM_SOURCE_BLOB,
            "expected_source_blob": audit.UPSTREAM_SOURCE_BLOB,
            "semantic_sha256": "a" * 64,
            "agent1_pr": audit.AGENT1_PR,
            "agent1_head": audit.AGENT1_HEAD,
        },
        "protocol": {
            "seed": audit.SEED,
            "spatial_steps": list(audit.SPATIAL_STEPS),
            "times": list(audit.TIMES),
            "radial_interval": list(audit.RADIAL_INTERVAL),
            "axial_interval": list(audit.AXIAL_INTERVAL),
            "points_per_time": audit.POINTS_PER_TIME,
            "integration_probe_count": len(audit.TIMES) * audit.POINTS_PER_TIME,
            "axis_axis_near_probe_count": 5,
            "operator": "centered Cartesian FD2 on public total and public leading velocity only",
            "divergence_gate": audit.DIVERGENCE_GATE,
            "stability_factor": audit.STABILITY_FACTOR,
            "stability_floor": audit.STABILITY_FLOOR,
            "nontrivial_total_speed_rms": audit.NONTRIVIAL_TOTAL_SPEED_RMS,
            "nontrivial_osc_speed_rms": audit.NONTRIVIAL_OSC_SPEED_RMS,
            "mutation_epsilon": audit.MUTATION_EPSILON,
            "mutation_detection_floor": audit.MUTATION_DETECTION_FLOOR,
            "volume_norm_note": "synthetic",
        },
        "resolutions": [
            _resolution(audit.SPATIAL_STEPS[0], 4.0e-6),
            _resolution(audit.SPATIAL_STEPS[1], 3.0e-6),
            _resolution(audit.SPATIAL_STEPS[2], 2.0e-6),
        ],
        "mutation": {
            "kind": "total u_x += 1e-3*x while leading reference is unchanged",
            "finest_total_sampled_max": 1.0e-3,
            "finest_total_weighted_rms": 1.0e-3,
        },
        "checks": {
            "upstream_source_blob_exact": True,
            "semantic_snapshot_replay": True,
            "upstream_truth_boundary_fail_closed": True,
            "velocity_mutation_detected": True,
            "offgrid_order_invariant": True,
        },
        "truth_boundary": copy.deepcopy(audit._TRUTH_BOUNDARY),
    }
    receipt["audit_pass"] = audit._derive_audit_pass(receipt)
    receipt["receipt_sha256"] = hashlib.sha256(
        audit._canonical_json(receipt).encode("utf-8")
    ).hexdigest()
    return receipt


def _rehash(receipt: dict[str, object]) -> None:
    payload = copy.deepcopy(receipt)
    payload.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = hashlib.sha256(
        audit._canonical_json(payload).encode("utf-8")
    ).hexdigest()


def test_synthetic_receipt_passes_frozen_enforcement():
    receipt = _synthetic_receipt()
    assert receipt["audit_pass"] is True
    audit.enforce_receipt(receipt)


def test_protocol_and_truth_mutations_fail_closed():
    protocol_drift = _synthetic_receipt()
    protocol_drift["protocol"]["divergence_gate"] = 2.0e-5
    _rehash(protocol_drift)
    with pytest.raises(ValueError, match="protocol drift"):
        audit.enforce_receipt(protocol_drift)

    truth_drift = _synthetic_receipt()
    truth_drift["truth_boundary"]["pde_validated"] = True
    _rehash(truth_drift)
    with pytest.raises(ValueError, match="truth-boundary drift"):
        audit.enforce_receipt(truth_drift)


def test_failed_scientific_gate_cannot_be_relabelled_pass():
    receipt = _synthetic_receipt()
    receipt["resolutions"][-1]["total"]["sampled_max"] = 2.0e-4
    receipt["audit_pass"] = True
    _rehash(receipt)
    with pytest.raises(ValueError, match="audit_pass"):
        audit.enforce_receipt(receipt)
