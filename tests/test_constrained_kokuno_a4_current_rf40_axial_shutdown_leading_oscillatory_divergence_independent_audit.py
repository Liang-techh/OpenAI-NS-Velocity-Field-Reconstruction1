from __future__ import annotations

import copy
import hashlib
import inspect
import math

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_a4_current_rf40_axial_shutdown_leading_oscillatory_divergence_independent_audit as a4


class LinearField:
    D = 1.0
    X_R = 2.0 * np.exp(5.0)
    X_1 = X_R * np.e
    T_d = 2.0
    X_2 = X_1 * np.exp(T_d)

    def __init__(self, matrix):
        self.matrix = np.asarray(matrix, dtype=float)

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, float), np.asarray(y, float), np.asarray(z, float), np.asarray(t, float)
        )
        p = np.stack((x, y, z), axis=-1)
        out = np.einsum("ij,...j->...i", self.matrix, p)
        out[..., 0] += 0.10 * t
        out[..., 1] -= 0.04 * t
        return out


LEAD = LinearField([[1.2, 0.3, 0.0], [-0.4, -0.7, 0.2], [0.2, 0.0, -0.5]])
OSC = LinearField([[0.2, -0.1, 0.0], [0.1, -0.3, 0.4], [0.0, -0.4, 0.1]])


class SumField:
    D = LEAD.D
    X_R = LEAD.X_R
    X_1 = LEAD.X_1
    T_d = LEAD.T_d
    X_2 = LEAD.X_2

    def velocity(self, x, y, z, t):
        return LEAD.velocity(x, y, z, t) + OSC.velocity(x, y, z, t)


def _seal(receipt):
    payload = copy.deepcopy(receipt)
    payload.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = hashlib.sha256(
        a4._canonical_json(payload).encode("utf-8")
    ).hexdigest()
    return receipt


def synthetic_receipt():
    stats = {
        "sampled_max": 2.0e-7,
        "rms": 1.0e-7,
        "normalized_sampled_max": 2.0e-7,
        "normalized_rms": 1.0e-7,
        "worst_divergence": 2.0e-7,
        "worst_normalized_divergence": 2.0e-7,
    }
    weighted = dict(
        stats,
        estimated_volume=1.0,
        weighted_rms=1.0e-7,
        volume_l2_estimate=1.0e-7,
        normalized_weighted_rms=1.0e-7,
    )

    def resolution(step):
        return {
            "step": step,
            "inner": {
                "total": copy.deepcopy(weighted),
                "leading": copy.deepcopy(weighted),
                "oscillatory_increment": copy.deepcopy(weighted),
                "total_speed_rms": 1.0,
                "oscillatory_speed_rms": 0.1,
                "oscillatory_speed_abs_max": 0.2,
            },
            "rf40_axial_shutdown": {
                "total": copy.deepcopy(weighted),
                "leading": copy.deepcopy(weighted),
                "oscillatory_increment": copy.deepcopy(weighted),
                "total_speed_rms": 1.0,
                "oscillatory_speed_rms": 0.0,
                "oscillatory_speed_abs_max": 0.0,
            },
            "x1_seam": {
                "total": copy.deepcopy(stats),
                "leading": copy.deepcopy(stats),
                "oscillatory_increment": copy.deepcopy(stats),
                "total_speed_rms": 1.0,
                "oscillatory_speed_rms": 0.0,
                "oscillatory_speed_abs_max": 0.0,
            },
            "axis_axis_near": {
                "total": copy.deepcopy(stats),
                "leading": copy.deepcopy(stats),
                "oscillatory_increment": copy.deepcopy(stats),
                "total_speed_rms": 1.0,
                "oscillatory_speed_rms": 0.0,
                "oscillatory_speed_abs_max": 0.0,
            },
            "axial_shutdown_by_time_zone": {},
            "worst_inner_witnesses": {},
            "worst_axial_shutdown_witnesses": {},
        }

    receipt = {
        "schema": a4.SCHEMA,
        "task": a4.TASK,
        "upstream": {
            "head": a4.UPSTREAM_HEAD,
            "expected_source_blob": a4.UPSTREAM_SOURCE_BLOB,
            "source_blob": a4.UPSTREAM_SOURCE_BLOB,
            "agent1_head": a4.AGENT1_HEAD,
            "agent1_source_blob": a4.AGENT1_SOURCE_BLOB,
        },
        "protocol": {
            "seed": a4.SEED,
            "steps": list(a4.STEPS),
            "times": list(a4.TIMES),
            "inner_radial_interval": list(a4.INNER_RADIAL_INTERVAL),
            "inner_axial_interval": list(a4.INNER_AXIAL_INTERVAL),
            "inner_points_per_time": a4.INNER_POINTS_PER_TIME,
            "eta_interval": list(a4.ETA_INTERVAL),
            "axial_shutdown_zones": [list(v) for v in a4.AXIAL_SHUTDOWN_ZONES],
            "axial_points_per_zone_time": a4.AXIAL_POINTS_PER_ZONE_TIME,
            "x1_seam_log_offsets": list(a4.X1_SEAM_LOG_OFFSETS),
            "divergence_gate": a4.DIVERGENCE_GATE,
            "stability_factor": a4.STABILITY_FACTOR,
            "stability_floor": a4.STABILITY_FLOOR,
            "total_speed_floor": a4.TOTAL_SPEED_FLOOR,
            "oscillatory_speed_floor": a4.OSCILLATORY_SPEED_FLOOR,
            "mutation_epsilon": a4.MUTATION_EPSILON,
            "mutation_detection_floor": a4.MUTATION_DETECTION_FLOOR,
            "replay_atol": a4.REPLAY_ATOL,
        },
        "resolutions": [resolution(step) for step in a4.STEPS],
        "checks": {
            "upstream_source_blob_exact": True,
            "save_reload_semantic_identity": True,
            "save_reload_velocity_replay": True,
            "independent_leading_identity_matches_saved_composite": True,
            "upstream_truth_boundary_fail_closed": True,
            "velocity_mutation_detected": True,
            "offgrid_order_invariant": True,
            "leading_parameter_perturbation_changes_semantic_identity": True,
            "oscillatory_runtime_parameter_mutation_rejected": True,
            "post_X2_fail_closed": True,
        },
        "truth_boundary": copy.deepcopy(a4.TRUTH_BOUNDARY),
    }
    receipt["audit_pass"] = a4._derive_audit_pass(receipt)
    return _seal(receipt)


def test_public_audit_entrypoints_do_not_accept_scientific_tuning_knobs():
    assert list(inspect.signature(a4.materialize_receipt).parameters) == []
    assert list(inspect.signature(a4.fd2_jacobian).parameters) == ["field", "points", "times", "step"]
    assert a4.DIVERGENCE_GATE == 1.0e-5
    assert a4.STEPS == (0.02, 0.01, 0.005)
    assert a4.SEED == 9173721


def test_manufactured_staged_fields_are_solenoidal_under_fd2():
    points = np.asarray([[0.2, -0.3, 0.4], [0.5, 0.1, -0.2], [-0.4, 0.7, 0.3]])
    times = np.asarray([0.31, 0.47, 0.71])
    total = SumField()
    for h in a4.STEPS:
        for field in (LEAD, OSC, total):
            jac = a4.fd2_jacobian(field, points, times, h)
            assert np.max(np.abs(np.trace(jac, axis1=1, axis2=2))) < 3.0e-13


def test_sampler_is_deterministic_strictly_in_axial_shutdown_and_contains_all_strata():
    f = SumField()
    p1, p2 = a4.make_probes(f), a4.make_probes(f)
    np.testing.assert_array_equal(p1.inner.points, p2.inner.points)
    np.testing.assert_array_equal(p1.axial_shutdown.points, p2.axial_shutdown.points)
    np.testing.assert_array_equal(p1.inner.weights, p2.inner.weights)
    np.testing.assert_array_equal(p1.axial_shutdown.weights, p2.axial_shutdown.weights)
    assert len(p1.inner.points) == len(a4.TIMES) * a4.INNER_POINTS_PER_TIME
    assert len(p1.axial_shutdown.points) == (
        len(a4.TIMES) * len(a4.AXIAL_SHUTDOWN_ZONES) * a4.AXIAL_POINTS_PER_ZONE_TIME
    )
    assert len(p1.seam_points) == len(a4.TIMES) * len(a4.X1_SEAM_LOG_OFFSETS)
    assert len(p1.axis_points) == 5
    assert np.all(p1.inner.weights > 0.0)
    assert np.all(p1.axial_shutdown.weights > 0.0)
    assert all(0.0 < r0 < r1 < 1.0 for _, r0, r1 in a4.AXIAL_SHUTDOWN_ZONES)
    assert math.isclose(math.log(f.X_2 / f.X_1), f.T_d, rel_tol=2.0e-15)


def test_staged_metrics_and_mutation_detection_on_manufactured_fields():
    f = SumField()
    probes = a4.make_probes(f)
    metrics = a4.metrics_for_step(f, LEAD, probes, a4.STEPS[-1])
    assert metrics["inner"]["total"]["sampled_max"] < 3.0e-11
    assert metrics["inner"]["leading"]["sampled_max"] < 3.0e-11
    assert metrics["inner"]["oscillatory_increment"]["sampled_max"] < 3.0e-11
    assert metrics["rf40_axial_shutdown"]["total"]["sampled_max"] < 3.0e-9
    assert len(metrics["axial_shutdown_by_time_zone"]["total"]) == (
        len(a4.TIMES) * len(a4.AXIAL_SHUTDOWN_ZONES)
    )

    mutated = a4.metrics_for_step(a4.XLinearMutation(f), LEAD, probes, a4.STEPS[-1])
    assert max(
        mutated["inner"]["total"]["sampled_max"],
        mutated["rf40_axial_shutdown"]["total"]["sampled_max"],
    ) >= 9.9e-4


def test_frozen_gate_accepts_good_receipt_and_rejects_laundering():
    receipt = synthetic_receipt()
    a4.enforce_receipt(receipt)
    bad = copy.deepcopy(receipt)
    bad["resolutions"][-1]["rf40_axial_shutdown"]["total"]["sampled_max"] = 2.0e-5
    bad["audit_pass"] = True
    _seal(bad)
    with pytest.raises(ValueError, match="laundering"):
        a4.enforce_receipt(bad)


def test_protocol_truth_and_checksum_mutations_fail_closed():
    receipt = synthetic_receipt()

    bad = copy.deepcopy(receipt)
    bad["protocol"]["divergence_gate"] = 2.0e-5
    _seal(bad)
    with pytest.raises(ValueError, match="protocol drift"):
        a4.enforce_receipt(bad)

    bad = copy.deepcopy(receipt)
    bad["truth_boundary"]["pde_validated"] = True
    _seal(bad)
    with pytest.raises(ValueError, match="truth-boundary drift"):
        a4.enforce_receipt(bad)

    bad = copy.deepcopy(receipt)
    bad["resolutions"][-1]["inner"]["total"]["sampled_max"] *= 1.01
    with pytest.raises(ValueError, match="checksum"):
        a4.enforce_receipt(bad)
