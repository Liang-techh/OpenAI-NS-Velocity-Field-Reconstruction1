import copy
import hashlib

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_current_cartesian_leading_divergence_independent_audit import (
    DIVERGENCE_GATE,
    ETA_INTERVAL,
    MUTATION_DETECTION_FLOOR,
    MUTATION_EPSILON,
    NONTRIVIAL_SPEED_RMS,
    POINTS_PER_TIME,
    SCHEMA,
    SEED,
    SPATIAL_STEPS,
    STABILITY_FACTOR,
    STABILITY_FLOOR,
    TIMES,
    UPSTREAM_HEAD,
    UPSTREAM_SOURCE_BLOB,
    X_INTERVAL,
    _TRUTH_BOUNDARY,
    _XLinearMutation,
    _canonical_json,
    _metrics_for_step,
    _physical_volume_jacobian,
    enforce_receipt,
    independent_fd2_jacobian,
    make_heldout_probes,
)


class _ManufacturedSolenoidalField:
    D = 0.37

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        # A nontrivial linear divergence-free field with time dependence.
        u = -2.0 * y + 0.3 * z + 0.2 * t
        v = 1.5 * x + y - 0.4 * z
        w = -z + 0.1 * x - 0.3 * t
        return np.stack((u, v, w), axis=-1)


@pytest.fixture(scope="module")
def manufactured():
    return _ManufacturedSolenoidalField()


def test_fd2_cartesian_jacobian_reconstructs_linear_solenoidal_field(manufactured):
    points = np.asarray(
        [[0.3, -0.2, 0.1], [-0.5, 0.4, -0.2], [0.0, 0.0, 0.0]], dtype=float
    )
    times = np.asarray([0.31, 0.47, 0.71], dtype=float)
    expected = np.asarray(
        [[0.0, -2.0, 0.3], [1.5, 1.0, -0.4], [0.1, 0.0, -1.0]], dtype=float
    )
    for step in SPATIAL_STEPS:
        got = independent_fd2_jacobian(manufactured, points, times, step)
        np.testing.assert_allclose(got, expected[None, :, :], atol=3e-13, rtol=3e-13)
        np.testing.assert_allclose(np.trace(got, axis1=1, axis2=2), 0.0, atol=3e-13)


def test_frozen_probe_generator_is_deterministic_positive_weight_and_offgrid(manufactured):
    a = make_heldout_probes(manufactured)
    b = make_heldout_probes(manufactured)
    np.testing.assert_array_equal(a.points, b.points)
    np.testing.assert_array_equal(a.times, b.times)
    np.testing.assert_array_equal(a.weights, b.weights)
    assert a.points.shape == (POINTS_PER_TIME * len(TIMES), 3)
    assert a.axis_points.shape == (5, 3)
    assert np.all(a.weights > 0.0)
    assert np.all(np.isfinite(a.weights))
    assert np.all(np.isfinite(a.points))
    assert np.any(np.abs(a.points[:, 0]) > 1e-6)
    assert np.any(np.abs(a.points[:, 1]) > 1e-6)


def test_public_velocity_mutation_is_detected_by_same_operator(manufactured):
    probes = make_heldout_probes(manufactured)
    clean = _metrics_for_step(manufactured, probes, SPATIAL_STEPS[-1])
    mutant = _metrics_for_step(_XLinearMutation(manufactured), probes, SPATIAL_STEPS[-1])
    assert clean["sampled_max"] < 1e-11
    assert mutant["sampled_max"] >= MUTATION_DETECTION_FLOOR
    assert mutant["pooled_weighted_rms"] >= MUTATION_DETECTION_FLOOR
    assert mutant["sampled_max"] == pytest.approx(MUTATION_EPSILON, rel=2e-10, abs=2e-12)


def test_axis_and_offgrid_order_invariance_on_manufactured_field(manufactured):
    probes = make_heldout_probes(manufactured)
    fine = _metrics_for_step(manufactured, probes, SPATIAL_STEPS[-1])
    reverse = np.arange(probes.points.shape[0] - 1, -1, -1)
    reversed_probes = type(probes)(
        points=probes.points[reverse],
        times=probes.times[reverse],
        weights=probes.weights[reverse],
        time_index=probes.time_index[reverse],
        axis_points=probes.axis_points,
        axis_times=probes.axis_times,
    )
    replay = _metrics_for_step(manufactured, reversed_probes, SPATIAL_STEPS[-1])
    assert fine["axis_axis_near_sampled_max"] < 1e-11
    assert replay["sampled_max"] == pytest.approx(fine["sampled_max"], abs=1e-13)
    assert replay["pooled_weighted_rms"] == pytest.approx(
        fine["pooled_weighted_rms"], abs=1e-13
    )


def test_physical_volume_jacobian_is_positive(manufactured):
    eta = np.asarray([-0.65, -0.2, 0.0, 0.2, 0.65])
    t = np.asarray([0.31, 0.47, 0.63, 0.71, 0.5])
    jac = _physical_volume_jacobian(manufactured, eta, t)
    assert np.all(jac > 0.0)
    assert np.all(np.isfinite(jac))


def _passing_receipt():
    resolution = {
        "step": SPATIAL_STEPS[-1],
        "sampled_max": 2e-7,
        "pooled_weighted_rms": 1e-7,
        "pooled_volume_l2_estimate": 2e-7,
        "normalized_sampled_max": 2e-7,
        "normalized_weighted_rms": 1e-7,
        "axis_axis_near_sampled_max": 2e-7,
        "speed_rms": 1.0,
        "per_time": [],
        "worst_witness": {},
    }
    receipt = {
        "schema": SCHEMA,
        "upstream": {
            "pr": 965,
            "head": UPSTREAM_HEAD,
            "source_blob": UPSTREAM_SOURCE_BLOB,
            "candidate_semantic_sha256": "1" * 64,
        },
        "protocol": {
            "seed": SEED,
            "spatial_steps": list(SPATIAL_STEPS),
            "times": list(TIMES),
            "X_interval": list(X_INTERVAL),
            "eta_interval": list(ETA_INTERVAL),
            "points_per_time": POINTS_PER_TIME,
            "integration_probe_count": 32,
            "axis_axis_near_probe_count": 5,
            "operator": "centered Cartesian FD2 on public velocity only",
            "divergence_gate": DIVERGENCE_GATE,
            "stability_factor": STABILITY_FACTOR,
            "stability_floor": STABILITY_FLOOR,
            "nontrivial_speed_rms": NONTRIVIAL_SPEED_RMS,
            "mutation_epsilon": MUTATION_EPSILON,
            "mutation_detection_floor": MUTATION_DETECTION_FLOOR,
            "volume_norm_note": "scoped",
        },
        "resolutions": [
            dict(resolution, step=SPATIAL_STEPS[0], sampled_max=4e-7, pooled_weighted_rms=3e-7),
            dict(resolution, step=SPATIAL_STEPS[1], sampled_max=3e-7, pooled_weighted_rms=2e-7),
            resolution,
        ],
        "mutation": {
            "kind": "u_x += 1e-3*x",
            "finest_sampled_max": 1e-3,
            "finest_weighted_rms": 1e-3,
        },
        "checks": {
            "save_reload_semantic_identity": True,
            "config_mutation_changes_semantic_identity": True,
            "velocity_mutation_detected": True,
            "offgrid_order_invariant": True,
        },
        "refinement": {},
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        "audit_pass": True,
    }
    receipt["receipt_sha256"] = hashlib.sha256(
        _canonical_json(receipt).encode("utf-8")
    ).hexdigest()
    return receipt


def _rehash(receipt):
    payload = copy.deepcopy(receipt)
    payload.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = hashlib.sha256(
        _canonical_json(payload).encode("utf-8")
    ).hexdigest()


def test_receipt_enforcement_fails_closed_on_threshold_or_truth_laundering():
    receipt = _passing_receipt()
    enforce_receipt(receipt)

    drift = copy.deepcopy(receipt)
    drift["protocol"]["divergence_gate"] = 1e-4
    _rehash(drift)
    with pytest.raises(ValueError, match="protocol drift"):
        enforce_receipt(drift)

    promoted = copy.deepcopy(receipt)
    promoted["truth_boundary"]["pde_validated"] = True
    _rehash(promoted)
    with pytest.raises(ValueError, match="truth-boundary drift"):
        enforce_receipt(promoted)

    false_pass = copy.deepcopy(receipt)
    false_pass["resolutions"][-1]["sampled_max"] = 2e-4
    _rehash(false_pass)
    with pytest.raises(ValueError, match="audit_pass"):
        enforce_receipt(false_pass)


def test_receipt_checksum_detects_posthoc_metric_mutation():
    receipt = _passing_receipt()
    receipt["resolutions"][-1]["normalized_sampled_max"] *= 2.0
    with pytest.raises(ValueError, match="checksum"):
        enforce_receipt(receipt)
