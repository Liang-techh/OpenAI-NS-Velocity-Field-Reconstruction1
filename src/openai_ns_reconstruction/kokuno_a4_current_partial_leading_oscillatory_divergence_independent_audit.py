"""Independent Cartesian divergence audit for current partial Kokuno leading+oscillation.

This Agent-4 increment consumes the exact Agent-2 PR #970 public composite

    velocity(x,y,z,t) = u_lead_current_through_Xh + u_osc_frozen_complete_curl

and differentiates only that public velocity surface (plus the separately public
Agent-1 leading velocity used for staged attribution).  No Agent-2 Jacobian,
divergence, vorticity, training tensor, pressure, forcing, residual, or source
coordinate derivative helper is used by the numerical reference path.

The result is deliberately scoped to the currently materialized Cartesian
domain through X_h.  It is not a complete Navier--Stokes momentum audit and is
not canonical whole-domain 24/48/96 admission.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
import hashlib
import inspect
import json
import math
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

from .kokuno_current_partial_leading_oscillatory_velocity import default_field as _upstream_default_field

SCHEMA = "kokuno-a4-current-partial-leading-oscillatory-divergence-audit-v1"
UPSTREAM_PR = 970
UPSTREAM_HEAD = "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4"
UPSTREAM_SOURCE_BLOB = "532435705bd5db341e82371180da8c04c3dd6c14"
AGENT1_PR = 965
AGENT1_HEAD = "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"

SEED = 9173671
SPATIAL_STEPS = (0.02, 0.01, 0.005)
TIMES = (0.31, 0.47, 0.63, 0.71)
RADIAL_INTERVAL = (0.50, 0.95)
AXIAL_INTERVAL = (-1.00, 1.00)
POINTS_PER_TIME = 8

DIVERGENCE_GATE = 1.0e-5
STABILITY_FACTOR = 1.25
STABILITY_FLOOR = 2.0e-8
NONTRIVIAL_TOTAL_SPEED_RMS = 1.0e-10
NONTRIVIAL_OSC_SPEED_RMS = 1.0e-12
MUTATION_EPSILON = 1.0e-3
MUTATION_DETECTION_FLOOR = 5.0e-4

_TRUTH_BOUNDARY = {
    "current_partial_leading_oscillatory_velocity_consumed": True,
    "semantic_artifact_replay_required": True,
    "upstream_composite_direct_save_load_available": False,
    "independent_cartesian_fd_operator_used": True,
    "leading_divergence_scoped_assessed": True,
    "leading_plus_oscillatory_divergence_scoped_assessed": True,
    "oscillatory_divergence_increment_scoped_assessed": True,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "velocity_beyond_Xh_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "agent3_correction_velocity_materialized": False,
    "canonical_24_48_96_volume_admission_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _git_blob_sha1(path: str | Path) -> str:
    raw = Path(path).read_bytes()
    header = f"blob {len(raw)}\0".encode("ascii")
    return hashlib.sha1(header + raw).hexdigest()


def _latin_hypercube(rng: np.random.Generator, count: int) -> np.ndarray:
    out = (np.arange(count, dtype=float) + rng.random(count)) / float(count)
    rng.shuffle(out)
    return out


@dataclass(frozen=True)
class HeldoutProbeSet:
    points: np.ndarray
    times: np.ndarray
    weights: np.ndarray
    time_index: np.ndarray
    axis_points: np.ndarray
    axis_times: np.ndarray


def make_heldout_probes() -> HeldoutProbeSet:
    """Frozen fresh physical-Cartesian probes, independent of candidate internals."""
    rng = np.random.default_rng(SEED)
    points_all: list[np.ndarray] = []
    times_all: list[np.ndarray] = []
    weights_all: list[np.ndarray] = []
    time_index_all: list[np.ndarray] = []

    r0, r1 = RADIAL_INTERVAL
    z0, z1 = AXIAL_INTERVAL
    param_measure = (r1 - r0) * (z1 - z0) * (2.0 * math.pi)

    for time_id, time in enumerate(TIMES):
        ur = _latin_hypercube(rng, POINTS_PER_TIME)
        uz = _latin_hypercube(rng, POINTS_PER_TIME)
        uth = _latin_hypercube(rng, POINTS_PER_TIME)
        radius = r0 + (r1 - r0) * ur
        z = z0 + (z1 - z0) * uz
        theta = 2.0 * math.pi * uth
        points = np.stack((radius * np.cos(theta), radius * np.sin(theta), z), axis=-1)
        times = np.full(POINTS_PER_TIME, time, dtype=float)

        # Sampling is uniform in (r,z,theta), so physical dV = r dr dz dtheta.
        weights = param_measure * radius / float(POINTS_PER_TIME)
        points_all.append(points)
        times_all.append(times)
        weights_all.append(weights)
        time_index_all.append(np.full(POINTS_PER_TIME, time_id, dtype=int))

    axis_radius = np.asarray([0.0, 0.0, 1.0e-10, 1.0e-8, 1.0e-6], dtype=float)
    axis_theta = np.asarray([0.0, 1.1, 0.7, 1.7, 2.3], dtype=float)
    axis_z = np.asarray([-0.60, 0.0, 0.55, -0.25, 0.30], dtype=float)
    axis_times = np.asarray([0.31, 0.47, 0.63, 0.71, 0.55], dtype=float)
    axis_points = np.stack(
        (axis_radius * np.cos(axis_theta), axis_radius * np.sin(axis_theta), axis_z),
        axis=-1,
    )

    return HeldoutProbeSet(
        points=np.concatenate(points_all, axis=0),
        times=np.concatenate(times_all, axis=0),
        weights=np.concatenate(weights_all, axis=0),
        time_index=np.concatenate(time_index_all, axis=0),
        axis_points=axis_points,
        axis_times=axis_times,
    )


def _velocity(field: Any, points: np.ndarray, times: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    times = np.asarray(times, dtype=float)
    out = np.asarray(field.velocity(points[:, 0], points[:, 1], points[:, 2], times), dtype=float)
    if out.shape != (points.shape[0], 3):
        raise RuntimeError(f"unexpected velocity shape {out.shape}")
    if np.any(~np.isfinite(out)):
        raise RuntimeError("public velocity returned non-finite values")
    return out


def independent_fd2_jacobian(
    field: Any,
    points: np.ndarray,
    times: np.ndarray,
    step: float,
) -> np.ndarray:
    """Independent centered Cartesian FD2 Jacobian from public velocity only."""
    points = np.asarray(points, dtype=float)
    times = np.asarray(times, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape [N,3]")
    if times.shape != (points.shape[0],):
        raise ValueError("times must have shape [N]")
    h = float(step)
    if not math.isfinite(h) or h <= 0.0:
        raise ValueError("step must be finite and positive")

    jac = np.empty((points.shape[0], 3, 3), dtype=float)
    for axis in range(3):
        plus = points.copy()
        minus = points.copy()
        plus[:, axis] += h
        minus[:, axis] -= h
        jac[:, :, axis] = (_velocity(field, plus, times) - _velocity(field, minus, times)) / (2.0 * h)
    if np.any(~np.isfinite(jac)):
        raise RuntimeError("independent Cartesian Jacobian became non-finite")
    return jac


def _field_divergence_metrics(
    jac: np.ndarray,
    weights: np.ndarray,
    time_index: np.ndarray,
) -> dict[str, Any]:
    div = np.trace(jac, axis1=1, axis2=2)
    frob = np.linalg.norm(jac, axis=(1, 2))
    normalized = np.abs(div) / np.maximum(1.0, frob)

    per_time: list[dict[str, Any]] = []
    for time_id, time in enumerate(TIMES):
        mask = time_index == time_id
        w = weights[mask]
        d = div[mask]
        volume = float(np.sum(w))
        per_time.append(
            {
                "time": float(time),
                "estimated_volume": volume,
                "volume_l2_estimate": float(np.sqrt(np.sum(w * d * d))),
                "weighted_rms": float(np.sqrt(np.sum(w * d * d) / volume)),
                "sampled_max": float(np.max(np.abs(d))),
            }
        )
    worst = int(np.argmax(np.abs(div)))
    return {
        "sampled_max": float(np.max(np.abs(div))),
        "pooled_weighted_rms": float(np.sqrt(np.sum(weights * div * div) / np.sum(weights))),
        "pooled_volume_l2_estimate": float(np.sqrt(np.sum(weights * div * div))),
        "normalized_sampled_max": float(np.max(normalized)),
        "normalized_weighted_rms": float(
            np.sqrt(np.sum(weights * normalized * normalized) / np.sum(weights))
        ),
        "per_time": per_time,
        "worst_index": worst,
        "worst_divergence": float(div[worst]),
        "worst_jacobian_frobenius": float(frob[worst]),
        "worst_normalized_divergence": float(normalized[worst]),
    }


def _metrics_for_step(total_field: Any, leading_field: Any, probes: HeldoutProbeSet, step: float) -> dict[str, Any]:
    total_jac = independent_fd2_jacobian(total_field, probes.points, probes.times, step)
    leading_jac = independent_fd2_jacobian(leading_field, probes.points, probes.times, step)
    osc_jac = total_jac - leading_jac

    total = _field_divergence_metrics(total_jac, probes.weights, probes.time_index)
    leading = _field_divergence_metrics(leading_jac, probes.weights, probes.time_index)
    oscillatory = _field_divergence_metrics(osc_jac, probes.weights, probes.time_index)

    total_axis_jac = independent_fd2_jacobian(total_field, probes.axis_points, probes.axis_times, step)
    lead_axis_jac = independent_fd2_jacobian(leading_field, probes.axis_points, probes.axis_times, step)
    total_axis_div = np.trace(total_axis_jac, axis1=1, axis2=2)
    lead_axis_div = np.trace(lead_axis_jac, axis1=1, axis2=2)
    osc_axis_div = total_axis_div - lead_axis_div

    total_velocity = _velocity(total_field, probes.points, probes.times)
    lead_velocity = _velocity(leading_field, probes.points, probes.times)
    osc_velocity = total_velocity - lead_velocity
    composition_delta = np.linalg.norm(osc_velocity, axis=1)

    idx_total = int(total["worst_index"])
    idx_osc = int(oscillatory["worst_index"])
    for metrics in (total, leading, oscillatory):
        metrics.pop("worst_index", None)

    return {
        "step": float(step),
        "total": total,
        "leading": leading,
        "oscillatory_increment": oscillatory,
        "axis_axis_near": {
            "total_sampled_max": float(np.max(np.abs(total_axis_div))),
            "leading_sampled_max": float(np.max(np.abs(lead_axis_div))),
            "oscillatory_increment_sampled_max": float(np.max(np.abs(osc_axis_div))),
        },
        "total_speed_rms": float(np.sqrt(np.mean(np.sum(total_velocity * total_velocity, axis=1)))),
        "oscillatory_speed_rms": float(np.sqrt(np.mean(composition_delta * composition_delta))),
        "worst_total_witness": {
            "point": [float(v) for v in probes.points[idx_total]],
            "time": float(probes.times[idx_total]),
            "divergence": float(total["worst_divergence"]),
            "normalized_divergence": float(total["worst_normalized_divergence"]),
        },
        "worst_oscillatory_witness": {
            "point": [float(v) for v in probes.points[idx_osc]],
            "time": float(probes.times[idx_osc]),
            "divergence_increment": float(oscillatory["worst_divergence"]),
            "normalized_divergence": float(oscillatory["worst_normalized_divergence"]),
        },
    }


class _XLinearMutation:
    def __init__(self, base: Any, epsilon: float = MUTATION_EPSILON):
        self.base = base
        self.epsilon = float(epsilon)

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        out = np.asarray(self.base.velocity(x, y, z, t), dtype=float).copy()
        out[..., 0] += self.epsilon * np.asarray(x, dtype=float)
        return out


def _semantic_snapshot_replay(field: Any) -> tuple[bool, str]:
    payload = field.semantic_payload()
    observed = str(field.semantic_sha256)
    if _sha256(payload) != observed:
        return False, observed
    with tempfile.TemporaryDirectory(prefix="kokuno-a4-partial-semantic-") as tmp:
        path = Path(tmp) / "partial_candidate_identity.json"
        path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
        rebound = json.loads(path.read_text(encoding="utf-8"))
    return bool(_sha256(rebound) == observed and rebound == payload), observed


def _stability_ok(medium: float, fine: float) -> bool:
    return bool(fine <= STABILITY_FACTOR * medium + STABILITY_FLOOR)


def _derive_audit_pass(receipt: Mapping[str, Any]) -> bool:
    resolutions = receipt["resolutions"]
    medium = resolutions[1]
    fine = resolutions[2]
    checks = receipt["checks"]

    scalar_metrics = (
        ("total", "sampled_max"),
        ("total", "pooled_weighted_rms"),
        ("leading", "sampled_max"),
        ("leading", "pooled_weighted_rms"),
        ("oscillatory_increment", "sampled_max"),
        ("oscillatory_increment", "pooled_weighted_rms"),
    )
    gates = all(float(fine[group][name]) <= DIVERGENCE_GATE for group, name in scalar_metrics)
    stable = all(
        _stability_ok(float(medium[group][name]), float(fine[group][name]))
        for group, name in scalar_metrics
    )
    axis = fine["axis_axis_near"]
    return bool(
        gates
        and stable
        and float(axis["total_sampled_max"]) <= DIVERGENCE_GATE
        and float(axis["leading_sampled_max"]) <= DIVERGENCE_GATE
        and float(axis["oscillatory_increment_sampled_max"]) <= DIVERGENCE_GATE
        and float(fine["total_speed_rms"]) >= NONTRIVIAL_TOTAL_SPEED_RMS
        and float(fine["oscillatory_speed_rms"]) >= NONTRIVIAL_OSC_SPEED_RMS
        and checks["upstream_source_blob_exact"]
        and checks["semantic_snapshot_replay"]
        and checks["upstream_truth_boundary_fail_closed"]
        and checks["velocity_mutation_detected"]
        and checks["offgrid_order_invariant"]
    )


def _upstream_truth_ok(field: Any) -> bool:
    truth = getattr(field, "truth_boundary", None)
    if not isinstance(truth, Mapping):
        return False
    required = {
        "leading_plus_oscillatory_partial_velocity_materialized": True,
        "velocity_beyond_Xh_materialized": False,
        "outer_global_leading_velocity_materialized": False,
        "agent3_correction_velocity_composed": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_velocity_pressure_forcing_api": False,
        "heldout_ns_residual_assessed": False,
        "pde_validated": False,
    }
    return all(truth.get(key) is expected for key, expected in required.items())


def materialize_receipt() -> dict[str, Any]:
    """Bind exact #970 identity, then audit only public Cartesian velocity surfaces."""
    field = _upstream_default_field()
    leading = field.leading_backend

    source_path = inspect.getsourcefile(type(field))
    if source_path is None:
        raise RuntimeError("cannot locate exact A2 #970 composition source")
    source_blob = _git_blob_sha1(source_path)
    semantic_ok, semantic_sha = _semantic_snapshot_replay(field)

    probes = make_heldout_probes()
    resolutions = [_metrics_for_step(field, leading, probes, step) for step in SPATIAL_STEPS]

    reverse = np.arange(probes.points.shape[0] - 1, -1, -1)
    reversed_probes = HeldoutProbeSet(
        points=probes.points[reverse],
        times=probes.times[reverse],
        weights=probes.weights[reverse],
        time_index=probes.time_index[reverse],
        axis_points=probes.axis_points,
        axis_times=probes.axis_times,
    )
    reversed_fine = _metrics_for_step(field, leading, reversed_probes, SPATIAL_STEPS[-1])
    fine = resolutions[-1]
    offgrid_order_invariant = bool(
        abs(reversed_fine["total"]["sampled_max"] - fine["total"]["sampled_max"]) <= 1.0e-13
        and abs(reversed_fine["total"]["pooled_weighted_rms"] - fine["total"]["pooled_weighted_rms"]) <= 1.0e-13
        and abs(
            reversed_fine["oscillatory_increment"]["pooled_weighted_rms"]
            - fine["oscillatory_increment"]["pooled_weighted_rms"]
        ) <= 1.0e-13
    )

    mutated = _XLinearMutation(field)
    mutated_fine = _metrics_for_step(mutated, leading, probes, SPATIAL_STEPS[-1])
    velocity_mutation_detected = bool(mutated_fine["total"]["sampled_max"] >= MUTATION_DETECTION_FLOOR)

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "upstream": {
            "pr": UPSTREAM_PR,
            "head": UPSTREAM_HEAD,
            "source_blob": source_blob,
            "expected_source_blob": UPSTREAM_SOURCE_BLOB,
            "semantic_sha256": semantic_sha,
            "agent1_pr": AGENT1_PR,
            "agent1_head": AGENT1_HEAD,
        },
        "protocol": {
            "seed": SEED,
            "spatial_steps": list(SPATIAL_STEPS),
            "times": list(TIMES),
            "radial_interval": list(RADIAL_INTERVAL),
            "axial_interval": list(AXIAL_INTERVAL),
            "points_per_time": POINTS_PER_TIME,
            "integration_probe_count": int(probes.points.shape[0]),
            "axis_axis_near_probe_count": int(probes.axis_points.shape[0]),
            "operator": "centered Cartesian FD2 on public total and public leading velocity only",
            "divergence_gate": DIVERGENCE_GATE,
            "stability_factor": STABILITY_FACTOR,
            "stability_floor": STABILITY_FLOOR,
            "nontrivial_total_speed_rms": NONTRIVIAL_TOTAL_SPEED_RMS,
            "nontrivial_osc_speed_rms": NONTRIVIAL_OSC_SPEED_RMS,
            "mutation_epsilon": MUTATION_EPSILON,
            "mutation_detection_floor": MUTATION_DETECTION_FLOOR,
            "volume_norm_note": (
                "physical cylindrical-shell Monte-Carlo integral over the frozen partial domain; "
                "not canonical 24/48/96 whole-domain admission"
            ),
        },
        "resolutions": resolutions,
        "mutation": {
            "kind": "total u_x += 1e-3*x while leading reference is unchanged",
            "finest_total_sampled_max": mutated_fine["total"]["sampled_max"],
            "finest_total_weighted_rms": mutated_fine["total"]["pooled_weighted_rms"],
        },
        "checks": {
            "upstream_source_blob_exact": source_blob == UPSTREAM_SOURCE_BLOB,
            "semantic_snapshot_replay": semantic_ok,
            "upstream_truth_boundary_fail_closed": _upstream_truth_ok(field),
            "velocity_mutation_detected": velocity_mutation_detected,
            "offgrid_order_invariant": offgrid_order_invariant,
        },
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
    }
    receipt["audit_pass"] = _derive_audit_pass(receipt)
    payload = copy.deepcopy(receipt)
    receipt["receipt_sha256"] = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
    return receipt


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA:
        raise ValueError("unexpected A4 leading+oscillatory divergence receipt schema")
    protocol = receipt.get("protocol", {})
    expected = {
        "seed": SEED,
        "spatial_steps": list(SPATIAL_STEPS),
        "times": list(TIMES),
        "radial_interval": list(RADIAL_INTERVAL),
        "axial_interval": list(AXIAL_INTERVAL),
        "points_per_time": POINTS_PER_TIME,
        "divergence_gate": DIVERGENCE_GATE,
        "stability_factor": STABILITY_FACTOR,
        "stability_floor": STABILITY_FLOOR,
        "nontrivial_total_speed_rms": NONTRIVIAL_TOTAL_SPEED_RMS,
        "nontrivial_osc_speed_rms": NONTRIVIAL_OSC_SPEED_RMS,
        "mutation_epsilon": MUTATION_EPSILON,
        "mutation_detection_floor": MUTATION_DETECTION_FLOOR,
    }
    for key, value in expected.items():
        if protocol.get(key) != value:
            raise ValueError(f"frozen protocol drift at {key}")

    upstream = receipt.get("upstream", {})
    if upstream.get("head") != UPSTREAM_HEAD or upstream.get("expected_source_blob") != UPSTREAM_SOURCE_BLOB:
        raise ValueError("upstream exact identity drift")
    if upstream.get("source_blob") != UPSTREAM_SOURCE_BLOB:
        raise ValueError("runtime upstream source blob drift")

    truth = receipt.get("truth_boundary", {})
    for key, expected_value in _TRUTH_BOUNDARY.items():
        if truth.get(key) is not expected_value:
            raise ValueError(f"truth-boundary drift at {key}")

    derived = _derive_audit_pass(receipt)
    if bool(receipt.get("audit_pass")) != derived:
        raise ValueError("audit_pass does not match frozen scientific gates")

    payload = copy.deepcopy(dict(receipt))
    observed_sha = payload.pop("receipt_sha256", None)
    expected_sha = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
    if observed_sha != expected_sha:
        raise ValueError("receipt checksum mismatch")

    if not derived:
        fine = receipt["resolutions"][-1]
        raise RuntimeError(
            "current partial leading+oscillatory divergence audit failed: "
            f"total_max={fine['total']['sampled_max']:.6e}, "
            f"total_rms={fine['total']['pooled_weighted_rms']:.6e}, "
            f"lead_max={fine['leading']['sampled_max']:.6e}, "
            f"osc_increment_max={fine['oscillatory_increment']['sampled_max']:.6e}, "
            f"axis_total={fine['axis_axis_near']['total_sampled_max']:.6e}"
        )


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify-receipt", type=Path)
    args = parser.parse_args()

    if args.verify_receipt is not None:
        payload = json.loads(args.verify_receipt.read_text(encoding="utf-8"))
        enforce_receipt(payload)
        print("A4 current leading+oscillatory divergence receipt passes frozen gates")
        return

    receipt = materialize_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False)
    if args.output is None:
        print(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
        print(text)


if __name__ == "__main__":
    _main()
