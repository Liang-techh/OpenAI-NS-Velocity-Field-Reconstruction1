"""Independent Cartesian divergence audit for the current Kokuno I1 leading field.

This Agent-4 validator is deliberately downstream of Agent-1 PR #1051.  It
serializes and reloads ``KokunoPA16CurrentCartesianI1FrozenGate`` and then
uses only its public ``velocity(x,y,z,t)`` surface for scientific numerical
differentiation.  No source-coordinate derivative helper, profile derivative,
construction tensor, pressure, forcing, or candidate-side residual routine is
used by the independent derivative path.

The audit is scoped to divergence consistency of the current *leading-only*
field through the reserved I1 interval.  It is not a complete Navier--Stokes
momentum audit and cannot promote ``pde_validated`` without the later global
candidate, matched pressure, preregistered restricted forcing, correction
velocity, and canonical whole-domain protocol.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
import json
import math
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_i1_frozen_gate import (
    I1_CLOSURE_TOLERANCE,
    KokunoPA16CurrentCartesianI1FrozenGate,
)

SCHEMA = "kokuno-a4-current-i1-leading-divergence-audit-v1"
UPSTREAM_PR = 1051
UPSTREAM_HEAD = "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5"
SEED = 9173801
SPATIAL_STEPS = (0.02, 0.01, 0.005)
TIMES = (0.31, 0.47, 0.63, 0.71)
POINTS_PER_TIME_PRE_I1 = 4
POINTS_PER_TIME_PER_I1_ZONE = 4
I1_ZONE_FRACTIONS = ((0.16, 0.30), (0.39, 0.55), (0.64, 0.78))
PRE_I1_LOG_OFFSETS = (0.08, 0.34)
ETA_INTERVAL = (-0.45, 0.45)
DIVERGENCE_GATE = 1.0e-5
STABILITY_FACTOR = 1.25
STABILITY_FLOOR = 2.0e-8
NONTRIVIAL_SPEED_RMS = 1.0e-10
MUTATION_EPSILON = 1.0e-3
MUTATION_DETECTION_FLOOR = 5.0e-4
SAVE_RELOAD_VELOCITY_TOLERANCE = 2.0e-12

_TRUTH_BOUNDARY = {
    "current_i1_leading_velocity_consumed": True,
    "candidate_save_reload_required": True,
    "independent_cartesian_fd_operator_used": True,
    "leading_only_divergence_scoped_assessed": True,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "canonical_24_48_96_volume_admission_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _source_to_cartesian(
    field: Any,
    X: np.ndarray,
    eta: np.ndarray,
    theta: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    """Independent public similarity-chart formula used only to place probes."""
    X, eta, theta, t = np.broadcast_arrays(
        np.asarray(X, dtype=float),
        np.asarray(eta, dtype=float),
        np.asarray(theta, dtype=float),
        np.asarray(t, dtype=float),
    )
    if np.any(~np.isfinite(X)) or np.any(X < 0.0):
        raise ValueError("X probes must be finite and nonnegative")
    if np.any(~np.isfinite(eta)) or np.any(np.abs(eta) >= 1.0):
        raise ValueError("eta probes must remain in the similarity chart")
    tau = 1.0 - t
    if np.any(~np.isfinite(tau)) or np.any(tau <= 0.0):
        raise ValueError("time probes must satisfy t<1")
    q = tau / (1.0 - eta * eta)
    r = np.sqrt(2.0 * q * X)
    z = np.power(q, float(field.D)) * eta
    return np.stack((r * np.cos(theta), r * np.sin(theta), z), axis=-1)


def _physical_volume_jacobian(field: Any, eta: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Independent dV/(dX d eta d theta) for the fixed-time chart."""
    eta, t = np.broadcast_arrays(np.asarray(eta, dtype=float), np.asarray(t, dtype=float))
    D = float(field.D)
    q = (1.0 - t) / (1.0 - eta * eta)
    factor = 1.0 + 2.0 * D * eta * eta / (1.0 - eta * eta)
    out = np.power(q, D + 1.0) * factor
    if np.any(~np.isfinite(out)) or np.any(out <= 0.0):
        raise RuntimeError("physical-volume Jacobian must be positive and finite")
    return out


def _latin_hypercube(rng: np.random.Generator, count: int) -> np.ndarray:
    values = (np.arange(count, dtype=float) + rng.random(count)) / float(count)
    rng.shuffle(values)
    return values


@dataclass(frozen=True)
class HeldoutProbeSet:
    points: np.ndarray
    times: np.ndarray
    weights: np.ndarray
    region: np.ndarray
    time_index: np.ndarray
    seam_points: np.ndarray
    seam_times: np.ndarray
    axis_points: np.ndarray
    axis_times: np.ndarray


def _append_log_region(
    *,
    field: Any,
    rng: np.random.Generator,
    time: float,
    time_index: int,
    label: str,
    log_lo: float,
    log_hi: float,
    count: int,
    points: list[np.ndarray],
    times: list[np.ndarray],
    weights: list[np.ndarray],
    regions: list[np.ndarray],
    time_indices: list[np.ndarray],
) -> None:
    if not log_hi > log_lo:
        raise ValueError("log-X region must have positive width")
    ux = _latin_hypercube(rng, count)
    ue = _latin_hypercube(rng, count)
    ut = _latin_hypercube(rng, count)
    logs = log_lo + (log_hi - log_lo) * ux
    X = np.exp(logs)
    e0, e1 = ETA_INTERVAL
    eta = e0 + (e1 - e0) * ue
    theta = 2.0 * math.pi * ut
    t = np.full(count, time, dtype=float)
    xyz = _source_to_cartesian(field, X, eta, theta, t)
    jac = _physical_volume_jacobian(field, eta, t)
    # dX = X d(log X).  Each Monte-Carlo weight estimates physical volume
    # over this disjoint log-X slab.
    cell_measure = (log_hi - log_lo) * (e1 - e0) * (2.0 * math.pi)
    w = cell_measure * X * jac / float(count)
    points.append(xyz)
    times.append(t)
    weights.append(w)
    regions.append(np.full(count, label, dtype=object))
    time_indices.append(np.full(count, time_index, dtype=int))


def make_heldout_probes(field: Any) -> HeldoutProbeSet:
    """Construct the preregistered fresh pre-I1/I1, seam and axis probes."""
    rng = np.random.default_rng(SEED)
    s0 = float(field.log_X_I1_start)
    s1 = float(field.log_X_I1_end)
    span = s1 - s0
    if not span > 0.0:
        raise RuntimeError("I1 log interval has nonpositive width")

    points: list[np.ndarray] = []
    times: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    regions: list[np.ndarray] = []
    time_indices: list[np.ndarray] = []

    for tid, time in enumerate(TIMES):
        pre_lo = s0 - PRE_I1_LOG_OFFSETS[1]
        pre_hi = s0 - PRE_I1_LOG_OFFSETS[0]
        _append_log_region(
            field=field,
            rng=rng,
            time=time,
            time_index=tid,
            label="pre_i1_control",
            log_lo=pre_lo,
            log_hi=pre_hi,
            count=POINTS_PER_TIME_PRE_I1,
            points=points,
            times=times,
            weights=weights,
            regions=regions,
            time_indices=time_indices,
        )
        for zid, (a, b) in enumerate(I1_ZONE_FRACTIONS):
            _append_log_region(
                field=field,
                rng=rng,
                time=time,
                time_index=tid,
                label=f"i1_zone_{zid}",
                log_lo=s0 + a * span,
                log_hi=s0 + b * span,
                count=POINTS_PER_TIME_PER_I1_ZONE,
                points=points,
                times=times,
                weights=weights,
                regions=regions,
                time_indices=time_indices,
            )

    # Probe both sides of the exact I1 handoff.  The offsets are fractions of
    # the reserved I1 log-width, not numerical differentiation steps.
    seam_xyz: list[np.ndarray] = []
    seam_t: list[np.ndarray] = []
    for time in (0.47, 0.63):
        offsets = np.asarray((-0.025, -0.008, 0.008, 0.025), dtype=float) * span
        for eta_value in (-0.22, 0.22):
            X = np.exp(s0 + offsets)
            eta = np.full(offsets.shape, eta_value, dtype=float)
            theta = np.asarray((0.37, 1.11, 2.03, 2.81), dtype=float)
            tt = np.full(offsets.shape, time, dtype=float)
            seam_xyz.append(_source_to_cartesian(field, X, eta, theta, tt))
            seam_t.append(tt)

    axis_X = np.asarray((0.0, 0.0, 1.0e-10, 1.0e-8, 1.0e-6), dtype=float)
    axis_eta = np.asarray((0.0, 0.40, -0.40, 0.22, 0.0), dtype=float)
    axis_theta = np.asarray((0.0, 0.9, 1.4, 2.0, 2.7), dtype=float)
    axis_times = np.asarray((0.50, 0.47, 0.63, 0.31, 0.71), dtype=float)
    axis_points = _source_to_cartesian(field, axis_X, axis_eta, axis_theta, axis_times)

    return HeldoutProbeSet(
        points=np.concatenate(points, axis=0),
        times=np.concatenate(times, axis=0),
        weights=np.concatenate(weights, axis=0),
        region=np.concatenate(regions, axis=0),
        time_index=np.concatenate(time_indices, axis=0),
        seam_points=np.concatenate(seam_xyz, axis=0),
        seam_times=np.concatenate(seam_t, axis=0),
        axis_points=axis_points,
        axis_times=axis_times,
    )


def _velocity(field: Any, points: np.ndarray, times: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    times = np.asarray(times, dtype=float)
    out = np.asarray(
        field.velocity(points[:, 0], points[:, 1], points[:, 2], times),
        dtype=float,
    )
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
    """Reconstruct the full Cartesian Jacobian from public velocity only."""
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
        jac[:, :, axis] = (
            _velocity(field, plus, times) - _velocity(field, minus, times)
        ) / (2.0 * h)
    if np.any(~np.isfinite(jac)):
        raise RuntimeError("independent Cartesian Jacobian became non-finite")
    return jac


def _divergence(field: Any, points: np.ndarray, times: np.ndarray, step: float) -> tuple[np.ndarray, np.ndarray]:
    jac = independent_fd2_jacobian(field, points, times, step)
    return np.trace(jac, axis1=1, axis2=2), jac


def _weighted_metrics(div: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    volume = float(np.sum(weights))
    if not volume > 0.0:
        raise RuntimeError("probe weights must have positive total volume")
    l2 = float(np.sqrt(np.sum(weights * div * div)))
    rms = float(l2 / math.sqrt(volume))
    return volume, l2, rms


def _metrics_for_step(field: Any, probes: HeldoutProbeSet, step: float) -> dict[str, Any]:
    div, jac = _divergence(field, probes.points, probes.times, step)
    frob = np.linalg.norm(jac, axis=(1, 2))
    normalized = np.abs(div) / np.maximum(1.0, frob)
    speed = np.linalg.norm(_velocity(field, probes.points, probes.times), axis=1)
    volume, l2, rms = _weighted_metrics(div, probes.weights)

    per_region: dict[str, Any] = {}
    for label in sorted({str(v) for v in probes.region.tolist()}):
        mask = probes.region == label
        rv, rl2, rrms = _weighted_metrics(div[mask], probes.weights[mask])
        per_region[label] = {
            "estimated_volume": rv,
            "sampled_max": float(np.max(np.abs(div[mask]))),
            "volume_l2_estimate": rl2,
            "weighted_rms": rrms,
        }

    per_time: list[dict[str, Any]] = []
    for tid, time in enumerate(TIMES):
        mask = probes.time_index == tid
        tv, tl2, trms = _weighted_metrics(div[mask], probes.weights[mask])
        per_time.append(
            {
                "time": float(time),
                "estimated_volume": tv,
                "sampled_max": float(np.max(np.abs(div[mask]))),
                "volume_l2_estimate": tl2,
                "weighted_rms": trms,
            }
        )

    seam_div, _ = _divergence(field, probes.seam_points, probes.seam_times, step)
    axis_div, _ = _divergence(field, probes.axis_points, probes.axis_times, step)
    worst = int(np.argmax(np.abs(div)))
    return {
        "step": float(step),
        "sampled_max": float(np.max(np.abs(div))),
        "estimated_volume": volume,
        "pooled_volume_l2_estimate": l2,
        "pooled_weighted_rms": rms,
        "normalized_sampled_max": float(np.max(normalized)),
        "normalized_weighted_rms": float(
            np.sqrt(np.sum(probes.weights * normalized * normalized) / np.sum(probes.weights))
        ),
        "i1_entry_seam_sampled_max": float(np.max(np.abs(seam_div))),
        "axis_axis_near_sampled_max": float(np.max(np.abs(axis_div))),
        "speed_rms": float(np.sqrt(np.mean(speed * speed))),
        "per_region": per_region,
        "per_time": per_time,
        "worst_witness": {
            "point": [float(v) for v in probes.points[worst]],
            "time": float(probes.times[worst]),
            "region": str(probes.region[worst]),
            "divergence": float(div[worst]),
            "jacobian_frobenius": float(frob[worst]),
            "normalized_divergence": float(normalized[worst]),
        },
    }


class _XLinearMutation:
    def __init__(self, base: Any, epsilon: float = MUTATION_EPSILON):
        self.base = base
        self.epsilon = float(epsilon)
        self.D = base.D

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        out = np.asarray(self.base.velocity(x, y, z, t), dtype=float).copy()
        out[..., 0] += self.epsilon * np.asarray(x, dtype=float)
        return out


def _stability_ok(medium: float, fine: float) -> bool:
    return bool(fine <= STABILITY_FACTOR * medium + STABILITY_FLOOR)


def _gate_mutation_rejected(field: KokunoPA16CurrentCartesianI1FrozenGate) -> bool:
    payload = copy.deepcopy(field.configuration())
    payload["frozen_i1_closure_gate"]["closure_tolerance"] = 1.0e-5
    try:
        KokunoPA16CurrentCartesianI1FrozenGate.from_configuration(payload)
    except ValueError:
        return True
    return False


def _nested_gate_mutation_rejected(field: KokunoPA16CurrentCartesianI1FrozenGate) -> bool:
    payload = copy.deepcopy(field.configuration())
    payload["parent_i1_candidate"]["i1_repair"]["closure_tolerance"] = 1.0e-5
    try:
        KokunoPA16CurrentCartesianI1FrozenGate.from_configuration(payload)
    except ValueError:
        return True
    return False


def _post_i1_fails_closed(field: Any) -> bool:
    X = np.asarray([1.001 * float(field.X_I1_end)], dtype=float)
    eta = np.asarray([0.0], dtype=float)
    theta = np.asarray([0.71], dtype=float)
    t = np.asarray([0.47], dtype=float)
    point = _source_to_cartesian(field, X, eta, theta, t)
    try:
        _velocity(field, point, t)
    except (ValueError, RuntimeError):
        return True
    return False


def _derive_audit_pass(receipt: Mapping[str, Any]) -> bool:
    resolutions = receipt["resolutions"]
    medium = resolutions[1]
    fine = resolutions[2]
    checks = receipt["checks"]
    return bool(
        fine["sampled_max"] <= DIVERGENCE_GATE
        and fine["pooled_weighted_rms"] <= DIVERGENCE_GATE
        and fine["i1_entry_seam_sampled_max"] <= DIVERGENCE_GATE
        and fine["axis_axis_near_sampled_max"] <= DIVERGENCE_GATE
        and fine["speed_rms"] >= NONTRIVIAL_SPEED_RMS
        and _stability_ok(medium["sampled_max"], fine["sampled_max"])
        and _stability_ok(medium["pooled_weighted_rms"], fine["pooled_weighted_rms"])
        and _stability_ok(
            medium["i1_entry_seam_sampled_max"], fine["i1_entry_seam_sampled_max"]
        )
        and _stability_ok(
            medium["axis_axis_near_sampled_max"], fine["axis_axis_near_sampled_max"]
        )
        and checks["save_reload_semantic_identity"]
        and checks["save_reload_velocity_replay_max_abs"] <= SAVE_RELOAD_VELOCITY_TOLERANCE
        and checks["velocity_mutation_detected"]
        and checks["offgrid_order_invariant"]
        and checks["frozen_gate_mutation_rejected"]
        and checks["nested_parent_gate_mutation_rejected"]
        and checks["post_i1_fails_closed"]
    )


def materialize_receipt() -> dict[str, Any]:
    """Save/reload exact #1051 and execute the preregistered independent audit."""
    original = KokunoPA16CurrentCartesianI1FrozenGate()
    with tempfile.TemporaryDirectory(prefix="kokuno-a4-current-i1-") as tmp:
        path = Path(tmp) / "candidate.json"
        original.save_configuration(path)
        field = KokunoPA16CurrentCartesianI1FrozenGate.load_configuration(path)

    probes = make_heldout_probes(field)
    resolutions = [_metrics_for_step(field, probes, h) for h in SPATIAL_STEPS]

    replay_count = min(8, probes.points.shape[0])
    replay_original = _velocity(original, probes.points[:replay_count], probes.times[:replay_count])
    replay_loaded = _velocity(field, probes.points[:replay_count], probes.times[:replay_count])
    replay_error = float(np.max(np.abs(replay_original - replay_loaded)))

    reverse = np.arange(probes.points.shape[0] - 1, -1, -1)
    reversed_probes = HeldoutProbeSet(
        points=probes.points[reverse],
        times=probes.times[reverse],
        weights=probes.weights[reverse],
        region=probes.region[reverse],
        time_index=probes.time_index[reverse],
        seam_points=probes.seam_points,
        seam_times=probes.seam_times,
        axis_points=probes.axis_points,
        axis_times=probes.axis_times,
    )
    reversed_fine = _metrics_for_step(field, reversed_probes, SPATIAL_STEPS[-1])
    fine = resolutions[-1]
    order_invariant = bool(
        abs(reversed_fine["sampled_max"] - fine["sampled_max"]) <= 1.0e-13
        and abs(reversed_fine["pooled_weighted_rms"] - fine["pooled_weighted_rms"]) <= 1.0e-13
        and abs(reversed_fine["pooled_volume_l2_estimate"] - fine["pooled_volume_l2_estimate"]) <= 1.0e-13
    )

    base_div, _ = _divergence(field, probes.points, probes.times, SPATIAL_STEPS[-1])
    mutation = _XLinearMutation(field)
    mutated_div, _ = _divergence(mutation, probes.points, probes.times, SPATIAL_STEPS[-1])
    mutation_signal = float(np.max(np.abs(mutated_div - base_div)))

    checks = {
        "save_reload_semantic_identity": field.semantic_sha256 == original.semantic_sha256,
        "save_reload_velocity_replay_max_abs": replay_error,
        "offgrid_order_invariant": order_invariant,
        "velocity_mutation_signal": mutation_signal,
        "velocity_mutation_detected": mutation_signal >= MUTATION_DETECTION_FLOOR,
        "frozen_gate_mutation_rejected": _gate_mutation_rejected(field),
        "nested_parent_gate_mutation_rejected": _nested_gate_mutation_rejected(field),
        "post_i1_fails_closed": _post_i1_fails_closed(field),
        "runtime_closure_gate_is_frozen": float(field.closure_tolerance) == I1_CLOSURE_TOLERANCE,
    }

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "upstream_pr": UPSTREAM_PR,
        "upstream_head": UPSTREAM_HEAD,
        "candidate_semantic_sha256": field.semantic_sha256,
        "protocol": {
            "seed": SEED,
            "spatial_steps": list(SPATIAL_STEPS),
            "times": list(TIMES),
            "pre_i1_log_offsets": list(PRE_I1_LOG_OFFSETS),
            "i1_zone_fractions": [list(v) for v in I1_ZONE_FRACTIONS],
            "eta_interval": list(ETA_INTERVAL),
            "integration_probe_count": int(probes.points.shape[0]),
            "i1_entry_seam_probe_count": int(probes.seam_points.shape[0]),
            "axis_axis_near_probe_count": int(probes.axis_points.shape[0]),
            "derivative_operator": "centered Cartesian FD2 from public velocity only",
            "physical_volume_weighting": "independent chart Jacobian with dX=X*dlogX",
            "frozen_i1_closure_tolerance": I1_CLOSURE_TOLERANCE,
        },
        "gates": {
            "divergence_finest_sampled_max": DIVERGENCE_GATE,
            "divergence_finest_weighted_rms": DIVERGENCE_GATE,
            "i1_entry_seam_sampled_max": DIVERGENCE_GATE,
            "axis_axis_near_sampled_max": DIVERGENCE_GATE,
            "nontrivial_speed_rms_min": NONTRIVIAL_SPEED_RMS,
            "stability_factor": STABILITY_FACTOR,
            "stability_floor": STABILITY_FLOOR,
            "mutation_detection_floor": MUTATION_DETECTION_FLOOR,
            "save_reload_velocity_tolerance": SAVE_RELOAD_VELOCITY_TOLERANCE,
            "final_project_momentum_gate_unchanged": 1.0e-3,
            "final_project_divergence_gate_unchanged": 1.0e-5,
        },
        "resolutions": resolutions,
        "checks": checks,
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        "limitations": [
            "scoped leading-only divergence evidence through current I1, not a complete NS residual",
            "the current A2 oscillatory runtime has not yet been recomposed on exact #1051 identity",
            "no current-lineage Cartesian correction velocity is consumed",
            "matched Cartesian pressure and preregistered restricted forcing remain unavailable",
            "canonical whole-domain 24/48/96 admission and same-protocol ST006 comparison remain pending",
        ],
    }
    receipt["audit_pass"] = _derive_audit_pass(receipt)
    return receipt


def enforce_preregistered_gates(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA:
        raise AssertionError("unexpected A4 current-I1 audit schema")
    if receipt.get("upstream_head") != UPSTREAM_HEAD:
        raise AssertionError("current-I1 audit is detached from exact #1051 head")
    if not _derive_audit_pass(receipt):
        fine = receipt["resolutions"][-1]
        raise AssertionError(
            "preregistered current-I1 divergence audit failed: "
            f"max={fine['sampled_max']:.6e}, "
            f"rms={fine['pooled_weighted_rms']:.6e}, "
            f"seam={fine['i1_entry_seam_sampled_max']:.6e}, "
            f"axis={fine['axis_axis_near_sampled_max']:.6e}"
        )


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--enforce-existing", action="store_true")
    args = parser.parse_args()
    path = Path(args.output)
    if args.enforce_existing:
        receipt = json.loads(path.read_text(encoding="utf-8"))
        enforce_preregistered_gates(receipt)
        return
    receipt = materialize_receipt()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    _main()
