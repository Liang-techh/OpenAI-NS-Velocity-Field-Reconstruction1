"""Independent Cartesian divergence audit for the current Kokuno leading field.

This Agent-4 validator consumes only the serialized public candidate artifact
introduced by Agent 1 PR #965.  The scientific derivative path treats
``velocity(x,y,z,t)`` as a black box and reconstructs its Cartesian Jacobian
with centered second-order finite differences on a frozen physical step
ladder.  It does not call candidate derivative helpers, source-coordinate
radial derivatives, training tensors, pressure, forcing, or residual code.

The result is intentionally scoped: it can assess divergence consistency of
the current *leading-only* Cartesian field on the currently materialized
region, but it cannot assess a complete Navier--Stokes momentum residual.
Matched pressure, preregistered restricted forcing, oscillatory composition,
correction velocity, canonical 24/48/96 whole-domain quadrature, and final PDE
admission remain separate downstream requirements.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_leading_velocity import (
    KokunoPA16CurrentCartesianLeadingVelocity,
)


SCHEMA = "kokuno-a4-current-cartesian-leading-divergence-audit-v1"
UPSTREAM_PR = 965
UPSTREAM_HEAD = "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"
UPSTREAM_SOURCE_BLOB = "ca8b80b0451be1a8f31deaf620f542d8c2e92c0c"
SEED = 9173661
SPATIAL_STEPS = (0.02, 0.01, 0.005)
TIMES = (0.31, 0.47, 0.63, 0.71)
X_INTERVAL = (0.02, 1.20)
ETA_INTERVAL = (-0.65, 0.65)
POINTS_PER_TIME = 8
DIVERGENCE_GATE = 1.0e-5
STABILITY_FACTOR = 1.25
STABILITY_FLOOR = 2.0e-8
NONTRIVIAL_SPEED_RMS = 1.0e-10
MUTATION_EPSILON = 1.0e-3
MUTATION_DETECTION_FLOOR = 5.0e-4

_TRUTH_BOUNDARY = {
    "current_cartesian_leading_velocity_consumed": True,
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


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _source_to_cartesian(
    field: Any,
    X: np.ndarray,
    eta: np.ndarray,
    theta: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    """Independent public coordinate formula used only to place audit probes."""
    X, eta, theta, t = np.broadcast_arrays(
        np.asarray(X, dtype=float),
        np.asarray(eta, dtype=float),
        np.asarray(theta, dtype=float),
        np.asarray(t, dtype=float),
    )
    if np.any(~np.isfinite(X)) or np.any(~np.isfinite(eta)):
        raise ValueError("probe coordinates must be finite")
    tau = 1.0 - t
    if np.any(tau <= 0.0) or np.any(np.abs(eta) >= 1.0):
        raise ValueError("probe coordinates leave the public similarity chart")
    q = tau / (1.0 - eta * eta)
    r = np.sqrt(2.0 * q * X)
    z = np.power(q, float(field.D)) * eta
    return np.stack((r * np.cos(theta), r * np.sin(theta), z), axis=-1)


def _physical_volume_jacobian(field: Any, eta: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Jacobian dV/(dX d eta d theta) for the public fixed-time map."""
    eta, t = np.broadcast_arrays(np.asarray(eta, dtype=float), np.asarray(t, dtype=float))
    D = float(field.D)
    q = (1.0 - t) / (1.0 - eta * eta)
    factor = 1.0 + 2.0 * D * eta * eta / (1.0 - eta * eta)
    out = np.power(q, D + 1.0) * factor
    if np.any(~np.isfinite(out)) or np.any(out <= 0.0):
        raise RuntimeError("physical volume Jacobian must be positive and finite")
    return out


@dataclass(frozen=True)
class HeldoutProbeSet:
    points: np.ndarray
    times: np.ndarray
    weights: np.ndarray
    time_index: np.ndarray
    axis_points: np.ndarray
    axis_times: np.ndarray


def _latin_hypercube(rng: np.random.Generator, count: int) -> np.ndarray:
    base = (np.arange(count, dtype=float) + rng.random(count)) / float(count)
    rng.shuffle(base)
    return base


def make_heldout_probes(field: Any) -> HeldoutProbeSet:
    """Create the frozen held-out off-grid and axis/near-axis probe set."""
    rng = np.random.default_rng(SEED)
    all_points: list[np.ndarray] = []
    all_times: list[np.ndarray] = []
    all_weights: list[np.ndarray] = []
    all_time_index: list[np.ndarray] = []
    x0, x1 = X_INTERVAL
    e0, e1 = ETA_INTERVAL
    cell_measure = (x1 - x0) * (e1 - e0) * (2.0 * math.pi)

    for time_id, time in enumerate(TIMES):
        ux = _latin_hypercube(rng, POINTS_PER_TIME)
        ue = _latin_hypercube(rng, POINTS_PER_TIME)
        ut = _latin_hypercube(rng, POINTS_PER_TIME)
        X = x0 + (x1 - x0) * ux
        eta = e0 + (e1 - e0) * ue
        theta = 2.0 * math.pi * ut
        t = np.full(POINTS_PER_TIME, time, dtype=float)
        points = _source_to_cartesian(field, X, eta, theta, t)
        jac = _physical_volume_jacobian(field, eta, t)
        weights = cell_measure * jac / float(POINTS_PER_TIME)
        all_points.append(points)
        all_times.append(t)
        all_weights.append(weights)
        all_time_index.append(np.full(POINTS_PER_TIME, time_id, dtype=int))

    axis_X = np.asarray([0.0, 0.0, 1.0e-10, 1.0e-8, 1.0e-6], dtype=float)
    axis_eta = np.asarray([0.0, 0.45, -0.45, 0.25, 0.0], dtype=float)
    axis_theta = np.asarray([0.0, 1.0, 0.7, 1.3, 2.1], dtype=float)
    axis_times = np.asarray([0.50, 0.47, 0.63, 0.31, 0.71], dtype=float)
    axis_points = _source_to_cartesian(field, axis_X, axis_eta, axis_theta, axis_times)

    return HeldoutProbeSet(
        points=np.concatenate(all_points, axis=0),
        times=np.concatenate(all_times, axis=0),
        weights=np.concatenate(all_weights, axis=0),
        time_index=np.concatenate(all_time_index, axis=0),
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
        v_plus = _velocity(field, plus, times)
        v_minus = _velocity(field, minus, times)
        jac[:, :, axis] = (v_plus - v_minus) / (2.0 * h)
    if np.any(~np.isfinite(jac)):
        raise RuntimeError("independent Cartesian Jacobian became non-finite")
    return jac


def _metrics_for_step(
    field: Any,
    probes: HeldoutProbeSet,
    step: float,
) -> dict[str, Any]:
    jac = independent_fd2_jacobian(field, probes.points, probes.times, step)
    div = np.trace(jac, axis1=1, axis2=2)
    frob = np.linalg.norm(jac, axis=(1, 2))
    normalized = np.abs(div) / np.maximum(1.0, frob)
    speed = np.linalg.norm(_velocity(field, probes.points, probes.times), axis=1)

    per_time: list[dict[str, Any]] = []
    for time_id, time in enumerate(TIMES):
        mask = probes.time_index == time_id
        w = probes.weights[mask]
        d = div[mask]
        volume = float(np.sum(w))
        integral_l2 = float(np.sqrt(np.sum(w * d * d)))
        weighted_rms = float(np.sqrt(np.sum(w * d * d) / volume))
        per_time.append(
            {
                "time": float(time),
                "estimated_volume": volume,
                "volume_l2_estimate": integral_l2,
                "weighted_rms": weighted_rms,
                "sampled_max": float(np.max(np.abs(d))),
            }
        )

    axis_jac = independent_fd2_jacobian(field, probes.axis_points, probes.axis_times, step)
    axis_div = np.trace(axis_jac, axis1=1, axis2=2)
    worst_index = int(np.argmax(np.abs(div)))
    return {
        "step": float(step),
        "sampled_max": float(np.max(np.abs(div))),
        "pooled_weighted_rms": float(
            np.sqrt(np.sum(probes.weights * div * div) / np.sum(probes.weights))
        ),
        "pooled_volume_l2_estimate": float(np.sqrt(np.sum(probes.weights * div * div))),
        "normalized_sampled_max": float(np.max(normalized)),
        "normalized_weighted_rms": float(
            np.sqrt(np.sum(probes.weights * normalized * normalized) / np.sum(probes.weights))
        ),
        "axis_axis_near_sampled_max": float(np.max(np.abs(axis_div))),
        "speed_rms": float(np.sqrt(np.mean(speed * speed))),
        "per_time": per_time,
        "worst_witness": {
            "point": [float(v) for v in probes.points[worst_index]],
            "time": float(probes.times[worst_index]),
            "divergence": float(div[worst_index]),
            "jacobian_frobenius": float(frob[worst_index]),
            "normalized_divergence": float(normalized[worst_index]),
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


def _identity_mutation_changes_semantic(field: KokunoPA16CurrentCartesianLeadingVelocity) -> bool:
    payload = copy.deepcopy(field.configuration())
    payload["eta_fd_step"] = float(payload["eta_fd_step"]) * 1.001
    mutated = KokunoPA16CurrentCartesianLeadingVelocity.from_configuration(payload)
    return mutated.semantic_sha256 != field.semantic_sha256


def _stability_ok(medium: float, fine: float) -> bool:
    return bool(fine <= STABILITY_FACTOR * medium + STABILITY_FLOOR)


def _derive_audit_pass(receipt: Mapping[str, Any]) -> bool:
    resolutions = receipt["resolutions"]
    medium = resolutions[1]
    fine = resolutions[2]
    checks = receipt["checks"]
    return bool(
        fine["sampled_max"] <= DIVERGENCE_GATE
        and fine["pooled_weighted_rms"] <= DIVERGENCE_GATE
        and fine["axis_axis_near_sampled_max"] <= DIVERGENCE_GATE
        and fine["speed_rms"] >= NONTRIVIAL_SPEED_RMS
        and _stability_ok(medium["sampled_max"], fine["sampled_max"])
        and _stability_ok(medium["pooled_weighted_rms"], fine["pooled_weighted_rms"])
        and _stability_ok(
            medium["axis_axis_near_sampled_max"], fine["axis_axis_near_sampled_max"]
        )
        and checks["save_reload_semantic_identity"]
        and checks["config_mutation_changes_semantic_identity"]
        and checks["velocity_mutation_detected"]
        and checks["offgrid_order_invariant"]
    )


def materialize_receipt() -> dict[str, Any]:
    """Save/reload the real #965 artifact and run the frozen independent audit."""
    original = KokunoPA16CurrentCartesianLeadingVelocity()
    with tempfile.TemporaryDirectory(prefix="kokuno-a4-cart-leading-") as tmp:
        path = Path(tmp) / "candidate.json"
        original.save_configuration(path)
        field = KokunoPA16CurrentCartesianLeadingVelocity.load_configuration(path)

    probes = make_heldout_probes(field)
    resolutions = [_metrics_for_step(field, probes, h) for h in SPATIAL_STEPS]

    # Order-invariance guard: the same finest calculation after a deterministic
    # reversal must reproduce the original pointwise-independent aggregate.
    reverse = np.arange(probes.points.shape[0] - 1, -1, -1)
    reversed_probes = HeldoutProbeSet(
        points=probes.points[reverse],
        times=probes.times[reverse],
        weights=probes.weights[reverse],
        time_index=probes.time_index[reverse],
        axis_points=probes.axis_points,
        axis_times=probes.axis_times,
    )
    reversed_fine = _metrics_for_step(field, reversed_probes, SPATIAL_STEPS[-1])
    offgrid_order_invariant = bool(
        abs(reversed_fine["sampled_max"] - resolutions[-1]["sampled_max"]) <= 1.0e-13
        and abs(
            reversed_fine["pooled_weighted_rms"] - resolutions[-1]["pooled_weighted_rms"]
        )
        <= 1.0e-13
    )

    mutated = _XLinearMutation(field)
    mutated_fine = _metrics_for_step(mutated, probes, SPATIAL_STEPS[-1])
    velocity_mutation_detected = bool(
        mutated_fine["sampled_max"] >= MUTATION_DETECTION_FLOOR
    )

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "upstream": {
            "pr": UPSTREAM_PR,
            "head": UPSTREAM_HEAD,
            "source_blob": UPSTREAM_SOURCE_BLOB,
            "candidate_semantic_sha256": field.semantic_sha256,
        },
        "protocol": {
            "seed": SEED,
            "spatial_steps": list(SPATIAL_STEPS),
            "times": list(TIMES),
            "X_interval": list(X_INTERVAL),
            "eta_interval": list(ETA_INTERVAL),
            "points_per_time": POINTS_PER_TIME,
            "integration_probe_count": int(probes.points.shape[0]),
            "axis_axis_near_probe_count": int(probes.axis_points.shape[0]),
            "operator": "centered Cartesian FD2 on public velocity only",
            "divergence_gate": DIVERGENCE_GATE,
            "stability_factor": STABILITY_FACTOR,
            "stability_floor": STABILITY_FLOOR,
            "nontrivial_speed_rms": NONTRIVIAL_SPEED_RMS,
            "mutation_epsilon": MUTATION_EPSILON,
            "mutation_detection_floor": MUTATION_DETECTION_FLOOR,
            "volume_norm_note": (
                "Monte-Carlo physical-volume estimate over the frozen source cell; "
                "not canonical 24/48/96 whole-domain admission"
            ),
        },
        "resolutions": resolutions,
        "mutation": {
            "kind": "u_x += 1e-3*x",
            "finest_sampled_max": mutated_fine["sampled_max"],
            "finest_weighted_rms": mutated_fine["pooled_weighted_rms"],
        },
        "checks": {
            "save_reload_semantic_identity": field.semantic_sha256 == original.semantic_sha256,
            "config_mutation_changes_semantic_identity": _identity_mutation_changes_semantic(field),
            "velocity_mutation_detected": velocity_mutation_detected,
            "offgrid_order_invariant": offgrid_order_invariant,
        },
        "refinement": {
            "coarse_to_medium_sampled_max_ratio": float(
                resolutions[0]["sampled_max"]
                / max(resolutions[1]["sampled_max"], STABILITY_FLOOR)
            ),
            "medium_to_fine_sampled_max_ratio": float(
                resolutions[1]["sampled_max"]
                / max(resolutions[2]["sampled_max"], STABILITY_FLOOR)
            ),
            "coarse_to_medium_weighted_rms_ratio": float(
                resolutions[0]["pooled_weighted_rms"]
                / max(resolutions[1]["pooled_weighted_rms"], STABILITY_FLOOR)
            ),
            "medium_to_fine_weighted_rms_ratio": float(
                resolutions[1]["pooled_weighted_rms"]
                / max(resolutions[2]["pooled_weighted_rms"], STABILITY_FLOOR)
            ),
        },
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
    }
    receipt["audit_pass"] = _derive_audit_pass(receipt)
    digest_payload = copy.deepcopy(receipt)
    receipt["receipt_sha256"] = hashlib.sha256(
        _canonical_json(digest_payload).encode("utf-8")
    ).hexdigest()
    return receipt


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    """Fail closed on protocol drift, receipt mutation, or scientific gate failure."""
    if receipt.get("schema") != SCHEMA:
        raise ValueError("unexpected A4 divergence receipt schema")
    protocol = receipt.get("protocol", {})
    expected_protocol = {
        "seed": SEED,
        "spatial_steps": list(SPATIAL_STEPS),
        "times": list(TIMES),
        "X_interval": list(X_INTERVAL),
        "eta_interval": list(ETA_INTERVAL),
        "points_per_time": POINTS_PER_TIME,
        "divergence_gate": DIVERGENCE_GATE,
        "stability_factor": STABILITY_FACTOR,
        "stability_floor": STABILITY_FLOOR,
        "nontrivial_speed_rms": NONTRIVIAL_SPEED_RMS,
        "mutation_epsilon": MUTATION_EPSILON,
        "mutation_detection_floor": MUTATION_DETECTION_FLOOR,
    }
    for key, value in expected_protocol.items():
        if protocol.get(key) != value:
            raise ValueError(f"frozen protocol drift at {key}")
    if receipt.get("upstream", {}).get("head") != UPSTREAM_HEAD:
        raise ValueError("upstream exact-head identity drift")
    if receipt.get("upstream", {}).get("source_blob") != UPSTREAM_SOURCE_BLOB:
        raise ValueError("upstream source-blob identity drift")
    truth = receipt.get("truth_boundary", {})
    for key, expected in _TRUTH_BOUNDARY.items():
        if truth.get(key) is not expected:
            raise ValueError(f"truth-boundary drift at {key}")
    derived = _derive_audit_pass(receipt)
    if bool(receipt.get("audit_pass")) != derived:
        raise ValueError("audit_pass does not match frozen scientific gates")
    digest_payload = copy.deepcopy(dict(receipt))
    observed_digest = digest_payload.pop("receipt_sha256", None)
    expected_digest = hashlib.sha256(
        _canonical_json(digest_payload).encode("utf-8")
    ).hexdigest()
    if observed_digest != expected_digest:
        raise ValueError("receipt checksum mismatch")
    if not derived:
        fine = receipt["resolutions"][-1]
        raise RuntimeError(
            "current Cartesian leading divergence audit failed: "
            f"sampled_max={fine['sampled_max']:.6e}, "
            f"weighted_rms={fine['pooled_weighted_rms']:.6e}, "
            f"axis_max={fine['axis_axis_near_sampled_max']:.6e}"
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
        print("A4 current Cartesian leading divergence receipt passes frozen gates")
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
