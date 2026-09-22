"""Independent Kokuno Agent-4 audit of A1 #1171's split relative-swirl field.

The upstream field preserves a mathematically nonzero, sub-binary64-relative-
epsilon correction as two public Cartesian channels

    velocity_split(x,y,z,t) -> (base_velocity, delta_velocity).

This audit deliberately does not use the upstream analytic/source-coordinate
profile derivatives.  After a save/load roundtrip it treats ``velocity_split``
as a black-box Cartesian evaluator and reconstructs the delta Jacobian with an
A4-owned centered finite-difference operator at three preregistered relative
physical resolutions.  The canonical project absolute FD ladder is checked
separately for binary64 representability and is never replaced for final PDE
admission.

This is scoped representation/divergence evidence only.  It is not a complete
NS residual and cannot promote ``pde_validated``.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_relative_swirl_split import (
    KokunoPA16CurrentCartesianRelativeSwirlSplit,
)

TASK = "K4-VAL-121"
SCHEMA = "kokuno-a4-relative-swirl-split-independent-audit-v1"
UPSTREAM_PR = 1171
UPSTREAM_HEAD = "657818dd83e119e6091bd2490d3804c04c3ef723"
SEED = 9173921
TIMES = (0.31, 0.47, 0.63, 0.71)
RELATIVE_FD_LADDER = (2.0e-6, 1.0e-6, 5.0e-7)
CANONICAL_ABSOLUTE_FD_LADDER = (0.02, 0.01, 0.005)

NORMALIZED_DIVERGENCE_GATE = 1.0e-5
LEAKAGE_RELATIVE_GATE = 5.0e-10
ZERO_ABSOLUTE_GATE = 1.0e-300
ACTIVE_RATIO_FLOOR = 1.0e-45
REFINEMENT_FLOOR = 2.0e-8
REFINEMENT_WORSEN_FACTOR = 1.25
SAVE_LOAD_REPLAY_GATE = 2.0e-12
FINAL_PROJECT_MOMENTUM_GATE = 1.0e-3
FINAL_PROJECT_DIVERGENCE_GATE = 1.0e-5

_TRUTH_BOUNDARY = {
    "split_relative_swirl_cartesian_channel_independently_audited": True,
    "canonical_absolute_fd_transfer_ready": False,
    "binary64_total_relative_swirl_sum_is_resolved": False,
    "current_cartesian_relative_swirl_composed": False,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
}


def default_candidate() -> KokunoPA16CurrentCartesianRelativeSwirlSplit:
    return KokunoPA16CurrentCartesianRelativeSwirlSplit()


def _as_vec(value: Any) -> np.ndarray:
    arr = np.asarray(value, dtype=float).reshape(-1)
    if arr.size != 3 or np.any(~np.isfinite(arr)):
        raise ValueError("velocity_split must return a finite Cartesian 3-vector")
    return arr


def _split(candidate: Any, point: Mapping[str, float]) -> tuple[np.ndarray, np.ndarray]:
    base, delta = candidate.velocity_split(
        float(point["x"]), float(point["y"]), float(point["z"]), float(point["t"])
    )
    return _as_vec(base), _as_vec(delta)


def _q_at_z0(candidate: Any, t: float) -> float:
    coords = candidate.similarity_coordinates_logX(1.0, 0.0, 0.0, float(t))
    q = float(np.asarray(coords["q"], dtype=float))
    if not (math.isfinite(q) and q > 0.0):
        raise RuntimeError("public similarity map returned invalid q at z=0")
    return q


def _physical_point(candidate: Any, *, s: float, theta: float, t: float, role: str) -> dict[str, float | str]:
    q = _q_at_z0(candidate, t)
    log_X = float(candidate.log_X_flatten_end + s)
    log_r = 0.5 * (math.log(2.0) + math.log(q) + log_X)
    if log_r >= math.log(np.finfo(float).max):
        raise RuntimeError("held-out relative-swirl radius is not finite in binary64")
    r = math.exp(log_r)
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    point: dict[str, float | str] = {
        "x": x,
        "y": y,
        "z": 0.0,
        "t": float(t),
        "r": r,
        "theta": float(theta),
        "q": q,
        "s_target": float(s),
        "log_X_target": log_X,
        "role": role,
    }
    coords = candidate.similarity_coordinates_logX(x, y, 0.0, float(t))
    got_log_X = float(np.asarray(coords["log_X"], dtype=float))
    got_eta = float(np.asarray(coords["eta"], dtype=float))
    point["log_X_reconstruction_error"] = abs(got_log_X - log_X)
    point["eta"] = got_eta
    if point["log_X_reconstruction_error"] > 5.0e-11 or abs(got_eta) > 5.0e-13:
        raise RuntimeError("held-out physical/source-coordinate reconstruction drift")
    return point


def _heldout_points(candidate: Any) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    rng = np.random.default_rng(SEED)
    width = float(candidate.compensator.width)
    y1 = float(candidate.compensator.y1)
    y2 = float(candidate.compensator.y2)
    active: list[dict[str, float | str]] = []
    controls: list[dict[str, float | str]] = []
    for t in TIMES:
        for label, center in (("bump1", y1), ("bump2", y2)):
            for k in range(2):
                # Strictly interior and off-grid.  The jitter is frozen before execution.
                offset = float(rng.uniform(-0.12 * width, 0.12 * width))
                theta = float(rng.uniform(0.31, 2.0 * math.pi - 0.31))
                active.append(
                    _physical_point(
                        candidate,
                        s=center + offset,
                        theta=theta,
                        t=t,
                        role=f"{label}_interior_{k}",
                    )
                )
        for role, s in (
            ("before_bump1", y1 - 2.0 * width),
            ("between_bumps", 0.5 * (y1 + y2)),
            ("after_bump2", y2 + 2.0 * width),
        ):
            theta = float(rng.uniform(0.37, 2.0 * math.pi - 0.37))
            controls.append(_physical_point(candidate, s=s, theta=theta, t=t, role=role))
    return active, controls


def _canonical_representability(points: list[Mapping[str, float | str]]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    all_ok = True
    for h in CANONICAL_ABSOLUTE_FD_LADDER:
        per_step = []
        for p in points:
            axes = []
            for key in ("x", "y", "z"):
                c = float(p[key])
                plus = np.float64(c) + np.float64(h)
                minus = np.float64(c) - np.float64(h)
                ok = bool(plus != c and minus != c and plus != minus)
                axes.append(ok)
                all_ok = all_ok and ok
            per_step.append(bool(all(axes)))
        rows.append(
            {
                "absolute_step": h,
                "all_points_all_axes_representable": bool(all(per_step)),
                "representable_point_fraction": float(np.mean(per_step)),
            }
        )
    return {"ladder": rows, "all_representable": bool(all_ok)}


def _jacobian_delta(candidate: Any, point: Mapping[str, float | str], rel_step: float) -> np.ndarray:
    r = float(point["r"])
    q = float(point["q"])
    h_xy = r * rel_step
    h_z = max(math.sqrt(q), 1.0e-6) * rel_step
    if not (math.isfinite(h_xy) and h_xy > 0.0 and math.isfinite(h_z) and h_z > 0.0):
        raise RuntimeError("invalid independent finite-difference step")
    jac = np.empty((3, 3), dtype=float)
    for axis, (key, h) in enumerate((("x", h_xy), ("y", h_xy), ("z", h_z))):
        plus = dict(point)
        minus = dict(point)
        plus[key] = float(point[key]) + h
        minus[key] = float(point[key]) - h
        if not (float(plus[key]) != float(point[key]) and float(minus[key]) != float(point[key])):
            raise RuntimeError("scoped relative finite-difference perturbation collapsed in binary64")
        _, d_plus = _split(candidate, plus)
        _, d_minus = _split(candidate, minus)
        jac[:, axis] = (d_plus - d_minus) / (2.0 * h)
    if np.any(~np.isfinite(jac)):
        raise RuntimeError("independent delta Jacobian became non-finite")
    return jac


def _point_observables(candidate: Any, point: Mapping[str, float | str]) -> dict[str, float]:
    base, delta = _split(candidate, point)
    theta = float(point["theta"])
    er = np.array([math.cos(theta), math.sin(theta), 0.0])
    et = np.array([-math.sin(theta), math.cos(theta), 0.0])
    dnorm = float(np.linalg.norm(delta))
    radial = float(np.dot(delta, er))
    tangential = float(np.dot(delta, et))
    axial = float(delta[2])
    base_tangential = float(np.dot(base, et))
    leakage = math.hypot(radial, axial) / max(dnorm, 1.0e-300)
    ratio = abs(tangential) / max(abs(base_tangential), 1.0e-300)
    return {
        "delta_norm": dnorm,
        "delta_radial": radial,
        "delta_tangential": tangential,
        "delta_axial": axial,
        "relative_leakage": leakage,
        "abs_delta_over_abs_base_tangential": ratio,
    }


def _resolution_metrics(candidate: Any, points: list[Mapping[str, float | str]], rel_step: float) -> dict[str, Any]:
    abs_div: list[float] = []
    normalized: list[float] = []
    leaks: list[float] = []
    ratios: list[float] = []
    witnesses: list[dict[str, Any]] = []
    for p in points:
        obs = _point_observables(candidate, p)
        jac = _jacobian_delta(candidate, p, rel_step)
        div = float(np.trace(jac))
        grad_norm = float(np.linalg.norm(jac))
        geometric_grad = obs["delta_norm"] / max(float(p["r"]), 1.0e-300)
        scale = max(grad_norm, geometric_grad, 1.0e-300)
        nd = abs(div) / scale
        abs_div.append(abs(div))
        normalized.append(nd)
        leaks.append(obs["relative_leakage"])
        ratios.append(obs["abs_delta_over_abs_base_tangential"])
        witnesses.append(
            {
                "role": str(p["role"]),
                "t": float(p["t"]),
                "x": float(p["x"]),
                "y": float(p["y"]),
                "z": float(p["z"]),
                "r": float(p["r"]),
                "s_target": float(p["s_target"]),
                "abs_divergence": abs(div),
                "normalized_divergence": nd,
                "gradient_frobenius": grad_norm,
                **obs,
            }
        )
    i = int(np.argmax(normalized))
    values = np.asarray(normalized, dtype=float)
    div_values = np.asarray(abs_div, dtype=float)
    return {
        "relative_step": rel_step,
        "sample_count": len(points),
        "divergence_sampled_max": float(np.max(div_values)),
        "divergence_sample_rms": float(np.sqrt(np.mean(div_values**2))),
        "normalized_divergence_sampled_max": float(np.max(values)),
        "normalized_divergence_sample_rms": float(np.sqrt(np.mean(values**2))),
        "relative_leakage_max": float(np.max(leaks)),
        "active_ratio_min": float(np.min(ratios)),
        "active_ratio_max": float(np.max(ratios)),
        "worst_witness": witnesses[i],
    }


def _replay_error(loaded: Any, original: Any, points: list[Mapping[str, float | str]]) -> float:
    worst = 0.0
    for p in points:
        b0, d0 = _split(original, p)
        b1, d1 = _split(loaded, p)
        for a, b in ((b0, b1), (d0, d1)):
            denom = max(float(np.linalg.norm(a)), float(np.linalg.norm(b)), 1.0e-300)
            worst = max(worst, float(np.linalg.norm(a - b)) / denom)
    return worst


def _zero_checks(candidate: Any, controls: list[Mapping[str, float | str]]) -> tuple[float, float]:
    outside = 0.0
    for p in controls:
        _, delta = _split(candidate, p)
        outside = max(outside, float(np.max(np.abs(delta))))
    axis = 0.0
    for t in TIMES:
        _, delta = candidate.velocity_split(0.0, 0.0, 0.0, t)
        axis = max(axis, float(np.max(np.abs(np.asarray(delta, dtype=float)))))
    return outside, axis


def _assess(report: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    fine = report["resolutions"][-1]
    medium = report["resolutions"][-2]
    if fine["normalized_divergence_sampled_max"] > NORMALIZED_DIVERGENCE_GATE:
        failures.append("fine_normalized_divergence_max")
    if fine["normalized_divergence_sample_rms"] > NORMALIZED_DIVERGENCE_GATE:
        failures.append("fine_normalized_divergence_rms")
    if fine["relative_leakage_max"] > LEAKAGE_RELATIVE_GATE:
        failures.append("radial_or_axial_leakage")
    if report["off_bump_delta_abs_max"] > ZERO_ABSOLUTE_GATE:
        failures.append("off_bump_not_zero")
    if report["exact_axis_delta_abs_max"] > ZERO_ABSOLUTE_GATE:
        failures.append("axis_not_zero")
    if fine["active_ratio_max"] < ACTIVE_RATIO_FLOOR:
        failures.append("delta_channel_collapsed")
    if fine["active_ratio_max"] >= np.finfo(float).eps:
        failures.append("delta_channel_not_sub_epsilon")
    if report["save_load_replay_relative_max"] > SAVE_LOAD_REPLAY_GATE:
        failures.append("save_load_replay")
    f = float(fine["normalized_divergence_sampled_max"])
    m = float(medium["normalized_divergence_sampled_max"])
    if f > REFINEMENT_FLOOR and f > REFINEMENT_WORSEN_FACTOR * max(m, 1.0e-300):
        failures.append("medium_to_fine_divergence_worsened")
    return failures


def audit_loaded_split(loaded: Any, pre_serialization_reference: Any) -> dict[str, Any]:
    """Run the fixed K4-VAL-121 audit; this API intentionally has no tuning knobs."""
    truth = dict(getattr(loaded, "truth_boundary"))
    required_false = (
        "binary64_total_relative_swirl_sum_is_resolved",
        "current_cartesian_relative_swirl_composed",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "pde_validated",
    )
    for key in required_false:
        if bool(truth.get(key, False)):
            raise ValueError(f"upstream truth promotion is not admissible: {key}")

    active, controls = _heldout_points(loaded)
    replay = _replay_error(loaded, pre_serialization_reference, active + controls)
    canonical = _canonical_representability(active)
    outside, axis = _zero_checks(loaded, controls)
    resolutions = [_resolution_metrics(loaded, active, h) for h in RELATIVE_FD_LADDER]

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "upstream_pr": UPSTREAM_PR,
        "upstream_head": UPSTREAM_HEAD,
        "seed": SEED,
        "times": list(TIMES),
        "relative_fd_ladder": list(RELATIVE_FD_LADDER),
        "canonical_absolute_fd_ladder": list(CANONICAL_ABSOLUTE_FD_LADDER),
        "active_points": active,
        "control_points": controls,
        "resolutions": resolutions,
        "canonical_absolute_fd_representability": canonical,
        "save_load_replay_relative_max": replay,
        "off_bump_delta_abs_max": outside,
        "exact_axis_delta_abs_max": axis,
        "gates": {
            "scoped_normalized_divergence_gate": NORMALIZED_DIVERGENCE_GATE,
            "leakage_relative_gate": LEAKAGE_RELATIVE_GATE,
            "zero_absolute_gate": ZERO_ABSOLUTE_GATE,
            "active_ratio_floor": ACTIVE_RATIO_FLOOR,
            "active_ratio_ceiling_float64_epsilon": float(np.finfo(float).eps),
            "refinement_floor": REFINEMENT_FLOOR,
            "refinement_worsen_factor": REFINEMENT_WORSEN_FACTOR,
            "save_load_replay_gate": SAVE_LOAD_REPLAY_GATE,
            "final_project_momentum_gate_unchanged": FINAL_PROJECT_MOMENTUM_GATE,
            "final_project_divergence_gate_unchanged": FINAL_PROJECT_DIVERGENCE_GATE,
        },
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
    }
    failures = _assess(report)
    report["failures"] = failures
    report["scoped_split_channel_gate_passed"] = len(failures) == 0
    # Even a scoped PASS cannot transfer to final project admission when the
    # canonical absolute stencil is not representable or the complete API is absent.
    report["truth_boundary"]["canonical_absolute_fd_transfer_ready"] = bool(
        canonical["all_representable"] and len(failures) == 0
    )
    report["final_project_admission_ready"] = False
    blockers = [
        "terminal_global_leading_missing",
        "matching_self_contained_oscillatory_composite_missing",
        "cartesian_correction_velocity_missing",
        "matched_pressure_missing",
        "preregistered_restricted_forcing_missing",
        "complete_ns_residual_not_assessed",
    ]
    if not canonical["all_representable"]:
        blockers.insert(0, "canonical_absolute_fd_ladder_not_binary64_representable_here")
    report["final_project_blockers"] = blockers
    return report


def enforce_preregistered_gates(report: Mapping[str, Any]) -> None:
    failures = _assess(report)
    if failures:
        raise AssertionError("K4-VAL-121 scoped gate failure: " + ", ".join(failures))
    if report.get("truth_boundary", {}).get("pde_validated") is not False:
        raise AssertionError("pde_validated must remain false")
    if report.get("final_project_admission_ready") is not False:
        raise AssertionError("scoped split audit cannot promote final project admission")
