"""Candidate-agnostic instantaneous streamline morphology diagnostics.

The diagnostic consumes only a frozen ``velocity(points, time) -> (..., 3)``
callable. Streamlines are integrated along normalized instantaneous velocity,
so uniform velocity-amplitude rescaling does not alter the geometry fingerprint.
This is visualization evidence only, not PDE validation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable

import numpy as np


VelocityFn = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class _Trace:
    status: str
    points: np.ndarray


class _SpeedFloor(RuntimeError):
    pass


def _as_seeds(seeds: np.ndarray) -> np.ndarray:
    arr = np.asarray(seeds, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3 or arr.shape[0] == 0:
        raise ValueError("seeds must have shape (n, 3) with n > 0")
    if not np.all(np.isfinite(arr)):
        raise ValueError("seeds must be finite")
    return arr


def _as_bounds(bounds: tuple[Iterable[float], Iterable[float]] | None):
    if bounds is None:
        return None
    if len(bounds) != 2:
        raise ValueError("bounds must be (lower, upper)")
    lo = np.asarray(tuple(bounds[0]), dtype=float)
    hi = np.asarray(tuple(bounds[1]), dtype=float)
    if lo.shape != (3,) or hi.shape != (3,):
        raise ValueError("bounds vectors must each have length 3")
    if not np.all(np.isfinite(lo)) or not np.all(np.isfinite(hi)):
        raise ValueError("bounds must be finite")
    if not np.all(hi > lo):
        raise ValueError("each upper bound must exceed its lower bound")
    return lo, hi


def _inside(point: np.ndarray, bounds) -> bool:
    if bounds is None:
        return True
    lo, hi = bounds
    return bool(np.all(point >= lo) and np.all(point <= hi))


def _velocity_vector(velocity: VelocityFn, point: np.ndarray, time: float) -> np.ndarray:
    raw = np.asarray(velocity(point.reshape(1, 3), float(time)), dtype=float)
    if raw.shape == (3,):
        raw = raw.reshape(1, 3)
    if raw.shape != (1, 3):
        raise ValueError("velocity must return shape (n, 3)")
    if not np.all(np.isfinite(raw)):
        raise ValueError("velocity returned non-finite values")
    return raw[0]


def _unit_tangent(
    velocity: VelocityFn,
    point: np.ndarray,
    time: float,
    speed_floor: float,
) -> np.ndarray:
    vec = _velocity_vector(velocity, point, time)
    speed = float(np.linalg.norm(vec))
    if not np.isfinite(speed):
        raise ValueError("velocity speed is non-finite")
    if speed <= speed_floor:
        raise _SpeedFloor
    return vec / speed


def _rk4_unit_step(
    velocity: VelocityFn,
    point: np.ndarray,
    time: float,
    ds: float,
    speed_floor: float,
) -> np.ndarray:
    k1 = _unit_tangent(velocity, point, time, speed_floor)
    k2 = _unit_tangent(velocity, point + 0.5 * ds * k1, time, speed_floor)
    k3 = _unit_tangent(velocity, point + 0.5 * ds * k2, time, speed_floor)
    k4 = _unit_tangent(velocity, point + ds * k3, time, speed_floor)
    result = point + (ds / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    if not np.all(np.isfinite(result)):
        raise ValueError("streamline integration produced non-finite coordinates")
    return result


def _trace_direction(
    velocity: VelocityFn,
    seed: np.ndarray,
    time: float,
    arc_length: float,
    step: float,
    direction: float,
    speed_floor: float,
    bounds,
) -> _Trace:
    point = np.asarray(seed, dtype=float)
    points = [point.copy()]
    remaining = float(arc_length)
    status = "complete"
    while remaining > 0.0:
        ds = min(float(step), remaining) * float(direction)
        try:
            new_point = _rk4_unit_step(velocity, point, time, ds, speed_floor)
        except _SpeedFloor:
            status = "speed_floor"
            break
        if not _inside(new_point, bounds):
            points.append(new_point.copy())
            status = "domain_exit"
            break
        points.append(new_point.copy())
        point = new_point
        remaining -= abs(ds)
    return _Trace(status=status, points=np.asarray(points, dtype=float))


def _combine_bidirectional(backward: _Trace, forward: _Trace) -> np.ndarray:
    back = backward.points[::-1]
    return np.concatenate([back[:-1], forward.points], axis=0)


def _quantiles(values: np.ndarray) -> dict[str, float]:
    if values.size == 0:
        return {"q10": float("nan"), "q50": float("nan"), "q90": float("nan")}
    q10, q50, q90 = np.quantile(values, [0.1, 0.5, 0.9])
    return {"q10": float(q10), "q50": float(q50), "q90": float(q90)}


def _path_metrics(points: np.ndarray, turn_floor: float) -> dict[str, float | None]:
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] == 0:
        raise ValueError("trace points must have shape (m, 3)")
    arclength = 0.0 if points.shape[0] == 1 else float(
        np.linalg.norm(np.diff(points, axis=0), axis=1).sum()
    )
    radius = np.hypot(points[:, 0], points[:, 1])
    z = points[:, 2]
    phi = np.unwrap(np.arctan2(points[:, 1], points[:, 0]))
    signed_turns = float((phi[-1] - phi[0]) / (2.0 * np.pi))
    abs_turns = abs(signed_turns)
    z_span = float(np.max(z) - np.min(z))
    pitch = None if abs_turns <= turn_floor else float(z_span / abs_turns)
    return {
        "arclength": arclength,
        "radial_min": float(np.min(radius)),
        "radial_max": float(np.max(radius)),
        "radial_drift": float(np.max(radius) - np.min(radius)),
        "z_min": float(np.min(z)),
        "z_max": float(np.max(z)),
        "z_span": z_span,
        "signed_turns": signed_turns,
        "abs_turns": abs_turns,
        "pitch_per_turn": pitch,
    }


def diagnose_streamline_fingerprint(
    velocity: VelocityFn,
    time: float,
    seeds: np.ndarray,
    *,
    half_arclength: float = 4.0,
    step: float = 0.05,
    speed_floor: float = 1e-12,
    turn_floor: float = 1e-8,
    bounds: tuple[Iterable[float], Iterable[float]] | None = None,
) -> dict[str, Any]:
    """Trace frozen-time streamlines and summarize observable geometry."""
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(half_arclength) or half_arclength <= 0.0:
        raise ValueError("half_arclength must be positive and finite")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    if step > half_arclength:
        raise ValueError("step must not exceed half_arclength")
    if not np.isfinite(speed_floor) or speed_floor < 0.0:
        raise ValueError("speed_floor must be finite and nonnegative")
    if not np.isfinite(turn_floor) or turn_floor < 0.0:
        raise ValueError("turn_floor must be finite and nonnegative")

    seed_array = _as_seeds(seeds)
    checked_bounds = _as_bounds(bounds)
    traces = []
    complete = domain_exit = speed_floor_count = 0

    for index, seed in enumerate(seed_array):
        if not _inside(seed, checked_bounds):
            raise ValueError("every seed must lie within bounds")
        backward = _trace_direction(
            velocity, seed, float(time), float(half_arclength), float(step),
            -1.0, float(speed_floor), checked_bounds,
        )
        forward = _trace_direction(
            velocity, seed, float(time), float(half_arclength), float(step),
            +1.0, float(speed_floor), checked_bounds,
        )
        points = _combine_bidirectional(backward, forward)
        metrics = _path_metrics(points, float(turn_floor))
        statuses = {backward.status, forward.status}
        if statuses == {"complete"}:
            status = "complete"
            complete += 1
        elif "domain_exit" in statuses:
            status = "domain_exit"
            domain_exit += 1
        else:
            status = "speed_floor"
            speed_floor_count += 1
        metrics.update({
            "seed_index": int(index),
            "status": status,
            "backward_status": backward.status,
            "forward_status": forward.status,
            "point_count": int(points.shape[0]),
        })
        traces.append(metrics)

    arclengths = np.asarray([row["arclength"] for row in traces], dtype=float)
    drifts = np.asarray([row["radial_drift"] for row in traces], dtype=float)
    z_spans = np.asarray([row["z_span"] for row in traces], dtype=float)
    signed_turns = np.asarray([row["signed_turns"] for row in traces], dtype=float)
    abs_turns = np.asarray([row["abs_turns"] for row in traces], dtype=float)
    pitches = np.asarray(
        [row["pitch_per_turn"] for row in traces if row["pitch_per_turn"] is not None],
        dtype=float,
    )
    denominator = float(np.sum(abs_turns))
    turn_coherence = None if denominator <= turn_floor else float(
        abs(np.sum(signed_turns)) / denominator
    )

    return {
        "claim_scope": "instantaneous_streamline_visualization_diagnostic_only",
        "time": float(time),
        "seed_count": int(seed_array.shape[0]),
        "half_arclength": float(half_arclength),
        "step": float(step),
        "complete_fraction": float(complete / seed_array.shape[0]),
        "domain_exit_fraction": float(domain_exit / seed_array.shape[0]),
        "speed_floor_fraction": float(speed_floor_count / seed_array.shape[0]),
        "radial_extent": [min(float(row["radial_min"]) for row in traces),
                          max(float(row["radial_max"]) for row in traces)],
        "axial_extent": [min(float(row["z_min"]) for row in traces),
                         max(float(row["z_max"]) for row in traces)],
        "arclength": _quantiles(arclengths),
        "radial_drift": _quantiles(drifts),
        "z_span": _quantiles(z_spans),
        "abs_turns": _quantiles(abs_turns),
        "pitch_per_turn": _quantiles(pitches),
        "signed_turn_coherence": turn_coherence,
        "traces": traces,
        "pde_validated": False,
        "paper_exact": False,
    }


def _metric_value(report: dict[str, Any], name: str) -> float | None:
    if name in {"complete_fraction", "signed_turn_coherence"}:
        value = report[name]
        return None if value is None else float(value)
    section, quantile = name.split(".", 1)
    value = float(report[section][quantile])
    return None if not np.isfinite(value) else value


def audit_streamline_step_resolution(
    velocity: VelocityFn,
    time: float,
    seeds: np.ndarray,
    *,
    steps: Iterable[float] = (0.1, 0.05, 0.025),
    half_arclength: float = 4.0,
    speed_floor: float = 1e-12,
    turn_floor: float = 1e-8,
    bounds: tuple[Iterable[float], Iterable[float]] | None = None,
) -> dict[str, Any]:
    """Compare morphology across at least three streamline integration steps."""
    step_values = tuple(float(value) for value in steps)
    if len(step_values) < 3:
        raise ValueError("at least three integration step sizes are required")
    if not all(np.isfinite(value) and value > 0.0 for value in step_values):
        raise ValueError("integration step sizes must be positive and finite")
    if not all(a > b for a, b in zip(step_values, step_values[1:])):
        raise ValueError("integration step sizes must be strictly decreasing")

    reports = [
        diagnose_streamline_fingerprint(
            velocity, time, seeds, half_arclength=half_arclength, step=value,
            speed_floor=speed_floor, turn_floor=turn_floor, bounds=bounds,
        )
        for value in step_values
    ]
    finest = reports[-1]
    metric_names = (
        "complete_fraction",
        "abs_turns.q50",
        "pitch_per_turn.q50",
        "z_span.q50",
        "radial_drift.q90",
        "signed_turn_coherence",
    )
    comparisons: dict[str, list[dict[str, float | None]]] = {}
    for metric_name in metric_names:
        finest_value = _metric_value(finest, metric_name)
        entries = []
        for step_value, report in zip(step_values, reports):
            value = _metric_value(report, metric_name)
            absolute = relative = None
            if value is not None and finest_value is not None:
                absolute = float(abs(value - finest_value))
                if abs(finest_value) > 1e-14:
                    relative = float(absolute / abs(finest_value))
            entries.append({
                "step": float(step_value),
                "value": value,
                "absolute_delta_to_finest": absolute,
                "relative_delta_to_finest": relative,
            })
        comparisons[metric_name] = entries

    return {
        "claim_scope": "streamline_integration_resolution_audit_only",
        "time": float(time),
        "steps": list(step_values),
        "reports": reports,
        "comparisons": comparisons,
        "pde_validated": False,
        "paper_exact": False,
    }
