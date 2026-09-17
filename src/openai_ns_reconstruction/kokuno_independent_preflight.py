"""Independent numerical preflight for the Kokuno-derived oscillatory correction.

This module is intentionally validation-only.  It does not import a training
loss, optimizer state, or derivative helper from the correction implementation.
It treats ``KokunoCompleteCurlCorrection`` only through its public velocity and
vector-potential evaluators, then reconstructs curl/divergence with a separate
fourth-order Cartesian finite-difference operator on held-out off-grid points.

This is *not* the registered Navier--Stokes acceptance gate.  The full
leading+oscillatory(+mean/correction) 3-D velocity plus pressure/fixed-forcing
contract is not yet present in the Kokuno lane, so the preregistered 1e-3 full
momentum threshold remains unassessed here.
"""
from __future__ import annotations

import json
from math import log2
from pathlib import Path
from typing import Callable

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection


HELD_OUT_SEED = 9172641
TIMES = (0.3125, 0.5, 0.6875)
STEPS = (0.02, 0.01, 0.005)
MUTATION_DIVERGENCE = 0.02
POINT_COUNT = 384
REGISTERED_PDE_THRESHOLD = 1.0e-3


def _velocity(correction: KokunoCompleteCurlCorrection, points: np.ndarray, time: float) -> np.ndarray:
    return np.asarray(correction.at_points(points, time), dtype=float)


def _potential(correction: KokunoCompleteCurlCorrection, points: np.ndarray, time: float) -> np.ndarray:
    return np.asarray(
        correction.vector_potential(points[:, 0], points[:, 1], points[:, 2], time),
        dtype=float,
    )


def _held_out_points(correction: KokunoCompleteCurlCorrection) -> tuple[np.ndarray, dict[str, slice]]:
    rng = np.random.default_rng(HELD_OUT_SEED)
    widths = np.asarray(correction.half_widths, dtype=float)
    interior = rng.uniform(-0.65, 0.65, (128, 3)) * widths
    collar = rng.uniform(-0.6, 0.6, (128, 3)) * widths
    axes = rng.integers(0, 3, 128)
    signs = rng.choice((-1.0, 1.0), 128)
    for index, axis in enumerate(axes):
        collar[index, axis] = signs[index] * rng.uniform(0.82, 0.94) * widths[axis]
    axis_near = np.zeros((64, 3), dtype=float)
    axis_near[:, :2] = rng.uniform(-2.0e-3, 2.0e-3, (64, 2))
    axis_near[:, 2] = rng.uniform(-0.65, 0.65, 64) * widths[2]
    corner = rng.uniform(0.75, 0.9, (64, 3)) * widths * rng.choice((-1.0, 1.0), (64, 3))
    points = np.vstack((interior, collar, axis_near, corner))
    return points, {
        "interior": slice(0, 128),
        "collar": slice(128, 256),
        "axis_near": slice(256, 320),
        "corner": slice(320, 384),
    }


def _fd4_first(function: Callable[[np.ndarray, float], np.ndarray], points: np.ndarray, time: float, axis: int, step: float) -> np.ndarray:
    direction = np.zeros(3, dtype=float)
    direction[axis] = step
    return (
        function(points - 2.0 * direction, time)
        - 8.0 * function(points - direction, time)
        + 8.0 * function(points + direction, time)
        - function(points + 2.0 * direction, time)
    ) / (12.0 * step)


def _curl_from_public_potential(correction: KokunoCompleteCurlCorrection, points: np.ndarray, time: float, step: float) -> np.ndarray:
    function = lambda p, t: _potential(correction, p, t)
    dx, dy, dz = [_fd4_first(function, points, time, axis, step) for axis in range(3)]
    return np.stack((dy[:, 2] - dz[:, 1], dz[:, 0] - dx[:, 2], dx[:, 1] - dy[:, 0]), axis=1)


def _divergence_from_public_velocity(correction: KokunoCompleteCurlCorrection, points: np.ndarray, time: float, step: float, *, inject_divergence: float = 0.0) -> np.ndarray:
    def function(p: np.ndarray, t: float) -> np.ndarray:
        values = _velocity(correction, p, t).copy()
        if inject_divergence:
            values[:, 0] += inject_divergence * p[:, 0]
        return values
    derivatives = [_fd4_first(function, points, time, axis, step) for axis in range(3)]
    return derivatives[0][:, 0] + derivatives[1][:, 1] + derivatives[2][:, 2]


def _norms(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    magnitudes = np.linalg.norm(values, axis=-1) if values.ndim > 1 else np.abs(values)
    return {"max": float(np.max(magnitudes)), "rms": float(np.sqrt(np.mean(magnitudes * magnitudes)))}


def _order(coarse: float, fine: float) -> float:
    if coarse <= 0.0 or fine <= 0.0:
        raise ValueError("positive errors are required for convergence order")
    return float(log2(coarse / fine))


def run_independent_preflight(correction: KokunoCompleteCurlCorrection | None = None) -> dict[str, object]:
    if correction is None:
        correction = KokunoCompleteCurlCorrection()
    points, regions = _held_out_points(correction)
    if points.shape != (POINT_COUNT, 3):
        raise RuntimeError("unexpected held-out point count")
    rows: list[dict[str, float]] = []
    fine_regions: dict[str, dict[str, dict[str, float]]] = {}
    velocity_activity: dict[str, dict[str, float]] = {}
    for time in TIMES:
        by_step: list[dict[str, float]] = []
        for step in STEPS:
            velocity = _velocity(correction, points, time)
            curl_error = _curl_from_public_potential(correction, points, time, step) - velocity
            divergence = _divergence_from_public_velocity(correction, points, time, step)
            mutated = _divergence_from_public_velocity(correction, points, time, step, inject_divergence=MUTATION_DIVERGENCE)
            curl_norms = _norms(curl_error)
            div_norms = _norms(divergence)
            mutation_norms = _norms(mutated)
            row = {
                "time": float(time), "step": float(step),
                "curl_error_max": curl_norms["max"], "curl_error_rms": curl_norms["rms"],
                "divergence_max": div_norms["max"], "divergence_rms": div_norms["rms"],
                "mutation_divergence_max": mutation_norms["max"], "mutation_divergence_rms": mutation_norms["rms"],
            }
            rows.append(row)
            by_step.append(row)
        for index in range(2):
            by_step[index]["curl_rms_order_to_next"] = _order(by_step[index]["curl_error_rms"], by_step[index + 1]["curl_error_rms"])
            by_step[index]["divergence_rms_order_to_next"] = _order(by_step[index]["divergence_rms"], by_step[index + 1]["divergence_rms"])
        finest = by_step[-1]
        fine_curl = _curl_from_public_potential(correction, points, time, STEPS[-1]) - _velocity(correction, points, time)
        fine_divergence = _divergence_from_public_velocity(correction, points, time, STEPS[-1])
        fine_regions[str(time)] = {name: {"curl_error": _norms(fine_curl[region]), "divergence": _norms(fine_divergence[region])} for name, region in regions.items()}
        velocity_activity[str(time)] = _norms(_velocity(correction, points, time))
        if finest["curl_error_rms"] > 1.0e-7 or finest["divergence_rms"] > 1.0e-7:
            raise AssertionError("finest oscillatory structural error exceeds preflight bound")
        if abs(finest["mutation_divergence_rms"] - MUTATION_DIVERGENCE) > 5.0e-5:
            raise AssertionError("divergence mutation calibration failed")
    orders = [value for row in rows for key, value in row.items() if key in {"curl_rms_order_to_next", "divergence_rms_order_to_next"}]
    if min(orders) < 3.5:
        raise AssertionError("fourth-order held-out convergence was not observed")
    widths = np.asarray(correction.half_widths, dtype=float)
    outside = []
    for axis in range(3):
        for sign in (-1.0, 1.0):
            point = np.zeros(3, dtype=float)
            point[axis] = sign * 1.02 * widths[axis]
            outside.append(point)
    outside_points = np.asarray(outside)
    outside_support_max = max(float(np.max(np.abs(_velocity(correction, outside_points, time)))) for time in TIMES)
    if outside_support_max != 0.0:
        raise AssertionError("public correction is not exactly zero outside its compact support")
    return {
        "task_id": "KOKUNO-VAL-OSC-PREFLIGHT-001",
        "candidate_family": correction.metadata()["family"],
        "operator": "independent fourth-order Cartesian centered finite differences",
        "public_interfaces_used": ["at_points/velocity", "vector_potential"],
        "held_out_seed": HELD_OUT_SEED,
        "point_count": POINT_COUNT,
        "times": list(TIMES),
        "steps": list(STEPS),
        "rows": rows,
        "finest_region_metrics": fine_regions,
        "velocity_activity": velocity_activity,
        "outside_support_max": outside_support_max,
        "structural_preflight_passed": True,
        "registered_full_ns_gate": {
            "normalized_residual_threshold": REGISTERED_PDE_THRESHOLD,
            "assessed": False,
            "reason": "No complete Kokuno leading+oscillatory(+correction) 3-D velocity with pressure and fixed/restricted forcing contract is available in this ancestry.",
        },
        "pde_validated": False,
        "truth_boundary": {
            "training_loss_used": False,
            "candidate_internal_derivative_helpers_used": False,
            "threshold_changed": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def main() -> None:
    report = run_independent_preflight()
    output = Path("artifacts/kokuno_independent_validation/oscillatory_preflight_report.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
