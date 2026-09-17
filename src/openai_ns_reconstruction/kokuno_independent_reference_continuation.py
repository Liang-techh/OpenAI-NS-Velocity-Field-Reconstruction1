"""Independent seam-stratified PDE preflight for Kokuno reference continuation.

Kokuno Agent 4 owns independent validation.  This module consumes only the
public ``velocity(x,y,z,t)`` and ``pressure(x,y,z,t)`` values of
``KokunoReferenceContinuationCandidate`` when evaluating derivatives.  It does
not call the candidate's profile derivatives, recurrence residuals, training
losses, or any Agent-2 analytic differential operator.

The numerical operator is deliberately different from the existing Agent-4 FD4
operator reused by Agent 2: ordinary second-order centered Cartesian finite
differences are evaluated on a three-step refinement ladder.  The purpose is to
locate persistent momentum/divergence structure across the inner core,
reference-continuation transition, post-core flat-reference region, and an
axis-near stratum.  This intermediate source stage is not the final global
Kokuno field, so the registered full-domain 1e-3 PDE gate remains unassessed.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_reference_continuation import KokunoReferenceContinuationCandidate


SCHEMA = "kokuno-independent-reference-continuation-fd2-v1"
TASK_ID = "KOKUNO-VAL-REFERENCE-CONTINUATION-FD2-005"
VALIDATION_SEED = 9_172_991
MUTATION_SEED = 9_172_997
DEFAULT_TIMES = (0.375, 0.5, 0.625)
DEFAULT_STEPS = (0.004, 0.002, 0.001)
NU = 0.01
REGISTERED_PDE_THRESHOLD = 1.0e-3
REGISTERED_DIVERGENCE_THRESHOLD = 1.0e-5

Velocity = Callable[[Any, Any, Any, Any], np.ndarray]
Pressure = Callable[[Any, Any, Any, Any], np.ndarray]


def _finite_points(points: Any) -> np.ndarray:
    array = np.asarray(points, dtype=float)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if array.shape[0] == 0 or not np.all(np.isfinite(array)):
        raise ValueError("points must be nonempty and finite")
    return array


def _velocity_call(
    velocity: Velocity, points: np.ndarray, time: float | np.ndarray
) -> np.ndarray:
    values = np.asarray(
        velocity(points[:, 0], points[:, 1], points[:, 2], time), dtype=float
    )
    if values.shape != points.shape or not np.all(np.isfinite(values)):
        raise ValueError("velocity must return a finite (n,3) array")
    return values


def _pressure_call(
    pressure: Pressure, points: np.ndarray, time: float | np.ndarray
) -> np.ndarray:
    values = np.asarray(
        pressure(points[:, 0], points[:, 1], points[:, 2], time), dtype=float
    )
    if values.shape != (points.shape[0],) or not np.all(np.isfinite(values)):
        raise ValueError("pressure must return a finite (n,) array")
    return values


def evaluate_fd2(
    velocity: Velocity,
    pressure: Pressure,
    points: Any,
    time: float,
    step: float,
    *,
    nu: float = NU,
) -> dict[str, np.ndarray]:
    """Reconstruct div and raw NS momentum with independent centered FD2.

    The raw diagnostic uses
        R = u_t + (u dot grad)u + grad p - nu Delta u
    with no fitted/free forcing.  ``step`` is used for both space and time so
    the refinement ladder changes every numerical derivative together.
    """

    points_array = _finite_points(points)
    time = float(time)
    step = float(step)
    nu = float(nu)
    if not math.isfinite(time):
        raise ValueError("time must be finite")
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    if not math.isfinite(nu) or nu < 0.0:
        raise ValueError("nu must be nonnegative and finite")

    u0 = _velocity_call(velocity, points_array, time)
    plus_t = _velocity_call(velocity, points_array, time + step)
    minus_t = _velocity_call(velocity, points_array, time - step)
    u_t = (plus_t - minus_t) / (2.0 * step)

    jacobian = np.empty((len(points_array), 3, 3), dtype=float)
    laplacian = np.zeros_like(u0)
    pressure_gradient = np.empty_like(u0)
    for axis in range(3):
        plus = points_array.copy()
        minus = points_array.copy()
        plus[:, axis] += step
        minus[:, axis] -= step

        u_plus = _velocity_call(velocity, plus, time)
        u_minus = _velocity_call(velocity, minus, time)
        jacobian[:, :, axis] = (u_plus - u_minus) / (2.0 * step)
        laplacian += (u_plus - 2.0 * u0 + u_minus) / (step * step)

        p_plus = _pressure_call(pressure, plus, time)
        p_minus = _pressure_call(pressure, minus, time)
        pressure_gradient[:, axis] = (p_plus - p_minus) / (2.0 * step)

    divergence = np.trace(jacobian, axis1=1, axis2=2)
    convection = np.einsum("nij,nj->ni", jacobian, u0)
    residual = u_t + convection + pressure_gradient - nu * laplacian
    return {
        "velocity": u0,
        "time_derivative": u_t,
        "jacobian": jacobian,
        "laplacian": laplacian,
        "pressure_gradient": pressure_gradient,
        "divergence": divergence,
        "residual": residual,
    }


def _solve_native_q(z: np.ndarray, time: float, h: float) -> np.ndarray:
    """Independent fixed-point solve of q-z^2 q^(2h)=1-t for sampling only."""

    z = np.asarray(z, dtype=float)
    tau = 1.0 - float(time)
    if not (0.0 < tau < 1.0):
        raise ValueError("sampling time must lie strictly inside (0,1)")
    q = np.full_like(z, tau, dtype=float) + z * z * tau ** (2.0 * h)
    for _ in range(80):
        updated = tau + z * z * q ** (2.0 * h)
        if float(np.max(np.abs(updated - q))) <= 2.0e-15:
            q = updated
            break
        q = updated
    residual = q - z * z * q ** (2.0 * h) - tau
    if float(np.max(np.abs(residual))) > 2.0e-13:
        raise ArithmeticError("independent native-q sampling solve did not converge")
    return q


def _points_from_X(
    rng: np.random.Generator,
    *,
    count: int,
    time: float,
    h: float,
    X_low: float,
    X_high: float,
    z_half: float,
) -> tuple[np.ndarray, np.ndarray]:
    X = rng.uniform(float(X_low), float(X_high), int(count))
    z = rng.uniform(-float(z_half), float(z_half), int(count))
    q = _solve_native_q(z, time, h)
    radius = np.sqrt(2.0 * q * X)
    angle = rng.uniform(-np.pi, np.pi, int(count))
    points = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
    return points, X


def sample_strata(
    candidate: KokunoReferenceContinuationCandidate,
    *,
    time: float,
    count: int,
    seed: int,
) -> dict[str, dict[str, np.ndarray]]:
    """Create deterministic off-grid points in four predeclared geometric strata."""

    if isinstance(count, bool) or not isinstance(count, (int, np.integer)):
        raise TypeError("count must be an integer")
    count = int(count)
    if count < 2:
        raise ValueError("count must be at least 2 per stratum")
    rng = np.random.default_rng(int(seed) + int(round(10_000.0 * float(time))))

    inner_points, inner_X = _points_from_X(
        rng,
        count=count,
        time=time,
        h=candidate.h,
        X_low=0.18,
        X_high=0.28,
        z_half=0.08,
    )
    transition_low = candidate.X1 + 0.20 * (candidate.X2 - candidate.X1)
    transition_high = candidate.X1 + 0.80 * (candidate.X2 - candidate.X1)
    transition_points, transition_X = _points_from_X(
        rng,
        count=count,
        time=time,
        h=candidate.h,
        X_low=transition_low,
        X_high=transition_high,
        z_half=0.08,
    )
    post_low = max(0.46, candidate.max_source_transition_X + 0.04)
    post_points, post_X = _points_from_X(
        rng,
        count=count,
        time=time,
        h=candidate.h,
        X_low=post_low,
        X_high=0.56,
        z_half=0.08,
    )

    radius = rng.uniform(0.002, 0.008, count)
    angle = rng.uniform(-np.pi, np.pi, count)
    z_abs = rng.uniform(0.12, 0.22, count)
    sign = rng.choice(np.array([-1.0, 1.0]), size=count)
    axis_z = sign * z_abs
    axis_points = np.column_stack(
        (radius * np.cos(angle), radius * np.sin(angle), axis_z)
    )
    axis_q = _solve_native_q(axis_z, time, candidate.h)
    axis_X = radius * radius / (2.0 * axis_q)

    return {
        "inner_core": {"points": inner_points, "X": inner_X},
        "reference_transition": {"points": transition_points, "X": transition_X},
        "post_core_reference": {"points": post_points, "X": post_X},
        "axis_near": {"points": axis_points, "X": axis_X},
    }


def _cylindrical_residual(
    residual: np.ndarray, points: np.ndarray
) -> np.ndarray:
    radius = np.hypot(points[:, 0], points[:, 1])
    if np.any(radius <= 0.0):
        raise ValueError("cylindrical projection requires nonzero radius")
    radial = (
        points[:, 0] * residual[:, 0] + points[:, 1] * residual[:, 1]
    ) / radius
    theta = (
        -points[:, 1] * residual[:, 0] + points[:, 0] * residual[:, 1]
    ) / radius
    return np.column_stack((radial, theta, residual[:, 2]))


def _rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _metrics(result: dict[str, np.ndarray], points: np.ndarray) -> dict[str, Any]:
    residual = np.asarray(result["residual"], dtype=float)
    divergence = np.asarray(result["divergence"], dtype=float)
    velocity = np.asarray(result["velocity"], dtype=float)
    norm = np.linalg.norm(residual, axis=1)
    speed = np.linalg.norm(velocity, axis=1)
    cylindrical = _cylindrical_residual(residual, points)
    worst = int(np.argmax(norm))
    component_rms = np.sqrt(np.mean(cylindrical * cylindrical, axis=0))
    sample_l2 = _rms(norm)
    return {
        "sample_l2": sample_l2,
        "normalized_sample_l2": sample_l2,
        "residual_max": float(norm[worst]),
        "divergence_rms": _rms(divergence),
        "divergence_max": float(np.max(np.abs(divergence))),
        "velocity_rms": _rms(speed),
        "velocity_max": float(np.max(speed)),
        "cylindrical_component_rms": {
            "radial": float(component_rms[0]),
            "theta": float(component_rms[1]),
            "axial": float(component_rms[2]),
        },
        "worst_point": points[worst].tolist(),
        "worst_residual_cartesian": residual[worst].tolist(),
        "worst_residual_cylindrical": cylindrical[worst].tolist(),
    }


def _relative_change(coarse: float, fine: float) -> float:
    scale = max(abs(float(fine)), 1.0e-30)
    return float(abs(float(coarse) - float(fine)) / scale)


def _mutation_calibration(
    candidate: KokunoReferenceContinuationCandidate,
    *,
    step: float,
) -> dict[str, float]:
    points = sample_strata(
        candidate, time=0.5, count=4, seed=MUTATION_SEED
    )["inner_core"]["points"]

    baseline = evaluate_fd2(
        candidate.velocity, candidate.pressure, points, 0.5, step, nu=NU
    )

    def velocity_mutated(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        values = np.asarray(candidate.velocity(x, y, z, t), dtype=float).copy()
        values[..., 0] += 0.02 * np.asarray(x, dtype=float)
        return values

    def pressure_mutated(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return np.asarray(candidate.pressure(x, y, z, t), dtype=float) + 0.02 * np.asarray(
            x, dtype=float
        )

    div_mutation = evaluate_fd2(
        velocity_mutated, candidate.pressure, points, 0.5, step, nu=NU
    )
    pressure_mutation = evaluate_fd2(
        candidate.velocity, pressure_mutated, points, 0.5, step, nu=NU
    )

    divergence_shift = np.asarray(div_mutation["divergence"]) - np.asarray(
        baseline["divergence"]
    )
    x_momentum_shift = np.asarray(pressure_mutation["residual"])[:, 0] - np.asarray(
        baseline["residual"]
    )[:, 0]
    return {
        "injected_slope": 0.02,
        "mean_divergence_shift": float(np.mean(divergence_shift)),
        "max_divergence_shift_error": float(np.max(np.abs(divergence_shift - 0.02))),
        "mean_pressure_x_momentum_shift": float(np.mean(x_momentum_shift)),
        "max_pressure_x_momentum_shift_error": float(
            np.max(np.abs(x_momentum_shift - 0.02))
        ),
    }


def run_audit(
    *,
    count: int = 4,
    seed: int = VALIDATION_SEED,
    times: tuple[float, ...] = DEFAULT_TIMES,
    steps: tuple[float, ...] = DEFAULT_STEPS,
) -> dict[str, Any]:
    """Run one fresh independent reference-continuation validation increment."""

    times = tuple(float(value) for value in times)
    steps = tuple(float(value) for value in steps)
    if len(steps) < 3:
        raise ValueError("at least three numerical resolutions are required")
    if any(not np.isfinite(value) for value in times + steps):
        raise ValueError("times and steps must be finite")
    if any(value <= 0.0 for value in steps):
        raise ValueError("steps must be positive")
    if sorted(steps, reverse=True) != list(steps):
        raise ValueError("steps must be listed coarse to fine")

    candidate = KokunoReferenceContinuationCandidate()
    rows: list[dict[str, Any]] = []
    samples_by_time: dict[float, dict[str, dict[str, np.ndarray]]] = {}
    for time in times:
        samples = sample_strata(candidate, time=time, count=count, seed=seed)
        samples_by_time[time] = samples
        for stratum, sample in samples.items():
            points = sample["points"]
            X = sample["X"]
            for step in steps:
                result = evaluate_fd2(
                    candidate.velocity,
                    candidate.pressure,
                    points,
                    time,
                    step,
                    nu=NU,
                )
                rows.append(
                    {
                        "stratum": stratum,
                        "time": time,
                        "step": step,
                        "X_min": float(np.min(X)),
                        "X_max": float(np.max(X)),
                        "metrics": _metrics(result, points),
                    }
                )

    finest = min(steps)
    medium = sorted(steps)[1]
    finest_rows = [row for row in rows if row["step"] == finest]
    strata = ("inner_core", "reference_transition", "post_core_reference", "axis_near")
    summaries: dict[str, dict[str, Any]] = {}
    for stratum in strata:
        selected = [row for row in finest_rows if row["stratum"] == stratum]
        l2_values = np.array([row["metrics"]["sample_l2"] for row in selected])
        component = {
            name: float(
                np.sqrt(
                    np.mean(
                        [
                            row["metrics"]["cylindrical_component_rms"][name] ** 2
                            for row in selected
                        ]
                    )
                )
            )
            for name in ("radial", "theta", "axial")
        }
        summaries[stratum] = {
            "equal_time_sample_l2": _rms(l2_values),
            "worst_residual_max": float(
                max(row["metrics"]["residual_max"] for row in selected)
            ),
            "worst_divergence_max": float(
                max(row["metrics"]["divergence_max"] for row in selected)
            ),
            "cylindrical_component_rms": component,
            "minimum_velocity_rms": float(
                min(row["metrics"]["velocity_rms"] for row in selected)
            ),
        }

    sensitivity_rows: list[dict[str, Any]] = []
    for stratum in strata:
        for time in times:
            med_row = next(
                row
                for row in rows
                if row["stratum"] == stratum
                and row["time"] == time
                and row["step"] == medium
            )
            fine_row = next(
                row
                for row in rows
                if row["stratum"] == stratum
                and row["time"] == time
                and row["step"] == finest
            )
            sensitivity_rows.append(
                {
                    "stratum": stratum,
                    "time": time,
                    "sample_l2_relative_change": _relative_change(
                        med_row["metrics"]["sample_l2"],
                        fine_row["metrics"]["sample_l2"],
                    ),
                    "residual_max_relative_change": _relative_change(
                        med_row["metrics"]["residual_max"],
                        fine_row["metrics"]["residual_max"],
                    ),
                    "divergence_rms_relative_change": _relative_change(
                        med_row["metrics"]["divergence_rms"],
                        fine_row["metrics"]["divergence_rms"],
                    ),
                }
            )

    global_sample_l2 = _rms(
        np.array([row["metrics"]["sample_l2"] for row in finest_rows], dtype=float)
    )
    worst_row = max(finest_rows, key=lambda row: row["metrics"]["residual_max"])
    largest_component = max(
        (
            (
                name,
                summary["cylindrical_component_rms"][name],
                stratum,
            )
            for stratum, summary in summaries.items()
            for name in ("radial", "theta", "axial")
        ),
        key=lambda item: item[1],
    )

    transition = candidate.to_payload()["transition"]
    transition_X = np.concatenate(
        [
            samples_by_time[time]["reference_transition"]["X"]
            for time in times
        ]
    )
    post_X = np.concatenate(
        [samples_by_time[time]["post_core_reference"]["X"] for time in times]
    )

    report = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "candidate": {
            "class": "KokunoReferenceContinuationCandidate",
            "sha256": candidate.sha256,
            "source_commit": candidate.to_payload()["source"]["commit"],
            "reference_continuation_stage_only": True,
        },
        "validation_contract": {
            "fresh_seed": int(seed),
            "points_per_stratum_per_time": int(count),
            "strata": list(strata),
            "times": list(times),
            "steps": list(steps),
            "nu": NU,
            "residual_normalization_scale": 1.0,
            "registered_pde_threshold": REGISTERED_PDE_THRESHOLD,
            "registered_divergence_threshold": REGISTERED_DIVERGENCE_THRESHOLD,
            "derivative_operator": "independent_second_order_centered_cartesian_fd",
            "candidate_profile_derivatives_used": False,
            "candidate_internal_residual_helpers_used": False,
            "training_loss_used": False,
            "forcing": "zero_raw_intermediate_stage_diagnostic_only_no_fit",
            "validation_used_for_parameter_selection": False,
        },
        "stage_geometry": {
            "X0": transition["X0"],
            "X1": transition["X1"],
            "X2": transition["X2"],
            "old_core_X_limit": transition["source_core_limit"],
            "transition_sample_X_min": float(np.min(transition_X)),
            "transition_sample_X_max": float(np.max(transition_X)),
            "post_core_sample_X_min": float(np.min(post_X)),
            "post_core_sample_X_max": float(np.max(post_X)),
            "transition_stratum_inside_declared_transition": bool(
                np.min(transition_X) > transition["X1"]
                and np.max(transition_X) < transition["X2"]
            ),
            "post_core_stratum_beyond_old_core_limit": bool(
                np.min(post_X) > transition["source_core_limit"]
            ),
        },
        "rows": rows,
        "finest_step": finest,
        "finest_summaries": summaries,
        "finest_equal_stratum_time_sample_l2": global_sample_l2,
        "finest_normalized_equal_stratum_time_sample_l2": global_sample_l2,
        "finest_worst_residual_max": float(worst_row["metrics"]["residual_max"]),
        "finest_worst_location": {
            "stratum": worst_row["stratum"],
            "time": worst_row["time"],
            "point": worst_row["metrics"]["worst_point"],
            "residual_cartesian": worst_row["metrics"]["worst_residual_cartesian"],
            "residual_cylindrical": worst_row["metrics"][
                "worst_residual_cylindrical"
            ],
        },
        "largest_finest_component_rms": {
            "component": largest_component[0],
            "rms": float(largest_component[1]),
            "stratum": largest_component[2],
        },
        "medium_to_fine_sensitivity": {
            "medium_step": medium,
            "fine_step": finest,
            "rows": sensitivity_rows,
            "worst_sample_l2_relative_change": float(
                max(row["sample_l2_relative_change"] for row in sensitivity_rows)
            ),
            "worst_residual_max_relative_change": float(
                max(row["residual_max_relative_change"] for row in sensitivity_rows)
            ),
        },
        "mutation_calibration": _mutation_calibration(candidate, step=finest),
        "gate": {
            "formal_full_domain_gate_assessed": False,
            "formal_full_domain_gate_passed": False,
            "reason": (
                "reference continuation still precedes cone modulation, five-moment "
                "repair, heat compensation, final core-to-heat assembly, and a "
                "compatible final fixed/restricted forcing contract"
            ),
        },
        "truth_boundary": {
            "reference_continuation_independently_audited": True,
            "global_leading_profile_reconstructed": False,
            "core_to_heat_matching_completed": False,
            "leading_plus_oscillatory_validated_by_this_increment": False,
            "after_correction_cycle_validated_by_this_increment": False,
            "global_support_validated": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    if not report["stage_geometry"]["transition_stratum_inside_declared_transition"]:
        raise AssertionError("transition sampling escaped the declared source stage")
    if not report["stage_geometry"]["post_core_stratum_beyond_old_core_limit"]:
        raise AssertionError("post-core sampling did not exercise the new continuation")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--count", type=int, default=4)
    parser.add_argument("--seed", type=int, default=VALIDATION_SEED)
    args = parser.parse_args(argv)
    report = run_audit(count=args.count, seed=args.seed)
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    if args.output is None:
        print(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
        print(
            json.dumps(
                {
                    "task_id": report["task_id"],
                    "finest_sample_l2": report[
                        "finest_equal_stratum_time_sample_l2"
                    ],
                    "finest_worst_max": report["finest_worst_residual_max"],
                    "largest_component": report["largest_finest_component_rms"],
                    "pde_validated": report["truth_boundary"]["pde_validated"],
                },
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
