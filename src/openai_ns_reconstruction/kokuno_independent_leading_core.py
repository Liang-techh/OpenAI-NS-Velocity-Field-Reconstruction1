"""Independent finite-difference preflight for the Kokuno leading-core candidate.

This module deliberately treats the Agent-1 candidate as a black-box public
``velocity``/``pressure`` provider.  It does not consume profile derivatives,
training tensors, ``PaperCoreSeries`` residual helpers, or any candidate-side
Jacobian/Laplacian implementation.

The check is *core-only*.  The public candidate is not yet the globally matched,
compactly supported Kokuno leading field, and the oscillatory/mean corrections
are not composed here.  Consequently the registered full-domain ``1e-3`` PDE
gate remains unassessed even though this module reports the unchanged numerical
threshold next to the local diagnostic.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate


VelocityCallable = Callable[[Any, Any, Any, Any], np.ndarray]
PressureCallable = Callable[[Any, Any, Any, Any], np.ndarray]

SCHEMA = "kokuno-independent-leading-core-preflight-v1"
VALIDATION_SEED = 9_172_841
DEFAULT_TIMES = (0.375, 0.5, 0.625)
DEFAULT_STEPS = (0.02, 0.01, 0.005)
NU = 0.01
REGISTERED_RESIDUAL_SCALE = 1.0
REGISTERED_PDE_THRESHOLD = 1.0e-3
REGISTERED_DIVERGENCE_THRESHOLD = 1.0e-5
CORE_BOX_HALF_WIDTH_XY = 0.18
CORE_BOX_HALF_HEIGHT_Z = 0.22
CORE_BOX_VOLUME = (2.0 * CORE_BOX_HALF_WIDTH_XY) ** 2 * (2.0 * CORE_BOX_HALF_HEIGHT_Z)


def _as_points(points_xyz: Any) -> np.ndarray:
    points = np.asarray(points_xyz, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points_xyz must have shape (n,3)")
    if points.shape[0] == 0 or not np.all(np.isfinite(points)):
        raise ValueError("points_xyz must be finite and nonempty")
    return points


def _velocity_values(
    velocity: VelocityCallable, points: np.ndarray, time: float
) -> np.ndarray:
    values = np.asarray(
        velocity(points[:, 0], points[:, 1], points[:, 2], float(time)),
        dtype=float,
    )
    if values.shape != (points.shape[0], 3) or not np.all(np.isfinite(values)):
        raise ValueError("velocity callable must return finite shape (n,3)")
    return values


def _pressure_values(
    pressure: PressureCallable, points: np.ndarray, time: float
) -> np.ndarray:
    values = np.asarray(
        pressure(points[:, 0], points[:, 1], points[:, 2], float(time)),
        dtype=float,
    )
    if values.shape != (points.shape[0],) or not np.all(np.isfinite(values)):
        raise ValueError("pressure callable must return finite shape (n,)")
    return values


def _shift(points: np.ndarray, axis: int, amount: float) -> np.ndarray:
    shifted = np.array(points, copy=True)
    shifted[:, axis] += float(amount)
    return shifted


def evaluate_fd4(
    velocity: VelocityCallable,
    pressure: PressureCallable,
    points_xyz: Any,
    time: float,
    step: float,
    *,
    nu: float = NU,
) -> dict[str, np.ndarray]:
    """Evaluate divergence and raw NS momentum with a separate FD4 operator.

    The diagnostic equation is

    ``R = u_t + (u dot grad)u + grad p - nu Delta u``.

    There is intentionally no fitted forcing in this routine.  ``f=0`` is only
    a raw leading-core diagnostic choice, not a claim that it is the compatible
    force for the final project candidate.
    """

    points = _as_points(points_xyz)
    time = float(time)
    step = float(step)
    nu = float(nu)
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    if not np.isfinite(nu) or nu < 0.0:
        raise ValueError("nu must be nonnegative and finite")

    u0 = _velocity_values(velocity, points, time)
    jacobian = np.empty((points.shape[0], 3, 3), dtype=float)
    laplacian = np.zeros_like(u0)
    grad_pressure = np.empty((points.shape[0], 3), dtype=float)

    for axis in range(3):
        xm2 = _shift(points, axis, -2.0 * step)
        xm1 = _shift(points, axis, -step)
        xp1 = _shift(points, axis, step)
        xp2 = _shift(points, axis, 2.0 * step)

        um2 = _velocity_values(velocity, xm2, time)
        um1 = _velocity_values(velocity, xm1, time)
        up1 = _velocity_values(velocity, xp1, time)
        up2 = _velocity_values(velocity, xp2, time)
        jacobian[:, :, axis] = (
            um2 - 8.0 * um1 + 8.0 * up1 - up2
        ) / (12.0 * step)
        laplacian += (
            -up2 + 16.0 * up1 - 30.0 * u0 + 16.0 * um1 - um2
        ) / (12.0 * step * step)

        pm2 = _pressure_values(pressure, xm2, time)
        pm1 = _pressure_values(pressure, xm1, time)
        pp1 = _pressure_values(pressure, xp1, time)
        pp2 = _pressure_values(pressure, xp2, time)
        grad_pressure[:, axis] = (
            pm2 - 8.0 * pm1 + 8.0 * pp1 - pp2
        ) / (12.0 * step)

    tm2 = _velocity_values(velocity, points, time - 2.0 * step)
    tm1 = _velocity_values(velocity, points, time - step)
    tp1 = _velocity_values(velocity, points, time + step)
    tp2 = _velocity_values(velocity, points, time + 2.0 * step)
    time_derivative = (tm2 - 8.0 * tm1 + 8.0 * tp1 - tp2) / (12.0 * step)

    divergence = np.trace(jacobian, axis1=1, axis2=2)
    convection = np.einsum("nij,nj->ni", jacobian, u0)
    residual = time_derivative + convection + grad_pressure - nu * laplacian
    return {
        "velocity": u0,
        "jacobian": jacobian,
        "laplacian": laplacian,
        "time_derivative": time_derivative,
        "grad_pressure": grad_pressure,
        "divergence": divergence,
        "residual": residual,
    }


def sample_held_out_core_points(
    *,
    interior_count: int = 48,
    axis_near_count: int = 24,
    seed: int = VALIDATION_SEED,
) -> dict[str, np.ndarray]:
    """Return fresh off-grid points safely inside the public core domain.

    The uniform interior box is used for a Monte-Carlo volume-L2 estimate.
    Axis-near points are a separate stress stratum and are never mixed into that
    volume quadrature.
    """

    if interior_count < 8 or axis_near_count < 4:
        raise ValueError("sample counts are too small for the governed preflight")
    rng = np.random.default_rng(int(seed))
    interior = np.empty((int(interior_count), 3), dtype=float)
    interior[:, 0] = rng.uniform(
        -CORE_BOX_HALF_WIDTH_XY, CORE_BOX_HALF_WIDTH_XY, interior_count
    )
    interior[:, 1] = rng.uniform(
        -CORE_BOX_HALF_WIDTH_XY, CORE_BOX_HALF_WIDTH_XY, interior_count
    )
    interior[:, 2] = rng.uniform(
        -CORE_BOX_HALF_HEIGHT_Z, CORE_BOX_HALF_HEIGHT_Z, interior_count
    )

    radius = rng.uniform(2.0e-4, 2.0e-2, axis_near_count)
    angle = rng.uniform(0.0, 2.0 * np.pi, axis_near_count)
    axis_near = np.column_stack(
        (
            radius * np.cos(angle),
            radius * np.sin(angle),
            rng.uniform(-CORE_BOX_HALF_HEIGHT_Z, CORE_BOX_HALF_HEIGHT_Z, axis_near_count),
        )
    )
    return {"interior": interior, "axis_near": axis_near}


def _metrics(values: dict[str, np.ndarray], *, volume: float | None) -> dict[str, float]:
    residual_norm = np.linalg.norm(values["residual"], axis=1)
    divergence_abs = np.abs(values["divergence"])
    speed = np.linalg.norm(values["velocity"], axis=1)
    result = {
        "residual_max": float(np.max(residual_norm)),
        "residual_rms": float(np.sqrt(np.mean(residual_norm * residual_norm))),
        "normalized_residual_max": float(np.max(residual_norm) / REGISTERED_RESIDUAL_SCALE),
        "normalized_residual_rms": float(
            np.sqrt(np.mean(residual_norm * residual_norm)) / REGISTERED_RESIDUAL_SCALE
        ),
        "divergence_max": float(np.max(divergence_abs)),
        "divergence_rms": float(np.sqrt(np.mean(values["divergence"] ** 2))),
        "velocity_rms": float(np.sqrt(np.mean(speed * speed))),
        "velocity_max": float(np.max(speed)),
    }
    if volume is not None:
        result["volume_l2_residual_estimate"] = float(
            np.sqrt(float(volume) * np.mean(residual_norm * residual_norm))
        )
        result["volume_l2_divergence_estimate"] = float(
            np.sqrt(float(volume) * np.mean(values["divergence"] ** 2))
        )
    return result


def _mutation_calibration(
    candidate: KokunoLeadingCoreSeriesCandidate,
    points: np.ndarray,
    *,
    time: float = 0.5,
    step: float = 0.005,
    slope: float = 0.02,
) -> dict[str, float]:
    baseline = evaluate_fd4(candidate.velocity, candidate.pressure, points, time, step)

    def mutated_velocity(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        value = np.array(candidate.velocity(x, y, z, t), copy=True)
        value[..., 0] += slope * np.asarray(x, dtype=float)
        return value

    mutated = evaluate_fd4(mutated_velocity, candidate.pressure, points, time, step)
    divergence_shift = mutated["divergence"] - baseline["divergence"]

    def mutated_pressure(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return candidate.pressure(x, y, z, t) + slope * np.asarray(x, dtype=float)

    pressure_mutated = evaluate_fd4(
        candidate.velocity, mutated_pressure, points, time, step
    )
    residual_shift = pressure_mutated["residual"] - baseline["residual"]
    return {
        "injected_linear_velocity_slope": float(slope),
        "divergence_shift_mean": float(np.mean(divergence_shift)),
        "divergence_shift_max_abs_error_from_expected": float(
            np.max(np.abs(divergence_shift - slope))
        ),
        "injected_linear_pressure_slope": float(slope),
        "pressure_gradient_x_shift_mean": float(np.mean(residual_shift[:, 0])),
        "pressure_gradient_x_shift_max_abs_error_from_expected": float(
            np.max(np.abs(residual_shift[:, 0] - slope))
        ),
        "pressure_gradient_transverse_shift_max": float(
            np.max(np.abs(residual_shift[:, 1:]))
        ),
    }


def run_leading_core_preflight(
    *,
    interior_count: int = 48,
    axis_near_count: int = 24,
    seed: int = VALIDATION_SEED,
    times: tuple[float, ...] = DEFAULT_TIMES,
    steps: tuple[float, ...] = DEFAULT_STEPS,
) -> dict[str, Any]:
    candidate = KokunoLeadingCoreSeriesCandidate()
    samples = sample_held_out_core_points(
        interior_count=interior_count,
        axis_near_count=axis_near_count,
        seed=seed,
    )
    all_points = np.vstack((samples["interior"], samples["axis_near"]))

    # Fail early if any governed FD stencil would leave the candidate's public
    # core domain or time domain.  This is a domain guard, not a support test.
    largest_step = max(float(value) for value in steps)
    for time in times:
        if time - 2.0 * largest_step <= 0.0 or time + 2.0 * largest_step >= 1.0:
            raise ValueError("time/step ladder leaves the candidate time domain")
        candidate.at_points(all_points, float(time))
        for axis in range(3):
            candidate.at_points(_shift(all_points, axis, -2.0 * largest_step), float(time))
            candidate.at_points(_shift(all_points, axis, 2.0 * largest_step), float(time))

    rows: list[dict[str, Any]] = []
    for time in times:
        for step in steps:
            interior_values = evaluate_fd4(
                candidate.velocity,
                candidate.pressure,
                samples["interior"],
                float(time),
                float(step),
            )
            axis_values = evaluate_fd4(
                candidate.velocity,
                candidate.pressure,
                samples["axis_near"],
                float(time),
                float(step),
            )
            combined_values = {
                key: np.concatenate((interior_values[key], axis_values[key]), axis=0)
                for key in (
                    "velocity",
                    "divergence",
                    "residual",
                )
            }
            rows.append(
                {
                    "time": float(time),
                    "step": float(step),
                    "interior": _metrics(interior_values, volume=CORE_BOX_VOLUME),
                    "axis_near": _metrics(axis_values, volume=None),
                    "all_points": _metrics(combined_values, volume=None),
                }
            )

    finest = float(min(steps))
    finest_rows = [row for row in rows if row["step"] == finest]
    worst_local_normalized_max = max(
        row["all_points"]["normalized_residual_max"] for row in finest_rows
    )
    worst_local_normalized_rms = max(
        row["all_points"]["normalized_residual_rms"] for row in finest_rows
    )
    worst_divergence_max = max(
        row["all_points"]["divergence_max"] for row in finest_rows
    )
    worst_divergence_rms = max(
        row["all_points"]["divergence_rms"] for row in finest_rows
    )

    mutation_points = samples["interior"][: min(12, samples["interior"].shape[0])]
    report: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": "KOKUNO-VAL-LEADING-CORE-MOMENTUM-PREFLIGHT-003",
        "candidate": {
            "schema": candidate.to_payload()["schema"],
            "sha256": candidate.sha256,
            "source_commit": candidate.to_payload()["source"]["commit"],
            "core_domain": "Lambda*X<=4.1",
            "global_leading_profile_reconstructed": False,
        },
        "operator": {
            "kind": "independent_fourth_order_centered_cartesian_finite_difference",
            "candidate_derivative_helpers_used": False,
            "training_tensors_used": False,
            "public_contract_only": ["velocity", "pressure"],
            "nu": NU,
            "forcing": "zero_only_for_raw_leading_core_diagnostic_not_a_fitted_or_compatible_force_claim",
            "steps": [float(value) for value in steps],
            "times": [float(value) for value in times],
        },
        "sampling": {
            "seed": int(seed),
            "interior_count": int(interior_count),
            "axis_near_count": int(axis_near_count),
            "off_grid": True,
            "interior_uniform_box": [
                [-CORE_BOX_HALF_WIDTH_XY, CORE_BOX_HALF_WIDTH_XY],
                [-CORE_BOX_HALF_WIDTH_XY, CORE_BOX_HALF_WIDTH_XY],
                [-CORE_BOX_HALF_HEIGHT_Z, CORE_BOX_HALF_HEIGHT_Z],
            ],
            "interior_box_volume": CORE_BOX_VOLUME,
            "axis_near_radius_range": [2.0e-4, 2.0e-2],
        },
        "rows": rows,
        "finest_summary": {
            "step": finest,
            "worst_local_normalized_residual_max": float(worst_local_normalized_max),
            "worst_local_normalized_residual_rms": float(worst_local_normalized_rms),
            "worst_divergence_max": float(worst_divergence_max),
            "worst_divergence_rms": float(worst_divergence_rms),
            "registered_residual_threshold_reference": REGISTERED_PDE_THRESHOLD,
            "registered_divergence_threshold_reference": REGISTERED_DIVERGENCE_THRESHOLD,
        },
        "mutation_calibration": _mutation_calibration(candidate, mutation_points),
        "formal_full_pde_gate": {
            "normalized_residual_threshold": REGISTERED_PDE_THRESHOLD,
            "divergence_threshold": REGISTERED_DIVERGENCE_THRESHOLD,
            "assessed": False,
            "passed": False,
            "reason": (
                "Agent-1 field is inner/core only: no globally matched compact leading field, "
                "no oscillatory composition, no applied Agent-3 mean correction, and no "
                "compatible fixed/restricted forcing contract for the final composite."
            ),
        },
        "truth_boundary": {
            "leading_core_preflight_only": True,
            "global_support_validated": False,
            "complete_kokuno_composite_velocity": False,
            "after_correction_cycle_validated": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    return report


def write_report(path: str | Path, **kwargs: Any) -> Path:
    report = run_leading_core_preflight(**kwargs)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return target


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_independent_validation/leading_core_preflight.json",
    )
    parser.add_argument("--interior-count", type=int, default=48)
    parser.add_argument("--axis-near-count", type=int, default=24)
    parser.add_argument("--seed", type=int, default=VALIDATION_SEED)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    path = write_report(
        args.output,
        interior_count=args.interior_count,
        axis_near_count=args.axis_near_count,
        seed=args.seed,
    )
    report = json.loads(path.read_text(encoding="utf-8"))
    print(json.dumps(report["finest_summary"], sort_keys=True))
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
