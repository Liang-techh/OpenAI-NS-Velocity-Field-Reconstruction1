"""Pressure-free curl obstruction comparison for the Eq. (4.5) deliveries.

This diagnostic asks a necessary, pressure-independent question for the frozen
velocity artifacts: can the curl of the momentum balance be matched by the
curl of the preregistered two-parameter force?  The force is fit only on a
deterministic training set with ``0 <= a,c <= 10``.  The fitted coefficients
are then held fixed while independent samples and finite-difference steps are
used for the reported comparison.

The calculation deliberately does not fit pressure or an arbitrary force.  A
small curl obstruction is necessary for a pressure to exist, but it is not
sufficient for the full Navier--Stokes equation.  This module is a bounded
diagnostic for the two delivery candidates and does not promote any PDE or
visualization truth flag.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np

from .constrained_force import RestrictedForce
from .constrained_validation import residual
from .eq45_supported_delivery import Eq45SupportedDeliveryField


SCHEMA = "eq45_bipolar_obstruction_v1"
DEFAULT_TRAIN_SEED = 20260916
DEFAULT_HOLDOUT_SEED = 914027
DEFAULT_TIMES = (0.3125, 0.5625, 0.6875)
DEFAULT_TRAIN_COUNTS = {"uniform": 128, "core": 64, "collar": 64}
DEFAULT_HOLDOUT_COUNTS = {"uniform": 256, "core": 64, "collar": 64}
DEFAULT_TRAIN_STEP = 0.01
DEFAULT_HOLDOUT_STEPS = (0.01, 0.005)
DEFAULT_NU = 0.01
FORCE_BOUNDS = (0.0, 10.0)
BOX = np.asarray(((-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0)), dtype=float)
BOX_VOLUME = float(np.prod(BOX[:, 1] - BOX[:, 0]))


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    digest = hashlib.sha256()
    digest.update(str(array.shape).encode("ascii"))
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _as_points(values: Any, *, name: str) -> np.ndarray:
    points = np.asarray(values, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3)")
    if len(points) == 0 or not np.all(np.isfinite(points)):
        raise ValueError(f"{name} must be nonempty and finite")
    return points


def _eval_velocity(field: Any, points: np.ndarray, time: float) -> np.ndarray:
    """Evaluate only the public point-batch velocity interface."""

    values = np.asarray(field.at_points(points, float(time)), dtype=float)
    expected = (len(points), 3)
    if values.shape != expected:
        raise ValueError(f"velocity returned {values.shape}, expected {expected}")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity returned nonfinite values")
    return values


def _zero_pressure(points: np.ndarray, _time: float) -> np.ndarray:
    """Pressure gauge used by the momentum oracle; curl removes any gauge."""

    return np.zeros(len(points), dtype=float)


def _zero_force(points: np.ndarray, _time: float) -> np.ndarray:
    return np.zeros((len(points), 3), dtype=float)


def _pressure_free_momentum(
    field: Any,
    points: np.ndarray,
    time: float,
    *,
    nu: float,
    step: float,
    time_bounds: tuple[float, float],
) -> np.ndarray:
    """Return ``u_t + (u dot grad)u - nu Delta u`` from the existing oracle.

    ``constrained_validation.residual`` supplies fourth-order Cartesian spatial
    stencils and a second-order time stencil.  Passing zero pressure and zero
    force makes the returned ``momentum`` independent of pressure and keeps
    this diagnostic separate from the force fit.
    """

    result = residual(
        field.at_points,
        _zero_pressure,
        _zero_force,
        _as_points(points, name="points"),
        float(time),
        nu=float(nu),
        step=float(step),
        time_bounds=time_bounds,
    )
    momentum = np.asarray(result["momentum"], dtype=float)
    if momentum.shape != (len(points), 3) or not np.all(np.isfinite(momentum)):
        raise RuntimeError("pressure-free momentum oracle returned malformed values")
    return momentum


def centered_curl(
    vector_fn: Callable[[np.ndarray, float], np.ndarray],
    points: Any,
    time: float,
    step: float,
) -> np.ndarray:
    """Compute a Cartesian curl with an independent centered second-order stencil."""

    points_arr = _as_points(points, name="points")
    step = float(step)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    if not np.isfinite(time):
        raise ValueError("time must be finite")

    gradient = np.empty((len(points_arr), 3, 3), dtype=float)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = step
        plus = np.asarray(vector_fn(points_arr + offset, float(time)), dtype=float)
        minus = np.asarray(vector_fn(points_arr - offset, float(time)), dtype=float)
        expected = (len(points_arr), 3)
        if plus.shape != expected or minus.shape != expected:
            raise ValueError("vector_fn must return shape (n, 3)")
        if not np.all(np.isfinite(plus)) or not np.all(np.isfinite(minus)):
            raise ValueError("vector_fn returned nonfinite values")
        gradient[:, :, axis] = (plus - minus) / (2.0 * step)

    curl = np.empty_like(points_arr)
    curl[:, 0] = gradient[:, 2, 1] - gradient[:, 1, 2]
    curl[:, 1] = gradient[:, 0, 2] - gradient[:, 2, 0]
    curl[:, 2] = gradient[:, 1, 0] - gradient[:, 0, 1]
    return curl


def curl_of_momentum(
    field: Any,
    points: Any,
    time: float,
    *,
    nu: float,
    step: float,
    time_bounds: tuple[float, float],
) -> np.ndarray:
    """Compute curl of the pressure-free momentum using nested FD operators."""

    momentum = lambda shifted, shifted_time: _pressure_free_momentum(
        field,
        shifted,
        shifted_time,
        nu=nu,
        step=step,
        time_bounds=time_bounds,
    )
    return centered_curl(momentum, points, float(time), float(step))


def curl_of_force(
    force: RestrictedForce,
    points: Any,
    time: float,
    *,
    step: float,
) -> np.ndarray:
    """Compute force curl through the same independent centered curl operator."""

    if not isinstance(force, RestrictedForce):
        raise TypeError("force must be RestrictedForce")
    return centered_curl(force, points, float(time), float(step))


def _sample_core(rng: np.random.Generator, count: int) -> np.ndarray:
    """Sample an axis-safe central cylinder inside the physical support."""

    theta = rng.uniform(0.0, 2.0 * np.pi, count)
    radius = 0.55 * np.sqrt(rng.uniform(0.0, 1.0, count))
    z = rng.uniform(-0.55, 0.55, count)
    return np.column_stack((radius * np.cos(theta), radius * np.sin(theta), z))


def _sample_collar(rng: np.random.Generator, count: int) -> np.ndarray:
    """Sample radial and axial support collars, away from the exact boundary."""

    radial_count = count // 2
    axial_count = count - radial_count

    theta_r = rng.uniform(0.0, 2.0 * np.pi, radial_count)
    radius_r = np.sqrt(rng.uniform(1.50**2, 1.95**2, radial_count))
    z_r = rng.uniform(-1.0, 1.0, radial_count)
    radial = np.column_stack((radius_r * np.cos(theta_r), radius_r * np.sin(theta_r), z_r))

    theta_z = rng.uniform(0.0, 2.0 * np.pi, axial_count)
    radius_z = 0.95 * np.sqrt(rng.uniform(0.0, 1.0, axial_count))
    axial_abs = np.sqrt(rng.uniform(1.50**2, 1.95**2, axial_count))
    axial_sign = np.where(rng.uniform(0.0, 1.0, axial_count) < 0.5, -1.0, 1.0)
    axial = np.column_stack(
        (radius_z * np.cos(theta_z), radius_z * np.sin(theta_z), axial_sign * axial_abs)
    )
    return np.concatenate((radial, axial), axis=0)


def make_samples(
    seed: int,
    counts: Mapping[str, int],
) -> dict[str, np.ndarray]:
    """Make deterministic, disjoint-by-construction regional sample pools."""

    rng = np.random.default_rng(int(seed))
    expected_regions = ("uniform", "core", "collar")
    if tuple(counts) != expected_regions:
        raise ValueError(f"counts must have keys {expected_regions}")
    samples = {
        "uniform": rng.uniform(BOX[:, 0], BOX[:, 1], (int(counts["uniform"]), 3)),
        "core": _sample_core(rng, int(counts["core"])),
        "collar": _sample_collar(rng, int(counts["collar"])),
    }
    for name, points in samples.items():
        _as_points(points, name=f"{name} samples")
    return samples


def _flatten_region_samples(samples: Mapping[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    names = ("uniform", "core", "collar")
    points = np.concatenate([_as_points(samples[name], name=name) for name in names], axis=0)
    labels = np.concatenate([np.full(len(samples[name]), name, dtype=object) for name in names])
    return points, labels


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("values must have shape (n, 3)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _stats(
    target: np.ndarray,
    predicted: np.ndarray,
    *,
    region: str,
) -> dict[str, float | int | str]:
    target = np.asarray(target, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if target.shape != predicted.shape or target.ndim != 2 or target.shape[1] != 3:
        raise ValueError("target and predicted must have matching shape (n, 3)")
    error = predicted - target
    error_norm = np.linalg.norm(error, axis=1)
    target_norm = np.linalg.norm(target, axis=1)
    predicted_norm = np.linalg.norm(predicted, axis=1)
    row: dict[str, float | int | str] = {
        "region": region,
        "sample_count": int(len(error)),
        "curl_momentum_max": float(np.max(target_norm)),
        "curl_force_max": float(np.max(predicted_norm)),
        "curl_obstruction_max": float(np.max(error_norm)),
        "curl_momentum_rms": _vector_rms(target),
        "curl_force_rms": _vector_rms(predicted),
        "curl_obstruction_rms": _vector_rms(error),
    }
    if region == "uniform":
        row["whole_volume_l2_estimate"] = float(np.sqrt(BOX_VOLUME * np.mean(error_norm**2)))
    return row


def _bounded_two_coefficient_fit(
    target: np.ndarray,
    columns: np.ndarray,
    *,
    bounds: tuple[float, float] = FORCE_BOUNDS,
) -> dict[str, Any]:
    """Solve the two-column least-squares problem exactly over a box.

    Enumerating the unconstrained point, the four one-face minimizers, and
    the four corners is sufficient for a convex quadratic in two variables.
    This keeps the [0, 10] restriction explicit and avoids a hidden optimizer
    or a candidate-dependent force basis.
    """

    target = np.asarray(target, dtype=float)
    columns = np.asarray(columns, dtype=float)
    if target.ndim != 2 or target.shape[1] != 3:
        raise ValueError("target must have shape (n, 3)")
    if columns.shape != (len(target), 2, 3):
        raise ValueError("columns must have shape (n, 2, 3)")
    matrix = np.concatenate((columns[:, 0, :], columns[:, 1, :]), axis=0)
    right = np.concatenate((target[:, 0], target[:, 1]))
    # The row-major flattening above is intentionally not used: each Cartesian
    # component is an observation with its matching coefficient columns.
    matrix = columns.transpose(0, 2, 1).reshape(-1, 2)
    right = target.reshape(-1)
    unconstrained, _, rank, singular_values = np.linalg.lstsq(matrix, right, rcond=None)
    lo, hi = (float(bounds[0]), float(bounds[1]))
    if not np.isfinite([lo, hi]).all() or lo > hi:
        raise ValueError("invalid coefficient bounds")

    candidates: list[np.ndarray] = []
    unconstrained = np.asarray(unconstrained, dtype=float)
    if np.all((unconstrained >= lo) & (unconstrained <= hi)):
        candidates.append(unconstrained)

    def append_face(index: int, fixed: float) -> None:
        free = 1 - index
        column_free = matrix[:, free]
        right_free = right - fixed * matrix[:, index]
        denominator = float(np.dot(column_free, column_free))
        value = float(np.dot(column_free, right_free) / denominator) if denominator > 0.0 else lo
        value = float(np.clip(value, lo, hi))
        candidate = np.empty(2, dtype=float)
        candidate[index] = fixed
        candidate[free] = value
        candidates.append(candidate)

    for index in (0, 1):
        append_face(index, lo)
        append_face(index, hi)
    for first in (lo, hi):
        for second in (lo, hi):
            candidates.append(np.asarray((first, second), dtype=float))

    def objective(coefficients: np.ndarray) -> float:
        error = matrix @ coefficients - right
        return float(np.dot(error, error))

    best = min(candidates, key=objective)
    singular_values = np.asarray(singular_values, dtype=float)
    if singular_values.size and singular_values[-1] > 0.0:
        condition = float(singular_values[0] / singular_values[-1])
    else:
        condition = float("inf")
    return {
        "unconstrained_coefficients": unconstrained.tolist(),
        "bounded_coefficients": best.tolist(),
        "bounds": [lo, hi],
        "rank": int(rank),
        "singular_values": singular_values.tolist(),
        "condition_number": condition,
        "objective_sum_squares": objective(best),
        "candidate_count": int(len(candidates)),
    }


def _fit_force_at_resolution(
    field: Any,
    samples: Mapping[str, np.ndarray],
    times: tuple[float, ...],
    *,
    nu: float,
    step: float,
    time_bounds: tuple[float, float],
    force_cache: Mapping[tuple[str, float], np.ndarray],
) -> tuple[dict[str, Any], np.ndarray, np.ndarray, np.ndarray]:
    points, _labels = _flatten_region_samples(samples)
    target_parts: list[np.ndarray] = []
    column_parts: list[np.ndarray] = []
    for time in times:
        target_parts.append(
            curl_of_momentum(
                field,
                points,
                time,
                nu=nu,
                step=step,
                time_bounds=time_bounds,
            )
        )
        basis = np.stack(
            (
                force_cache[("all", float(time))][:, 0],
                force_cache[("all", float(time))][:, 1],
            ),
            axis=1,
        )
        column_parts.append(basis)
    target = np.concatenate(target_parts, axis=0)
    columns = np.concatenate(column_parts, axis=0)
    fit = _bounded_two_coefficient_fit(target, columns)
    coefficient = np.asarray(fit["bounded_coefficients"], dtype=float)
    predicted = np.einsum("nac,a->nc", columns, coefficient)
    fit["training_sample_count"] = int(len(target))
    fit["training_times"] = [float(value) for value in times]
    fit["training_step"] = float(step)
    fit["training_target_rms"] = _vector_rms(target)
    fit["training_fit_rms"] = _vector_rms(predicted - target)
    fit["training_fit_max"] = float(np.max(np.linalg.norm(predicted - target, axis=1)))
    return fit, points, target, predicted


def _force_cache(
    samples: Mapping[str, np.ndarray],
    times: tuple[float, ...],
    *,
    step: float,
) -> dict[tuple[str, float], np.ndarray]:
    points, _labels = _flatten_region_samples(samples)
    cache: dict[tuple[str, float], np.ndarray] = {}
    for time in times:
        cache[("all", float(time))] = np.stack(
            (
                curl_of_force(RestrictedForce(a=1.0, c=0.0), points, time, step=step),
                curl_of_force(RestrictedForce(a=0.0, c=1.0), points, time, step=step),
            ),
            axis=1,
        )
    return cache


def _holdout_rows(
    field: Any,
    coefficients: np.ndarray,
    samples: Mapping[str, np.ndarray],
    times: tuple[float, ...],
    *,
    nu: float,
    step: float,
    time_bounds: tuple[float, float],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for region in ("uniform", "core", "collar"):
        points = _as_points(samples[region], name=region)
        for time in times:
            target = curl_of_momentum(
                field,
                points,
                time,
                nu=nu,
                step=step,
                time_bounds=time_bounds,
            )
            force = RestrictedForce(a=float(coefficients[0]), c=float(coefficients[1]))
            force_curl = curl_of_force(force, points, time, step=step)
            row = _stats(target, force_curl, region=region)
            row.update({"time": float(time), "step": float(step)})
            rows.append(row)
    return rows


def _sample_manifest(
    samples: Mapping[str, np.ndarray],
    *,
    seed: int,
) -> dict[str, Any]:
    return {
        "seed": int(seed),
        "regions": {
            name: {
                "count": int(len(values)),
                "sha256": _sha256_array(values),
            }
            for name, values in samples.items()
        },
    }


def _relative_change(new: float, old: float) -> float:
    return float(abs(float(new) - float(old)) / max(abs(float(old)), 1e-14))


def _resolution_summary(rows: list[dict[str, Any]], steps: tuple[float, ...]) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    times = sorted({float(row["time"]) for row in rows})
    for time in times:
        by_step = {
            float(step): next(
                row
                for row in rows
                if row["region"] == "uniform"
                and float(row["time"]) == float(time)
                and float(row["step"]) == float(step)
            )
            for step in steps
        }
        if len(steps) >= 2:
            coarse = by_step[float(steps[0])]
            fine = by_step[float(steps[-1])]
            summaries.append(
                {
                    "time": float(time),
                    "coarse_step": float(steps[0]),
                    "fine_step": float(steps[-1]),
                    "uniform_obstruction_max_relative_change": _relative_change(
                        fine["curl_obstruction_max"], coarse["curl_obstruction_max"]
                    ),
                    "uniform_whole_volume_l2_relative_change": _relative_change(
                        fine["whole_volume_l2_estimate"], coarse["whole_volume_l2_estimate"]
                    ),
                }
            )
    return summaries


def _resolve_path(path: str | Path, repo_root: Path) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else repo_root / candidate


def _load_candidate(path: Path) -> tuple[Eq45SupportedDeliveryField, dict[str, Any]]:
    field = Eq45SupportedDeliveryField.load_candidate(path)
    metadata = field.metadata()
    return field, metadata


def _candidate_specs(repo_root: Path) -> list[dict[str, Any]]:
    return [
        {
            "name": "default",
            "path": repo_root / "artifacts/delivery/eq45_supported/candidate.json",
            "comparison_role": "required_baseline",
        },
        {
            "name": "bipolar",
            "path": repo_root / "artifacts/delivery/eq45_bipolar/candidate.json",
            "comparison_role": "required_qualitative_variant",
        },
        {
            "name": "bipolar_energy_normalized",
            "path": repo_root / "artifacts/bipolar_energy/normalized_candidate.json",
            "comparison_role": "optional_auxiliary_rescaled_variant",
        },
    ]


def run(
    *,
    repo_root: str | Path | None = None,
    output: str | Path | None = None,
    train_seed: int = DEFAULT_TRAIN_SEED,
    holdout_seed: int = DEFAULT_HOLDOUT_SEED,
    times: tuple[float, ...] = DEFAULT_TIMES,
    train_counts: Mapping[str, int] = DEFAULT_TRAIN_COUNTS,
    holdout_counts: Mapping[str, int] = DEFAULT_HOLDOUT_COUNTS,
    train_step: float = DEFAULT_TRAIN_STEP,
    holdout_steps: tuple[float, ...] = DEFAULT_HOLDOUT_STEPS,
) -> dict[str, Any]:
    """Run the bounded train/holdout comparison and write one JSON report."""

    module_path = Path(__file__).resolve()
    root = Path(repo_root).resolve() if repo_root is not None else module_path.parents[2]
    config_path = root / "configs/constraints.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    nu = float(config["nu"])
    if nu != DEFAULT_NU:
        raise ValueError(f"fixed viscosity must remain {DEFAULT_NU}, got {nu}")
    domain_interval = tuple(float(value) for value in config["domain"]["time_interval"])
    if len(domain_interval) != 2 or not domain_interval[0] < domain_interval[1]:
        raise ValueError("invalid configured time interval")
    time_values = tuple(float(value) for value in times)
    if not time_values or not all(domain_interval[0] < value < domain_interval[1] for value in time_values):
        raise ValueError("obstruction times must lie strictly inside the configured interval")
    step_values = tuple(float(value) for value in holdout_steps)
    if not step_values or not all(np.isfinite(value) and value > 0.0 for value in step_values):
        raise ValueError("holdout_steps must be positive and finite")
    train_step = float(train_step)
    if not np.isfinite(train_step) or train_step <= 0.0:
        raise ValueError("train_step must be positive and finite")

    train_samples = make_samples(int(train_seed), train_counts)
    holdout_samples = make_samples(int(holdout_seed), holdout_counts)
    if _sha256_array(_flatten_region_samples(train_samples)[0]) == _sha256_array(
        _flatten_region_samples(holdout_samples)[0]
    ):
        raise RuntimeError("train and holdout sample pools unexpectedly match")

    train_force_cache = _force_cache(train_samples, time_values, step=train_step)
    candidate_reports: dict[str, Any] = {}
    errors: list[dict[str, str]] = []

    required_names = {"default", "bipolar"}
    for spec in _candidate_specs(root):
        name = str(spec["name"])
        path = Path(spec["path"])
        if not path.is_file():
            if name in required_names:
                raise FileNotFoundError(path)
            errors.append({"candidate": name, "error": f"optional artifact missing: {path}"})
            continue
        try:
            field, metadata = _load_candidate(path)
            fit, _train_points, train_target, train_predicted = _fit_force_at_resolution(
                field,
                train_samples,
                time_values,
                nu=nu,
                step=train_step,
                time_bounds=domain_interval,
                force_cache=train_force_cache,
            )
            coefficients = np.asarray(fit["bounded_coefficients"], dtype=float)
            holdout_rows = []
            for step in step_values:
                holdout_rows.extend(
                    _holdout_rows(
                        field,
                        coefficients,
                        holdout_samples,
                        time_values,
                        nu=nu,
                        step=step,
                        time_bounds=domain_interval,
                    )
                )
            candidate_reports[name] = {
                "path": str(path.relative_to(root)).replace("\\", "/"),
                "comparison_role": spec["comparison_role"],
                "file_sha256": _sha256_file(path),
                "candidate_sha256": field.sha256,
                "parent_sha256": metadata.get("parent_sha256"),
                "metadata_truth_boundary": metadata.get("truth_boundary"),
                "fitted_force": {
                    "family": "RestrictedForce",
                    "coefficients": {"a": float(coefficients[0]), "c": float(coefficients[1])},
                    **fit,
                },
                "training": {
                    "sample_manifest": _sample_manifest(train_samples, seed=int(train_seed)),
                    "rows": {
                        "all_regions": _stats(
                            train_target,
                            train_predicted,
                            region="all_regions",
                        )
                    },
                },
                "holdout": {
                    "sample_manifest": _sample_manifest(holdout_samples, seed=int(holdout_seed)),
                    "rows": holdout_rows,
                    "resolution_summaries": _resolution_summary(holdout_rows, step_values),
                },
            }
        except Exception as exc:  # preserve per-candidate evidence while identifying failures
            errors.append({"candidate": name, "error": f"{type(exc).__name__}: {exc}"})

    # The report is comparison-oriented; rows are retained even when a candidate
    # fit fails so that reruns expose the exact failure rather than hiding it.
    fine_step = float(step_values[-1])
    fine_uniform: dict[str, float] = {}
    for name, report in candidate_reports.items():
        rows = report["holdout"]["rows"]
        values = [
            float(row["whole_volume_l2_estimate"])
            for row in rows
            if row["region"] == "uniform" and float(row["step"]) == fine_step
        ]
        if values:
            fine_uniform[name] = float(max(values))
    ordering = sorted(fine_uniform, key=fine_uniform.get)
    comparison = {
        "fine_step": fine_step,
        "uniform_whole_volume_l2_max_by_candidate": fine_uniform,
        "fine_step_uniform_l2_ordering": ordering,
        "default_vs_bipolar_fine_uniform_l2_relative_change": (
            _relative_change(fine_uniform["bipolar"], fine_uniform["default"])
            if "default" in fine_uniform and "bipolar" in fine_uniform
            else None
        ),
    }
    if "default" in fine_uniform and "bipolar" in fine_uniform:
        lead = "bipolar" if fine_uniform["bipolar"] < fine_uniform["default"] else "default"
        comparison["lower_fine_uniform_l2_candidate"] = lead

    report = {
        "schema": SCHEMA,
        "claim_scope": "pressure_free_necessary_condition_diagnostic_only",
        "equation": "curl(u_t + (u dot grad)u - nu*laplacian(u)) = curl(f)",
        "interpretation": (
            "A nonzero curl obstruction is incompatible with any pressure gradient "
            "at the sampled points and finite-difference resolutions. A small curl "
            "obstruction would remain necessary but would not establish the full NS "
            "equation, divergence, boundary regularity, or a continuum bound."
        ),
        "fixed_constraints": {
            "viscosity_nu": nu,
            "domain": config["domain"],
            "force_family": "RestrictedForce",
            "force_definition": config["forcing"],
            "force_coefficient_bounds": list(FORCE_BOUNDS),
            "pressure": "excluded by taking curl; no pressure coefficients fitted",
            "arbitrary_force": "excluded; only bounded a,c in the preregistered curl-bump family",
        },
        "operator": {
            "momentum": "constrained_validation.residual with pressure=0 and force=0",
            "spatial_momentum_stencil": "fourth_order_cartesian",
            "time_momentum_stencil": "second_order_centered_cartesian",
            "outer_curl_stencil": "second_order_centered_cartesian",
            "train_step": train_step,
            "holdout_steps": list(step_values),
            "time_values": list(time_values),
            "note": "Nested finite differences make this a numerical obstruction estimate; no interval certification is claimed.",
        },
        "sampling": {
            "box": BOX.tolist(),
            "box_volume": BOX_VOLUME,
            "train": _sample_manifest(train_samples, seed=int(train_seed)),
            "holdout": _sample_manifest(holdout_samples, seed=int(holdout_seed)),
            "region_definitions": {
                "uniform": "uniform Cartesian samples in [-2,2]^3; only this region receives the whole-volume Monte Carlo L2 estimate",
                "core": "r <= 0.55 and |z| <= 0.55, sampled in an axis-safe central cylinder",
                "collar": "equal radial collar 1.50 <= r <= 1.95, |z| <= 1 and axial collar 1.50 <= |z| <= 1.95, r <= 0.95",
            },
            "train_holdout_separation": "independent seeds and independently generated point pools",
        },
        "inputs": {
            "config_path": str(config_path.relative_to(root)).replace("\\", "/"),
            "config_sha256": _sha256_file(config_path),
            "candidates": {
                name: {
                    "path": details["path"],
                    "file_sha256": details["file_sha256"],
                    "candidate_sha256": details["candidate_sha256"],
                    "parent_sha256": details["parent_sha256"],
                }
                for name, details in candidate_reports.items()
            },
        },
        "candidates": candidate_reports,
        "comparison": comparison,
        "errors": errors,
        "conclusion": {
            "status": "diagnostic_complete" if not errors else "diagnostic_complete_with_errors",
            "statement": (
                "The fitted force comparison is evidence about pressure-free curl compatibility "
                "for the frozen sampled candidates only. It is not a full NS validation or an "
                "impossibility result for other velocity families."
            ),
            "next_action": (
                "Use the regional and resolution-stability pattern to choose the next bounded "
                "velocity-shape experiment; preserve the force family and independently validate "
                "the resulting full momentum residual."
            ),
        },
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "module_sha256": _sha256_file(module_path),
        },
    }
    target = _resolve_path(output or "artifacts/bipolar_obstruction/report.json", root)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/bipolar_obstruction/report.json")
    parser.add_argument("--repo-root", default=None)
    args = parser.parse_args(argv)
    result = run(repo_root=args.repo_root, output=args.output)
    print(json.dumps(result["comparison"], sort_keys=True))
    print(result["conclusion"]["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
