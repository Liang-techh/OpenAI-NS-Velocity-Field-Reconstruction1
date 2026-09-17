"""Independent public-output time-derivative audit for the cubic Phi(1,0) trial.

The cubic-localized temporal candidate is attractive for visualization because it
returns exactly to the static supported field at several key frames.  Snapshot
identity, however, does not imply equality of the time derivative.  This module
checks that distinction without reusing any training loss, pressure fit, force
fit, or internal analytic velocity derivative.

All candidate derivatives are reconstructed only from serialized/reloaded
public ``at_points(...)->[u,v,w]`` evaluations.  A three-level second-order time
finite-difference ladder is compared with a separate fourth-order finer
reference.  The report is a convergence/generalization diagnostic only; it is
not the preregistered full momentum PDE gate.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import tempfile
from typing import Callable, Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_cubic_temporal_mode import (
    Eq45SupportedPhi10CubicLocalizedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_early_localized_temporal_mode import (
    Eq45SupportedPhi10EarlyLocalizedTemporalCandidate,
)


SCHEMA = "eq45_supported_phi10_cubic_temporal_derivative_audit_v1"
PROBE_SEED = 914509
PROBE_COUNT = 24
TIME_LEVELS = (0.01, 0.005, 0.0025)
REFERENCE_STEP = 0.00125
AUDIT_TIMES = (0.25, 0.375, 0.50, 0.625, 0.6875, 0.75)

_TRUTH_BOUNDARY = {
    "serialized_candidate_reloaded": True,
    "public_velocity_only": True,
    "training_loss_reused": False,
    "pressure_fitted": False,
    "forcing_fitted": False,
    "formal_pde_gate_assessed": False,
    "visualization_candidate_only": True,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


@dataclass(frozen=True)
class ErrorMetrics:
    sampled_max_vector: float
    sampled_rms: float
    normalized_rms: float

    def to_dict(self) -> dict[str, float]:
        return {
            "sampled_max_vector": float(self.sampled_max_vector),
            "sampled_rms": float(self.sampled_rms),
            "normalized_rms": float(self.normalized_rms),
        }


def _fresh_probe_points(seed: int = PROBE_SEED) -> tuple[np.ndarray, tuple[str, ...]]:
    """Create deterministic off-grid probes across plateau and support collars."""
    rng = np.random.default_rng(int(seed))
    points: list[list[float]] = []
    regions: list[str] = []

    def add_polar(radius: float, angle: float, z: float, region: str) -> None:
        points.append([radius * np.cos(angle), radius * np.sin(angle), z])
        regions.append(region)

    for _ in range(8):
        r = np.sqrt(rng.uniform(0.04, 1.0)) * 1.30
        add_polar(r, rng.uniform(-np.pi, np.pi), rng.uniform(-1.30, 1.30), "plateau")
    for _ in range(6):
        add_polar(
            rng.uniform(1.65, 1.85),
            rng.uniform(-np.pi, np.pi),
            rng.uniform(-1.25, 1.25),
            "radial_collar",
        )
    for _ in range(6):
        z = rng.choice((-1.0, 1.0)) * rng.uniform(1.65, 1.85)
        add_polar(
            rng.uniform(0.20, 1.30),
            rng.uniform(-np.pi, np.pi),
            z,
            "axial_collar",
        )
    for _ in range(4):
        z = rng.choice((-1.0, 1.0)) * rng.uniform(1.65, 1.85)
        add_polar(
            rng.uniform(1.65, 1.85),
            rng.uniform(-np.pi, np.pi),
            z,
            "corner_collar",
        )

    array = np.asarray(points, dtype=float)
    if array.shape != (PROBE_COUNT, 3) or not np.all(np.isfinite(array)):
        raise RuntimeError("fresh temporal-derivative probe construction failed")
    return array, tuple(regions)


def _evaluate_public(candidate, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(candidate.at_points(points, float(time)), dtype=float)
    if values.shape != points.shape:
        raise ValueError("public velocity evaluator must return points.shape")
    if not np.all(np.isfinite(values)):
        raise ValueError("public velocity evaluator returned nonfinite values")
    return values


def _second_order_time_derivative(candidate, points: np.ndarray, time: float, step: float) -> np.ndarray:
    """Second-order derivative, using centered or endpoint one-sided formulas."""
    h = float(step)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("time step must be positive and finite")
    start = float(candidate.time_start)
    end = float(candidate.time_end)
    t = float(time)
    if t < start or t > end:
        raise ValueError("audit time lies outside candidate interval")

    if np.isclose(t, start, rtol=0.0, atol=1e-14):
        if t + 2.0 * h > end:
            raise ValueError("forward stencil leaves candidate interval")
        f0 = _evaluate_public(candidate, points, t)
        f1 = _evaluate_public(candidate, points, t + h)
        f2 = _evaluate_public(candidate, points, t + 2.0 * h)
        return (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * h)
    if np.isclose(t, end, rtol=0.0, atol=1e-14):
        if t - 2.0 * h < start:
            raise ValueError("backward stencil leaves candidate interval")
        f0 = _evaluate_public(candidate, points, t)
        f1 = _evaluate_public(candidate, points, t - h)
        f2 = _evaluate_public(candidate, points, t - 2.0 * h)
        return (3.0 * f0 - 4.0 * f1 + f2) / (2.0 * h)
    if t - h < start or t + h > end:
        raise ValueError("centered stencil leaves candidate interval")
    return (
        _evaluate_public(candidate, points, t + h)
        - _evaluate_public(candidate, points, t - h)
    ) / (2.0 * h)


def _fourth_order_time_derivative(candidate, points: np.ndarray, time: float, step: float) -> np.ndarray:
    """Independent finer fourth-order reference derivative from public velocity."""
    h = float(step)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("reference time step must be positive and finite")
    start = float(candidate.time_start)
    end = float(candidate.time_end)
    t = float(time)
    if t < start or t > end:
        raise ValueError("audit time lies outside candidate interval")

    if np.isclose(t, start, rtol=0.0, atol=1e-14):
        if t + 4.0 * h > end:
            raise ValueError("forward reference stencil leaves candidate interval")
        samples = [_evaluate_public(candidate, points, t + j * h) for j in range(5)]
        return (-25.0 * samples[0] + 48.0 * samples[1] - 36.0 * samples[2] + 16.0 * samples[3] - 3.0 * samples[4]) / (12.0 * h)
    if np.isclose(t, end, rtol=0.0, atol=1e-14):
        if t - 4.0 * h < start:
            raise ValueError("backward reference stencil leaves candidate interval")
        samples = [_evaluate_public(candidate, points, t - j * h) for j in range(5)]
        return (25.0 * samples[0] - 48.0 * samples[1] + 36.0 * samples[2] - 16.0 * samples[3] + 3.0 * samples[4]) / (12.0 * h)
    if t - 2.0 * h < start or t + 2.0 * h > end:
        raise ValueError("centered reference stencil leaves candidate interval")
    fm2 = _evaluate_public(candidate, points, t - 2.0 * h)
    fm1 = _evaluate_public(candidate, points, t - h)
    fp1 = _evaluate_public(candidate, points, t + h)
    fp2 = _evaluate_public(candidate, points, t + 2.0 * h)
    return (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h)


def _metrics(error: np.ndarray, reference: np.ndarray) -> ErrorMetrics:
    error = np.asarray(error, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if error.shape != reference.shape or error.shape[-1] != 3:
        raise ValueError("error/reference arrays must have matching (...,3) shape")
    rms = float(np.sqrt(np.mean(error * error)))
    reference_rms = float(np.sqrt(np.mean(reference * reference)))
    normalized = rms / max(reference_rms, np.finfo(float).tiny)
    point_norm = np.linalg.norm(error.reshape((-1, 3)), axis=1)
    return ErrorMetrics(
        sampled_max_vector=float(np.max(point_norm)),
        sampled_rms=rms,
        normalized_rms=float(normalized),
    )


def _derivative_stack(candidate, points: np.ndarray, times: Iterable[float], step: float, *, reference: bool) -> np.ndarray:
    derivative = _fourth_order_time_derivative if reference else _second_order_time_derivative
    return np.stack([derivative(candidate, points, float(time), step) for time in times], axis=0)


def _convergence_report(candidate, points: np.ndarray) -> dict[str, object]:
    reference = _derivative_stack(candidate, points, AUDIT_TIMES, REFERENCE_STEP, reference=True)
    levels: list[dict[str, object]] = []
    errors: list[float] = []
    for step in TIME_LEVELS:
        estimate = _derivative_stack(candidate, points, AUDIT_TIMES, step, reference=False)
        metrics = _metrics(estimate - reference, reference)
        levels.append({"step": float(step), **metrics.to_dict()})
        errors.append(metrics.sampled_rms)

    observed_orders = []
    for left, right in zip(errors[:-1], errors[1:]):
        observed_orders.append(float(np.log(left / right) / np.log(2.0)))

    per_time_finest: dict[str, object] = {}
    finest = float(TIME_LEVELS[-1])
    for index, time in enumerate(AUDIT_TIMES):
        estimate = _second_order_time_derivative(candidate, points, time, finest)
        metrics = _metrics(estimate - reference[index], reference[index])
        per_time_finest[f"{time:.6f}"] = metrics.to_dict()

    return {
        "levels": levels,
        "observed_orders": observed_orders,
        "per_time_finest": per_time_finest,
        "reference_step": float(REFERENCE_STEP),
        "reference_operator": "fourth_order_public_velocity_time_difference",
    }


def _rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def _snapshot_dynamic_comparison(cubic, quadratic, base, points: np.ndarray) -> list[dict[str, object]]:
    """Compare fields that are exactly equal at a frame but can have different u_t."""
    comparisons = [
        (0.25, "cubic", cubic, "quadratic", quadratic),
        (0.50, "cubic", cubic, "base", base),
        (0.50, "quadratic", quadratic, "base", base),
        (0.625, "cubic", cubic, "base", base),
        (0.75, "cubic", cubic, "base", base),
        (0.75, "quadratic", quadratic, "base", base),
    ]
    rows: list[dict[str, object]] = []
    for time, left_name, left, right_name, right in comparisons:
        instantaneous = _evaluate_public(left, points, time) - _evaluate_public(right, points, time)
        left_dt = _fourth_order_time_derivative(left, points, time, REFERENCE_STEP)
        right_dt = _fourth_order_time_derivative(right, points, time, REFERENCE_STEP)
        derivative_difference = left_dt - right_dt
        denom = max(_rms(right_dt), np.finfo(float).tiny)
        rows.append(
            {
                "time": float(time),
                "left": left_name,
                "right": right_name,
                "instantaneous_velocity_rms_difference": _rms(instantaneous),
                "instantaneous_velocity_max_vector_difference": float(
                    np.max(np.linalg.norm(instantaneous, axis=1))
                ),
                "time_derivative_rms_difference": _rms(derivative_difference),
                "time_derivative_difference_over_right_rms": float(
                    _rms(derivative_difference) / denom
                ),
            }
        )
    return rows


def _manufactured_sign_mutation_calibration() -> dict[str, float]:
    """Calibrate sign sensitivity on an analytic public-output-style evaluator."""
    points = np.array([[0.2, -0.3, 0.4], [-0.7, 0.5, -0.1]], dtype=float)
    t = 0.4
    h = 0.005

    def good(time: float) -> np.ndarray:
        return np.column_stack(
            (
                points[:, 0] + 2.0 * time,
                points[:, 1] - 3.0 * time,
                points[:, 2] + 0.5 * time,
            )
        )

    def mutated(time: float) -> np.ndarray:
        return np.column_stack(
            (
                points[:, 0] + 2.0 * time,
                points[:, 1] + 3.0 * time,
                points[:, 2] + 0.5 * time,
            )
        )

    exact = np.tile(np.array([2.0, -3.0, 0.5], dtype=float), (len(points), 1))
    good_dt = (good(t + h) - good(t - h)) / (2.0 * h)
    mutant_dt = (mutated(t + h) - mutated(t - h)) / (2.0 * h)
    good_error = _metrics(good_dt - exact, exact)
    mutant_error = _metrics(mutant_dt - exact, exact)
    return {
        "good_sampled_rms_error": good_error.sampled_rms,
        "mutated_sampled_rms_error": mutant_error.sampled_rms,
        "mutated_sampled_max_vector_error": mutant_error.sampled_max_vector,
    }


def audit_reloaded_candidate(candidate: Eq45SupportedPhi10CubicLocalizedTemporalCandidate) -> dict[str, object]:
    if not isinstance(candidate, Eq45SupportedPhi10CubicLocalizedTemporalCandidate):
        raise TypeError("candidate must be a reloaded cubic-localized temporal candidate")
    points, regions = _fresh_probe_points()
    base = candidate.base
    quadratic = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(
        base=base, early_delta=float(candidate.early_delta)
    )

    report = {
        "schema": SCHEMA,
        "candidate_sha256": candidate.sha256,
        "base_sha256": base.sha256,
        "quadratic_sha256": quadratic.sha256,
        "probe_seed": PROBE_SEED,
        "probe_count": int(points.shape[0]),
        "probe_region_counts": {
            name: int(sum(region == name for region in regions))
            for name in sorted(set(regions))
        },
        "audit_times": [float(value) for value in AUDIT_TIMES],
        "time_levels": [float(value) for value in TIME_LEVELS],
        "cubic_convergence": _convergence_report(candidate, points),
        "quadratic_convergence": _convergence_report(quadratic, points),
        "base_convergence": _convergence_report(base, points),
        "same_snapshot_dynamic_comparisons": _snapshot_dynamic_comparison(
            candidate, quadratic, base, points
        ),
        "manufactured_sign_mutation": _manufactured_sign_mutation_calibration(),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
    return report


def build_default_report() -> dict[str, object]:
    """Serialize/reload the production-shaped cubic candidate, then audit it."""
    base = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    candidate = Eq45SupportedPhi10CubicLocalizedTemporalCandidate(base=base)
    with tempfile.TemporaryDirectory(prefix="eq45-cubic-ut-") as directory:
        path = Path(directory) / "candidate.json"
        candidate.save_json(path)
        loaded = Eq45SupportedPhi10CubicLocalizedTemporalCandidate.load_json(path)
        if loaded.sha256 != candidate.sha256 or loaded.to_dict() != candidate.to_dict():
            raise RuntimeError("serialized/reloaded cubic candidate identity drifted")
        return audit_reloaded_candidate(loaded)


def main() -> None:
    print(json.dumps(build_default_report(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
