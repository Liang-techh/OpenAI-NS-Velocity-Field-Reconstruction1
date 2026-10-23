"""Independent bounded-parameter sensitivity audit for the frozen Eq. (4.5) candidate.

This module is deliberately downstream of the serialized public velocity candidate.
It does not inspect optimizer state or profile jets while measuring the response:
every baseline and perturbed value is obtained through ``at_points`` after
loading (or reloading) a candidate artifact.

The audit is a CR009 generalization/sensitivity diagnostic. It is not a
Navier--Stokes validation and does not identify OpenAI's hidden profiles.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import random
from pathlib import Path
import tempfile
from typing import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_profile_basis import Eq45CompactProfileBasis


DEFAULT_PROBE_SEED = 914131
DEFAULT_PERTURBATION_SEEDS = (271828, 314159, 161803)
DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_SCALES = (0.005, 0.010)
DEFAULT_PROBE_COUNT = 128


@dataclass(frozen=True)
class PerturbationResponse:
    seed: int
    scale: float
    coefficient_relative_l2: float
    time: float
    rms_relative_velocity_change: float
    max_change_over_base_rms: float
    response_gain: float


def _coefficient_vector(candidate: Eq45VelocityCandidate) -> np.ndarray:
    basis = candidate.profile_basis
    return np.asarray(basis.phi_coefficients + basis.swirl_coefficients, dtype=float)


def _perturbed_candidate(
    candidate: Eq45VelocityCandidate,
    *,
    seed: int,
    scale: float,
) -> tuple[Eq45VelocityCandidate, float]:
    if not isinstance(seed, (int, np.integer)):
        raise ValueError("perturbation seed must be an integer")
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("perturbation scale must be finite and positive")

    basis = candidate.profile_basis
    base = _coefficient_vector(candidate)
    rng = random.Random(int(seed))
    direction = np.asarray(
        [rng.uniform(-1.0, 1.0) for _ in range(base.size)], dtype=float
    )
    direction /= np.linalg.norm(direction)

    target_l2 = float(scale) * float(basis.coefficient_limit)
    delta = target_l2 * direction
    trial = base + delta
    if np.any(np.abs(trial) > basis.coefficient_limit):
        raise ValueError(
            "requested perturbation leaves the registered coefficient bounds; "
            "reduce scale rather than clipping the direction"
        )

    n = basis.mode_count
    perturbed_basis = Eq45CompactProfileBasis(
        radial_degree=basis.radial_degree,
        eta_degree=basis.eta_degree,
        phi_coefficients=tuple(float(value) for value in trial[:n]),
        swirl_coefficients=tuple(float(value) for value in trial[n:]),
        x_cut=basis.x_cut,
        eta_cut=basis.eta_cut,
        cutoff_power=basis.cutoff_power,
        coefficient_limit=basis.coefficient_limit,
    )
    perturbed = Eq45VelocityCandidate(
        profile_basis=perturbed_basis,
        h=candidate.h,
        time_start=candidate.time_start,
        time_end=candidate.time_end,
    )
    relative_l2 = float(np.linalg.norm(delta) / np.linalg.norm(base))
    return perturbed, relative_l2


def _held_out_points(seed: int, count: int) -> np.ndarray:
    if not isinstance(seed, (int, np.integer)):
        raise ValueError("probe seed must be an integer")
    if not isinstance(count, (int, np.integer)) or count < 16:
        raise ValueError("probe_count must be an integer >= 16")
    rng = random.Random(int(seed))
    return np.column_stack(
        (
            [rng.uniform(-0.9, 0.9) for _ in range(int(count))],
            [rng.uniform(-0.9, 0.9) for _ in range(int(count))],
            [rng.uniform(-0.45, 0.45) for _ in range(int(count))],
        )
    )


def _validate_times(candidate: Eq45VelocityCandidate, times: Iterable[float]) -> tuple[float, ...]:
    result = tuple(float(value) for value in times)
    if len(result) < 3 or not np.all(np.isfinite(result)):
        raise ValueError("times must contain at least three finite values")
    if any(value < candidate.time_start or value > candidate.time_end for value in result):
        raise ValueError("audit times must stay inside the candidate delivery interval")
    return result


def run_eq45_parameter_perturbation_audit(
    artifact_path: str | Path,
    *,
    probe_seed: int = DEFAULT_PROBE_SEED,
    perturbation_seeds: Iterable[int] = DEFAULT_PERTURBATION_SEEDS,
    times: Iterable[float] = DEFAULT_TIMES,
    scales: Iterable[float] = DEFAULT_SCALES,
    probe_count: int = DEFAULT_PROBE_COUNT,
) -> dict:
    """Audit public-velocity sensitivity to unseen bounded coefficient perturbations."""

    source = Path(artifact_path)
    candidate = Eq45VelocityCandidate.load_json(source)
    checked_times = _validate_times(candidate, times)
    checked_scales = tuple(float(value) for value in scales)
    if len(checked_scales) < 2 or any(
        (not np.isfinite(value)) or value <= 0.0 for value in checked_scales
    ):
        raise ValueError("scales must contain at least two finite positive values")
    if tuple(sorted(checked_scales)) != checked_scales:
        raise ValueError("scales must be strictly nondecreasing")
    seeds = tuple(int(value) for value in perturbation_seeds)
    if len(seeds) < 3 or len(set(seeds)) != len(seeds):
        raise ValueError("use at least three distinct perturbation seeds")

    points = _held_out_points(probe_seed, probe_count)
    baseline = {time: candidate.at_points(points, time) for time in checked_times}
    baseline_rms = {
        time: float(np.sqrt(np.mean(np.sum(values * values, axis=-1))))
        for time, values in baseline.items()
    }
    if any(value <= 0.0 or not np.isfinite(value) for value in baseline_rms.values()):
        raise RuntimeError("baseline public velocity must have finite nonzero RMS")

    responses: list[PerturbationResponse] = []
    hashes: dict[str, str] = {}
    with tempfile.TemporaryDirectory(prefix="eq45-perturb-audit-") as directory:
        tempdir = Path(directory)
        for seed in seeds:
            for scale in checked_scales:
                perturbed, coefficient_relative_l2 = _perturbed_candidate(
                    candidate, seed=seed, scale=scale
                )
                target = tempdir / f"candidate-seed{seed}-scale{scale:.6g}.json"
                perturbed.save_json(target)
                reloaded = Eq45VelocityCandidate.load_json(target)
                if reloaded.sha256 != perturbed.sha256:
                    raise RuntimeError("perturbed candidate changed across JSON replay")
                if reloaded.sha256 == candidate.sha256:
                    raise RuntimeError("nonzero perturbation did not change candidate identity")
                hashes[f"{seed}:{scale:.12g}"] = reloaded.sha256

                for time in checked_times:
                    changed = reloaded.at_points(points, time)
                    difference = changed - baseline[time]
                    diff_norm = np.linalg.norm(difference, axis=-1)
                    rms_change = float(np.sqrt(np.mean(diff_norm * diff_norm)))
                    rms_relative = rms_change / baseline_rms[time]
                    max_over_base_rms = float(np.max(diff_norm) / baseline_rms[time])
                    responses.append(
                        PerturbationResponse(
                            seed=seed,
                            scale=scale,
                            coefficient_relative_l2=coefficient_relative_l2,
                            time=time,
                            rms_relative_velocity_change=rms_relative,
                            max_change_over_base_rms=max_over_base_rms,
                            response_gain=rms_relative / coefficient_relative_l2,
                        )
                    )

    by_seed: dict[str, dict[str, float]] = {}
    low_scale, high_scale = checked_scales[0], checked_scales[-1]
    for seed in seeds:
        low = np.asarray(
            [
                row.rms_relative_velocity_change
                for row in responses
                if row.seed == seed and row.scale == low_scale
            ]
        )
        high = np.asarray(
            [
                row.rms_relative_velocity_change
                for row in responses
                if row.seed == seed and row.scale == high_scale
            ]
        )
        ratios = high / low
        by_seed[str(seed)] = {
            "low_scale_mean_rms_relative_change": float(np.mean(low)),
            "high_scale_mean_rms_relative_change": float(np.mean(high)),
            "high_over_low_response_ratio_min": float(np.min(ratios)),
            "high_over_low_response_ratio_max": float(np.max(ratios)),
        }

    low_rows = [row for row in responses if row.scale == low_scale]
    report = {
        "schema": "eq45_parameter_perturbation_audit_v1",
        "candidate_sha256": candidate.sha256,
        "artifact_path": str(source),
        "probe_seed": int(probe_seed),
        "probe_count": int(probe_count),
        "times": list(checked_times),
        "perturbation_seeds": list(seeds),
        "scales_times_coefficient_limit": list(checked_scales),
        "baseline_velocity_rms": {str(time): baseline_rms[time] for time in checked_times},
        "responses": [row.__dict__ for row in responses],
        "low_scale_summary": {
            "rms_relative_change_min": float(
                min(row.rms_relative_velocity_change for row in low_rows)
            ),
            "rms_relative_change_max": float(
                max(row.rms_relative_velocity_change for row in low_rows)
            ),
            "rms_relative_change_mean": float(
                np.mean([row.rms_relative_velocity_change for row in low_rows])
            ),
            "response_gain_min": float(min(row.response_gain for row in low_rows)),
            "response_gain_max": float(max(row.response_gain for row in low_rows)),
            "response_gain_mean": float(np.mean([row.response_gain for row in low_rows])),
        },
        "scale_calibration": by_seed,
        "perturbed_candidate_sha256": hashes,
        "truth_boundary": {
            "velocity_changed_in_repository": False,
            "generalization_diagnostic_only": True,
            "physical_support_validated": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }
    return report


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Audit frozen Eq45 public velocity under bounded coefficient perturbations."
    )
    parser.add_argument("artifact", type=Path)
    args = parser.parse_args()
    print(json.dumps(run_eq45_parameter_perturbation_audit(args.artifact), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
