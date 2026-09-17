"""Bounded Agent-2 phase-offset screen on the callable Kokuno leading core.

This is a deliberately narrow engineering screen.  The corrected Kokuno reader
motivates oscillatory complete curls and a structured phase, but the current
Agent-2 correction still uses an autonomous Cartesian affine phase.  This module
therefore varies only that already-declared phase offset inside its existing
[-pi, pi] bound; wave vector, polarization, support, frequency and amplitude are
frozen.

Training and holdout point sets are disjoint.  The screen uses Agent-1's public
leading-core pressure unchanged and f=0 only as a raw local diagnostic.  The
core is not globally matched, so neither phase selection nor a lower raw
momentum norm assesses the registered full-domain 1e-3 PDE gate.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from math import pi
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_core_composite_checkpoint import KokunoCoreCompositeCandidate
from .kokuno_independent_leading_core import evaluate_fd4
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate


TASK = "K2-OSC-004"
NU = 0.01
AMPLITUDE = 0.125
TRAIN_SEED = 9_172_881
HOLDOUT_SEED = 9_172_891
TRAIN_POINTS = 32
HOLDOUT_POINTS = 64
TIMES = (0.375, 0.5, 0.625)
STEPS = (0.02, 0.01, 0.005)
PHASES = tuple(float(value) for value in np.linspace(-pi, pi, 8, endpoint=False))


def _sample_core_points(seed: int, count: int) -> np.ndarray:
    """Sample a small box whose FD4 stencil stays inside the Agent-1 core."""
    if count <= 0:
        raise ValueError("count must be positive")
    rng = np.random.default_rng(int(seed))
    points = np.empty((int(count), 3), dtype=float)
    points[:, 0] = rng.uniform(-0.11, 0.11, count)
    points[:, 1] = rng.uniform(-0.11, 0.11, count)
    points[:, 2] = rng.uniform(-0.15, 0.15, count)
    return points


def _rms(result: dict[str, np.ndarray]) -> float:
    residual = np.asarray(result["residual"], dtype=float)
    norms2 = np.sum(residual * residual, axis=-1)
    return float(np.sqrt(np.mean(norms2)))


def _maxnorm(result: dict[str, np.ndarray]) -> float:
    residual = np.asarray(result["residual"], dtype=float)
    return float(np.max(np.linalg.norm(residual, axis=-1)))


def _divergence_max(result: dict[str, np.ndarray]) -> float:
    return float(np.max(np.abs(np.asarray(result["divergence"], dtype=float))))


def _composite_velocity(
    leading: KokunoLeadingCoreSeriesCandidate,
    correction: KokunoCompleteCurlCorrection,
):
    def velocity(x: Any, y: Any, z: Any, time: Any) -> np.ndarray:
        base = np.asarray(leading.velocity(x, y, z, time), dtype=float)
        wave = np.asarray(correction.velocity(x, y, z, time), dtype=float)
        if base.shape != wave.shape:
            raise ValueError("leading and correction velocity shapes differ")
        out = base + wave
        if not np.all(np.isfinite(out)):
            raise RuntimeError("phase-screen composite velocity became nonfinite")
        return out

    return velocity


def evaluate_phase_grid(
    leading: KokunoLeadingCoreSeriesCandidate,
    points: np.ndarray,
    *,
    phases: Sequence[float] = PHASES,
    times: Sequence[float] = TIMES,
    steps: Sequence[float] = STEPS,
    template: KokunoCompleteCurlCorrection | None = None,
) -> list[dict[str, float]]:
    """Evaluate a predeclared one-dimensional phase grid without fitting it."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) == 0:
        raise ValueError("points must have shape (N,3)")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must be finite")

    phase_values = tuple(float(value) for value in phases)
    time_values = tuple(float(value) for value in times)
    step_values = tuple(float(value) for value in steps)
    if not phase_values or not time_values or not step_values:
        raise ValueError("phases, times and steps must be nonempty")
    if any(not np.isfinite(value) or value < -pi or value > pi for value in phase_values):
        raise ValueError("phase values must lie in [-pi, pi]")

    base_template = template or KokunoCompleteCurlCorrection(amplitude=AMPLITUDE)
    if abs(base_template.amplitude - AMPLITUDE) > 1.0e-15:
        raise ValueError("this governed screen freezes amplitude at 0.125")

    rows: list[dict[str, float]] = []
    for phase in phase_values:
        correction = replace(base_template, phase=phase)
        velocity = _composite_velocity(leading, correction)
        for time in time_values:
            for step in step_values:
                result = evaluate_fd4(
                    velocity,
                    leading.pressure,
                    points,
                    time,
                    step,
                    nu=NU,
                )
                rows.append(
                    {
                        "phase": float(phase),
                        "time": float(time),
                        "step": float(step),
                        "momentum_rms": _rms(result),
                        "momentum_max": _maxnorm(result),
                        "divergence_max": _divergence_max(result),
                    }
                )
    return rows


def _aggregate_finest(rows: list[dict[str, float]]) -> dict[str, dict[str, Any]]:
    finest = min(STEPS)
    result: dict[str, dict[str, Any]] = {}
    for phase in PHASES:
        selected = [
            row for row in rows
            if row["step"] == finest and abs(row["phase"] - phase) < 1.0e-15
        ]
        if len(selected) != len(TIMES):
            raise RuntimeError("phase screen did not produce every governed time")
        result[repr(float(phase))] = {
            "phase": float(phase),
            "mean_momentum_rms": float(np.mean([row["momentum_rms"] for row in selected])),
            "max_momentum": float(np.max([row["momentum_max"] for row in selected])),
            "max_divergence": float(np.max([row["divergence_max"] for row in selected])),
            "per_time_momentum_rms": {
                repr(float(row["time"])): float(row["momentum_rms"]) for row in selected
            },
        }
    return result


def _select_training_phase(aggregates: dict[str, dict[str, Any]]) -> float:
    """Finite deterministic selection: minimum training mean, then |phase|."""
    choices = [entry for entry in aggregates.values()]
    best = min(
        choices,
        key=lambda entry: (
            float(entry["mean_momentum_rms"]),
            abs(float(entry["phase"])),
            float(entry["phase"]),
        ),
    )
    return float(best["phase"])


def generate_core_phase_report(
    *,
    output_dir: str | Path = "artifacts/kokuno_agent2/core_phase_screen_004",
    train_seed: int = TRAIN_SEED,
    holdout_seed: int = HOLDOUT_SEED,
    train_count: int = TRAIN_POINTS,
    holdout_count: int = HOLDOUT_POINTS,
) -> dict[str, Any]:
    """Select phase on training points, then evaluate it on disjoint holdout."""
    if int(train_seed) == int(holdout_seed):
        raise ValueError("training and holdout seeds must differ")

    leading = KokunoLeadingCoreSeriesCandidate()
    template = KokunoCompleteCurlCorrection(amplitude=AMPLITUDE)
    train_points = _sample_core_points(train_seed, train_count)
    holdout_points = _sample_core_points(holdout_seed, holdout_count)

    train_rows = evaluate_phase_grid(leading, train_points, template=template)
    train_aggregate = _aggregate_finest(train_rows)
    selected_phase = _select_training_phase(train_aggregate)

    selected = replace(template, phase=selected_phase)
    baseline_phase = replace(template, phase=0.0)
    selected_velocity = _composite_velocity(leading, selected)
    baseline_velocity = _composite_velocity(leading, baseline_phase)

    holdout_rows: list[dict[str, Any]] = []
    for time in TIMES:
        for step in STEPS:
            leading_result = evaluate_fd4(
                leading.velocity, leading.pressure, holdout_points, time, step, nu=NU
            )
            phase0_result = evaluate_fd4(
                baseline_velocity, leading.pressure, holdout_points, time, step, nu=NU
            )
            selected_result = evaluate_fd4(
                selected_velocity, leading.pressure, holdout_points, time, step, nu=NU
            )
            holdout_rows.append(
                {
                    "time": float(time),
                    "step": float(step),
                    "leading_only_momentum_rms": _rms(leading_result),
                    "phase0_momentum_rms": _rms(phase0_result),
                    "selected_phase_momentum_rms": _rms(selected_result),
                    "selected_over_phase0": _rms(selected_result) / _rms(phase0_result),
                    "selected_over_leading": _rms(selected_result) / _rms(leading_result),
                    "selected_divergence_max": _divergence_max(selected_result),
                }
            )

    finest = min(STEPS)
    finest_rows = [row for row in holdout_rows if row["step"] == finest]
    phase0_mean = float(np.mean([row["phase0_momentum_rms"] for row in finest_rows]))
    selected_mean = float(
        np.mean([row["selected_phase_momentum_rms"] for row in finest_rows])
    )
    leading_mean = float(
        np.mean([row["leading_only_momentum_rms"] for row in finest_rows])
    )
    selected_generalizes_vs_phase0 = bool(selected_mean < phase0_mean)

    candidate = KokunoCoreCompositeCandidate(leading=leading, oscillatory=selected)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_path = candidate.save_json(output_dir / "training_selected_phase_core_candidate.json")
    replay = KokunoCoreCompositeCandidate.load_json(candidate_path)
    if replay.sha256 != candidate.sha256:
        raise RuntimeError("training-selected phase candidate did not round-trip")

    report = {
        "task": TASK,
        "source": {
            "reader": "KokunoYumeto corrected 208-page reconstruction",
            "edition": "2026.09.09-consolidated",
            "doi": "10.5281/zenodo.22678406",
            "source_structure_used": [
                "oscillatory velocity represented as a complete curl",
                "potential-level localization retains cutoff-gradient curl remainder",
            ],
            "autonomous_choices_retained": [
                "Cartesian affine phase/frame",
                "C4 Cartesian box support",
                "wave vector, polarization, frequency and support widths",
                "amplitude=0.125 diagnostic probe",
                "eight-point phase-offset grid",
            ],
            "paper_exact": False,
        },
        "screen": {
            "parameter": "phase_offset",
            "bound": [-pi, pi],
            "grid": list(PHASES),
            "amplitude_frozen": AMPLITUDE,
            "training_seed": int(train_seed),
            "training_points": int(train_count),
            "holdout_seed": int(holdout_seed),
            "holdout_points": int(holdout_count),
            "times": list(TIMES),
            "steps": list(STEPS),
            "nu": NU,
            "pressure": "Agent-1 leading-core pressure unchanged",
            "forcing": "zero raw diagnostic; no force fitted",
            "selection_rule": "minimum training mean momentum RMS at h=0.005; tie by |phase| then phase",
            "holdout_acceptance_rule": "training-selected phase must improve disjoint holdout mean momentum RMS relative to frozen phase=0",
        },
        "training": {
            "finest_step": finest,
            "aggregate_by_phase": train_aggregate,
            "selected_phase": selected_phase,
        },
        "holdout": {
            "rows": holdout_rows,
            "finest_step": finest,
            "mean_leading_only_momentum_rms": leading_mean,
            "mean_phase0_momentum_rms": phase0_mean,
            "mean_selected_phase_momentum_rms": selected_mean,
            "selected_over_phase0": selected_mean / max(phase0_mean, np.finfo(float).tiny),
            "selected_over_leading": selected_mean / max(leading_mean, np.finfo(float).tiny),
            "selected_generalizes_vs_phase0": selected_generalizes_vs_phase0,
            "max_selected_divergence": float(
                max(row["selected_divergence_max"] for row in finest_rows)
            ),
        },
        "candidate": {
            "path": str(candidate_path),
            "sha256": candidate.sha256,
            "velocity_api": "velocity(x,y,z,t)->[...,3] Cartesian [u,v,w]",
            "training_selected_diagnostic_only": True,
            "accepted_for_next_cycle": selected_generalizes_vs_phase0,
            "global_candidate_promoted": False,
        },
        "routing": {
            "retain_phase_zero_if_holdout_rejects_training_winner": True,
            "recommended_phase_offset": selected_phase if selected_generalizes_vs_phase0 else 0.0,
            "reason": (
                "training-selected phase improved the disjoint holdout relative to phase=0"
                if selected_generalizes_vs_phase0
                else "training-selected phase did not improve the disjoint holdout relative to phase=0"
            ),
        },
        "truth_boundary": {
            "training_holdout_separated": True,
            "training_phase_selected": True,
            "phase_accepted_for_next_cycle": selected_generalizes_vs_phase0,
            "global_outer_matching_complete": False,
            "mean_correction_applied": False,
            "compatible_final_forcing_attached": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    (output_dir / "phase_screen_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Screen bounded Agent-2 phase offset on the Agent-1 leading core")
    parser.add_argument(
        "--output-dir", default="artifacts/kokuno_agent2/core_phase_screen_004"
    )
    parser.add_argument("--train-points", type=int, default=TRAIN_POINTS)
    parser.add_argument("--holdout-points", type=int, default=HOLDOUT_POINTS)
    args = parser.parse_args()
    report = generate_core_phase_report(
        output_dir=args.output_dir,
        train_count=args.train_points,
        holdout_count=args.holdout_points,
    )
    print(
        json.dumps(
            {
                "selected_phase": report["training"]["selected_phase"],
                "holdout": report["holdout"],
                "candidate": report["candidate"],
                "routing": report["routing"],
                "truth_boundary": report["truth_boundary"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
