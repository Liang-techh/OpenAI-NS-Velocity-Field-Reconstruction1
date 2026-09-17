"""Localize the independent leading-core residual in cylindrical components.

This is a reporting companion to :mod:`kokuno_independent_leading_core`.  It
reuses only that independent FD operator and the Agent-1 public velocity and
pressure callables.  No candidate-side derivative helper is introduced here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_independent_leading_core import (
    DEFAULT_TIMES,
    NU,
    VALIDATION_SEED,
    evaluate_fd4,
    sample_held_out_core_points,
)
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate

SCHEMA = "kokuno-independent-leading-core-hotspots-v1"
FINEST_STEP = 0.005


def _cylindrical_residual(points: np.ndarray, residual: np.ndarray) -> np.ndarray:
    radius = np.hypot(points[:, 0], points[:, 1])
    if np.any(radius <= 0.0):
        raise ValueError("hotspot points must stay off the exact axis")
    radial = (points[:, 0] * residual[:, 0] + points[:, 1] * residual[:, 1]) / radius
    theta = (-points[:, 1] * residual[:, 0] + points[:, 0] * residual[:, 1]) / radius
    return np.column_stack((radial, theta, residual[:, 2]))


def _component_metrics(points: np.ndarray, residual: np.ndarray) -> dict[str, Any]:
    cylindrical = _cylindrical_residual(points, residual)
    labels = ("radial", "theta", "axial")
    metrics: dict[str, Any] = {}
    for index, label in enumerate(labels):
        values = cylindrical[:, index]
        absolute = np.abs(values)
        worst = int(np.argmax(absolute))
        metrics[label] = {
            "max_abs": float(absolute[worst]),
            "rms": float(np.sqrt(np.mean(values * values))),
            "worst_point_xyz": [float(value) for value in points[worst]],
            "worst_signed_value": float(values[worst]),
        }
    norm = np.linalg.norm(residual, axis=1)
    worst_norm = int(np.argmax(norm))
    metrics["vector"] = {
        "max_norm": float(norm[worst_norm]),
        "rms_norm": float(np.sqrt(np.mean(norm * norm))),
        "worst_point_xyz": [float(value) for value in points[worst_norm]],
        "worst_cylindrical_components": [
            float(value) for value in cylindrical[worst_norm]
        ],
    }
    return metrics


def run_hotspot_report(
    *,
    interior_count: int = 48,
    axis_near_count: int = 24,
    seed: int = VALIDATION_SEED,
    times: tuple[float, ...] = DEFAULT_TIMES,
    step: float = FINEST_STEP,
) -> dict[str, Any]:
    candidate = KokunoLeadingCoreSeriesCandidate()
    samples = sample_held_out_core_points(
        interior_count=interior_count,
        axis_near_count=axis_near_count,
        seed=seed,
    )
    rows: list[dict[str, Any]] = []
    for time in times:
        stratum_reports: dict[str, Any] = {}
        for name, points in samples.items():
            values = evaluate_fd4(
                candidate.velocity,
                candidate.pressure,
                points,
                float(time),
                float(step),
                nu=NU,
            )
            stratum_reports[name] = {
                "residual": _component_metrics(points, values["residual"]),
                "divergence": {
                    "max_abs": float(np.max(np.abs(values["divergence"]))),
                    "rms": float(np.sqrt(np.mean(values["divergence"] ** 2))),
                },
            }
        rows.append({"time": float(time), "strata": stratum_reports})

    all_component_rms = {
        label: max(
            row["strata"][stratum]["residual"][label]["rms"]
            for row in rows
            for stratum in ("interior", "axis_near")
        )
        for label in ("radial", "theta", "axial")
    }
    dominant = max(all_component_rms, key=all_component_rms.get)
    return {
        "schema": SCHEMA,
        "task_id": "KOKUNO-VAL-LEADING-CORE-MOMENTUM-PREFLIGHT-003",
        "candidate_sha256": candidate.sha256,
        "operator": "same_independent_fd4_public_velocity_pressure_only",
        "step": float(step),
        "seed": int(seed),
        "rows": rows,
        "worst_component_rms_over_time_and_strata": {
            key: float(value) for key, value in all_component_rms.items()
        },
        "dominant_component_by_worst_rms": dominant,
        "truth_boundary": {
            "local_core_hotspot_diagnostic_only": True,
            "formal_full_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }


def write_report(path: str | Path, **kwargs: Any) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(run_hotspot_report(**kwargs), indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_independent_validation/leading_core_hotspots.json",
    )
    args = parser.parse_args(argv)
    path = write_report(args.output)
    report = json.loads(path.read_text(encoding="utf-8"))
    print(
        json.dumps(
            {
                "dominant_component_by_worst_rms": report[
                    "dominant_component_by_worst_rms"
                ],
                "worst_component_rms_over_time_and_strata": report[
                    "worst_component_rms_over_time_and_strata"
                ],
            },
            sort_keys=True,
        )
    )
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
