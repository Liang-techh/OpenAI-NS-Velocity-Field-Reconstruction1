"""Measure the source-required radial compact-stress side effect across spacetime.

Agent-3 #357 established at the reference cell ``(t=.5,z=.08)`` that the
corrected-reader compact stress has a nonzero radial divergence channel

    (div T)_r = d_z sigma_1,

in addition to the theta/e=2 and axial/e=1 rows used by the covariance inverse.
This module makes one narrower routing fact executable: the same radial side
effect is measured on the frozen finite-cycle held-out times and core-safe
axial slices already used by Agent 3's 3x3 second-column screen.

Every cell reuses #357's real Agent-1-leading + routed Agent-2 oscillatory
phase-mean defect, the same compact radial inverse, and the unchanged centered
z-derivative ladder ``(.02,.01,.005)``.  The 3x3 product grid is autonomous
repository engineering; it is not a Kokuno source constant or a PDE acceptance
threshold.  No second oscillatory column, correction velocity, pressure/force
fit, or finite correction cycle is introduced here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .kokuno_multislice_covariance_preflight import HELD_OUT_CYCLE_TIMES
from .kokuno_radial_force_completeness_audit import (
    FINE_PAIR_RELATIVE_STABILITY_TOLERANCE,
    Z_DERIVATIVE_STEP_LADDER,
    generate_actual_core_report as generate_radial_force_cell_report,
)
from .kokuno_spacetime_covariance_preflight import SPATIAL_SCREEN_Z

TASK = "KOKUNO-A3-RADIAL-FORCE-SPACETIME-ENVELOPE-017"


def _validate_axis(values: Iterable[float], *, name: str) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if len(result) < 2 or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain at least two finite values")
    if any(right <= left for left, right in zip(result, result[1:])):
        raise ValueError(f"{name} must be strictly increasing")
    return result


def _cell_key(report: dict[str, Any]) -> tuple[float, float]:
    try:
        inputs = report["inputs"]
        return float(inputs["time"]), float(inputs["z"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("cell report is missing finite time/z inputs") from exc


def summarize_radial_force_spacetime(
    reports: Iterable[dict[str, Any]],
    *,
    times: Iterable[float] = HELD_OUT_CYCLE_TIMES,
    z_values: Iterable[float] = SPATIAL_SCREEN_Z,
) -> dict[str, Any]:
    """Aggregate one #357-style real radial-force report for every grid cell.

    Scientific failures are retained as false routing states rather than raised.
    Structural/provenance inconsistencies fail closed with ``ValueError``.
    """
    frozen_times = _validate_axis(times, name="times")
    frozen_z = _validate_axis(z_values, name="z_values")
    expected = {(time, z) for z in frozen_z for time in frozen_times}

    checked: list[dict[str, Any]] = []
    seen: set[tuple[float, float]] = set()
    leading_sha: str | None = None
    source_identity: tuple[str, str, str, str] | None = None

    for report in reports:
        if not isinstance(report, dict):
            raise TypeError("every cell report must be a dictionary")
        key = _cell_key(report)
        if key in seen:
            raise ValueError("cell time/z pairs must be unique")
        seen.add(key)
        if key not in expected:
            raise ValueError("cell lies outside the frozen spacetime product grid")

        inputs = report.get("inputs", {})
        truth = report.get("truth_boundary", {})
        derivative = report.get("radial_force_derivative_audit", {})
        force = report.get("full_tensor_force_completeness", {})
        source = report.get("source", {})

        if inputs.get("surrogate_defect_used") is not False:
            raise ValueError("spacetime envelope accepts only real candidate defect cells")
        if truth.get("real_candidate_defect_consumed") is not True:
            raise ValueError("cell must certify real candidate defect consumption")
        if truth.get("surrogate_defect_used") is not False:
            raise ValueError("cell truth boundary must reject surrogate defect use")
        if tuple(inputs.get("z_derivative_steps", ())) != tuple(Z_DERIVATIVE_STEP_LADDER):
            raise ValueError("cell changed the frozen z-derivative ladder")
        if float(derivative.get("fine_pair_relative_stability_tolerance", np.nan)) != float(
            FINE_PAIR_RELATIVE_STABILITY_TOLERANCE
        ):
            raise ValueError("cell changed the frozen derivative-stability tolerance")

        cell_sha = str(inputs.get("leading_candidate_sha256", ""))
        if not cell_sha:
            raise ValueError("cell is missing the leading candidate SHA")
        if leading_sha is None:
            leading_sha = cell_sha
        elif cell_sha != leading_sha:
            raise ValueError("leading candidate changed across spacetime cells")

        identity = (
            str(source.get("repository", "")),
            str(source.get("commit", "")),
            str(source.get("path", "")),
            str(source.get("formulas", "")),
        )
        if not all(identity):
            raise ValueError("cell is missing source provenance")
        if source_identity is None:
            source_identity = identity
        elif identity != source_identity:
            raise ValueError("source provenance changed across spacetime cells")

        required_force_keys = (
            "finest_radial_force_rms",
            "finest_radial_force_max_abs",
            "radial_force_over_tangential_force_rms",
            "radial_force_over_raw_ring_radial_defect_rms",
            "radial_force_over_gated_ring_radial_defect_rms",
        )
        values = {name: float(force.get(name, np.nan)) for name in required_force_keys}
        if not np.isfinite(list(values.values())).all() or any(
            value < 0.0 for value in values.values()
        ):
            raise ValueError("cell radial-force metrics must be finite and nonnegative")

        checked.append(
            {
                "time": key[0],
                "z": key[1],
                **values,
                "radial_force_machine_nonzero": bool(
                    force.get("radial_force_machine_nonzero", False)
                ),
                "finest_pair_relative_rms_difference": float(
                    derivative.get("finest_pair_relative_rms_difference", np.inf)
                ),
                "radial_force_derivative_stability_preflight_passed": bool(
                    derivative.get(
                        "radial_force_derivative_stability_preflight_passed", False
                    )
                ),
            }
        )

    if seen != expected:
        missing = sorted(expected - seen)
        raise ValueError(f"missing frozen spacetime cells: {missing}")

    checked.sort(key=lambda item: (item["z"], item["time"]))
    rms_values = np.array([item["finest_radial_force_rms"] for item in checked])
    max_values = np.array([item["finest_radial_force_max_abs"] for item in checked])
    tangential_ratios = np.array(
        [item["radial_force_over_tangential_force_rms"] for item in checked]
    )
    raw_ratios = np.array(
        [item["radial_force_over_raw_ring_radial_defect_rms"] for item in checked]
    )
    gated_ratios = np.array(
        [item["radial_force_over_gated_ring_radial_defect_rms"] for item in checked]
    )
    derivative_differences = np.array(
        [item["finest_pair_relative_rms_difference"] for item in checked]
    )

    all_stable = all(
        item["radial_force_derivative_stability_preflight_passed"] for item in checked
    )
    all_nonzero = all(item["radial_force_machine_nonzero"] for item in checked)
    strongest_index = int(np.argmax(rms_values))

    return {
        "cell_count": len(checked),
        "times": list(frozen_times),
        "z_values": list(frozen_z),
        "leading_candidate_sha256": leading_sha,
        "source_identity": {
            "repository": source_identity[0] if source_identity else None,
            "commit": source_identity[1] if source_identity else None,
            "path": source_identity[2] if source_identity else None,
            "formulas": source_identity[3] if source_identity else None,
        },
        "z_derivative_steps": list(Z_DERIVATIVE_STEP_LADDER),
        "fine_pair_relative_stability_tolerance": FINE_PAIR_RELATIVE_STABILITY_TOLERANCE,
        "all_cells_derivative_stable": bool(all_stable),
        "all_cells_radial_force_machine_nonzero": bool(all_nonzero),
        "finest_pair_relative_rms_difference_max": float(np.max(derivative_differences)),
        "radial_force_rms_min": float(np.min(rms_values)),
        "radial_force_rms_max": float(np.max(rms_values)),
        "radial_force_rms_max_over_min": float(
            np.max(rms_values) / max(float(np.min(rms_values)), np.finfo(float).tiny)
        ),
        "radial_force_max_abs_across_cells": float(np.max(max_values)),
        "radial_force_over_tangential_force_rms_min": float(np.min(tangential_ratios)),
        "radial_force_over_tangential_force_rms_max": float(np.max(tangential_ratios)),
        "radial_force_over_raw_ring_radial_defect_rms_max": float(np.max(raw_ratios)),
        "radial_force_over_gated_ring_radial_defect_rms_max": float(
            np.max(gated_ratios)
        ),
        "strongest_radial_force_cell": {
            "time": checked[strongest_index]["time"],
            "z": checked[strongest_index]["z"],
            "finest_radial_force_rms": checked[strongest_index][
                "finest_radial_force_rms"
            ],
        },
        "cells": checked,
        "interpretation": (
            "This is the conservative spacetime envelope of the source-required radial "
            "d_z sigma_1 side effect generated by the real compact-stress target. It is "
            "not a residual contraction estimate. A future materialized two-column "
            "correction must retain this channel and then be judged by the complete "
            "held-in/held-out Navier-Stokes residual."
        ),
    }


def generate_actual_core_report(
    *,
    output_dir: str | Path = "artifacts/kokuno_agent3/radial_force_spacetime_envelope_v17",
    times: Iterable[float] = HELD_OUT_CYCLE_TIMES,
    z_values: Iterable[float] = SPATIAL_SCREEN_Z,
) -> dict[str, Any]:
    """Run the real #357 radial-force audit on every frozen spacetime cell."""
    frozen_times = _validate_axis(times, name="times")
    frozen_z = _validate_axis(z_values, name="z_values")
    root = Path(output_dir)
    cells_dir = root / "cells"
    cells_dir.mkdir(parents=True, exist_ok=True)

    reports: list[dict[str, Any]] = []
    for z in frozen_z:
        for time in frozen_times:
            label = f"t{time:.3f}_z{z:.3f}".replace(".", "p")
            report = generate_radial_force_cell_report(
                output=cells_dir / f"{label}.json",
                time=time,
                z=z,
                z_steps=Z_DERIVATIVE_STEP_LADDER,
            )
            reports.append(report)

    envelope = summarize_radial_force_spacetime(
        reports,
        times=frozen_times,
        z_values=frozen_z,
    )
    first = reports[0]
    summary: dict[str, Any] = {
        "task": TASK,
        "source": first["source"],
        "inputs": {
            "times": list(frozen_times),
            "z_values": list(frozen_z),
            "z_derivative_steps": list(Z_DERIVATIVE_STEP_LADDER),
            "fine_pair_relative_stability_tolerance": FINE_PAIR_RELATIVE_STABILITY_TOLERANCE,
            "surrogate_defect_used": False,
            "spacetime_grid_origin": (
                "same autonomous held-out-time/core-safe-z product already frozen by "
                "Agent 3's second-column spacetime preflight"
            ),
        },
        "radial_force_spacetime_envelope": envelope,
        "routing": {
            "second_public_covariance_column_available": False,
            "radial_force_spacetime_envelope_ready": True,
            "radial_force_spacetime_envelope_numerically_stable": envelope[
                "all_cells_derivative_stable"
            ],
            "finite_correction_cycle_rerun_allowed": False,
            "next_required": (
                "Agent 2 must still supply a genuinely independent source-motivated "
                "public velocity/covariance column. If it passes the frozen theta/z "
                "bounded-inverse spacetime screen, Agent 3 must materialize the full "
                "correction while retaining the measured radial d_z sigma_1 channel "
                "before any held-in/held-out finite-cycle residual verdict."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_targets_consumed": True,
            "surrogate_defect_used": False,
            "source_radial_force_channel_measured_across_spacetime": True,
            "spacetime_grid_is_autonomous_engineering": True,
            "new_oscillatory_column_constructed": False,
            "agent2_complete_curl_reimplemented": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
            "residual_reduction_claimed": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }

    summary_path = root / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="artifacts/kokuno_agent3/radial_force_spacetime_envelope_v17",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output_dir=args.output_dir)
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
