"""Preflight a future second Kokuno covariance column across time and axial slices.

Agent-3 #329 screens the real compact theta/axial mean-stress target at the
frozen finite-cycle held-out times ``t=(.375,.5,.625)`` but only on one axial
slice ``z=.08``.  This module adds one fail-closed repository-engineering
screen: the same physical second velocity column must span the real target on a
small preregistered core-safe axial set ``z=(.06,.08,.10)`` at every frozen
time before another finite correction cycle may be run.

The signed-coefficient budget is *not* re-fit on the new slices.  It remains the
existing Agent-3 #302/#329 budget measured at the reference cell
``(t=.5,z=.08)``.  Every target is reconstructed from the actual Agent-1
leading candidate plus the routed Agent-2 oscillatory field through the existing
real phase-mean defect -> compact radial-stress pipeline; no surrogate defect is
accepted.

This is only a routing/preflight layer.  It does not construct Agent-2's missing
independent source pair, reimplement complete curl/pulse sensitivities, fit
pressure/forcing, change a PDE threshold, or claim residual contraction.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_covariance_tangent_preflight import (
    complete_curl_phase_evaluator,
    evaluate_covariance_tangent_preflight,
)
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    DEFAULT_RADIAL_COUNT,
    PROFILE_Z,
    MissingCovarianceColumnTarget,
)
from .kokuno_multislice_covariance_preflight import (
    HELD_OUT_CYCLE_TIMES,
    REFERENCE_BUDGET_TIME,
    build_real_covariance_slice,
)
from .kokuno_second_column_bounded_inverse import (
    current_signed_coefficient_budget,
    required_second_column_envelope,
)
from .kokuno_signed_covariance_inverse import (
    PROFILE_ANNULUS,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)
from .kokuno_velocity_column_preflight import (
    PhaseVelocityEvaluator,
    measure_multiplicative_covariance_response,
)

TASK = "KOKUNO-A3-SPACETIME-COVARIANCE-PREFLIGHT-014"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_SECTION = "Signed covariance / compact stress / finite correction iteration"

# Autonomous repository sampling choices, frozen before this report is run.
# They are not source constants and are not PDE acceptance thresholds.
SPATIAL_SCREEN_Z = (0.06, float(PROFILE_Z), 0.10)
REFERENCE_BUDGET_Z = float(PROFILE_Z)


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("values must have shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _validate_axis_values(values: Iterable[float], *, name: str) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if len(result) < 2 or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain at least two finite values")
    if any(right <= left for left, right in zip(result, result[1:])):
        raise ValueError(f"{name} must be strictly increasing")
    return result


def _validate_cells(
    cells: Iterable[tuple[float, float, MissingCovarianceColumnTarget]],
) -> tuple[tuple[float, float, MissingCovarianceColumnTarget], ...]:
    result = tuple((float(time), float(z), receipt) for time, z, receipt in cells)
    if not result:
        raise ValueError("cells must not be empty")
    pairs: set[tuple[float, float]] = set()
    reference_radii: np.ndarray | None = None
    for time, z, receipt in result:
        if not np.isfinite([time, z]).all():
            raise ValueError("cell time/z values must be finite")
        if (time, z) in pairs:
            raise ValueError("cell time/z pairs must be unique")
        pairs.add((time, z))
        if not isinstance(receipt, MissingCovarianceColumnTarget):
            raise TypeError("each cell receipt must be MissingCovarianceColumnTarget")
        if reference_radii is None:
            reference_radii = np.asarray(receipt.radii, dtype=float)
        elif not np.array_equal(reference_radii, receipt.radii):
            raise ValueError("all cell receipts must use the same radial grid")
        if receipt.metrics()["nodes_requiring_second_direction"] <= 0:
            raise ValueError("every cell must require a genuine second covariance direction")
    return result


def evaluate_spacetime_velocity_column_preflight(
    cells: Iterable[tuple[float, float, MissingCovarianceColumnTarget]],
    unit_velocity: PhaseVelocityEvaluator,
    *,
    amplitude: float,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
    coefficient_budget: float,
) -> dict[str, Any]:
    """Require one physical velocity column to pass every frozen spacetime cell."""
    checked = _validate_cells(cells)
    if not callable(unit_velocity):
        raise TypeError("unit_velocity must be callable")
    amplitude = float(amplitude)
    coefficient_budget = float(coefficient_budget)
    if not np.isfinite(amplitude) or amplitude <= 0.0:
        raise ValueError("amplitude must be finite and positive")
    if not np.isfinite(coefficient_budget) or coefficient_budget <= 0.0:
        raise ValueError("coefficient_budget must be finite and positive")
    if int(angular_count) < 1 or int(phase_count) < 1:
        raise ValueError("angular_count and phase_count must be positive")

    reports: list[dict[str, Any]] = []
    for time, z, receipt in checked:
        response = measure_multiplicative_covariance_response(
            unit_velocity,
            receipt.radii,
            amplitude=amplitude,
            time=time,
            z=z,
            angular_count=int(angular_count),
            phase_count=int(phase_count),
        )
        routed = evaluate_covariance_tangent_preflight(
            receipt,
            response,
            coefficient_budget=coefficient_budget,
        )
        bounded = routed["bounded_inverse"]
        reports.append(
            {
                "time": time,
                "z": z,
                "target_missing_relative_vector_rms": receipt.metrics()[
                    "missing_relative_vector_rms"
                ],
                "candidate_response_vector_rms": _vector_rms(
                    response[receipt.active_target_mask]
                ),
                "rank2_required_nodes": bounded["rank2_required_nodes"],
                "nodes_requiring_second_direction": bounded[
                    "nodes_requiring_second_direction"
                ],
                "two_column_relative_stress_residual_rms": bounded[
                    "two_column_relative_stress_residual_rms"
                ],
                "bounded_inverse_preflight_passed": bounded[
                    "bounded_inverse_preflight_passed"
                ],
                "finite_cycle_rerun_allowed_for_cell": routed[
                    "finite_cycle_rerun_allowed"
                ],
                "bounded_preflight": routed,
            }
        )

    all_pass = all(item["finite_cycle_rerun_allowed_for_cell"] for item in reports)
    required_total = sum(item["nodes_requiring_second_direction"] for item in reports)
    rank2_total = sum(item["rank2_required_nodes"] for item in reports)
    return {
        "cell_count": len(reports),
        "amplitude": amplitude,
        "coefficient_budget": coefficient_budget,
        "required_nodes_total": int(required_total),
        "rank2_required_nodes_total": int(rank2_total),
        "all_cells_bounded_inverse_preflight_passed": bool(all_pass),
        "finite_cycle_rerun_allowed": bool(all_pass),
        "cells": reports,
        "interpretation": (
            "A future physical second covariance column must pass the unchanged bounded "
            "inverse on every preregistered time/z cell before Agent 3 may materialize "
            "a two-column correction. Passing is necessary, not sufficient, for an "
            "actual Navier-Stokes residual contraction."
        ),
    }


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/spacetime_covariance_preflight_report.json",
    times: Iterable[float] = HELD_OUT_CYCLE_TIMES,
    z_values: Iterable[float] = SPATIAL_SCREEN_Z,
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    """Reconstruct real targets over the frozen time/z product and reject duplicates."""
    frozen_times = _validate_axis_values(times, name="times")
    frozen_z = _validate_axis_values(z_values, name="z_values")
    if REFERENCE_BUDGET_TIME not in frozen_times:
        raise ValueError("times must include the existing t=.5 budget reference")
    if REFERENCE_BUDGET_Z not in frozen_z:
        raise ValueError("z_values must include the existing z=.08 budget reference")

    cells: list[tuple[float, float, MissingCovarianceColumnTarget]] = []
    summaries: list[dict[str, Any]] = []
    reference_receipt: MissingCovarianceColumnTarget | None = None
    leading_sha: str | None = None

    for z in frozen_z:
        for time in frozen_times:
            receipt, metadata = build_real_covariance_slice(
                time=time,
                z=z,
                radial_count=radial_count,
                angular_count=angular_count,
                phase_count=phase_count,
            )
            cells.append((time, z, receipt))
            if time == REFERENCE_BUDGET_TIME and z == REFERENCE_BUDGET_Z:
                reference_receipt = receipt
            if leading_sha is None:
                leading_sha = str(metadata["leading_candidate_sha256"])
            elif leading_sha != metadata["leading_candidate_sha256"]:
                raise RuntimeError("leading candidate changed across spacetime reconstruction")
            summaries.append(
                {
                    "time": time,
                    "z": z,
                    "real_defect_projection": metadata["real_defect_projection"],
                    "target_metrics": metadata["metrics"],
                }
            )

    if reference_receipt is None:
        raise RuntimeError("reference budget receipt was not constructed")
    frozen_budget = current_signed_coefficient_budget(reference_receipt)
    for summary, (_, _, receipt) in zip(summaries, cells):
        summary["required_second_column_envelope"] = required_second_column_envelope(
            receipt,
            coefficient_budget=frozen_budget,
        )

    duplicate_unit = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=1.0,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    duplicate_control = evaluate_spacetime_velocity_column_preflight(
        cells,
        complete_curl_phase_evaluator(duplicate_unit),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        angular_count=angular_count,
        phase_count=phase_count,
        coefficient_budget=frozen_budget,
    )

    missing_relative = [
        float(item["target_metrics"]["missing_relative_vector_rms"])
        for item in summaries
    ]
    transverse_ratios = [
        float(
            item["required_second_column_envelope"][
                "required_transverse_over_current_max"
            ]
        )
        for item in summaries
    ]
    transverse_response_max = [
        float(
            item["required_second_column_envelope"][
                "required_transverse_response_max"
            ]
        )
        for item in summaries
    ]

    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "section": SOURCE_SECTION,
            "source_structure": (
                "The corrected source uses compact mean-stress reconstruction, an "
                "invertible two-column signed covariance map, and finite correction "
                "iteration. The extra axial-slice screen is autonomous repository "
                "engineering and is not a source theorem or threshold."
            ),
        },
        "inputs": {
            "leading_candidate_sha256": leading_sha,
            "oscillatory_amplitude": ROUTED_OSCILLATORY_AMPLITUDE,
            "oscillatory_phase": ROUTED_OSCILLATORY_PHASE,
            "times": list(frozen_times),
            "z_values": list(frozen_z),
            "profile_annulus": [float(PROFILE_ANNULUS[0]), float(PROFILE_ANNULUS[1])],
            "radial_count": int(radial_count),
            "angular_count": int(angular_count),
            "phase_count": int(phase_count),
            "coefficient_budget": frozen_budget,
            "coefficient_budget_reference": {
                "time": REFERENCE_BUDGET_TIME,
                "z": REFERENCE_BUDGET_Z,
                "origin": (
                    "unchanged max |delta a_1| from the pre-existing real rank-one "
                    "Agent-3 #302/#329 reference target"
                ),
            },
            "spatial_sampling_origin": (
                "autonomous core-safe robustness screen frozen before measurement; "
                "not a Kokuno/OpenAI hidden numerical parameter"
            ),
            "surrogate_defect_used": False,
        },
        "real_spacetime_target": {
            "cell_count": len(summaries),
            "cells": summaries,
            "missing_relative_vector_rms_min": float(min(missing_relative)),
            "missing_relative_vector_rms_max": float(max(missing_relative)),
            "required_transverse_over_current_max_across_cells": float(max(transverse_ratios)),
            "required_transverse_response_max_across_cells": float(max(transverse_response_max)),
        },
        "duplicate_existing_column_negative_control": duplicate_control,
        "routing": {
            "agent2_latest_dependency": (
                "Agent 2 #336 exposes homogeneous-pulse D_r/D_z sensitivities but still "
                "does not provide the actual independent public source velocity column"
            ),
            "second_public_covariance_column_available": False,
            "finite_correction_cycle_rerun_allowed": bool(
                duplicate_control["finite_cycle_rerun_allowed"]
            ),
            "next_required": (
                "Once Agent 2 materializes a genuinely independent source-motivated "
                "real-pair public velocity column, route exactly that column through "
                "this frozen spacetime guard; only then may Agent 3 construct a "
                "two-column signed correction and rerun held-in/held-out residuals."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_targets_consumed": True,
            "surrogate_defect_used": False,
            "held_out_cycle_times_screened": True,
            "multiple_axial_slices_screened": True,
            "spatial_screen_is_autonomous_engineering": True,
            "coefficient_budget_changed": False,
            "new_oscillatory_column_constructed": False,
            "agent2_complete_curl_reimplemented": False,
            "agent2_pulse_sensitivity_reimplemented": False,
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

    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Screen a future Kokuno covariance column across frozen time/z cells"
    )
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/spacetime_covariance_preflight_report.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output=args.output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
