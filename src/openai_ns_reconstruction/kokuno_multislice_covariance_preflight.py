"""Preflight a future second Kokuno covariance column across held-out cycle times.

The existing Agent-3 bounded inverse (#302) and public velocity-column adapter
(#322) screen one real defect/stress slice at ``t=.5, z=.08``.  A finite
correction cycle, however, is judged on frozen held-out times as well.  This
module extends the *same* real-defect compact-stress contract to
``t=(.375,.5,.625)`` while keeping one coefficient budget frozen from the
existing ``t=.5`` receipt.

For every slice we reconstruct the theta/e=2 and axial/e=1 compact radial
stresses from the real Agent-1 + Agent-2 phase-mean Navier--Stokes defect,
measure the current covariance response, and retain the pointwise rank-one
remainder.  A caller-supplied public phase-dependent velocity column is then
measured through the already-tested multiplicative covariance adapter and must
pass the existing rank/conditioning/signed-budget guard on *every* slice before
another finite correction cycle is authorized.

This is a routing/preflight layer only.  It does not construct Agent-2's second
oscillatory pair, fit pressure/forcing, change any PDE threshold, or claim a
residual reduction.
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
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_mean_defect import PhaseMeanDefectContract
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    DEFAULT_RADIAL_COUNT,
    PROFILE_Z,
    MissingCovarianceColumnTarget,
    build_missing_covariance_column_target,
)
from .kokuno_phase_orbit_covariance_rank import measure_phase_mean_covariance_vector
from .kokuno_radial_stress import (
    CompactRadialStressInverse,
    sample_real_phase_mean_radial_profiles,
)
from .kokuno_second_column_bounded_inverse import (
    current_signed_coefficient_budget,
    required_second_column_envelope,
)
from .kokuno_signed_covariance_inverse import (
    PROFILE_ANNULUS,
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)
from .kokuno_velocity_column_preflight import (
    PhaseVelocityEvaluator,
    measure_multiplicative_covariance_response,
)

TASK = "KOKUNO-A3-MULTISLICE-COVARIANCE-PREFLIGHT-013"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_SECTION = "Signed covariance / compact stress / finite correction iteration"

# These are the frozen held-out times already used by the rejected one-column
# Agent-3 correction-cycle replay.  They are not new acceptance thresholds.
HELD_OUT_CYCLE_TIMES = (0.375, 0.5, 0.625)
REFERENCE_BUDGET_TIME = PROFILE_TIME


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("values must have shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _validate_times(times: Iterable[float]) -> tuple[float, ...]:
    values = tuple(float(value) for value in times)
    if len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("times must contain at least two finite values")
    if any(right <= left for left, right in zip(values, values[1:])):
        raise ValueError("times must be strictly increasing")
    return values


def build_real_covariance_slice(
    *,
    time: float,
    z: float = PROFILE_Z,
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> tuple[MissingCovarianceColumnTarget, dict[str, Any]]:
    """Build one real defect/stress/covariance receipt at a fixed ``(time,z)``."""
    if not np.isfinite([time, z]).all():
        raise ValueError("time and z must be finite")
    if radial_count < 9:
        raise ValueError("radial_count must be at least 9")

    leading = KokunoLeadingCoreSeriesCandidate()
    correction = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    contract = PhaseMeanDefectContract(phase_count=int(phase_count))
    theta_profile, axial_profile, projection = sample_real_phase_mean_radial_profiles(
        leading.at_points,
        correction,
        contract,
        time=float(time),
        z=float(z),
        r_inner=PROFILE_ANNULUS[0],
        r_outer=PROFILE_ANNULUS[1],
        radial_count=int(radial_count),
        angular_count=int(angular_count),
    )
    theta_inverse = CompactRadialStressInverse(theta_profile)
    axial_inverse = CompactRadialStressInverse(axial_profile)
    radii = np.asarray(theta_profile.radii, dtype=float)
    target = np.stack(
        (theta_inverse.stress(radii), axial_inverse.stress(radii)), axis=-1
    )
    unit_covariance = measure_phase_mean_covariance_vector(
        correction,
        radii,
        phase_offset=0.0,
        time=float(time),
        z=float(z),
        angular_count=int(angular_count),
        phase_count=int(phase_count),
    )
    current_response = 2.0 * float(correction.amplitude) * unit_covariance
    receipt = build_missing_covariance_column_target(radii, target, current_response)
    return receipt, {
        "time": float(time),
        "z": float(z),
        "leading_candidate_sha256": leading.sha256,
        "oscillatory_amplitude": float(correction.amplitude),
        "oscillatory_phase": float(correction.phase),
        "real_defect_projection": projection,
        "metrics": receipt.metrics(),
    }


def evaluate_multislice_velocity_column_preflight(
    slices: Iterable[tuple[float, MissingCovarianceColumnTarget]],
    unit_velocity: PhaseVelocityEvaluator,
    *,
    amplitude: float,
    z: float = PROFILE_Z,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
    coefficient_budget: float,
) -> dict[str, Any]:
    """Require one public velocity column to pass the frozen guard on every slice."""
    pairs = tuple((float(time), receipt) for time, receipt in slices)
    times = _validate_times(time for time, _ in pairs)
    if len(pairs) != len(times):
        raise AssertionError("unreachable slice/time mismatch")
    if not callable(unit_velocity):
        raise TypeError("unit_velocity must be callable")
    amplitude = float(amplitude)
    coefficient_budget = float(coefficient_budget)
    if not np.isfinite(amplitude) or amplitude <= 0.0:
        raise ValueError("amplitude must be finite and positive")
    if not np.isfinite(coefficient_budget) or coefficient_budget <= 0.0:
        raise ValueError("coefficient_budget must be finite and positive")
    if not np.isfinite(z):
        raise ValueError("z must be finite")

    reference_radii: np.ndarray | None = None
    reports: list[dict[str, Any]] = []
    for time, receipt in pairs:
        if not isinstance(receipt, MissingCovarianceColumnTarget):
            raise TypeError("each slice receipt must be MissingCovarianceColumnTarget")
        if reference_radii is None:
            reference_radii = np.asarray(receipt.radii, dtype=float)
        elif not np.array_equal(reference_radii, receipt.radii):
            raise ValueError("all slice receipts must use the same radial grid")
        if receipt.metrics()["nodes_requiring_second_direction"] <= 0:
            raise ValueError("each slice must contain at least one node requiring a second direction")

        response = measure_multiplicative_covariance_response(
            unit_velocity,
            receipt.radii,
            amplitude=amplitude,
            time=time,
            z=float(z),
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
                "finite_cycle_rerun_allowed_for_slice": routed[
                    "finite_cycle_rerun_allowed"
                ],
                "bounded_preflight": routed,
            }
        )

    all_pass = all(item["finite_cycle_rerun_allowed_for_slice"] for item in reports)
    required_total = sum(item["nodes_requiring_second_direction"] for item in reports)
    rank2_total = sum(item["rank2_required_nodes"] for item in reports)
    return {
        "times": list(times),
        "z": float(z),
        "amplitude": amplitude,
        "coefficient_budget": coefficient_budget,
        "slice_count": len(reports),
        "required_nodes_total": int(required_total),
        "rank2_required_nodes_total": int(rank2_total),
        "all_slices_bounded_inverse_preflight_passed": bool(all_pass),
        "finite_cycle_rerun_allowed": bool(all_pass),
        "slices": reports,
        "interpretation": (
            "A future physical second covariance column must pass the unchanged bounded "
            "inverse on every frozen held-out cycle time before Agent 3 reruns a finite "
            "velocity correction cycle. Passing remains necessary, not sufficient, for "
            "Navier-Stokes residual contraction."
        ),
    }


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/multislice_covariance_preflight_report.json",
    times: Iterable[float] = HELD_OUT_CYCLE_TIMES,
    z: float = PROFILE_Z,
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    """Build real held-out slice targets and run the duplicate-column negative control."""
    frozen_times = _validate_times(times)
    if REFERENCE_BUDGET_TIME not in frozen_times:
        raise ValueError("times must include the existing t=.5 budget reference slice")

    receipts: list[tuple[float, MissingCovarianceColumnTarget]] = []
    slice_summaries: list[dict[str, Any]] = []
    reference_receipt: MissingCovarianceColumnTarget | None = None
    leading_sha: str | None = None

    for time in frozen_times:
        receipt, metadata = build_real_covariance_slice(
            time=time,
            z=z,
            radial_count=radial_count,
            angular_count=angular_count,
            phase_count=phase_count,
        )
        receipts.append((time, receipt))
        if time == REFERENCE_BUDGET_TIME:
            reference_receipt = receipt
        if leading_sha is None:
            leading_sha = str(metadata["leading_candidate_sha256"])
        elif leading_sha != metadata["leading_candidate_sha256"]:
            raise RuntimeError("leading candidate changed across slice construction")

        metrics = metadata["metrics"]
        envelope = required_second_column_envelope(
            receipt,
            coefficient_budget=(
                current_signed_coefficient_budget(receipt)
                if time == REFERENCE_BUDGET_TIME
                else 1.0
            ),
        )
        slice_summaries.append(
            {
                "time": time,
                "z": float(z),
                "real_defect_projection": metadata["real_defect_projection"],
                "target_metrics": metrics,
                # Replaced below for every slice after the reference budget is frozen.
                "required_second_column_envelope": envelope,
            }
        )

    if reference_receipt is None:
        raise RuntimeError("reference receipt was not constructed")
    frozen_budget = current_signed_coefficient_budget(reference_receipt)

    # Recompute every envelope with the one pre-existing t=.5 budget.  The
    # provisional value above avoids a second expensive defect reconstruction.
    for summary, (_, receipt) in zip(slice_summaries, receipts):
        summary["required_second_column_envelope"] = required_second_column_envelope(
            receipt,
            coefficient_budget=frozen_budget,
        )

    duplicate_unit = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=1.0,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    duplicate_control = evaluate_multislice_velocity_column_preflight(
        receipts,
        complete_curl_phase_evaluator(duplicate_unit),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        z=z,
        angular_count=angular_count,
        phase_count=phase_count,
        coefficient_budget=frozen_budget,
    )

    missing_relative = [
        float(summary["target_metrics"]["missing_relative_vector_rms"])
        for summary in slice_summaries
    ]
    required_transverse_over_current_max = [
        float(
            summary["required_second_column_envelope"][
                "required_transverse_over_current_max"
            ]
        )
        for summary in slice_summaries
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
                "iteration. The multislice held-out guard is repository engineering "
                "that prevents a one-slice covariance fit from authorizing a cycle."
            ),
        },
        "inputs": {
            "leading_candidate_sha256": leading_sha,
            "oscillatory_amplitude": ROUTED_OSCILLATORY_AMPLITUDE,
            "oscillatory_phase": ROUTED_OSCILLATORY_PHASE,
            "times": list(frozen_times),
            "profile_z": float(z),
            "profile_annulus": [float(PROFILE_ANNULUS[0]), float(PROFILE_ANNULUS[1])],
            "radial_count": int(radial_count),
            "angular_count": int(angular_count),
            "phase_count": int(phase_count),
            "coefficient_budget": frozen_budget,
            "coefficient_budget_source": (
                "unchanged max |delta a_1| measured from the existing t=.5 rank-one "
                "real-defect target used by Agent-3 #302"
            ),
            "surrogate_defect_used": False,
        },
        "real_multislice_target": {
            "slice_count": len(slice_summaries),
            "slices": slice_summaries,
            "missing_relative_vector_rms_min": float(min(missing_relative)),
            "missing_relative_vector_rms_max": float(max(missing_relative)),
            "required_transverse_over_current_max_across_slices": float(
                max(required_transverse_over_current_max)
            ),
        },
        "duplicate_existing_column_negative_control": duplicate_control,
        "routing": {
            "second_public_covariance_column_available": False,
            "finite_correction_cycle_rerun_allowed": bool(
                duplicate_control["finite_cycle_rerun_allowed"]
            ),
            "next_required": (
                "Agent 2 must provide a genuinely independent source-motivated public "
                "real-pair velocity column. Agent 3 should run this same multislice "
                "preflight first and only materialize a two-column correction if every "
                "held-out time passes the frozen rank/conditioning/budget guard."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_targets_consumed": True,
            "surrogate_defect_used": False,
            "held_out_cycle_times_screened": True,
            "coefficient_budget_changed": False,
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

    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Screen a future Kokuno covariance velocity column across held-out cycle times"
    )
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/multislice_covariance_preflight_report.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output=args.output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
