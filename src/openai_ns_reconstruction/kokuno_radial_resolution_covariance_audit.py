"""Audit radial-resolution stability of the real Kokuno mean-stress target.

Agent-3 #338 screens a future second covariance column on a frozen 3x3
spacetime set, but each real compact-stress target currently uses 33 radial
nodes.  Before treating the measured transverse-response requirement as a
stable engineering handoff, this module repeats the same real
Agent-1-leading + Agent-2-oscillatory defect -> compact radial stress -> missing
covariance construction on one preregistered nested radial ladder.

The ladder is a repository numerical-audit choice, not a Kokuno/OpenAI source
constant and not a Navier--Stokes acceptance threshold.  The signed coefficient
budget is frozen from the pre-existing 33-node reference receipt and is never
refit on the finer grids.  No new oscillatory column is constructed and no
finite correction cycle is rerun here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    MissingCovarianceColumnTarget,
    PROFILE_Z,
)
from .kokuno_multislice_covariance_preflight import build_real_covariance_slice
from .kokuno_second_column_bounded_inverse import (
    current_signed_coefficient_budget,
    evaluate_bounded_second_column,
    required_second_column_envelope,
)
from .kokuno_signed_covariance_inverse import PROFILE_ANNULUS, PROFILE_TIME

TASK = "KOKUNO-A3-RADIAL-RESOLUTION-COVARIANCE-AUDIT-015"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_SECTION = "Signed covariance / compact stress / finite correction iteration"

# Autonomous numerical-audit choices frozen before the real report is run.
RADIAL_RESOLUTION_LADDER = (33, 65, 129)
REFERENCE_RADIAL_COUNT = 33
FINE_PAIR_RELATIVE_STABILITY_TOLERANCE = 2.0e-2


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("values must have shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _flat_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _relative_change(a: float, b: float) -> float:
    scale = max(abs(float(b)), np.finfo(float).tiny)
    return abs(float(a) - float(b)) / scale


def _relative_vector_difference(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 2:
        raise ValueError("vector arrays must have equal shape (N,2)")
    return _vector_rms(a - b) / max(_vector_rms(b), np.finfo(float).tiny)


def _validate_radial_counts(counts: Iterable[int]) -> tuple[int, ...]:
    values = tuple(int(value) for value in counts)
    if len(values) < 2:
        raise ValueError("radial_counts must contain at least two levels")
    if any(value < 9 for value in values):
        raise ValueError("every radial count must be at least 9")
    if any(right <= left for left, right in zip(values, values[1:])):
        raise ValueError("radial_counts must be strictly increasing")
    for coarse, fine in zip(values, values[1:]):
        if (fine - 1) % (coarse - 1) != 0:
            raise ValueError("radial_counts must form a nested endpoint-preserving ladder")
    return values


def _nested_indices(coarse_count: int, fine_count: int) -> np.ndarray:
    coarse_count = int(coarse_count)
    fine_count = int(fine_count)
    if coarse_count < 2 or fine_count <= coarse_count:
        raise ValueError("expected fine_count > coarse_count >= 2")
    if (fine_count - 1) % (coarse_count - 1) != 0:
        raise ValueError("radial grids are not nested")
    stride = (fine_count - 1) // (coarse_count - 1)
    return np.arange(coarse_count, dtype=int) * stride


def compare_nested_receipts(
    coarse: MissingCovarianceColumnTarget,
    fine: MissingCovarianceColumnTarget,
    *,
    coefficient_budget: float,
) -> dict[str, Any]:
    """Compare a finer real target on the exact nodes of a nested coarse grid."""
    if not isinstance(coarse, MissingCovarianceColumnTarget) or not isinstance(
        fine, MissingCovarianceColumnTarget
    ):
        raise TypeError("coarse and fine must be MissingCovarianceColumnTarget")
    if not np.isclose(coarse.radii[0], fine.radii[0], rtol=0.0, atol=1e-14) or not np.isclose(
        coarse.radii[-1], fine.radii[-1], rtol=0.0, atol=1e-14
    ):
        raise ValueError("coarse and fine receipts must share radial endpoints")

    indices = _nested_indices(len(coarse.radii), len(fine.radii))
    if not np.allclose(coarse.radii, fine.radii[indices], rtol=0.0, atol=1e-14):
        raise ValueError("receipt radial grids do not realize the declared nested ladder")

    target_difference = _relative_vector_difference(
        coarse.target_stress, fine.target_stress[indices]
    )
    current_difference = _relative_vector_difference(
        coarse.current_response, fine.current_response[indices]
    )
    missing_difference = _relative_vector_difference(
        coarse.missing_response, fine.missing_response[indices]
    )

    common_safe = coarse.safe_active_mask & fine.safe_active_mask[indices]
    if not np.any(common_safe):
        raise ValueError("nested receipts have no common safe active nodes")
    coefficient_difference = _flat_rms(
        coarse.current_coefficient[common_safe]
        - fine.current_coefficient[indices][common_safe]
    ) / max(
        _flat_rms(fine.current_coefficient[indices][common_safe]),
        np.finfo(float).tiny,
    )

    coarse_metrics = coarse.metrics()
    fine_metrics = fine.metrics()
    coarse_envelope = required_second_column_envelope(
        coarse, coefficient_budget=coefficient_budget
    )
    fine_envelope = required_second_column_envelope(
        fine, coefficient_budget=coefficient_budget
    )

    guarded_changes = {
        "target_stress_relative_vector_rms_difference": target_difference,
        "current_response_relative_vector_rms_difference": current_difference,
        "missing_response_relative_vector_rms_difference": missing_difference,
        "missing_relative_vector_rms_relative_change": _relative_change(
            coarse_metrics["missing_relative_vector_rms"],
            fine_metrics["missing_relative_vector_rms"],
        ),
        "required_transverse_response_max_relative_change": _relative_change(
            coarse_envelope["required_transverse_response_max"],
            fine_envelope["required_transverse_response_max"],
        ),
        "required_transverse_over_current_max_relative_change": _relative_change(
            coarse_envelope["required_transverse_over_current_max"],
            fine_envelope["required_transverse_over_current_max"],
        ),
    }
    return {
        "coarse_radial_count": len(coarse.radii),
        "fine_radial_count": len(fine.radii),
        "common_safe_active_nodes": int(np.count_nonzero(common_safe)),
        "current_coefficient_relative_rms_difference": float(coefficient_difference),
        **guarded_changes,
        "max_guarded_relative_change": float(max(guarded_changes.values())),
        "coarse_nodes_requiring_second_direction": coarse_metrics[
            "nodes_requiring_second_direction"
        ],
        "fine_nodes_requiring_second_direction": fine_metrics[
            "nodes_requiring_second_direction"
        ],
    }


def audit_radial_resolution_ladder(
    receipts: Iterable[MissingCovarianceColumnTarget],
    *,
    coefficient_budget: float,
    fine_pair_relative_stability_tolerance: float = FINE_PAIR_RELATIVE_STABILITY_TOLERANCE,
) -> dict[str, Any]:
    """Audit convergence while keeping one pre-existing signed-update budget frozen."""
    checked = tuple(receipts)
    counts = _validate_radial_counts(len(receipt.radii) for receipt in checked)
    if len(checked) != len(counts):
        raise AssertionError("unreachable radial-level mismatch")
    if not np.isfinite(coefficient_budget) or coefficient_budget <= 0.0:
        raise ValueError("coefficient_budget must be positive and finite")
    tolerance = float(fine_pair_relative_stability_tolerance)
    if not np.isfinite(tolerance) or not 0.0 < tolerance < 0.2:
        raise ValueError("fine-pair stability tolerance must lie in (0,0.2)")

    levels: list[dict[str, Any]] = []
    for receipt in checked:
        levels.append(
            {
                "radial_count": len(receipt.radii),
                "target_metrics": receipt.metrics(),
                "required_second_column_envelope": required_second_column_envelope(
                    receipt, coefficient_budget=coefficient_budget
                ),
            }
        )

    comparisons = [
        compare_nested_receipts(
            coarse,
            fine,
            coefficient_budget=coefficient_budget,
        )
        for coarse, fine in zip(checked, checked[1:])
    ]
    finest = comparisons[-1]
    stable = finest["max_guarded_relative_change"] <= tolerance
    return {
        "radial_counts": list(counts),
        "coefficient_budget": float(coefficient_budget),
        "coefficient_budget_refit_on_fine_levels": False,
        "fine_pair_relative_stability_tolerance": tolerance,
        "levels": levels,
        "comparisons": comparisons,
        "finest_pair_max_guarded_relative_change": finest[
            "max_guarded_relative_change"
        ],
        "radial_resolution_stability_preflight_passed": bool(stable),
        "interpretation": (
            "This is a numerical-stability guard on the real compact-stress/covariance "
            "handoff, not a PDE acceptance threshold. A stable target still requires "
            "a genuinely independent physical second column and an actual held-out "
            "finite correction cycle."
        ),
    }


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/radial_resolution_covariance_audit_report.json",
    radial_counts: Iterable[int] = RADIAL_RESOLUTION_LADDER,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    """Run the nested radial audit on the real reference defect/stress target."""
    counts = _validate_radial_counts(radial_counts)
    if counts[0] != REFERENCE_RADIAL_COUNT:
        raise ValueError("radial ladder must begin at the pre-existing 33-node reference")
    if int(angular_count) < 1 or int(phase_count) < 1:
        raise ValueError("angular_count and phase_count must be positive")

    receipts: list[MissingCovarianceColumnTarget] = []
    raw_levels: list[dict[str, Any]] = []
    leading_sha: str | None = None
    for count in counts:
        receipt, metadata = build_real_covariance_slice(
            time=PROFILE_TIME,
            z=PROFILE_Z,
            radial_count=count,
            angular_count=int(angular_count),
            phase_count=int(phase_count),
        )
        receipts.append(receipt)
        if leading_sha is None:
            leading_sha = str(metadata["leading_candidate_sha256"])
        elif leading_sha != metadata["leading_candidate_sha256"]:
            raise RuntimeError("leading candidate changed across radial-resolution audit")
        raw_levels.append(
            {
                "radial_count": count,
                "real_defect_projection": metadata["real_defect_projection"],
                "target_metrics": metadata["metrics"],
            }
        )

    frozen_budget = current_signed_coefficient_budget(receipts[0])
    audit = audit_radial_resolution_ladder(
        receipts,
        coefficient_budget=frozen_budget,
    )
    finest_duplicate = evaluate_bounded_second_column(
        receipts[-1],
        receipts[-1].current_response,
        coefficient_budget=frozen_budget,
    )

    for raw, audited in zip(raw_levels, audit["levels"]):
        raw["required_second_column_envelope"] = audited[
            "required_second_column_envelope"
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
                "The corrected source uses compact mean-stress reconstruction and an "
                "invertible signed covariance map. The radial-resolution ladder and "
                "2% fine-pair stability guard are autonomous repository numerical-audit "
                "choices, not source constants or PDE thresholds."
            ),
        },
        "inputs": {
            "leading_candidate_sha256": leading_sha,
            "profile_time": float(PROFILE_TIME),
            "profile_z": float(PROFILE_Z),
            "profile_annulus": [float(PROFILE_ANNULUS[0]), float(PROFILE_ANNULUS[1])],
            "radial_counts": list(counts),
            "angular_count": int(angular_count),
            "phase_count": int(phase_count),
            "coefficient_budget": frozen_budget,
            "coefficient_budget_origin": (
                "unchanged max |delta a_1| from the pre-existing 33-node real rank-one "
                "reference target; never refit on 65/129 nodes"
            ),
            "surrogate_defect_used": False,
        },
        "real_resolution_levels": raw_levels,
        "resolution_audit": audit,
        "finest_duplicate_existing_column_negative_control": finest_duplicate,
        "routing": {
            "second_public_covariance_column_available": False,
            "radial_target_stable_enough_for_future_column_screen": audit[
                "radial_resolution_stability_preflight_passed"
            ],
            "finite_correction_cycle_rerun_allowed": False,
            "next_required_if_stable": (
                "Consume Agent 2's first genuinely independent source-motivated public "
                "velocity column, screen it with the frozen #338 spacetime guard and "
                "this resolution-audited target before materializing a two-column correction."
            ),
            "next_required_if_unstable": (
                "Refine the compact radial-stress/covariance target numerically before "
                "using its transverse-response envelope as a second-column gate."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_targets_consumed": True,
            "surrogate_defect_used": False,
            "coefficient_budget_changed": False,
            "radial_resolution_audit_completed": True,
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
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/radial_resolution_covariance_audit_report.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output=args.output)
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
