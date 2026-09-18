"""Expose the real covariance direction still missing from Agent 3's rank-one mean correction.

The corrected Kokuno reader realizes a two-component compact stress through two
*genuinely distinct* signed covariance columns.  PR #275 established that phase
translations of Agent 2's current complete-curl column collapse to one quadratic
covariance direction, so rerunning the rejected one-column finite correction is
not justified.

This module turns that negative rank result into an executable handoff contract.
It consumes the same real Agent-1 + Agent-2 phase-mean Navier--Stokes defect,
reconstructs the compact theta/e=2 and axial/e=1 radial stresses, measures the
current public Agent-2 covariance response ``(<w_r w_theta>, <w_r w_z>)``, and
removes the best pointwise rank-one projection.  The orthogonal remainder is the
minimum local two-component response that is *not* representable by the current
column.

The remainder is a target for a future independent source-motivated Agent-2
column; it is not itself a velocity, pulse, polarization, or proof that such a
column exists.  No residual-selected parameter, pressure, force, PDE threshold,
or finite correction cycle is introduced here.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_mean_defect import PhaseMeanDefectContract
from .kokuno_phase_orbit_covariance_rank import measure_phase_mean_covariance_vector
from .kokuno_radial_stress import CompactRadialStressInverse, sample_real_phase_mean_radial_profiles
from .kokuno_signed_covariance_inverse import (
    PROFILE_ANNULUS,
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)


TASK = "KOKUNO-A3-MISSING-COVARIANCE-COLUMN-TARGET-008"
SOURCE_READER = "KokunoYumeto corrected 208-page reconstruction"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_SECTION = "Signed covariance and exact signed differential inverse"
PROFILE_Z = 0.08
DEFAULT_RADIAL_COUNT = 33
DEFAULT_ANGULAR_COUNT = 8
DEFAULT_PHASE_COUNT = 8
TARGET_FLOOR_FRACTION = 1.0e-3
RESPONSE_FLOOR_FRACTION = 1.0e-3
RANK_RELATIVE_TOLERANCE = 1.0e-10


def _flat_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("values must have shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


@dataclass(frozen=True)
class MissingCovarianceColumnTarget:
    """Pointwise remainder after the best current one-column covariance fit."""

    radii: np.ndarray
    target_stress: np.ndarray
    current_response: np.ndarray
    current_coefficient: np.ndarray
    current_projection: np.ndarray
    missing_response: np.ndarray
    active_target_mask: np.ndarray
    safe_response_mask: np.ndarray
    target_floor_fraction: float = TARGET_FLOOR_FRACTION
    response_floor_fraction: float = RESPONSE_FLOOR_FRACTION

    def __post_init__(self) -> None:
        radii = np.asarray(self.radii, dtype=float)
        target = np.asarray(self.target_stress, dtype=float)
        current = np.asarray(self.current_response, dtype=float)
        coefficient = np.asarray(self.current_coefficient, dtype=float)
        projection = np.asarray(self.current_projection, dtype=float)
        missing = np.asarray(self.missing_response, dtype=float)
        active = np.asarray(self.active_target_mask, dtype=bool)
        safe = np.asarray(self.safe_response_mask, dtype=bool)
        if radii.ndim != 1 or len(radii) < 9 or np.any(np.diff(radii) <= 0.0):
            raise ValueError("radii must be a strictly increasing one-dimensional grid")
        if target.shape != (len(radii), 2) or current.shape != target.shape:
            raise ValueError("target_stress and current_response must have shape (N,2)")
        if projection.shape != target.shape or missing.shape != target.shape:
            raise ValueError("projection and missing_response must have shape (N,2)")
        if coefficient.shape != radii.shape or active.shape != radii.shape or safe.shape != radii.shape:
            raise ValueError("coefficient and masks must match radii")
        numeric = (radii, target, current, coefficient, projection, missing)
        if not all(np.isfinite(value).all() for value in numeric):
            raise ValueError("all numeric contract arrays must be finite")
        if not 0.0 < self.target_floor_fraction < 0.1:
            raise ValueError("target_floor_fraction must lie in (0,0.1)")
        if not 0.0 < self.response_floor_fraction < 0.1:
            raise ValueError("response_floor_fraction must lie in (0,0.1)")
        for name, value in (
            ("radii", radii),
            ("target_stress", target),
            ("current_response", current),
            ("current_coefficient", coefficient),
            ("current_projection", projection),
            ("missing_response", missing),
            ("active_target_mask", active),
            ("safe_response_mask", safe),
        ):
            copied = value.copy()
            copied.setflags(write=False)
            object.__setattr__(self, name, copied)

    @property
    def safe_active_mask(self) -> np.ndarray:
        return self.active_target_mask & self.safe_response_mask

    @property
    def unresolved_active_mask(self) -> np.ndarray:
        return self.active_target_mask & ~self.safe_response_mask

    def metrics(self) -> dict[str, Any]:
        active = self.active_target_mask
        safe_active = self.safe_active_mask
        target_scale = max(_vector_rms(self.target_stress[active]), np.finfo(float).tiny)
        missing_scale = _vector_rms(self.missing_response[active])
        dot = np.sum(self.current_response * self.missing_response, axis=1)
        denom = np.linalg.norm(self.current_response, axis=1) * np.linalg.norm(
            self.missing_response, axis=1
        )
        orthogonal = np.zeros_like(dot)
        valid = safe_active & (denom > np.finfo(float).tiny)
        orthogonal[valid] = np.abs(dot[valid]) / denom[valid]
        theta_missing = self.missing_response[active, 0]
        axial_missing = self.missing_response[active, 1]
        missing_norm = np.linalg.norm(self.missing_response, axis=1)
        target_max = max(float(np.max(np.linalg.norm(self.target_stress, axis=1))), np.finfo(float).tiny)
        second_needed = active & (missing_norm > self.target_floor_fraction * target_max)
        return {
            "active_target_nodes": int(np.count_nonzero(active)),
            "safe_active_nodes": int(np.count_nonzero(safe_active)),
            "unresolved_active_nodes": int(np.count_nonzero(self.unresolved_active_mask)),
            "nodes_requiring_second_direction": int(np.count_nonzero(second_needed)),
            "target_vector_rms": target_scale,
            "current_projection_vector_rms": _vector_rms(self.current_projection[active]),
            "missing_vector_rms": missing_scale,
            "missing_relative_vector_rms": missing_scale / target_scale,
            "missing_theta_rms": _flat_rms(theta_missing),
            "missing_axial_rms": _flat_rms(axial_missing),
            "missing_max_vector_norm": float(np.max(missing_norm[active])),
            "max_pointwise_current_missing_cosine_abs": float(np.max(orthogonal[valid])) if np.any(valid) else 0.0,
            "current_coefficient_rms": _flat_rms(self.current_coefficient[safe_active]) if np.any(safe_active) else 0.0,
            "current_coefficient_max_abs": float(np.max(np.abs(self.current_coefficient[safe_active]))) if np.any(safe_active) else 0.0,
        }


def build_missing_covariance_column_target(
    radii: np.ndarray,
    target_stress: np.ndarray,
    current_response: np.ndarray,
    *,
    target_floor_fraction: float = TARGET_FLOOR_FRACTION,
    response_floor_fraction: float = RESPONSE_FLOOR_FRACTION,
) -> MissingCovarianceColumnTarget:
    """Remove the best pointwise rank-one response from a real two-component target."""
    radii = np.asarray(radii, dtype=float)
    target = np.asarray(target_stress, dtype=float)
    current = np.asarray(current_response, dtype=float)
    if radii.ndim != 1 or target.shape != (len(radii), 2) or current.shape != target.shape:
        raise ValueError("expected radii (N,) and target/current arrays (N,2)")
    if len(radii) < 9 or not np.isfinite(radii).all() or np.any(np.diff(radii) <= 0.0):
        raise ValueError("radii must be finite, increasing, and contain at least 9 nodes")
    if not np.isfinite(target).all() or not np.isfinite(current).all():
        raise ValueError("target and current responses must be finite")
    if not 0.0 < target_floor_fraction < 0.1 or not 0.0 < response_floor_fraction < 0.1:
        raise ValueError("floor fractions must lie in (0,0.1)")

    target_norm = np.linalg.norm(target, axis=1)
    response_norm = np.linalg.norm(current, axis=1)
    target_max = float(np.max(target_norm))
    response_max = float(np.max(response_norm))
    if target_max <= np.finfo(float).tiny:
        raise ValueError("target stress is numerically zero")
    if response_max <= np.finfo(float).tiny:
        raise ValueError("current covariance response is numerically zero")

    active = target_norm > target_floor_fraction * target_max
    safe = response_norm > response_floor_fraction * response_max
    coefficient = np.zeros(len(radii), dtype=float)
    safe_denominator = np.sum(current[safe] * current[safe], axis=1)
    coefficient[safe] = (
        np.sum(current[safe] * target[safe], axis=1) / safe_denominator
    )
    projection = coefficient[:, None] * current
    missing = target - projection
    # A response too small to invert is deliberately left wholly unresolved.
    projection[~safe] = 0.0
    missing[~safe] = target[~safe]
    coefficient[~safe] = 0.0
    return MissingCovarianceColumnTarget(
        radii=radii,
        target_stress=target,
        current_response=current,
        current_coefficient=coefficient,
        current_projection=projection,
        missing_response=missing,
        active_target_mask=active,
        safe_response_mask=safe,
        target_floor_fraction=float(target_floor_fraction),
        response_floor_fraction=float(response_floor_fraction),
    )


def evaluate_second_column_candidate(
    receipt: MissingCovarianceColumnTarget,
    second_response: np.ndarray,
    *,
    relative_rank_tolerance: float = RANK_RELATIVE_TOLERANCE,
) -> dict[str, Any]:
    """Evaluate a future public second covariance column without selecting one here."""
    if not isinstance(receipt, MissingCovarianceColumnTarget):
        raise TypeError("receipt must be MissingCovarianceColumnTarget")
    second = np.asarray(second_response, dtype=float)
    if second.shape != receipt.target_stress.shape or not np.isfinite(second).all():
        raise ValueError("second_response must be finite and match target shape")
    if not 0.0 < relative_rank_tolerance < 1.0e-2:
        raise ValueError("relative_rank_tolerance must lie in (0,1e-2)")

    active = receipt.active_target_mask
    residual1 = receipt.missing_response[active]
    target = receipt.target_stress[active]
    residual2 = np.zeros_like(target)
    ranks = []
    conditions = []
    for out_index, source_index in enumerate(np.flatnonzero(active)):
        matrix = np.column_stack(
            (receipt.current_response[source_index], second[source_index])
        )
        singular = np.linalg.svd(matrix, compute_uv=False)
        if singular[0] <= np.finfo(float).tiny:
            rank = 0
            condition = float("inf")
            residual2[out_index] = target[out_index]
        else:
            rank = int(np.count_nonzero(singular > relative_rank_tolerance * singular[0]))
            condition = float(singular[0] / singular[-1]) if rank == 2 else float("inf")
            coefficients, _, _, _ = np.linalg.lstsq(matrix, target[out_index], rcond=None)
            residual2[out_index] = matrix @ coefficients - target[out_index]
        ranks.append(rank)
        conditions.append(condition)

    rank_array = np.asarray(ranks, dtype=int)
    finite_condition = np.asarray([v for v in conditions if np.isfinite(v)], dtype=float)
    target_scale = max(_vector_rms(target), np.finfo(float).tiny)
    rank1_relative = _vector_rms(residual1) / target_scale
    rank2_relative = _vector_rms(residual2) / target_scale
    return {
        "active_target_nodes": int(len(rank_array)),
        "rank2_active_nodes": int(np.count_nonzero(rank_array == 2)),
        "rank_deficient_active_nodes": int(np.count_nonzero(rank_array < 2)),
        "rank1_relative_residual_rms": rank1_relative,
        "two_column_relative_residual_rms": rank2_relative,
        "incremental_relative_residual_reduction": (
            1.0 - rank2_relative / rank1_relative if rank1_relative > np.finfo(float).tiny else 0.0
        ),
        "finite_condition_max": float(np.max(finite_condition)) if len(finite_condition) else None,
        "finite_condition_median": float(np.median(finite_condition)) if len(finite_condition) else None,
        "all_active_nodes_rank2": bool(np.all(rank_array == 2)),
    }


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/missing_covariance_column_target_report.json",
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    """Build the real compact-stress remainder that a future second column must span."""
    leading = KokunoLeadingCoreSeriesCandidate()
    correction = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    contract = PhaseMeanDefectContract(phase_count=phase_count)
    theta_profile, axial_profile, projection = sample_real_phase_mean_radial_profiles(
        leading.at_points,
        correction,
        contract,
        time=PROFILE_TIME,
        z=PROFILE_Z,
        r_inner=PROFILE_ANNULUS[0],
        r_outer=PROFILE_ANNULUS[1],
        radial_count=radial_count,
        angular_count=angular_count,
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
        time=PROFILE_TIME,
        z=PROFILE_Z,
        angular_count=angular_count,
        phase_count=phase_count,
    )
    current_response = 2.0 * float(correction.amplitude) * unit_covariance
    receipt = build_missing_covariance_column_target(radii, target, current_response)
    metrics = receipt.metrics()

    duplicate_audit = evaluate_second_column_candidate(receipt, current_response)
    ideal_handoff_audit = evaluate_second_column_candidate(receipt, receipt.missing_response)

    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "reader": SOURCE_READER,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "section": SOURCE_SECTION,
            "source_requirement": (
                "two distinct signed covariance columns are inverted through H; "
                "a phase translation of one column is not a replacement for the second column"
            ),
        },
        "inputs": {
            "leading_candidate_sha256": leading.sha256,
            "oscillatory_amplitude": float(correction.amplitude),
            "oscillatory_phase": float(correction.phase),
            "profile_time": float(PROFILE_TIME),
            "profile_z": float(PROFILE_Z),
            "profile_annulus": [float(PROFILE_ANNULUS[0]), float(PROFILE_ANNULUS[1])],
            "radial_count": int(radial_count),
            "angular_count": int(angular_count),
            "phase_count": int(phase_count),
            "target_floor_fraction": TARGET_FLOOR_FRACTION,
            "response_floor_fraction": RESPONSE_FLOOR_FRACTION,
            "surrogate_defect_used": False,
        },
        "real_defect_projection": projection,
        "missing_second_column_target": {
            **metrics,
            "radii": [float(v) for v in receipt.radii],
            "target_stress_theta_axial": receipt.target_stress.tolist(),
            "current_response_theta_axial": receipt.current_response.tolist(),
            "current_rank1_coefficient": receipt.current_coefficient.tolist(),
            "missing_response_theta_axial": receipt.missing_response.tolist(),
            "active_target_mask": receipt.active_target_mask.tolist(),
            "safe_response_mask": receipt.safe_response_mask.tolist(),
            "interpretation": (
                "pointwise minimum orthogonal remainder after the best existing one-column fit; "
                "this is a measured covariance-response target, not a velocity or source pulse"
            ),
        },
        "calibration": {
            "duplicate_existing_column": duplicate_audit,
            "ideal_missing_response_as_abstract_second_column": ideal_handoff_audit,
            "ideal_handoff_used_as_candidate_data": False,
        },
        "routing": {
            "finite_correction_cycle_rerun_allowed": False,
            "second_public_covariance_column_available": False,
            "next_required": (
                "Agent 2 should materialize a genuinely independent source-motivated public velocity column; "
                "Agent 3 can then measure its two-component covariance response with this contract and only "
                "rerun the finite correction cycle if the measured two-column map is rank two and reduces "
                "the real target projection without changing registered PDE thresholds"
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_consumed": True,
            "theta_and_axial_radial_stress_reconstructed": True,
            "missing_second_column_target_measured": True,
            "missing_target_is_public_velocity": False,
            "new_oscillatory_column_constructed": False,
            "agent2_complete_curl_reimplemented": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
            "residual_reduction_claimed": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/missing_covariance_column_target_report.json",
    )
    parser.add_argument("--radial-count", type=int, default=DEFAULT_RADIAL_COUNT)
    parser.add_argument("--angular-count", type=int, default=DEFAULT_ANGULAR_COUNT)
    parser.add_argument("--phase-count", type=int, default=DEFAULT_PHASE_COUNT)
    return parser


def main() -> None:
    args = _parser().parse_args()
    report = generate_actual_core_report(
        output=args.output,
        radial_count=args.radial_count,
        angular_count=args.angular_count,
        phase_count=args.phase_count,
    )
    summary = report["missing_second_column_target"]
    print(
        json.dumps(
            {
                "task": report["task"],
                "missing_relative_vector_rms": summary["missing_relative_vector_rms"],
                "missing_theta_rms": summary["missing_theta_rms"],
                "missing_axial_rms": summary["missing_axial_rms"],
                "nodes_requiring_second_direction": summary["nodes_requiring_second_direction"],
                "duplicate_rank2_nodes": report["calibration"]["duplicate_existing_column"]["rank2_active_nodes"],
                "ideal_handoff_relative_residual_rms": report["calibration"]["ideal_missing_response_as_abstract_second_column"]["two_column_relative_residual_rms"],
            },
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
