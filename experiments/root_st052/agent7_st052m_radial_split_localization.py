"""Localize the frozen ST052-M 24/24 radial path split by seed metadata.

Preregistered in issue #691 before execution and stacked directly on exact
Agent-7 PR #683.  This diagnostic changes no velocity formula or coefficient.
It reuses the exact 48-path DOP853 protocol and asks whether the inward/outward
label is explained by existing z-sign, band, or radius metadata before any new
basis is proposed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_public_source_inward_spiral_persistence as parent

TASK_ID = "CR003-ST052M-RADIAL-SPLIT-LOCALIZATION-097"
PREREG_ISSUE = 691
SOURCE_PARENT_PR = 683
SOURCE_PARENT_HEAD = "c8046f42814a45ecaa9dbe44d5a0521ee0289948"
SOURCE_TEMPORAL_PR = 587
SOURCE_WITNESS_PR = 652
PUBLIC_SOURCE_URL = "https://openai.com/index/navier-stokes-solution/"

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "witness_retuned": False,
    "parameter_grid_scan_performed": False,
    "optimization_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "visual_acceptance_threshold_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _radial_velocity(velocity_fn, points: np.ndarray, time: float) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    vel = np.asarray(velocity_fn(pts, float(time)), dtype=float)
    radius = np.hypot(pts[:, 0], pts[:, 1])
    if vel.shape != pts.shape or np.any(radius <= 0.0) or not np.isfinite(vel).all():
        raise ValueError("invalid off-axis radial-velocity sample")
    return (pts[:, 0] * vel[:, 0] + pts[:, 1] * vel[:, 1]) / radius


def _group_summary(
    delta_r: np.ndarray,
    start_ur: np.ndarray,
    metadata: list[dict[str, Any]],
    key: str,
) -> dict[str, dict[str, Any]]:
    values = []
    for meta in metadata:
        value = meta[key]
        if isinstance(value, float):
            value = f"{value:.12g}"
        else:
            value = str(value)
        values.append(value)
    result: dict[str, dict[str, Any]] = {}
    for value in sorted(set(values)):
        idx = np.asarray([i for i, v in enumerate(values) if v == value], dtype=int)
        dr = np.asarray(delta_r[idx], dtype=float)
        ur = np.asarray(start_ur[idx], dtype=float)
        inward = dr < 0.0
        count = int(np.sum(inward))
        result[value] = {
            "path_count": int(len(idx)),
            "inward_path_count": count,
            "outward_or_zero_path_count": int(len(idx) - count),
            "mean_delta_r": float(np.mean(dr)),
            "min_delta_r": float(np.min(dr)),
            "max_delta_r": float(np.max(dr)),
            "mean_start_radial_velocity": float(np.mean(ur)),
            "min_start_radial_velocity": float(np.min(ur)),
            "max_start_radial_velocity": float(np.max(ur)),
            "uniform_inward_label": bool(count == 0 or count == len(idx)),
            "all_inward": bool(count == len(idx)),
            "all_outward_or_zero": bool(count == 0),
        }
    return result


def _interval_record(
    velocity_fn,
    times: np.ndarray,
    positions: np.ndarray,
    metadata: list[dict[str, Any]],
    i0: int,
    i1: int,
) -> dict[str, Any]:
    radii = np.hypot(positions[:, :, 0], positions[:, :, 1])
    delta_r = np.asarray(radii[i1] - radii[i0], dtype=float)
    start_ur = _radial_velocity(velocity_fn, positions[i0], float(times[i0]))
    return {
        "time_start": float(times[i0]),
        "time_end": float(times[i1]),
        "mean_delta_r": float(np.mean(delta_r)),
        "inward_path_count": int(np.sum(delta_r < 0.0)),
        "groups": {
            "sign": _group_summary(delta_r, start_ur, metadata, "sign"),
            "band": _group_summary(delta_r, start_ur, metadata, "band"),
            "radius": _group_summary(delta_r, start_ur, metadata, "radius"),
            "angle_index": _group_summary(delta_r, start_ur, metadata, "angle_index"),
        },
    }


def field_records(velocity_fn) -> dict[str, Any]:
    times, positions, metadata = parent._integrate_positions(velocity_fn)
    segments = {}
    for i0, i1 in parent.SEGMENTS:
        key = f"{times[i0]:.3f}-{times[i1]:.3f}"
        segments[key] = _interval_record(velocity_fn, times, positions, metadata, i0, i1)
    whole = _interval_record(
        velocity_fn,
        times,
        positions,
        metadata,
        parent.CHECKPOINT_INDICES[0],
        parent.CHECKPOINT_INDICES[-1],
    )
    return {"segments": segments, "whole_interval": whole}


def _stable_uniform_pattern(
    records_by_field: dict[str, dict[str, Any]], dimension: str
) -> dict[str, Any]:
    expected_pattern: dict[str, bool] | None = None
    failures = []
    for field_name in ("linear", "child"):
        for segment, row in records_by_field[field_name]["segments"].items():
            groups = row["groups"][dimension]
            pattern: dict[str, bool] = {}
            for value, summary in groups.items():
                if not summary["uniform_inward_label"]:
                    failures.append(f"{field_name}:{segment}:{value}:mixed")
                    continue
                pattern[value] = bool(summary["all_inward"])
            if len(pattern) != len(groups):
                continue
            if expected_pattern is None:
                expected_pattern = pattern
            elif pattern != expected_pattern:
                failures.append(f"{field_name}:{segment}:pattern_drift")
    nontrivial = bool(expected_pattern and any(expected_pattern.values()) and not all(expected_pattern.values()))
    return {
        "stable": bool(expected_pattern is not None and not failures and nontrivial),
        "pattern_all_inward": expected_pattern,
        "failures": failures,
        "nontrivial": nontrivial,
    }


def localization_decision(records_by_field: dict[str, dict[str, Any]]) -> dict[str, Any]:
    sign = _stable_uniform_pattern(records_by_field, "sign")
    band = _stable_uniform_pattern(records_by_field, "band")
    radius = _stable_uniform_pattern(records_by_field, "radius")

    z_sign_exact = bool(
        sign["stable"]
        and sign["pattern_all_inward"] is not None
        and set(sign["pattern_all_inward"].keys()) == {"-1", "1"}
        and sum(bool(v) for v in sign["pattern_all_inward"].values()) == 1
    )

    if z_sign_exact:
        inward_sign = next(k for k, v in sign["pattern_all_inward"].items() if v)
        route = (
            "24/24 split is exactly axial-sign/parity locked across #587/#652 and all four segments; "
            "next representation experiment should be one minimal poloidal/radial axial-parity correction, "
            "not additional generic swirl or time capacity"
        )
        classification = "z_sign_parity"
    elif band["stable"]:
        inward_sign = None
        route = (
            "radial split is stably band-local rather than z-sign locked; next representation experiment should "
            "be one compact poloidal correction confined to the implicated shoulder/tip window"
        )
        classification = "band_local"
    elif radius["stable"]:
        inward_sign = None
        route = (
            "radial split is stably radius-local; next representation experiment should be one compact poloidal "
            "correction confined to the implicated radial window"
        )
        classification = "radius_local"
    else:
        inward_sign = None
        route = (
            "radial split is mixed under the frozen sign/band/radius metadata; do not grow the basis yet and "
            "localize with a finer radial/poloidal Jacobian or sensitivity audit"
        )
        classification = "mixed"

    return {
        "z_sign_exact_split": z_sign_exact,
        "inward_z_sign": inward_sign,
        "band_uniform_pattern": band,
        "radius_uniform_pattern": radius,
        "sign_uniform_pattern": sign,
        "classification": classification,
        "routing": route,
        "basis_growth_justified_by_this_audit": False,
    }


def _build_fields():
    grid = parent.grid
    witness = parent.witness
    field, raw = grid.witness.base.morph.prior.replay_st052.reconstruct()
    redistribution_scale, _, _ = grid.witness.base.morph.prior.base.child_scale(field, raw)
    energy_solve = grid.witness.base.morph.prior.comp.solve_energy_beta(field, raw, redistribution_scale)
    if not energy_solve["root_exists"]:
        raise RuntimeError("#559 frozen energy root disappeared")
    beta = float(energy_solve["beta"])
    if abs(beta - grid.witness.base.EXPECTED_BETA) > 2.0e-10:
        raise RuntimeError("#559 beta identity drift")

    control_fn = lambda p, t: grid.witness.base.control_velocity(field, raw, p, t, redistribution_scale)
    linear_fn = lambda p, t: grid.witness.base.ramp_velocity(field, raw, p, t, redistribution_scale, beta)
    epsilon, normalization = grid.witness.swirl.derive_epsilon(linear_fn)
    child_fn = lambda p, t: grid.witness.combined_velocity(
        field, raw, p, t, redistribution_scale, beta, epsilon
    )
    return control_fn, linear_fn, child_fn, float(epsilon), normalization


def run(out: Path) -> dict[str, Any]:
    control_fn, linear_fn, child_fn, epsilon, normalization = _build_fields()
    records = {
        "control": field_records(control_fn),
        "linear": field_records(linear_fn),
        "child": field_records(child_fn),
    }
    decision = localization_decision(records)

    # Preserve the parent-observed 24/48 split as an implementation-consistency guard.
    parent_split_guard = bool(
        all(row["inward_path_count"] == 24 for name in ("linear", "child") for row in records[name]["segments"].values())
    )
    if not parent_split_guard:
        raise RuntimeError("frozen #683 24/48 segment split drifted")

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "candidate_pair": {"baseline_pr": SOURCE_TEMPORAL_PR, "challenger_pr": SOURCE_WITNESS_PR},
        "public_source": {"url": PUBLIC_SOURCE_URL, "qualitative_only": True},
        "frozen_witness": {
            "swirl_a": parent.witness.SWIRL_A,
            "shoulder_lambda": parent.witness.SHOULDER_LAMBDA,
            "compact_swirl_epsilon": epsilon,
            "parameter_grid_scan_performed": False,
            "optimization_performed": False,
        },
        "trajectory_protocol": {
            "seed_radii": list(parent.base.SEED_RADII),
            "seed_bands": [[name, z] for name, z in parent.base.SEED_BANDS],
            "seed_angles": parent.base.SEED_ANGLES,
            "path_count": 48,
            "checkpoint_times": list(parent.CHECKPOINT_TIMES),
            "checkpoint_indices": list(parent.CHECKPOINT_INDICES),
            "solver": parent.base.SOLVER_METHOD,
        },
        "normalization": normalization,
        "records": records,
        "parent_24_of_48_split_replayed": parent_split_guard,
        "localization_decision": decision,
        "source_652_target_free_pareto_rejection_remains_binding": True,
        "source_683_temporal_inward_spiral_failure_remains_binding": True,
        **TRUTH,
    }
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    print("classification=", report["localization_decision"]["classification"])
    print("z_sign_exact_split=", report["localization_decision"]["z_sign_exact_split"])
    print("inward_z_sign=", report["localization_decision"]["inward_z_sign"])
    print("routing=", report["localization_decision"]["routing"])


if __name__ == "__main__":
    main()
