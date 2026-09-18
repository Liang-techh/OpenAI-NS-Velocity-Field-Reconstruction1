"""Truth-bounded material-path observables for the retained ST006 field.

This module measures candidate-side kinematics only.  It does not infer OpenAI
seed locations, camera geometry, hidden frame times, or PDE validity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from . import CANDIDATE_SHA256, load_best

SCHEMA = "st006_material_path_observables_v1"
REGISTERED_TIME_INTERVAL = (0.25, 0.75)
EVALUATION_BOX = (-2.0, 2.0)
SEED_RADII = (0.6, 0.9, 1.2)
SEED_Z = (-0.3, 0.3)
SEED_ANGLES = 8
OUTPUT_SAMPLES = 33
SOLVER_METHOD = "DOP853"
SOLVER_RTOL = 1.0e-9
SOLVER_ATOL = 1.0e-11
SOLVER_MAX_STEP = 0.01

VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


def _canonical_hash(payload: dict) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def _fixed_seed_table() -> tuple[np.ndarray, list[dict]]:
    rows: list[list[float]] = []
    metadata: list[dict] = []
    for radius in SEED_RADII:
        for angle_index in range(SEED_ANGLES):
            angle = 2.0 * np.pi * angle_index / SEED_ANGLES
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            for z in SEED_Z:
                rows.append([x, y, z])
                metadata.append(
                    {
                        "radius": float(radius),
                        "angle_index": int(angle_index),
                        "angle_radians": float(angle),
                        "z": float(z),
                    }
                )
    seeds = np.asarray(rows, dtype=float)
    seeds.setflags(write=False)
    return seeds, metadata


def _evaluate_velocity(velocity: VelocityCallable, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != points.shape:
        raise ValueError("velocity must return shape (n,3) for an (n,3) point array")
    if not np.isfinite(values).all():
        raise ValueError("velocity returned non-finite values")
    return values


def measure_material_paths(
    velocity: VelocityCallable,
    *,
    candidate: str,
    candidate_sha256: str,
    provenance: str,
) -> dict:
    """Measure a frozen, target-free material-path contract for one velocity field."""
    if not callable(velocity):
        raise TypeError("velocity must be callable")
    if not candidate or not isinstance(candidate, str):
        raise ValueError("candidate must be a nonempty string")
    if not provenance or not isinstance(provenance, str):
        raise ValueError("provenance must be a nonempty string")
    if not isinstance(candidate_sha256, str) or len(candidate_sha256) != 64:
        raise ValueError("candidate_sha256 must be a 64-character SHA-256 hex digest")
    try:
        int(candidate_sha256, 16)
    except ValueError as exc:
        raise ValueError("candidate_sha256 must be hexadecimal") from exc

    seeds, seed_metadata = _fixed_seed_table()
    t0, t1 = REGISTERED_TIME_INTERVAL
    times = np.linspace(t0, t1, OUTPUT_SAMPLES, dtype=float)
    n_paths = len(seeds)

    # Preflight catches malformed/zero fields before the adaptive solver obscures the cause.
    initial_velocity = _evaluate_velocity(velocity, seeds, t0)
    initial_speed_rms = float(np.sqrt(np.mean(np.sum(initial_velocity * initial_velocity, axis=1))))
    if initial_speed_rms <= 128.0 * np.finfo(float).eps:
        raise ValueError("sampled initial velocity is numerically zero")

    def rhs(time: float, flat_positions: np.ndarray) -> np.ndarray:
        points = np.asarray(flat_positions, dtype=float).reshape(n_paths, 3)
        return _evaluate_velocity(velocity, points, time).reshape(-1)

    solution = solve_ivp(
        rhs,
        (t0, t1),
        seeds.reshape(-1),
        method=SOLVER_METHOD,
        t_eval=times,
        rtol=SOLVER_RTOL,
        atol=SOLVER_ATOL,
        max_step=SOLVER_MAX_STEP,
    )
    if not solution.success:
        raise RuntimeError(f"material-path integration failed: {solution.message}")
    positions = np.asarray(solution.y.T, dtype=float).reshape(len(times), n_paths, 3)
    if not np.isfinite(positions).all():
        raise ValueError("material-path integration produced non-finite positions")

    all_inside_box = bool(np.all((positions >= EVALUATION_BOX[0]) & (positions <= EVALUATION_BOX[1])))
    if not all_inside_box:
        raise ValueError("sampled material path left the registered evaluation box")

    velocities = np.empty_like(positions)
    for time_index, time in enumerate(times):
        velocities[time_index] = _evaluate_velocity(velocity, positions[time_index], float(time))
    speeds = np.linalg.norm(velocities, axis=2)

    radii = np.hypot(positions[:, :, 0], positions[:, :, 1])
    angles = np.unwrap(np.arctan2(positions[:, :, 1], positions[:, :, 0]), axis=0)
    abs_z = np.abs(positions[:, :, 2])
    radius_change = radii[-1] - radii[0]
    axial_abs_change = abs_z[-1] - abs_z[0]
    angular_change = angles[-1] - angles[0]
    turns = angular_change / (2.0 * np.pi)
    speed_change = speeds[-1] - speeds[0]
    path_lengths = np.sum(np.linalg.norm(np.diff(positions, axis=0), axis=2), axis=0)

    per_path: list[dict] = []
    for index, meta in enumerate(seed_metadata):
        per_path.append(
            {
                "path_index": index,
                "seed": meta,
                "radius_start": float(radii[0, index]),
                "radius_end": float(radii[-1, index]),
                "radius_change": float(radius_change[index]),
                "abs_z_start": float(abs_z[0, index]),
                "abs_z_end": float(abs_z[-1, index]),
                "abs_z_change": float(axial_abs_change[index]),
                "angular_change_radians": float(angular_change[index]),
                "signed_turns": float(turns[index]),
                "speed_start": float(speeds[0, index]),
                "speed_end": float(speeds[-1, index]),
                "speed_change": float(speed_change[index]),
                "sampled_path_length": float(path_lengths[index]),
            }
        )

    # Seeds are ordered as radius, angle, z=-0.3, z=+0.3.  Pair separation is a
    # direct material-line diagnostic and does not require a visual target value.
    pair_rows: list[dict] = []
    for radius_index, radius in enumerate(SEED_RADII):
        for angle_index in range(SEED_ANGLES):
            base = (radius_index * SEED_ANGLES + angle_index) * 2
            minus_index, plus_index = base, base + 1
            separation = positions[:, plus_index, 2] - positions[:, minus_index, 2]
            pair_rows.append(
                {
                    "radius": float(radius),
                    "angle_index": int(angle_index),
                    "initial_axial_separation": float(separation[0]),
                    "final_axial_separation": float(separation[-1]),
                    "axial_separation_change": float(separation[-1] - separation[0]),
                    "axial_separation_ratio": float(separation[-1] / separation[0]),
                    "mean_pair_radius_start": float(0.5 * (radii[0, minus_index] + radii[0, plus_index])),
                    "mean_pair_radius_end": float(0.5 * (radii[-1, minus_index] + radii[-1, plus_index])),
                }
            )

    pair_changes = np.asarray([row["axial_separation_change"] for row in pair_rows], dtype=float)
    pair_ratios = np.asarray([row["axial_separation_ratio"] for row in pair_rows], dtype=float)
    abs_turns = np.abs(turns)

    summary = {
        "path_count": int(n_paths),
        "paired_material_line_count": int(len(pair_rows)),
        "initial_velocity_rms": initial_speed_rms,
        "mean_radius_change": float(np.mean(radius_change)),
        "median_radius_change": float(np.median(radius_change)),
        "inward_path_count": int(np.sum(radius_change < 0.0)),
        "outward_path_count": int(np.sum(radius_change > 0.0)),
        "mean_absolute_turns": float(np.mean(abs_turns)),
        "median_absolute_turns": float(np.median(abs_turns)),
        "minimum_absolute_turns": float(np.min(abs_turns)),
        "maximum_absolute_turns": float(np.max(abs_turns)),
        "mean_abs_z_change": float(np.mean(axial_abs_change)),
        "axially_outward_path_count": int(np.sum(axial_abs_change > 0.0)),
        "axially_inward_path_count": int(np.sum(axial_abs_change < 0.0)),
        "mean_pair_axial_separation_change": float(np.mean(pair_changes)),
        "median_pair_axial_separation_change": float(np.median(pair_changes)),
        "mean_pair_axial_separation_ratio": float(np.mean(pair_ratios)),
        "pair_axial_separation_growth_count": int(np.sum(pair_changes > 0.0)),
        "pair_axial_separation_shrink_count": int(np.sum(pair_changes < 0.0)),
        "mean_speed_change": float(np.mean(speed_change)),
        "mean_sampled_path_length": float(np.mean(path_lengths)),
    }

    receipt = {
        "schema": SCHEMA,
        "task_id": "CR-A9-040",
        "candidate": candidate,
        "candidate_sha256": candidate_sha256,
        "provenance": provenance,
        "registered_contract": {
            "physical_domain": "R^3",
            "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
            "time_interval": [t0, t1],
            "support": "r < 2 and abs(z) < 2 with smooth zero extension",
            "seed_radii": list(SEED_RADII),
            "seed_z": list(SEED_Z),
            "seed_angles": SEED_ANGLES,
            "seed_count": int(n_paths),
            "output_samples": OUTPUT_SAMPLES,
            "solver": {
                "method": SOLVER_METHOD,
                "rtol": SOLVER_RTOL,
                "atol": SOLVER_ATOL,
                "max_step": SOLVER_MAX_STEP,
                "status": "autonomous_numerical_choice_not_visual_acceptance_threshold",
            },
        },
        "summary": summary,
        "per_path": per_path,
        "paired_axial_material_lines": pair_rows,
        "public_observable_evidence": {
            "inward_spiraling": {
                "status": "measured_candidate_side_only",
                "evidence": ["radius_change", "signed_turns"],
                "no_openai_seed_matching": True,
            },
            "axial_stretching": {
                "status": "measured_candidate_side_only",
                "evidence": ["paired_axial_separation_change", "abs_z_change"],
                "no_openai_seed_matching": True,
            },
        },
        "interpretation_boundary": {
            "claim_scope": "candidate_material_kinematics_only",
            "all_sampled_points_inside_registered_box": all_inside_box,
            "openai_seed_locations_used": False,
            "hidden_time_alignment_used": False,
            "camera_registration_used": False,
            "image_fit_used": False,
            "visual_acceptance_threshold_defined": False,
            "pde_derivatives_evaluated": False,
            "pde_acceptance_evidence": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    receipt["report_sha256"] = _canonical_hash(receipt)
    return receipt


def measure_st006_material_paths() -> dict:
    field = load_best()
    if field.sha256 != CANDIDATE_SHA256 or field.pde_validated is not False:
        raise ValueError("retained ST006 identity/truth state changed")
    return measure_material_paths(
        field.at_points,
        candidate="ST006",
        candidate_sha256=CANDIDATE_SHA256,
        provenance="research_baseline.load_best().at_points",
    )


def write_st006_material_path_receipt(path: str | Path) -> dict:
    output = Path(path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    report = measure_st006_material_paths()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return report


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = write_st006_material_path_receipt(args.output)
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    _main()
