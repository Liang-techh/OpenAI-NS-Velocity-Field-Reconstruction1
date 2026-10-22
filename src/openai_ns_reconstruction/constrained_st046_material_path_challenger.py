"""Truth-bounded material-path replay for the frozen ST046-A challenger.

The contract is inherited unchanged from CR-A9-040/ST006: fixed seeds, time
interval, and ODE tolerances. This module compares candidate-side kinematics
only; it never upgrades visual evidence into PDE acceptance.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

TASK_ID = "CR-A9-042"
SCHEMA = "st046a_material_path_challenger_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"
ST046_PR = 366
ST046_HEAD_SHA = "04d2fa64798d9639bddadb289cc9c23f116eec56"
ST046_ORIGINAL_RAW_SHA256 = "94f5eeb0d94568587c9f3d88e69876a8051632830ded796f7c8b0db8e6707619"
ST046_RECIPE_ID = "ST046-A"
ST006_REFERENCE_HEAD = "3911718884f021f0c0e3ba41ed694b9688764084"
ST006_REFERENCE_REPORT_SHA256 = "2595c8dafd0716b3b07cb9dc761b112e2178bf44200b9c50eff3c000c6166772"
ST006_REFERENCE = {
    "candidate": "ST006",
    "path_count": 48,
    "inward_path_count": 48,
    "mean_radius_change": -0.08601065505093726,
    "mean_absolute_turns": 0.024485286231482325,
    "minimum_absolute_turns": 0.003880821959757115,
    "maximum_absolute_turns": 0.0436185396886357,
    "paired_material_line_count": 24,
    "pair_growth_count": 16,
    "pair_shrink_count": 8,
    "mean_pair_separation_change": 0.07130825328092755,
    "mean_pair_separation_ratio": 1.1188470888015458,
}

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
    return np.asarray(rows, dtype=float), metadata


def _evaluate_velocity(velocity: VelocityCallable, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != points.shape:
        raise ValueError("velocity must return shape (n,3) for an (n,3) point array")
    if not np.isfinite(values).all():
        raise ValueError("velocity returned non-finite values")
    return values


def measure_material_paths(velocity: VelocityCallable) -> dict:
    """Replay the frozen CR-A9-040 material-path contract on one velocity field."""
    if not callable(velocity):
        raise TypeError("velocity must be callable")
    seeds, seed_metadata = _fixed_seed_table()
    t0, t1 = REGISTERED_TIME_INTERVAL
    times = np.linspace(t0, t1, OUTPUT_SAMPLES, dtype=float)
    n_paths = len(seeds)

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
    if not np.all((positions >= EVALUATION_BOX[0]) & (positions <= EVALUATION_BOX[1])):
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
    turns = (angles[-1] - angles[0]) / (2.0 * np.pi)
    speed_change = speeds[-1] - speeds[0]
    path_lengths = np.sum(np.linalg.norm(np.diff(positions, axis=0), axis=2), axis=0)

    per_path: list[dict] = []
    for index, meta in enumerate(seed_metadata):
        per_path.append(
            {
                "path_index": index,
                "seed": meta,
                "radius_change": float(radius_change[index]),
                "abs_z_change": float(axial_abs_change[index]),
                "signed_turns": float(turns[index]),
                "speed_change": float(speed_change[index]),
                "sampled_path_length": float(path_lengths[index]),
            }
        )

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
                }
            )

    pair_changes = np.asarray([row["axial_separation_change"] for row in pair_rows], dtype=float)
    pair_ratios = np.asarray([row["axial_separation_ratio"] for row in pair_rows], dtype=float)
    abs_turns = np.abs(turns)
    return {
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
        "per_path": per_path,
        "paired_axial_material_lines": pair_rows,
    }


def _git_head(repository_root: Path) -> str:
    completed = subprocess.run(
        ["git", "-C", str(repository_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def load_st046a_velocity(candidate_root: str | Path) -> tuple[VelocityCallable, dict]:
    """Load the frozen ST046-A mathematical field from the pinned PR head."""
    root = Path(candidate_root).resolve()
    if _git_head(root) != ST046_HEAD_SHA:
        raise ValueError("ST046 checkout is not the pinned audited head")
    experiment = root / "experiments" / "root_st046"
    records = json.loads((experiment / "recipes.json").read_text())
    record = records.get(ST046_RECIPE_ID)
    if not isinstance(record, dict):
        raise ValueError("missing ST046-A recipe")
    if record.get("parent_id") != "ST045-H" or record.get("pde_validated") is not False:
        raise ValueError("ST046-A recipe parent/truth state drift")

    previous_path = list(sys.path)
    try:
        sys.path.insert(0, str(experiment))
        for name in ("replay_st046", "acceleration_fit", "mechanism_audit"):
            sys.modules.pop(name, None)
        replay = importlib.import_module("replay_st046")
        family, raw = replay.reconstruct(ST046_RECIPE_ID)
    finally:
        sys.path[:] = previous_path

    def velocity(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        times = np.full(len(points), float(time), dtype=float)
        values, _pressure = family.fields(raw, points, times)
        return np.asarray(values, dtype=float)

    metadata = {
        "source_pr": ST046_PR,
        "source_head_sha": ST046_HEAD_SHA,
        "recipe_id": ST046_RECIPE_ID,
        "original_raw_candidate_sha256": ST046_ORIGINAL_RAW_SHA256,
        "recipe_modifiers_sha256": record.get("modifiers_sha256"),
        "reconstruction_note": (
            "PR #366 states the readable recipe reproduces the mathematical field and "
            "reference values, while regenerated JSON metadata is not byte-identical to the original raw candidate."
        ),
        "pde_validated": False,
    }
    return velocity, metadata


def build_receipt(summary: dict, candidate_metadata: dict) -> dict:
    """Bind one measured challenger summary to ST006 without defining a pass score."""
    required = (
        "path_count",
        "inward_path_count",
        "mean_radius_change",
        "mean_absolute_turns",
        "mean_pair_axial_separation_ratio",
        "pair_axial_separation_growth_count",
        "pair_axial_separation_shrink_count",
    )
    if any(key not in summary for key in required):
        raise ValueError("incomplete challenger summary")
    if summary["path_count"] != 48:
        raise ValueError("frozen path population drift")
    if candidate_metadata.get("pde_validated") is not False:
        raise ValueError("challenger PDE truth state must remain false")

    st006 = ST006_REFERENCE
    comparison = {
        "mean_radius_change_delta_vs_st006": float(summary["mean_radius_change"] - st006["mean_radius_change"]),
        "absolute_radial_contraction_ratio_vs_st006": float(
            abs(summary["mean_radius_change"]) / abs(st006["mean_radius_change"])
        ),
        "mean_absolute_turns_delta_vs_st006": float(summary["mean_absolute_turns"] - st006["mean_absolute_turns"]),
        "mean_absolute_turns_ratio_vs_st006": float(summary["mean_absolute_turns"] / st006["mean_absolute_turns"]),
        "mean_pair_separation_ratio_delta_vs_st006": float(
            summary["mean_pair_axial_separation_ratio"] - st006["mean_pair_separation_ratio"]
        ),
        "inward_path_count_delta_vs_st006": int(summary["inward_path_count"] - st006["inward_path_count"]),
    }
    payload = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_main_sha": BASE_MAIN_SHA,
        "challenger": {
            "candidate": ST046_RECIPE_ID,
            **candidate_metadata,
        },
        "registered_contract": {
            "physical_domain": "R^3",
            "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
            "time_interval": list(REGISTERED_TIME_INTERVAL),
            "support": "r < 2 and abs(z) < 2 with smooth zero extension",
            "seed_radii": list(SEED_RADII),
            "seed_z": list(SEED_Z),
            "seed_angles": SEED_ANGLES,
            "seed_count": 48,
            "output_samples": OUTPUT_SAMPLES,
            "solver": {
                "method": SOLVER_METHOD,
                "rtol": SOLVER_RTOL,
                "atol": SOLVER_ATOL,
                "max_step": SOLVER_MAX_STEP,
                "status": "frozen CR-A9-040 numerical contract; not a visual/PDE acceptance threshold",
            },
        },
        "challenger_measurement": summary,
        "st006_reference": {
            "source_task": "CR-A9-040",
            "source_head_sha": ST006_REFERENCE_HEAD,
            "source_report_sha256": ST006_REFERENCE_REPORT_SHA256,
            "summary": ST006_REFERENCE,
        },
        "descriptive_comparison": comparison,
        "external_method": {
            "classification": "direct_migration_public_api_only",
            "source_repo": "scipy/scipy",
            "source_commit": "eff78058ed0d4cb8e1f0b5c5585d62d712e948d2",
            "license": "BSD-3-Clause",
            "scope": "scipy.integrate.solve_ivp with DOP853 for the frozen material-path IVP",
            "copied_source_code": False,
        },
        "truth_boundary": {
            "comparison_is_descriptive_not_acceptance": True,
            "visual_acceptance_threshold_defined": False,
            "openai_numerical_velocity_used": False,
            "openai_seed_locations_used": False,
            "hidden_time_alignment_used": False,
            "camera_registration_used": False,
            "free_residual_force_used": False,
            "collapsed_velocity_accepted": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    payload["report_sha256"] = _canonical_hash(payload)
    return payload


def make_st046a_receipt(candidate_root: str | Path) -> dict:
    velocity, metadata = load_st046a_velocity(candidate_root)
    summary = measure_material_paths(velocity)
    return build_receipt(summary, metadata)


def write_receipt(candidate_root: str | Path, output: str | Path) -> dict:
    path = Path(output)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {path}")
    report = make_st046a_receipt(candidate_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return report


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = write_receipt(args.candidate_root, args.output)
    print(json.dumps({
        "task_id": report["task_id"],
        "report_sha256": report["report_sha256"],
        "challenger_summary": {
            key: report["challenger_measurement"][key]
            for key in (
                "inward_path_count",
                "mean_radius_change",
                "mean_absolute_turns",
                "maximum_absolute_turns",
                "pair_axial_separation_growth_count",
                "pair_axial_separation_shrink_count",
                "mean_pair_axial_separation_ratio",
                "mean_speed_change",
            )
        },
        "comparison": report["descriptive_comparison"],
    }, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    _main()
