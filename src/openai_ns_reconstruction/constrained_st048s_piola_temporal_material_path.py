"""Replay frozen Agent-9 material paths on the first ST048-S temporal Piola crossing.

This is a visualization-side governance experiment only.  It reconstructs the
frozen ST048-S field from PR #390, applies Agent 7's exact axial Piola map with
beta(t)=.075+.025*(4t-2)^2 from PR #418, restores E(.25)=1 by one positive
common scale, and compares the resulting material paths with the fixed beta=.075
warp under the unchanged CR-A9 material-path contract.

No production child, pressure/forcing transfer, PDE promotion, visual target
threshold, hidden OpenAI time/camera/seed, or recovered numerical field is used.
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
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

TASK_ID = "CR-A9-047"
SCHEMA = "st048s_piola_temporal_material_path_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"
ST048_SOURCE_PR = 390
ST048_SOURCE_HEAD = "97695a86f85ce68fb4ae70c41fc81c904d655183"
ST048_PARENT_ID = "ST047-E"
ST048_PARENT_RAW_SHA256 = "dfe6e51af93c9c42f1322798b82f89523e6a193c1b6495894cb60f63b4eb07e0"
ST048S_RAW_SHA256 = "6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09"
AGENT7_SOURCE_PR = 418
AGENT7_SOURCE_HEAD = "ae10c3ccacd16e1f5e0054502f1dc5a6f68320b5"
AGENT7_TASK_ID = "CR003-ST048S-PIOLA-TEMPORAL-COEFFICIENT-062"
BASE_BETA = 0.075
TEMPORAL_GAMMA = 0.025
AGENT7_REPORTED_TEMPORAL_SCALE = 0.998791257501046
PHYSICAL_Z_CUT = 2.0
REFERENCE_ENERGY = 1.0
ENERGY_ORDER = 64

ST048S_REFERENCE_HEAD = "50da502ef03dae05fad7327b3371a774f18e7d23"
ST048S_REFERENCE_REPORT_SHA256 = "6c5aa6e4472c7feac96a01fc07e0a5705d0a9f1c9ccdec943ee49f78b8769be0"
ST048S_REFERENCE = {
    "candidate": "ST048-S",
    "path_count": 48,
    "inward_path_count": 48,
    "mean_radius_change": -0.0853545661,
    "mean_absolute_turns": 0.0221757297,
    "minimum_absolute_turns": 0.0061777407,
    "maximum_absolute_turns": 0.0346689104,
    "paired_material_line_count": 24,
    "pair_growth_count": 16,
    "pair_shrink_count": 8,
    "mean_pair_separation_change": 0.0778727503,
    "mean_pair_separation_ratio": 1.1297879172,
    "mean_speed_change": 0.0674821286,
}
ST006_REFERENCE = {
    "candidate": "ST006",
    "mean_radius_change": -0.08601065505093726,
    "mean_absolute_turns": 0.024485286231482325,
    "maximum_absolute_turns": 0.0436185396886357,
    "mean_pair_separation_change": 0.07130825328092755,
    "mean_pair_separation_ratio": 1.1188470888015458,
    "pair_growth_count": 16,
    "pair_shrink_count": 8,
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

SCIPY_SOURCE = {
    "repository": "scipy/scipy",
    "commit": "eff78058ed0d4cb8e1f0b5c5585d62d712e948d2",
    "api": "scipy.integrate.solve_ivp / DOP853",
    "license": "BSD-3-Clause",
    "classification": "direct migration / public API only",
    "copied_upstream_implementation": False,
}

TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "production_candidate_selected": False,
    "production_gamma_selected": False,
    "pressure_or_force_transferred_from_parent": False,
    "held_out_pde_residual_evaluated": False,
    "held_out_pde_residual_transferred_from_parent": False,
    "visual_acceptance_threshold_defined": False,
    "comparison_is_descriptive_not_acceptance": True,
    "hidden_openai_time_camera_seed_or_velocity_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


def _canonical_hash(payload: dict) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def _git_head(repository_root: Path) -> str:
    completed = subprocess.run(
        ["git", "-C", str(repository_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _validate_gamma(gamma: float) -> float:
    gamma = float(gamma)
    if not np.isfinite(gamma) or not 0.0 <= gamma <= 0.125:
        raise ValueError("gamma must lie in [0,.125]")
    return gamma


def beta_at_time(time: float, gamma: float = TEMPORAL_GAMMA, *, base_beta: float = BASE_BETA) -> float:
    time = float(time)
    gamma = _validate_gamma(gamma)
    base_beta = float(base_beta)
    if not np.isfinite(time) or not REGISTERED_TIME_INTERVAL[0] <= time <= REGISTERED_TIME_INTERVAL[1]:
        raise ValueError("time must lie in [0.25,0.75]")
    if not np.isfinite(base_beta) or not 0.0 <= base_beta <= 0.20:
        raise ValueError("base_beta must lie in [0,.20]")
    beta = base_beta + gamma * (4.0 * time - 2.0) ** 2
    if not 0.0 <= beta <= 0.20 + 1.0e-15:
        raise ValueError("temporal schedule exceeds the screened beta range")
    return float(beta)


def _warp_z_and_jacobian(z, beta: float):
    beta = float(beta)
    if not np.isfinite(beta) or not 0.0 <= beta <= 0.20 + 1.0e-15:
        raise ValueError("beta must lie in [0,.20]")
    z = np.asarray(z, dtype=float)
    if not np.all(np.isfinite(z)):
        raise ValueError("z must be finite")
    s = z / PHYSICAL_Z_CUT
    inside = np.abs(s) < 1.0
    one_minus = np.maximum(1.0 - s * s, 0.0)
    bump = one_minus**4
    mapped = np.where(inside, z * (1.0 - beta * bump), z)
    jac_inside = 1.0 - beta * bump + 8.0 * beta * s * s * one_minus**3
    jac = np.where(inside, jac_inside, 1.0)
    if np.any(jac <= 0.0) or not np.all(np.isfinite(mapped)):
        raise RuntimeError("Piola coordinate warp lost orientation or became nonfinite")
    return mapped, jac


class TemporalPiolaWarpVelocity:
    """Time-dependent contravariant Piola pullback followed by a common scale."""

    def __init__(self, parent: VelocityCallable, *, gamma: float, scale: float):
        if not callable(parent):
            raise TypeError("parent must be callable")
        self.parent = parent
        self.gamma = _validate_gamma(gamma)
        self.scale = float(scale)
        if not np.isfinite(self.scale) or self.scale <= 0.0:
            raise ValueError("scale must be positive and finite")
        beta_at_time(0.25, self.gamma)
        beta_at_time(0.75, self.gamma)

    def __call__(self, points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
            raise ValueError("points must be finite shape (n,3)")
        beta = beta_at_time(float(time), self.gamma)
        mapped_z, jac = _warp_z_and_jacobian(points[:, 2], beta)
        mapped = np.array(points, copy=True)
        mapped[:, 2] = mapped_z
        velocity = np.asarray(self.parent(mapped, float(time)), dtype=float)
        if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
            raise ValueError("parent velocity returned invalid values")
        out = self.scale * velocity
        out[:, 0] *= jac
        out[:, 1] *= jac
        return out


def load_st048s_velocity(candidate_root: str | Path) -> tuple[VelocityCallable, dict]:
    root = Path(candidate_root).resolve()
    if _git_head(root) != ST048_SOURCE_HEAD:
        raise ValueError("ST048 checkout is not the pinned audited head")
    experiment = root / "experiments" / "root_st048"
    records = json.loads((experiment / "recipes.json").read_text())
    record = records.get("ST048-S")
    if not isinstance(record, dict):
        raise ValueError("missing ST048-S recipe")
    if (
        record.get("parent_id") != ST048_PARENT_ID
        or record.get("parent_original_raw_sha256") != ST048_PARENT_RAW_SHA256
        or record.get("original_local_raw_sha256") != ST048S_RAW_SHA256
        or record.get("pde_validated") is not False
        or record.get("source_correspondence_verified") is not False
    ):
        raise ValueError("ST048-S provenance/truth state drift")

    previous_path = list(sys.path)
    stale_modules = (
        "replay_st048", "boundary_shear", "replay_st047", "continuation",
        "structure_and_particles", "replay_st046", "acceleration_fit",
        "mechanism_audit", "localized_model", "spacetime", "validate",
    )
    try:
        sys.path.insert(0, str(experiment))
        for name in stale_modules:
            sys.modules.pop(name, None)
        replay = importlib.import_module("replay_st048")
        family, raw = replay.reconstruct("ST048-S")
    finally:
        sys.path[:] = previous_path

    def velocity(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        times = np.full(len(points), float(time), dtype=float)
        values, _pressure = family.fields(raw, points, times)
        return np.asarray(values, dtype=float)

    metadata = {
        "source_pr": ST048_SOURCE_PR,
        "source_head_sha": ST048_SOURCE_HEAD,
        "candidate_id": "ST048-S",
        "parent_id": ST048_PARENT_ID,
        "parent_original_raw_sha256": ST048_PARENT_RAW_SHA256,
        "original_raw_candidate_sha256": ST048S_RAW_SHA256,
        "recipe_modifiers_sha256": record.get("modifiers_sha256"),
        "pde_validated": False,
        "source_correspondence_verified": False,
    }
    return velocity, metadata


def axisymmetric_energy(velocity: VelocityCallable, *, time: float = 0.25, order: int = ENERGY_ORDER) -> float:
    if int(order) < 12:
        raise ValueError("energy quadrature order must be >=12")
    node, weight = leggauss(int(order))
    radius = node + 1.0
    z = 2.0 * node
    rr, zz = np.meshgrid(radius, z, indexing="ij")
    points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
    values = np.asarray(velocity(points, float(time)), dtype=float)
    speed_sq = np.sum(values * values, axis=1).reshape(rr.shape)
    return float(np.sum((weight[:, None] * (2.0 * weight)[None, :]) * (np.pi * rr * speed_sq)))


def normalized_temporal_piola_velocity(parent: VelocityCallable, gamma: float) -> tuple[TemporalPiolaWarpVelocity, dict]:
    gamma = _validate_gamma(gamma)
    raw = TemporalPiolaWarpVelocity(parent, gamma=gamma, scale=1.0)
    raw_energy = axisymmetric_energy(raw, time=0.25, order=ENERGY_ORDER)
    if not np.isfinite(raw_energy) or raw_energy <= 0.0:
        raise ValueError("warped reference energy must be positive and finite")
    scale = float(np.sqrt(REFERENCE_ENERGY / raw_energy))
    warped = TemporalPiolaWarpVelocity(parent, gamma=gamma, scale=scale)
    normalized_energy = axisymmetric_energy(warped, time=0.25, order=ENERGY_ORDER)
    return warped, {
        "base_beta": BASE_BETA,
        "gamma": gamma,
        "beta_at_025": beta_at_time(0.25, gamma),
        "beta_at_050": beta_at_time(0.50, gamma),
        "beta_at_075": beta_at_time(0.75, gamma),
        "energy_quadrature_order": ENERGY_ORDER,
        "raw_reference_energy": raw_energy,
        "normalization_scale": scale,
        "normalized_reference_energy": normalized_energy,
    }


def _fixed_seed_table() -> tuple[np.ndarray, list[dict]]:
    rows = []
    metadata = []
    for radius in SEED_RADII:
        for angle_index in range(SEED_ANGLES):
            angle = 2.0 * np.pi * angle_index / SEED_ANGLES
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            for z in SEED_Z:
                rows.append([x, y, z])
                metadata.append({"radius": radius, "angle_index": angle_index, "angle_radians": angle, "z": z})
    return np.asarray(rows, dtype=float), metadata


def _evaluate_velocity(velocity: VelocityCallable, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != points.shape or not np.isfinite(values).all():
        raise ValueError("velocity must return finite shape (n,3)")
    return values


def measure_material_paths(velocity: VelocityCallable) -> dict:
    if not callable(velocity):
        raise TypeError("velocity must be callable")
    seeds, seed_metadata = _fixed_seed_table()
    t0, t1 = REGISTERED_TIME_INTERVAL
    times = np.linspace(t0, t1, OUTPUT_SAMPLES, dtype=float)
    n_paths = len(seeds)
    initial = _evaluate_velocity(velocity, seeds, t0)
    initial_velocity_rms = float(np.sqrt(np.mean(np.sum(initial * initial, axis=1))))
    if initial_velocity_rms <= 128.0 * np.finfo(float).eps:
        raise ValueError("sampled initial velocity is numerically zero")

    def rhs(time: float, flat_positions: np.ndarray) -> np.ndarray:
        points = np.asarray(flat_positions, dtype=float).reshape(n_paths, 3)
        return _evaluate_velocity(velocity, points, time).reshape(-1)

    solution = solve_ivp(
        rhs, (t0, t1), seeds.reshape(-1), method=SOLVER_METHOD, t_eval=times,
        rtol=SOLVER_RTOL, atol=SOLVER_ATOL, max_step=SOLVER_MAX_STEP,
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

    pair_rows = []
    for radius_index, radius in enumerate(SEED_RADII):
        for angle_index in range(SEED_ANGLES):
            base = (radius_index * SEED_ANGLES + angle_index) * 2
            separation = positions[:, base + 1, 2] - positions[:, base, 2]
            pair_rows.append({
                "radius": float(radius),
                "angle_index": int(angle_index),
                "initial_axial_separation": float(separation[0]),
                "final_axial_separation": float(separation[-1]),
                "axial_separation_change": float(separation[-1] - separation[0]),
                "axial_separation_ratio": float(separation[-1] / separation[0]),
            })

    pair_changes = np.asarray([row["axial_separation_change"] for row in pair_rows], dtype=float)
    pair_ratios = np.asarray([row["axial_separation_ratio"] for row in pair_rows], dtype=float)
    abs_turns = np.abs(turns)
    return {
        "path_count": int(n_paths),
        "paired_material_line_count": int(len(pair_rows)),
        "initial_velocity_rms": initial_velocity_rms,
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
        "per_path": [
            {
                "path_index": int(i), "seed": seed_metadata[i],
                "radius_change": float(radius_change[i]), "abs_z_change": float(axial_abs_change[i]),
                "signed_turns": float(turns[i]), "speed_change": float(speed_change[i]),
                "sampled_path_length": float(path_lengths[i]),
            }
            for i in range(n_paths)
        ],
        "paired_axial_material_lines": pair_rows,
    }


def _comparison(summary: dict, reference: dict) -> dict:
    return {
        "mean_radius_change_delta": float(summary["mean_radius_change"] - reference["mean_radius_change"]),
        "absolute_radial_contraction_ratio": float(abs(summary["mean_radius_change"]) / abs(reference["mean_radius_change"])),
        "mean_absolute_turns_delta": float(summary["mean_absolute_turns"] - reference["mean_absolute_turns"]),
        "mean_absolute_turns_ratio": float(summary["mean_absolute_turns"] / reference["mean_absolute_turns"]),
        "maximum_absolute_turns_delta": float(summary["maximum_absolute_turns"] - reference["maximum_absolute_turns"]),
        "maximum_absolute_turns_ratio": float(summary["maximum_absolute_turns"] / reference["maximum_absolute_turns"]),
        "mean_pair_separation_change_delta": float(summary["mean_pair_axial_separation_change"] - reference["mean_pair_separation_change"]),
        "mean_pair_separation_change_ratio": float(summary["mean_pair_axial_separation_change"] / reference["mean_pair_separation_change"]),
        "mean_pair_separation_ratio_delta": float(summary["mean_pair_axial_separation_ratio"] - reference["mean_pair_separation_ratio"]),
        "pair_growth_count_delta": int(summary["pair_axial_separation_growth_count"] - reference["pair_growth_count"]),
        "pair_shrink_count_delta": int(summary["pair_axial_separation_shrink_count"] - reference["pair_shrink_count"]),
    }


def build_report(candidate_root: str | Path) -> dict:
    parent, candidate = load_st048s_velocity(candidate_root)
    fixed, fixed_norm = normalized_temporal_piola_velocity(parent, 0.0)
    temporal, temporal_norm = normalized_temporal_piola_velocity(parent, TEMPORAL_GAMMA)
    fixed_measurement = measure_material_paths(fixed)
    temporal_measurement = measure_material_paths(temporal)
    report = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_main_sha": BASE_MAIN_SHA,
        "candidate_source": candidate,
        "agent7_temporal_piola_source": {
            "source_pr": AGENT7_SOURCE_PR,
            "source_head_sha": AGENT7_SOURCE_HEAD,
            "source_task_id": AGENT7_TASK_ID,
            "migration_classification": "direct internal method reuse / independent reimplementation",
            "migration_scope": "beta(t)=.075+.025*(4t-2)^2, exact Piola map, separate positive E(.25)=1 normalization",
            "production_gamma_selected": False,
        },
        "external_method": SCIPY_SOURCE,
        "contract": {
            "time_interval": list(REGISTERED_TIME_INTERVAL),
            "seed_radii": list(SEED_RADII), "seed_z": list(SEED_Z),
            "seed_angles": SEED_ANGLES, "path_count": 48, "paired_material_line_count": 24,
            "output_samples": OUTPUT_SAMPLES, "solver_method": SOLVER_METHOD,
            "solver_rtol": SOLVER_RTOL, "solver_atol": SOLVER_ATOL, "solver_max_step": SOLVER_MAX_STEP,
        },
        "normalization": {"fixed_beta_0075": fixed_norm, "temporal_gamma_0025": temporal_norm},
        "measurements": {"fixed_beta_0075": fixed_measurement, "temporal_gamma_0025": temporal_measurement},
        "comparisons": {
            "temporal_vs_same_run_fixed_beta": _comparison(temporal_measurement, {
                "mean_radius_change": fixed_measurement["mean_radius_change"],
                "mean_absolute_turns": fixed_measurement["mean_absolute_turns"],
                "maximum_absolute_turns": fixed_measurement["maximum_absolute_turns"],
                "mean_pair_separation_change": fixed_measurement["mean_pair_axial_separation_change"],
                "mean_pair_separation_ratio": fixed_measurement["mean_pair_axial_separation_ratio"],
                "pair_growth_count": fixed_measurement["pair_axial_separation_growth_count"],
                "pair_shrink_count": fixed_measurement["pair_axial_separation_shrink_count"],
            }),
            "temporal_vs_unwarped_st048s": _comparison(temporal_measurement, ST048S_REFERENCE),
            "temporal_vs_st006": _comparison(temporal_measurement, ST006_REFERENCE),
        },
        "references": {
            "unwarped_st048s": {"source_head": ST048S_REFERENCE_HEAD, "source_report_sha256": ST048S_REFERENCE_REPORT_SHA256, "measurement": ST048S_REFERENCE},
            "st006": {"measurement": ST006_REFERENCE},
        },
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    report["agent7_scale_replay"] = {
        "reported_temporal_scale": AGENT7_REPORTED_TEMPORAL_SCALE,
        "computed_temporal_scale": temporal_norm["normalization_scale"],
        "absolute_error": abs(temporal_norm["normalization_scale"] - AGENT7_REPORTED_TEMPORAL_SCALE),
    }
    report["report_sha256"] = _canonical_hash(report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    report = build_report(args.candidate_root)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    summary = report["measurements"]["temporal_gamma_0025"]
    print(json.dumps({
        "report_sha256": report["report_sha256"],
        "normalization_scale": report["normalization"]["temporal_gamma_0025"]["normalization_scale"],
        "mean_radius_change": summary["mean_radius_change"],
        "mean_absolute_turns": summary["mean_absolute_turns"],
        "maximum_absolute_turns": summary["maximum_absolute_turns"],
        "mean_pair_axial_separation_change": summary["mean_pair_axial_separation_change"],
        "pair_growth_count": summary["pair_axial_separation_growth_count"],
        "pair_shrink_count": summary["pair_axial_separation_shrink_count"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
