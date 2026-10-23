"""Replay Agent-7's first swirl-gain crossing on the frozen Agent-9 material paths.

This is a visualization-side governance experiment. It reconstructs frozen ST048-S
from PR #390, applies the exact temporal Piola schedule used by Agent 7 and one
axisymmetric azimuthal gain kappa=.05 from PR #428, restores E(.25)=1 with one
positive common scale, then compares real material trajectories with the same
temporal-Piola field at kappa=0.

No pressure/forcing/PDE evidence is transferred. No public-image target, hidden
OpenAI time/camera/seed, visual pass threshold, or production coefficient is used.
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

TASK_ID = "CR-A9-048"
SCHEMA = "st048s_piola_swirl_material_path_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"

ST048_SOURCE_PR = 390
ST048_SOURCE_HEAD = "97695a86f85ce68fb4ae70c41fc81c904d655183"
ST048_PARENT_ID = "ST047-E"
ST048_PARENT_RAW_SHA256 = "dfe6e51af93c9c42f1322798b82f89523e6a193c1b6495894cb60f63b4eb07e0"
ST048S_RAW_SHA256 = "6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09"

AGENT7_SOURCE_PR = 428
AGENT7_SOURCE_HEAD = "b2b4888ece83a9f862435ce1adb847c94ad3ad7d"
AGENT7_TASK_ID = "CR003-ST048S-TEMPORAL-PIOLA-SWIRL-GAIN-063"
BASE_BETA = 0.075
BASE_GAMMA = 0.025
SWIRL_KAPPA = 0.05
AGENT7_REPORTED_KAPPA_SCALE = 0.9726627718

PHYSICAL_Z_CUT = 2.0
REFERENCE_ENERGY = 1.0
ENERGY_ORDER = 64

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

ST048S_REFERENCE = {
    "candidate": "ST048-S",
    "mean_radius_change": -0.0853545661,
    "mean_absolute_turns": 0.0221757297,
    "maximum_absolute_turns": 0.0346689104,
    "mean_pair_separation_change": 0.0778727503,
    "mean_pair_separation_ratio": 1.1297879172,
    "pair_growth_count": 16,
    "pair_shrink_count": 8,
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

SCIPY_SOURCE = {
    "repository": "scipy/scipy",
    "commit": "eff78058ed0d4cb8e1f0b5c5585d62d712e948d2",
    "api": "scipy.integrate.solve_ivp / DOP853",
    "license": "BSD-3-Clause",
    "classification": "direct migration / public API only",
    "copied_upstream_implementation": False,
}
INTERNAL_SOURCE = {
    "repository": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
    "pr": AGENT7_SOURCE_PR,
    "commit": AGENT7_SOURCE_HEAD,
    "task_id": AGENT7_TASK_ID,
    "classification": "direct internal method reuse / independent reimplementation",
    "migrated_scope": "beta(t), Piola map, cylindrical swirl multiplier, common energy normalization",
}
TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "production_candidate_selected": False,
    "production_kappa_selected": False,
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


def _git_head(repository_root: str | Path) -> str:
    completed = subprocess.run(
        ["git", "-C", str(Path(repository_root).resolve()), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    )
    return completed.stdout.strip()


def beta_at_time(time: float, gamma: float = BASE_GAMMA, *, base_beta: float = BASE_BETA) -> float:
    time = float(time)
    gamma = float(gamma)
    base_beta = float(base_beta)
    if not REGISTERED_TIME_INTERVAL[0] <= time <= REGISTERED_TIME_INTERVAL[1]:
        raise ValueError("time must lie in [0.25,0.75]")
    if not 0.0 <= gamma <= 0.125 or not 0.0 <= base_beta <= 0.20:
        raise ValueError("Piola coefficients outside screened bounds")
    beta = base_beta + gamma * (4.0 * time - 2.0) ** 2
    if not 0.0 <= beta <= 0.20 + 1.0e-15:
        raise ValueError("temporal schedule exceeds screened beta range")
    return float(beta)


def _warp_z_and_jacobian(z, beta: float):
    z = np.asarray(z, dtype=float)
    beta = float(beta)
    if not np.all(np.isfinite(z)) or not 0.0 <= beta <= 0.20 + 1.0e-15:
        raise ValueError("invalid Piola input")
    s = z / PHYSICAL_Z_CUT
    inside = np.abs(s) < 1.0
    one_minus = np.maximum(1.0 - s * s, 0.0)
    bump = one_minus**4
    mapped = np.where(inside, z * (1.0 - beta * bump), z)
    jac_inside = 1.0 - beta * bump + 8.0 * beta * s * s * one_minus**3
    jac = np.where(inside, jac_inside, 1.0)
    if np.any(jac <= 0.0) or not np.all(np.isfinite(mapped)):
        raise RuntimeError("Piola map lost orientation")
    return mapped, jac


def _apply_swirl_gain_xy(points: np.ndarray, velocity: np.ndarray, kappa: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    kappa = float(kappa)
    if points.shape != velocity.shape or points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points and velocity must both have shape (n,3)")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(velocity)):
        raise ValueError("non-finite swirl input")
    if not 0.0 <= kappa <= 0.15:
        raise ValueError("kappa outside preregistered Agent-7 screen")

    x = points[:, 0]
    y = points[:, 1]
    radius = np.hypot(x, y)
    out = np.array(velocity, copy=True)
    mask = radius > 1.0e-14
    if np.any(mask):
        rx = x[mask] / radius[mask]
        ry = y[mask] / radius[mask]
        ux = velocity[mask, 0]
        uy = velocity[mask, 1]
        radial = rx * ux + ry * uy
        azimuthal = -ry * ux + rx * uy
        azimuthal *= 1.0 + kappa
        out[mask, 0] = rx * radial - ry * azimuthal
        out[mask, 1] = ry * radial + rx * azimuthal
    return out


class TemporalPiolaSwirlVelocity:
    """Temporal contravariant Piola field plus one constant cylindrical swirl gain."""

    def __init__(self, parent: VelocityCallable, *, kappa: float, scale: float = 1.0):
        if not callable(parent):
            raise TypeError("parent must be callable")
        self.parent = parent
        self.kappa = float(kappa)
        self.scale = float(scale)
        if not 0.0 <= self.kappa <= 0.15:
            raise ValueError("kappa outside preregistered Agent-7 screen")
        if not np.isfinite(self.scale) or self.scale <= 0.0:
            raise ValueError("scale must be positive and finite")

    def __call__(self, points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
            raise ValueError("points must be finite shape (n,3)")
        beta = beta_at_time(float(time))
        mapped_z, jac = _warp_z_and_jacobian(points[:, 2], beta)
        mapped = np.array(points, copy=True)
        mapped[:, 2] = mapped_z
        parent_velocity = np.asarray(self.parent(mapped, float(time)), dtype=float)
        if parent_velocity.shape != points.shape or not np.all(np.isfinite(parent_velocity)):
            raise ValueError("parent velocity returned invalid values")
        piola = np.array(parent_velocity, copy=True)
        piola[:, 0] *= jac
        piola[:, 1] *= jac
        gained = _apply_swirl_gain_xy(points, piola, self.kappa)
        return self.scale * gained


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

    return velocity, {
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


def normalized_field(parent: VelocityCallable, *, kappa: float) -> tuple[TemporalPiolaSwirlVelocity, dict]:
    raw = TemporalPiolaSwirlVelocity(parent, kappa=kappa, scale=1.0)
    raw_energy = axisymmetric_energy(raw, time=0.25, order=ENERGY_ORDER)
    if not np.isfinite(raw_energy) or raw_energy <= 0.0:
        raise ValueError("raw reference energy must be positive")
    scale = float(np.sqrt(REFERENCE_ENERGY / raw_energy))
    field = TemporalPiolaSwirlVelocity(parent, kappa=kappa, scale=scale)
    normalized_energy = axisymmetric_energy(field, time=0.25, order=ENERGY_ORDER)
    return field, {
        "base_beta": BASE_BETA,
        "gamma": BASE_GAMMA,
        "kappa": float(kappa),
        "beta_at_025": beta_at_time(0.25),
        "beta_at_050": beta_at_time(0.50),
        "beta_at_075": beta_at_time(0.75),
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
                metadata.append({
                    "radius": float(radius), "angle_index": int(angle_index),
                    "angle_radians": float(angle), "z": float(z),
                })
    return np.asarray(rows, dtype=float), metadata


def _evaluate_velocity(velocity: VelocityCallable, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != points.shape or not np.isfinite(values).all():
        raise ValueError("velocity must return finite shape (n,3)")
    return values


def measure_material_paths(velocity: VelocityCallable) -> dict:
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
        raise ValueError("sampled material path left registered evaluation box")

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
    abs_turns = np.abs(turns)
    speed_change = speeds[-1] - speeds[0]
    segment_lengths = np.linalg.norm(np.diff(positions, axis=0), axis=2)
    path_lengths = np.sum(segment_lengths, axis=0)

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

    return {
        "path_count": int(n_paths),
        "paired_material_line_count": int(len(pair_rows)),
        "initial_velocity_rms": initial_velocity_rms,
        "mean_radius_change": float(np.mean(radius_change)),
        "median_radius_change": float(np.median(radius_change)),
        "inward_path_count": int(np.sum(radius_change < 0.0)),
        "outward_path_count": int(np.sum(radius_change > 0.0)),
        "mean_absolute_turns": float(np.mean(abs_turns)),
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
        "solver": {
            "method": SOLVER_METHOD, "rtol": SOLVER_RTOL, "atol": SOLVER_ATOL,
            "max_step": SOLVER_MAX_STEP, "output_samples": OUTPUT_SAMPLES,
        },
        "per_path": [
            {
                "path_index": int(i), "seed": seed_metadata[i],
                "radius_change": float(radius_change[i]),
                "absolute_turns": float(abs_turns[i]),
                "abs_z_change": float(axial_abs_change[i]),
                "speed_change": float(speed_change[i]),
                "sampled_path_length": float(path_lengths[i]),
            }
            for i in range(n_paths)
        ],
        "pair_rows": pair_rows,
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


def crosscheck_agent7_formula(agent7_root: str | Path) -> dict:
    root = Path(agent7_root).resolve()
    if _git_head(root) != AGENT7_SOURCE_HEAD:
        raise ValueError("Agent-7 checkout is not the pinned #428 head")
    experiment = root / "experiments" / "root_st048"
    previous_path = list(sys.path)
    module_names = (
        "agent7_st048s_piola_swirl_gain_screen",
        "agent7_st048s_piola_temporal_screen",
        "agent7_st048s_piola_warp_screen",
        "agent7_st048_morphology_transfer",
    )
    try:
        sys.path.insert(0, str(experiment))
        for name in module_names:
            sys.modules.pop(name, None)
        screen = importlib.import_module("agent7_st048s_piola_swirl_gain_screen")
    finally:
        sys.path[:] = previous_path
    if (
        getattr(screen, "TASK_ID", None) != AGENT7_TASK_ID
        or abs(float(screen.BASE_BETA) - BASE_BETA) > 0.0
        or abs(float(screen.BASE_GAMMA) - BASE_GAMMA) > 0.0
        or SWIRL_KAPPA not in tuple(float(x) for x in screen.DEFAULT_KAPPAS)
    ):
        raise ValueError("Agent-7 source constants drift")

    class SyntheticParent:
        def at_points(self, points, time):
            points = np.asarray(points, dtype=float)
            t = float(time)
            return np.column_stack((
                0.11 + 0.07 * points[:, 0] - 0.03 * points[:, 2] + 0.02 * t,
                -0.09 + 0.05 * points[:, 1] + 0.04 * points[:, 2] - 0.01 * t,
                0.06 + 0.02 * points[:, 0] - 0.01 * points[:, 1] + 0.03 * t,
            ))

    parent_obj = SyntheticParent()
    parent_callable = lambda pts, t: parent_obj.at_points(pts, t)
    ours = TemporalPiolaSwirlVelocity(parent_callable, kappa=SWIRL_KAPPA, scale=1.0)
    theirs = screen._raw_field(parent_obj, base_beta=BASE_BETA, gamma=BASE_GAMMA, kappa=SWIRL_KAPPA)
    points = np.asarray([
        [0.6, 0.0, -0.3], [0.0, 0.9, 0.3], [-0.8, 0.8, 1.1],
        [1.2, -0.4, -1.4], [0.0, 0.0, 0.2],
    ], dtype=float)
    errors = []
    for time in (0.25, 0.5, 0.75):
        expected = np.asarray(theirs.at_points(points, time), dtype=float)
        actual = ours(points, time)
        errors.append(float(np.max(np.abs(actual - expected))))
    worst = max(errors)
    if worst > 2.0e-15:
        raise RuntimeError(f"independent Agent-7 formula cross-check failed: {worst}")
    return {
        "source_pr": AGENT7_SOURCE_PR,
        "source_head_sha": AGENT7_SOURCE_HEAD,
        "task_id": AGENT7_TASK_ID,
        "kappa": SWIRL_KAPPA,
        "tested_times": [0.25, 0.5, 0.75],
        "max_abs_formula_difference": worst,
        "crosscheck_passed": True,
    }


def build_report(candidate_root: str | Path, *, agent7_root: str | Path | None = None) -> dict:
    parent, candidate_metadata = load_st048s_velocity(candidate_root)
    baseline, baseline_norm = normalized_field(parent, kappa=0.0)
    child, child_norm = normalized_field(parent, kappa=SWIRL_KAPPA)
    baseline_measurement = measure_material_paths(baseline)
    child_measurement = measure_material_paths(child)
    agent7_crosscheck = crosscheck_agent7_formula(agent7_root) if agent7_root is not None else {
        "source_pr": AGENT7_SOURCE_PR,
        "source_head_sha": AGENT7_SOURCE_HEAD,
        "crosscheck_passed": False,
        "reason": "agent7_root_not_supplied",
    }

    report = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_main_sha": BASE_MAIN_SHA,
        "candidate": candidate_metadata,
        "internal_source": INTERNAL_SOURCE,
        "external_source": SCIPY_SOURCE,
        "frozen_contract": {
            "time_interval": list(REGISTERED_TIME_INTERVAL),
            "seed_radii": list(SEED_RADII),
            "seed_z": list(SEED_Z),
            "seed_angles": SEED_ANGLES,
            "path_count": 48,
            "paired_material_line_count": 24,
            "output_samples": OUTPUT_SAMPLES,
            "solver_method": SOLVER_METHOD,
            "solver_rtol": SOLVER_RTOL,
            "solver_atol": SOLVER_ATOL,
            "solver_max_step": SOLVER_MAX_STEP,
        },
        "normalization": {
            "temporal_piola_kappa_0": baseline_norm,
            "temporal_piola_kappa_005": child_norm,
            "agent7_reported_kappa_005_scale_rounded": AGENT7_REPORTED_KAPPA_SCALE,
        },
        "agent7_formula_crosscheck": agent7_crosscheck,
        "measurements": {
            "temporal_piola_kappa_0": baseline_measurement,
            "temporal_piola_kappa_005": child_measurement,
        },
        "comparisons": {
            "kappa_005_vs_same_run_kappa_0": _comparison(child_measurement, {
                "mean_radius_change": baseline_measurement["mean_radius_change"],
                "mean_absolute_turns": baseline_measurement["mean_absolute_turns"],
                "maximum_absolute_turns": baseline_measurement["maximum_absolute_turns"],
                "mean_pair_separation_change": baseline_measurement["mean_pair_axial_separation_change"],
                "mean_pair_separation_ratio": baseline_measurement["mean_pair_axial_separation_ratio"],
                "pair_growth_count": baseline_measurement["pair_axial_separation_growth_count"],
                "pair_shrink_count": baseline_measurement["pair_axial_separation_shrink_count"],
            }),
            "kappa_005_vs_unwarped_st048s": _comparison(child_measurement, ST048S_REFERENCE),
            "kappa_005_vs_st006": _comparison(child_measurement, ST006_REFERENCE),
        },
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    report["report_sha256"] = _canonical_hash(report)
    return report


def audit_truth_boundary(report: dict) -> None:
    if report.get("task_id") != TASK_ID:
        raise ValueError("task identity drift")
    truth = report.get("truth_boundary")
    if truth != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drift")
    if report["candidate"].get("pde_validated") is not False:
        raise ValueError("candidate PDE state changed")
    if report["agent7_formula_crosscheck"].get("crosscheck_passed") is not True:
        raise ValueError("exact Agent-7 formula cross-check is required")
    for key in (
        "canonical_velocity_changed", "production_candidate_selected",
        "production_kappa_selected", "pressure_or_force_transferred_from_parent",
        "held_out_pde_residual_evaluated", "held_out_pde_residual_transferred_from_parent",
        "visual_acceptance_threshold_defined", "hidden_openai_time_camera_seed_or_velocity_used",
        "visualization_ready", "visual_correspondence_verified", "pde_validated",
        "source_correspondence_verified", "paper_exact", "openai_field_identified", "blowup_proved",
    ):
        if truth[key] is not False:
            raise ValueError(f"forbidden truth promotion: {key}")
    if truth["comparison_is_descriptive_not_acceptance"] is not True:
        raise ValueError("comparison must stay descriptive")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-root", required=True)
    parser.add_argument("--agent7-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    report = build_report(args.candidate_root, agent7_root=args.agent7_root)
    audit_truth_boundary(report)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    summary = report["measurements"]["temporal_piola_kappa_005"]
    baseline = report["measurements"]["temporal_piola_kappa_0"]
    print(json.dumps({
        "task_id": TASK_ID,
        "report_sha256": report["report_sha256"],
        "normalization_scale": report["normalization"]["temporal_piola_kappa_005"]["normalization_scale"],
        "baseline_mean_absolute_turns": baseline["mean_absolute_turns"],
        "child_mean_absolute_turns": summary["mean_absolute_turns"],
        "baseline_maximum_absolute_turns": baseline["maximum_absolute_turns"],
        "child_maximum_absolute_turns": summary["maximum_absolute_turns"],
        "child_mean_radius_change": summary["mean_radius_change"],
        "child_mean_pair_axial_separation_change": summary["mean_pair_axial_separation_change"],
        "pair_growth_count": summary["pair_axial_separation_growth_count"],
        "pair_shrink_count": summary["pair_axial_separation_shrink_count"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
