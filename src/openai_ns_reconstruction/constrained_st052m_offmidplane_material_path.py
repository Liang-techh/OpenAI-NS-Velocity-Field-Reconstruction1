"""Replay the frozen ST052-M local-swirl energy child on active off-midplane paths.

CR-A9-063 closes the localization gap left by CR-A9-062. The earlier central
seed set at z=+/-0.3 never entered the shoulder compensation or localized tip
windows. This increment keeps the exact Agent-7 PR #559 child frozen and
introduces one target-free material-path protocol whose seeds begin strictly
inside those active windows.

No coefficient/window scan, image fit, pressure/forcing refit, PDE promotion,
production selection, or visual/path acceptance threshold is introduced.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

TASK_ID = "CR-A9-063"
SCHEMA = "st052m_offmidplane_material_path_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"

AGENT7_PR = 559
AGENT7_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
AGENT7_TASK_ID = "CR003-ST052M-LOCAL-SWIRL-ENERGY-082"
CONTROL_ID = "ST052-M+frozen-redistribution-.05"
CONTROL_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"
PARENT_ID = "ST052-M"
PARENT_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
SOURCE_TAPER_HEAD = "093c7171cd61c6bd439afa30b2da69598a02d182"
SOURCE_FAILED_COMP_HEAD = "62e8c170427d5d830d7f897ba31768e0fc4ce56a"
MORPHOLOGY_SIBLING_PR = 568
MORPHOLOGY_SIBLING_HEAD = "85f2f0940d80757313651d2ab50c6ef89dc22914"

REDISTRIBUTION_ALPHA = 2.520520814687742
REDISTRIBUTION_GAIN = 0.05
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
TAPER_TAU = 0.05
TIP_WINDOW = (0.50, 0.82)
SHOULDER_WINDOW = (0.36, 0.49)
BETA_BRACKET = (0.0, 0.9)
ENERGY_ORDER = 64
EXPECTED_BETA = 0.08837490297155456
BETA_IDENTITY_TOL = 5.0e-10

PATH_ORIGIN_PR = 435
PATH_ORIGIN_HEAD = "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"
REGISTERED_TIME_INTERVAL = (0.25, 0.75)
EVALUATION_BOX = (-2.0, 2.0)
SEED_RADII = (0.6, 0.9, 1.2)
SEED_BANDS = (("shoulder", 0.85), ("tip", 1.15))
SEED_ANGLES = 4
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
    "production_taper_selected": False,
    "production_compensation_selected": False,
    "pressure_or_force_refit_in_this_increment": False,
    "held_out_pde_residual_recomputed_in_this_increment": False,
    "parent_pde_receipt_transferred": False,
    "upstream_eulerian_receipt_promoted_to_pde_truth": False,
    "public_image_used_as_numeric_target": False,
    "hidden_openai_time_camera_seed_or_velocity_used": False,
    "visual_acceptance_threshold_defined": False,
    "trajectory_acceptance_threshold_defined": False,
    "comparison_is_descriptive_not_acceptance": True,
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


def _git_head(root: str | Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(Path(root).resolve()), "rev-parse", "HEAD"], text=True
    ).strip()


def load_agent7_source(root: str | Path):
    root = Path(root).resolve()
    if _git_head(root) != AGENT7_HEAD:
        raise ValueError("Agent-7 checkout is not pinned PR #559 exact head")
    path = root / "experiments" / "root_st052" / "agent7_st052m_local_swirl_energy.py"
    spec = importlib.util.spec_from_file_location("agent7_st052m_local_swirl_energy", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen Agent-7 local swirl energy source")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    expected = {
        "TASK_ID": AGENT7_TASK_ID,
        "PARENT_ID": CONTROL_ID,
        "PARENT_HEAD": CONTROL_HEAD,
        "SOURCE_TAPER_HEAD": SOURCE_TAPER_HEAD,
        "SOURCE_FAILED_COMP_HEAD": SOURCE_FAILED_COMP_HEAD,
        "TAPER_TAU": TAPER_TAU,
        "SHOULDER_WINDOW": SHOULDER_WINDOW,
        "BETA_BRACKET": BETA_BRACKET,
        "ENERGY_ORDER": ENERGY_ORDER,
    }
    for name, value in expected.items():
        if getattr(mod, name) != value:
            raise ValueError(f"frozen Agent-7 identity drift: {name}")
    if mod.prior.WINDOW != TIP_WINDOW:
        raise ValueError("frozen tip-window drift")
    if mod.base.SOURCE_ALPHA != REDISTRIBUTION_ALPHA or mod.base.GAIN != REDISTRIBUTION_GAIN:
        raise ValueError("frozen redistribution coefficient drift")
    if mod.base.INNER_WINDOW != INNER_WINDOW or mod.base.OUTER_WINDOW != OUTER_WINDOW:
        raise ValueError("frozen redistribution radial-window drift")
    if mod.TRUTH["material_paths_integrated"] is not False or mod.TRUTH["pde_validated"] is not False:
        raise ValueError("upstream truth-boundary drift")
    return mod


def verify_path_origin(root: str | Path) -> dict:
    root = Path(root).resolve()
    if _git_head(root) != PATH_ORIGIN_HEAD:
        raise ValueError("Agent-9 path-origin checkout is not pinned PR #435 exact head")
    path = root / "src" / "openai_ns_reconstruction" / "constrained_st048s_piola_swirl_material_path.py"
    spec = importlib.util.spec_from_file_location("agent9_path_origin", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen Agent-9 path-origin module")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    expected = {
        "REGISTERED_TIME_INTERVAL": REGISTERED_TIME_INTERVAL,
        "OUTPUT_SAMPLES": OUTPUT_SAMPLES,
        "SOLVER_METHOD": SOLVER_METHOD,
        "SOLVER_RTOL": SOLVER_RTOL,
        "SOLVER_ATOL": SOLVER_ATOL,
        "SOLVER_MAX_STEP": SOLVER_MAX_STEP,
    }
    for name, value in expected.items():
        if getattr(mod, name) != value:
            raise ValueError(f"frozen Agent-9 solver contract drift: {name}")
    return {
        "source_pr": PATH_ORIGIN_PR,
        "source_head": PATH_ORIGIN_HEAD,
        "reused_exactly": ["time_interval", "output_samples", "solver_method", "rtol", "atol", "max_step"],
        "changed_autonomously": ["seed_z_bands", "azimuth_count"],
    }


def _offmidplane_seed_table() -> tuple[np.ndarray, list[dict]]:
    rows: list[list[float]] = []
    metadata: list[dict] = []
    for radius in SEED_RADII:
        for angle_index in range(SEED_ANGLES):
            angle = 2.0 * np.pi * angle_index / SEED_ANGLES
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            for band, z_abs in SEED_BANDS:
                for sign in (-1, 1):
                    z = sign * z_abs
                    rows.append([x, y, z])
                    metadata.append({
                        "radius": float(radius),
                        "angle_index": int(angle_index),
                        "angle_radians": float(angle),
                        "band": band,
                        "z_abs": float(z_abs),
                        "sign": int(sign),
                        "z": float(z),
                    })
    seeds = np.asarray(rows, dtype=float)
    if seeds.shape != (48, 3) or len(metadata) != 48:
        raise RuntimeError("off-midplane seed construction drift")
    return seeds, metadata


def _evaluate_velocity(velocity: VelocityCallable, points: np.ndarray, time: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != points.shape or points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("velocity must return shape (n,3)")
    if not np.isfinite(values).all():
        raise ValueError("velocity returned non-finite values")
    return values


def _window_fraction(abs_z_over_2: np.ndarray, window: tuple[float, float]) -> np.ndarray:
    return np.mean((abs_z_over_2 > window[0]) & (abs_z_over_2 < window[1]), axis=0)


def _summarize_band(
    band: str,
    metadata: list[dict],
    radius_change: np.ndarray,
    abs_turns: np.ndarray,
    axial_abs_change: np.ndarray,
    path_lengths: np.ndarray,
    shoulder_fraction: np.ndarray,
    tip_fraction: np.ndarray,
    pair_rows: list[dict],
) -> dict:
    indices = [i for i, row in enumerate(metadata) if row["band"] == band]
    pairs = [row for row in pair_rows if row["band"] == band]
    pair_changes = np.asarray([row["axial_separation_change"] for row in pairs], dtype=float)
    return {
        "path_count": len(indices),
        "paired_material_line_count": len(pairs),
        "mean_radius_change": float(np.mean(radius_change[indices])),
        "inward_path_count": int(np.sum(radius_change[indices] < 0.0)),
        "mean_absolute_turns": float(np.mean(abs_turns[indices])),
        "maximum_absolute_turns": float(np.max(abs_turns[indices])),
        "mean_abs_z_change": float(np.mean(axial_abs_change[indices])),
        "mean_pair_axial_separation_change": float(np.mean(pair_changes)),
        "pair_axial_separation_growth_count": int(np.sum(pair_changes > 0.0)),
        "pair_axial_separation_shrink_count": int(np.sum(pair_changes < 0.0)),
        "mean_sampled_path_length": float(np.mean(path_lengths[indices])),
        "mean_shoulder_window_sample_fraction": float(np.mean(shoulder_fraction[indices])),
        "mean_tip_window_sample_fraction": float(np.mean(tip_fraction[indices])),
    }


def measure_offmidplane_paths(velocity: VelocityCallable) -> dict:
    seeds, metadata = _offmidplane_seed_table()
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
        raise ValueError("sampled off-midplane material path left registered evaluation box")

    velocities = np.empty_like(positions)
    for time_index, time in enumerate(times):
        velocities[time_index] = _evaluate_velocity(velocity, positions[time_index], float(time))
    speeds = np.linalg.norm(velocities, axis=2)
    radii = np.hypot(positions[:, :, 0], positions[:, :, 1])
    angles = np.unwrap(np.arctan2(positions[:, :, 1], positions[:, :, 0]), axis=0)
    abs_z = np.abs(positions[:, :, 2])
    abs_z_over_2 = abs_z / 2.0
    radius_change = radii[-1] - radii[0]
    axial_abs_change = abs_z[-1] - abs_z[0]
    abs_turns = np.abs((angles[-1] - angles[0]) / (2.0 * np.pi))
    speed_change = speeds[-1] - speeds[0]
    path_lengths = np.sum(np.linalg.norm(np.diff(positions, axis=0), axis=2), axis=0)
    shoulder_fraction = _window_fraction(abs_z_over_2, SHOULDER_WINDOW)
    tip_fraction = _window_fraction(abs_z_over_2, TIP_WINDOW)

    lookup: dict[tuple[float, int, str], dict[int, int]] = {}
    for index, row in enumerate(metadata):
        key = (row["radius"], row["angle_index"], row["band"])
        lookup.setdefault(key, {})[row["sign"]] = index
    pair_rows: list[dict] = []
    for (radius, angle_index, band), signs in lookup.items():
        if set(signs) != {-1, 1}:
            raise RuntimeError("off-midplane axial pair construction drift")
        neg = signs[-1]
        pos = signs[1]
        separation = positions[:, pos, 2] - positions[:, neg, 2]
        pair_rows.append({
            "radius": float(radius),
            "angle_index": int(angle_index),
            "band": band,
            "initial_axial_separation": float(separation[0]),
            "final_axial_separation": float(separation[-1]),
            "axial_separation_change": float(separation[-1] - separation[0]),
            "axial_separation_ratio": float(separation[-1] / separation[0]),
        })
    if len(pair_rows) != 24:
        raise RuntimeError("off-midplane pair population drift")
    pair_changes = np.asarray([row["axial_separation_change"] for row in pair_rows], dtype=float)

    result = {
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
        "mean_pair_axial_separation_change": float(np.mean(pair_changes)),
        "pair_axial_separation_growth_count": int(np.sum(pair_changes > 0.0)),
        "pair_axial_separation_shrink_count": int(np.sum(pair_changes < 0.0)),
        "mean_speed_change": float(np.mean(speed_change)),
        "mean_sampled_path_length": float(np.mean(path_lengths)),
        "mean_shoulder_window_sample_fraction": float(np.mean(shoulder_fraction)),
        "mean_tip_window_sample_fraction": float(np.mean(tip_fraction)),
        "solver": {
            "method": SOLVER_METHOD, "rtol": SOLVER_RTOL, "atol": SOLVER_ATOL,
            "max_step": SOLVER_MAX_STEP, "output_samples": OUTPUT_SAMPLES,
        },
        "per_path": [
            {
                "path_index": int(i), "seed": metadata[i],
                "radius_change": float(radius_change[i]),
                "absolute_turns": float(abs_turns[i]),
                "abs_z_change": float(axial_abs_change[i]),
                "speed_change": float(speed_change[i]),
                "sampled_path_length": float(path_lengths[i]),
                "shoulder_window_sample_fraction": float(shoulder_fraction[i]),
                "tip_window_sample_fraction": float(tip_fraction[i]),
            }
            for i in range(n_paths)
        ],
        "pair_rows": pair_rows,
    }
    result["by_seed_band"] = {
        band: _summarize_band(
            band, metadata, radius_change, abs_turns, axial_abs_change, path_lengths,
            shoulder_fraction, tip_fraction, pair_rows,
        )
        for band, _z in SEED_BANDS
    }
    return result


def _relative(new: float, old: float) -> float | None:
    if abs(old) <= 1.0e-14:
        return None
    return float(new / old - 1.0)


def compare_measurements(child: dict, control: dict) -> dict:
    return {
        "mean_absolute_turns_relative": _relative(child["mean_absolute_turns"], control["mean_absolute_turns"]),
        "maximum_absolute_turns_relative": _relative(child["maximum_absolute_turns"], control["maximum_absolute_turns"]),
        "radial_contraction_magnitude_relative": _relative(abs(child["mean_radius_change"]), abs(control["mean_radius_change"])),
        "mean_pair_axial_separation_change_relative": _relative(
            child["mean_pair_axial_separation_change"], control["mean_pair_axial_separation_change"]
        ),
        "mean_abs_z_change_delta": float(child["mean_abs_z_change"] - control["mean_abs_z_change"]),
        "inward_path_count_delta": int(child["inward_path_count"] - control["inward_path_count"]),
        "pair_growth_count_delta": int(
            child["pair_axial_separation_growth_count"] - control["pair_axial_separation_growth_count"]
        ),
        "pair_shrink_count_delta": int(
            child["pair_axial_separation_shrink_count"] - control["pair_axial_separation_shrink_count"]
        ),
    }


def _seed_response_by_band(control_velocity: VelocityCallable, child_velocity: VelocityCallable) -> dict:
    seeds, metadata = _offmidplane_seed_table()
    control = _evaluate_velocity(control_velocity, seeds, REGISTERED_TIME_INTERVAL[0])
    child = _evaluate_velocity(child_velocity, seeds, REGISTERED_TIME_INTERVAL[0])
    delta = child - control
    result = {}
    for band, _z in SEED_BANDS:
        idx = [i for i, row in enumerate(metadata) if row["band"] == band]
        c = control[idx]
        d = delta[idx]
        control_rms = float(np.sqrt(np.mean(np.sum(c * c, axis=1))))
        delta_rms = float(np.sqrt(np.mean(np.sum(d * d, axis=1))))
        result[band] = {
            "path_count": len(idx),
            "control_velocity_rms": control_rms,
            "child_minus_control_velocity_rms": delta_rms,
            "relative_delta_rms": float(delta_rms / control_rms),
            "max_pointwise_delta_norm": float(np.max(np.linalg.norm(d, axis=1))),
        }
    return result


def build_report(agent7_root: str | Path, path_origin_root: str | Path) -> dict:
    source = load_agent7_source(agent7_root)
    path_origin = verify_path_origin(path_origin_root)
    family, raw = source.replay_st052.reconstruct()
    redistribution_scale, parent_energy, redistribution_raw_energy = source.base.child_scale(family, raw)
    energy_solve = source.solve_energy_beta(family, raw, redistribution_scale)
    if not energy_solve["root_exists"]:
        raise ValueError("frozen local-swirl energy root disappeared")
    beta = float(energy_solve["beta"])
    if abs(beta - EXPECTED_BETA) > BETA_IDENTITY_TOL:
        raise ValueError("Agent-7 local-swirl beta identity drift")
    if abs(float(energy_solve["post_transform_common_scale"]) - 1.0) > 1.0e-15:
        raise ValueError("post-transform common scale is forbidden by frozen contract")

    def control_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return source.prior.control_velocity(
            family, raw, np.asarray(points, dtype=float), float(time), redistribution_scale
        )

    def child_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return source.compensated_velocity(
            family, raw, np.asarray(points, dtype=float), float(time),
            redistribution_scale=redistribution_scale, beta=beta,
        )

    control = measure_offmidplane_paths(control_velocity)
    child = measure_offmidplane_paths(child_velocity)
    by_band_comparison = {
        band: compare_measurements(child["by_seed_band"][band], control["by_seed_band"][band])
        for band, _z in SEED_BANDS
    }

    report = {
        "task_id": TASK_ID,
        "schema": SCHEMA,
        "base_main_sha": BASE_MAIN_SHA,
        "lineage": {
            "st052_parent_candidate": PARENT_ID,
            "st052_parent_head": PARENT_HEAD,
            "control_candidate": CONTROL_ID,
            "control_head": CONTROL_HEAD,
            "agent7_source_pr": AGENT7_PR,
            "agent7_source_head": AGENT7_HEAD,
            "agent7_task_id": AGENT7_TASK_ID,
            "source_taper_head": SOURCE_TAPER_HEAD,
            "source_failed_compensation_head": SOURCE_FAILED_COMP_HEAD,
            "morphology_sibling_pr": MORPHOLOGY_SIBLING_PR,
            "morphology_sibling_head": MORPHOLOGY_SIBLING_HEAD,
            "path_origin_pr": PATH_ORIGIN_PR,
            "path_origin_head": PATH_ORIGIN_HEAD,
        },
        "frozen_transform": {
            "redistribution_alpha": REDISTRIBUTION_ALPHA,
            "redistribution_gain": REDISTRIBUTION_GAIN,
            "redistribution_inner_window": list(INNER_WINDOW),
            "redistribution_outer_window": list(OUTER_WINDOW),
            "redistribution_normalization": float(redistribution_scale),
            "redistribution_parent_reference_energy": float(parent_energy),
            "redistribution_unscaled_reference_energy": float(redistribution_raw_energy),
            "taper_tau": TAPER_TAU,
            "tip_window_abs_z_over_2": list(TIP_WINDOW),
            "shoulder_window_abs_z_over_2": list(SHOULDER_WINDOW),
            "beta_bracket": list(BETA_BRACKET),
            "beta": beta,
            "energy_order": ENERGY_ORDER,
            "control_reference_energy": float(energy_solve["control_energy"]),
            "child_reference_energy": float(energy_solve["child_energy"]),
            "relative_reference_energy_error": float(energy_solve["relative_energy_error"]),
            "post_transform_common_scale": float(energy_solve["post_transform_common_scale"]),
        },
        "frozen_offmidplane_path_contract": {
            "time_interval": list(REGISTERED_TIME_INTERVAL),
            "seed_radii": list(SEED_RADII),
            "seed_bands_abs_z": {band: z for band, z in SEED_BANDS},
            "seed_angles": SEED_ANGLES,
            "path_count": 48,
            "paired_material_line_count": 24,
            "output_samples": OUTPUT_SAMPLES,
            "solver": SOLVER_METHOD,
            "rtol": SOLVER_RTOL,
            "atol": SOLVER_ATOL,
            "max_step": SOLVER_MAX_STEP,
            "shoulder_seed_abs_z_over_2": 0.85 / 2.0,
            "tip_seed_abs_z_over_2": 1.15 / 2.0,
            "scientific_acceptance_threshold_defined": False,
        },
        "path_origin_verification": path_origin,
        "seed_response_by_band": _seed_response_by_band(control_velocity, child_velocity),
        "redistributed_control_material_paths": control,
        "local_swirl_energy_child_material_paths": child,
        "child_vs_control": compare_measurements(child, control),
        "child_vs_control_by_seed_band": by_band_comparison,
        "external_method": dict(SCIPY_SOURCE),
        "internal_method": {
            "repository": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
            "source_pr": AGENT7_PR,
            "source_commit": AGENT7_HEAD,
            "license": "same-repository internal source; no external code copied",
            "classification": "direct internal method reuse / independent off-midplane material-path replay",
            "migrated_scope": "exact PR #559 compensated velocity child; solver/time controls inherited from Agent-9 PR #435",
            "difference": "new autonomous seeds begin in the frozen shoulder and tip active regions; no source coefficient or PDE state is changed",
        },
        "truth_boundary": dict(TRUTH_BOUNDARY),
        "direct_contribution_to_final_velocity": (
            "tests whether the exact energy-neutral localized morphology child changes real trajectories where its shoulder/tip transforms are actually active, closing the blind spot of the prior central-path replay"
        ),
    }
    report["report_sha256"] = _canonical_hash(report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--agent7-root", required=True)
    parser.add_argument("--path-origin-root", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(args.agent7_root, args.path_origin_root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
