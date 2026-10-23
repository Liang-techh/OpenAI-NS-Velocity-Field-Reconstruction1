"""Replay Agent-7's callable ST051-B Piola child on frozen Agent-9 paths.

CR-A9-057 consumes exactly the callable ST051-B `.025` redistribution +
`beta=.075` axial Piola field from Agent 7 PR #510 and compares it with the
same callable `.025` child before Piola under Agent 9's unchanged 48-path
material-trajectory protocol.

This is visualization-side routing evidence only. No pressure/forcing/PDE
receipt is transferred and no OpenAI hidden numerical target is used.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from typing import Callable

import numpy as np

TASK_ID = "CR-A9-057"
SCHEMA = "st051b_callable_piola_material_path_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"

PARENT_ID = "ST051-B"
PARENT_PR = 460
PARENT_HEAD = "4b784f1b8457af2ead49295631d834d4e882000b"
AGENT7_PR = 510
AGENT7_HEAD = "c87ffa3f1a798df0429f16f5ee12dfb41d6c8ded"
AGENT7_TASK_ID = "CR003-ST051B-CALLABLE-PIOLA-TRANSFER-076"
AGENT9_PATH_PR = 435
AGENT9_PATH_HEAD = "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"

REDISTRIBUTION_GAIN = 0.025
REDISTRIBUTION_ALPHA = 2.520520814687742
PIOLA_BETA = 0.075
SEED_RADII = (0.6, 0.9, 1.2)

# CR-A9-055 same-backbone `.025` aggregate control. This is replayed using
# the previously audited redistribution scale, not the slightly different
# PR #510 quadrature replay scale. The original tight identity tolerance is
# intentionally preserved rather than relaxed after seeing a new result.
BASE_REFERENCE = {
    "mean_absolute_turns": 0.02229297191872165,
    "maximum_absolute_turns": 0.027371569482050228,
    "radial_contraction_magnitude": 0.07049762300296604,
    "mean_pair_axial_separation_change": 0.06509150801226499,
}
BASE_REFERENCE_ATOL = 2.5e-8

SCIPY_SOURCE = {
    "repository": "scipy/scipy",
    "commit": "eff78058ed0d4cb8e1f0b5c5585d62d712e948d2",
    "api": "scipy.integrate.solve_ivp / DOP853",
    "license": "BSD-3-Clause",
    "classification": "direct migration / public API only through frozen Agent-9 path engine",
    "copied_upstream_implementation": False,
}

TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "production_beta_selected": False,
    "production_candidate_selected": False,
    "pressure_or_force_transferred_from_parent": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_used_as_numeric_target": False,
    "hidden_openai_time_camera_seed_or_velocity_used": False,
    "visual_acceptance_threshold_defined": False,
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


def _git_head(root: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()


def _canonical_hash(payload: dict) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def load_agent7(root: str | Path):
    root = Path(root).resolve()
    if _git_head(root) != AGENT7_HEAD:
        raise ValueError("Agent-7 checkout is not pinned PR #510 exact head")
    exp = root / "experiments" / "root_st051"
    old_path = list(sys.path)
    names = (
        "agent7_st051b_callable_piola_transfer",
        "replay_st051",
        "aligned_continuation",
        "replay_recovery",
    )
    try:
        sys.path.insert(0, str(exp))
        for name in names:
            sys.modules.pop(name, None)
        mod = importlib.import_module("agent7_st051b_callable_piola_transfer")
    finally:
        sys.path[:] = old_path
    if getattr(mod, "TASK_ID", None) != AGENT7_TASK_ID:
        raise ValueError("Agent-7 task identity drift")
    if abs(float(mod.SOURCE_REDISTRIBUTION_GAIN) - REDISTRIBUTION_GAIN) > 0.0:
        raise ValueError("redistribution gain drift")
    if abs(float(mod.SOURCE_REDISTRIBUTION_ALPHA) - REDISTRIBUTION_ALPHA) > 0.0:
        raise ValueError("redistribution alpha drift")
    if abs(float(mod.PIOLA_BETA) - PIOLA_BETA) > 0.0:
        raise ValueError("Piola beta drift")
    return mod


def load_path_engine(root: str | Path):
    root = Path(root).resolve()
    if _git_head(root) != AGENT9_PATH_HEAD:
        raise ValueError("Agent-9 path checkout is not pinned PR #435 exact head")
    path = root / "src" / "openai_ns_reconstruction" / "constrained_st048s_piola_swirl_material_path.py"
    spec = importlib.util.spec_from_file_location("agent9_frozen_path_engine", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen Agent-9 path engine")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    expected = {
        "REGISTERED_TIME_INTERVAL": (0.25, 0.75),
        "SEED_RADII": SEED_RADII,
        "SEED_Z": (-0.3, 0.3),
        "SEED_ANGLES": 8,
        "OUTPUT_SAMPLES": 33,
        "SOLVER_METHOD": "DOP853",
        "SOLVER_RTOL": 1.0e-9,
        "SOLVER_ATOL": 1.0e-11,
        "SOLVER_MAX_STEP": 0.01,
    }
    for name, value in expected.items():
        actual = getattr(mod, name)
        if actual != value:
            raise ValueError(f"frozen material-path contract drift: {name}={actual!r}")
    return mod


def by_seed_radius(measurement: dict) -> dict[str, dict[str, float]]:
    rows = measurement["per_path"]
    result = {}
    for radius in SEED_RADII:
        selected = [row for row in rows if abs(float(row["seed"]["radius"]) - radius) < 1.0e-14]
        if len(selected) != 16:
            raise ValueError("unexpected frozen seed population")
        result[f"{radius:.1f}"] = {
            "path_count": len(selected),
            "mean_absolute_turns": float(np.mean([row["absolute_turns"] for row in selected])),
            "mean_radius_change": float(np.mean([row["radius_change"] for row in selected])),
            "mean_abs_z_change": float(np.mean([row["abs_z_change"] for row in selected])),
        }
    return result


def comparison(child: dict, base: dict) -> dict:
    def ratio(key: str) -> float:
        return float(child[key] / base[key])
    return {
        "mean_absolute_turns_relative": ratio("mean_absolute_turns") - 1.0,
        "maximum_absolute_turns_relative": ratio("maximum_absolute_turns") - 1.0,
        "radial_contraction_magnitude_relative": (
            abs(float(child["mean_radius_change"])) / abs(float(base["mean_radius_change"])) - 1.0
        ),
        "mean_pair_axial_separation_change_relative": ratio("mean_pair_axial_separation_change") - 1.0,
        "inward_path_count_delta": int(child["inward_path_count"] - base["inward_path_count"]),
        "pair_growth_count_delta": int(
            child["pair_axial_separation_growth_count"] - base["pair_axial_separation_growth_count"]
        ),
        "pair_shrink_count_delta": int(
            child["pair_axial_separation_shrink_count"] - base["pair_axial_separation_shrink_count"]
        ),
    }


def _check_base_reference(base: dict) -> None:
    observed = {
        "mean_absolute_turns": base["mean_absolute_turns"],
        "maximum_absolute_turns": base["maximum_absolute_turns"],
        "radial_contraction_magnitude": abs(base["mean_radius_change"]),
        "mean_pair_axial_separation_change": base["mean_pair_axial_separation_change"],
    }
    for key, ref in BASE_REFERENCE.items():
        if abs(float(observed[key]) - ref) > BASE_REFERENCE_ATOL:
            raise RuntimeError(f"frozen `.025` path control drift for {key}: {observed[key]} vs {ref}")


def build_report(agent7_root: str | Path, agent9_root: str | Path) -> dict:
    a7 = load_agent7(agent7_root)
    path_engine = load_path_engine(agent9_root)
    family, raw = a7.replay_st051.reconstruct(PARENT_ID)
    energy = a7.energy_scales(family, raw)
    finest = energy[-1]
    red_scale = float(finest["redistribution_scale"])
    prior_red_scale = float(a7.SOURCE_REDISTRIBUTION_SCALE)
    piola_scale = float(finest["piola_scale"])

    def frozen_reference_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return a7.redistributed_velocity(
            family, raw, np.asarray(points, dtype=float), float(time),
            gain=REDISTRIBUTION_GAIN, scale=prior_red_scale,
        )

    def base_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return a7.redistributed_velocity(
            family, raw, np.asarray(points, dtype=float), float(time),
            gain=REDISTRIBUTION_GAIN, scale=red_scale,
        )

    def child_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return a7.piola_velocity(
            family, raw, np.asarray(points, dtype=float), float(time),
            redistribution_gain=REDISTRIBUTION_GAIN,
            redistribution_scale=red_scale,
            beta=PIOLA_BETA,
            piola_scale=piola_scale,
        )

    # Preserve the original CR-A9-055 replay check at its original scale and
    # tolerance. The current Agent-7 callable screen independently recomputes
    # the common normalization by quadrature; that current-scale control is
    # compared apples-to-apples with its Piola child below.
    frozen_reference = path_engine.measure_material_paths(frozen_reference_velocity)
    _check_base_reference(frozen_reference)
    base = path_engine.measure_material_paths(base_velocity)
    child = path_engine.measure_material_paths(child_velocity)
    base_by_radius = by_seed_radius(base)
    child_by_radius = by_seed_radius(child)
    radius_cmp = {}
    for key in base_by_radius:
        b = base_by_radius[key]
        c = child_by_radius[key]
        radius_cmp[key] = {
            "mean_absolute_turns_relative": float(c["mean_absolute_turns"] / b["mean_absolute_turns"] - 1.0),
            "radial_contraction_magnitude_relative": float(abs(c["mean_radius_change"]) / abs(b["mean_radius_change"]) - 1.0),
            "mean_abs_z_change_delta": float(c["mean_abs_z_change"] - b["mean_abs_z_change"]),
        }

    report = {
        "task_id": TASK_ID,
        "schema": SCHEMA,
        "base_main_sha": BASE_MAIN_SHA,
        "lineage": {
            "parent_candidate": PARENT_ID,
            "parent_pr": PARENT_PR,
            "parent_head": PARENT_HEAD,
            "agent7_callable_piola_pr": AGENT7_PR,
            "agent7_callable_piola_head": AGENT7_HEAD,
            "agent7_task_id": AGENT7_TASK_ID,
            "agent9_path_engine_pr": AGENT9_PATH_PR,
            "agent9_path_engine_head": AGENT9_PATH_HEAD,
        },
        "transform": {
            "redistribution_gain": REDISTRIBUTION_GAIN,
            "redistribution_alpha": REDISTRIBUTION_ALPHA,
            "piola_beta": PIOLA_BETA,
            "gain_scan_performed": False,
            "beta_scan_performed": False,
            "callable_parent_used": True,
            "sampled_grid_interpolation_used": False,
        },
        "normalization": {
            "energy_order": int(finest["order"]),
            "parent_energy": float(finest["parent_energy"]),
            "previously_audited_redistribution_scale": prior_red_scale,
            "redistribution_scale": red_scale,
            "redistribution_scale_shift_from_prior": float(red_scale - prior_red_scale),
            "redistribution_energy": float(finest["redistribution_energy"]),
            "piola_raw_energy": float(finest["piola_raw_energy"]),
            "piola_scale": piola_scale,
            "final_energy": float(finest["final_energy"]),
            "final_relative_to_redistribution": float(finest["final_relative_to_redistribution"]),
        },
        "frozen_material_path_contract": {
            "time_interval": [0.25, 0.75],
            "seed_radii": list(SEED_RADII),
            "seed_z": [-0.3, 0.3],
            "seed_angles": 8,
            "path_count": 48,
            "paired_material_line_count": 24,
            "output_samples": 33,
            "solver": "DOP853",
            "rtol": 1.0e-9,
            "atol": 1.0e-11,
            "max_step": 0.01,
        },
        "frozen_cr_a9_055_reference_material_paths": frozen_reference,
        "redistribution_only_material_paths": base,
        "piola_child_material_paths": child,
        "callable_control_vs_frozen_reference": comparison(base, frozen_reference),
        "redistribution_only_by_seed_radius": base_by_radius,
        "piola_child_by_seed_radius": child_by_radius,
        "piola_vs_redistribution": comparison(child, base),
        "piola_vs_redistribution_by_seed_radius": radius_cmp,
        "external_method": SCIPY_SOURCE,
        "internal_method": {
            "classification": "direct internal method reuse / independent material-path replay",
            "migrated_scope": "exact callable ST051-B `.025` redistribution and beta=.075 Piola field from PR #510; exact frozen material-path measurement from PR #435",
            "difference": "Agent-7 Eulerian morphology/capacity evidence is replaced by cumulative real-time material trajectories",
        },
        "truth_boundary": dict(TRUTH_BOUNDARY),
        "direct_contribution_to_final_velocity": "tests whether the independent callable axial Piola coordinate improves or harms cumulative trajectory geometry before any combined child is materialized",
    }
    report["report_sha256"] = _canonical_hash(report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--agent7-root", required=True)
    parser.add_argument("--agent9-root", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(args.agent7_root, args.agent9_root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
