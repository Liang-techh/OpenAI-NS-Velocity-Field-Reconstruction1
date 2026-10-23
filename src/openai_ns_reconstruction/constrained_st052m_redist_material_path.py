"""Replay frozen ST052-M + kappa=.05 redistribution on Agent-9 material paths.

CR-A9-059 consumes Agent 7 PR #528's exact frozen one-degree swirl
redistribution on ST052-M and compares it against the exact ST052-M parent
under Agent 9's unchanged 48-path / 24-pair material-trajectory contract.

This module performs no gain scan, pressure/forcing refit, PDE promotion,
OpenAI image fitting, or visualization acceptance test.
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

TASK_ID = "CR-A9-059"
SCHEMA = "st052m_frozen_redistribution_material_path_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"

AGENT7_PR = 528
AGENT7_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"
AGENT7_TASK_ID = "CR003-ST052M-FROZEN-REDISTRIBUTION-078"
PARENT_ID = "ST052-M"
PARENT_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
REDISTRIBUTION_ALPHA = 2.520520814687742
REDISTRIBUTION_GAIN = 0.05
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
EXPECTED_AGENT7_NORMALIZATION = 1.0032534663681094

AGENT9_PATH_PR = 435
AGENT9_PATH_HEAD = "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"
SEED_RADII = (0.6, 0.9, 1.2)

SCIPY_SOURCE = {
    "repository": "scipy/scipy",
    "commit": "eff78058ed0d4cb8e1f0b5c5585d62d712e948d2",
    "api": "scipy.integrate.solve_ivp / DOP853",
    "license": "BSD-3-Clause",
    "classification": "direct migration / public API only through frozen Agent-9 path engine",
    "copied_upstream_implementation": False,
}

UPSTREAM_AGENT7_RECEIPT = {
    "classification": "upstream Eulerian expression-capacity receipt only; independently replayed here on material paths",
    "clean_frozen_profile_transfer": True,
    "normalization": EXPECTED_AGENT7_NORMALIZATION,
    "angular_gain_r06": 0.05136897375553384,
    "angular_gain_r09": 0.03183529606752522,
    "angular_gain_r12": -0.09559449216522342,
    "divergence_fd_max": 2.496960971321016e-09,
    "response_rank": 2,
    "response_condition_number": 1.0740026848082393,
    "response_cosine": -0.07127144948660563,
    "pde_validated": False,
    "material_paths_integrated": False,
}

TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "production_candidate_selected": False,
    "production_redistribution_gain_selected": False,
    "pressure_or_force_refit_in_this_increment": False,
    "held_out_pde_residual_recomputed_in_this_increment": False,
    "parent_pde_receipt_transferred": False,
    "upstream_eulerian_receipt_promoted_to_pde_truth": False,
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


def load_agent7_source(root: str | Path):
    root = Path(root).resolve()
    if _git_head(root) != AGENT7_HEAD:
        raise ValueError("Agent-7 checkout is not pinned PR #528 exact head")
    path = root / "experiments" / "root_st052" / "agent7_st052m_frozen_redistribution.py"
    spec = importlib.util.spec_from_file_location("agent7_st052m_redist", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen Agent-7 ST052-M redistribution source")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    expected = {
        "TASK_ID": AGENT7_TASK_ID,
        "PARENT_ID": PARENT_ID,
        "PARENT_HEAD": PARENT_HEAD,
        "SOURCE_ALPHA": REDISTRIBUTION_ALPHA,
        "GAIN": REDISTRIBUTION_GAIN,
        "INNER_WINDOW": INNER_WINDOW,
        "OUTER_WINDOW": OUTER_WINDOW,
    }
    for name, value in expected.items():
        if getattr(mod, name) != value:
            raise ValueError(f"frozen Agent-7 identity drift: {name}")
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
        if getattr(mod, name) != value:
            raise ValueError(f"frozen material-path contract drift: {name}")
    return mod


def by_seed_radius(measurement: dict) -> dict[str, dict[str, float]]:
    rows = measurement["per_path"]
    result: dict[str, dict[str, float]] = {}
    for radius in SEED_RADII:
        selected = [
            row for row in rows
            if abs(float(row["seed"]["radius"]) - radius) < 1.0e-14
        ]
        if len(selected) != 16:
            raise ValueError("unexpected frozen seed population")
        result[f"{radius:.1f}"] = {
            "path_count": len(selected),
            "mean_absolute_turns": float(np.mean([row["absolute_turns"] for row in selected])),
            "mean_radius_change": float(np.mean([row["radius_change"] for row in selected])),
            "mean_abs_z_change": float(np.mean([row["abs_z_change"] for row in selected])),
        }
    return result


def comparison(child: dict, parent: dict) -> dict:
    return {
        "mean_absolute_turns_relative": float(child["mean_absolute_turns"] / parent["mean_absolute_turns"] - 1.0),
        "maximum_absolute_turns_relative": float(child["maximum_absolute_turns"] / parent["maximum_absolute_turns"] - 1.0),
        "radial_contraction_magnitude_relative": float(
            abs(child["mean_radius_change"]) / abs(parent["mean_radius_change"]) - 1.0
        ),
        "mean_pair_axial_separation_change_relative": float(
            child["mean_pair_axial_separation_change"] / parent["mean_pair_axial_separation_change"] - 1.0
        ),
        "inward_path_count_delta": int(child["inward_path_count"] - parent["inward_path_count"]),
        "pair_growth_count_delta": int(
            child["pair_axial_separation_growth_count"] - parent["pair_axial_separation_growth_count"]
        ),
        "pair_shrink_count_delta": int(
            child["pair_axial_separation_shrink_count"] - parent["pair_axial_separation_shrink_count"]
        ),
    }


def build_report(agent7_root: str | Path, agent9_root: str | Path) -> dict:
    source = load_agent7_source(agent7_root)
    path_engine = load_path_engine(agent9_root)

    family, raw = source.replay_st052.reconstruct()
    scale, parent_energy, unscaled_child_energy = source.child_scale(family, raw)
    if abs(scale - EXPECTED_AGENT7_NORMALIZATION) > 2.0e-12:
        raise ValueError("Agent-7 normalization drift")

    def parent_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return source.base_velocity(family, raw, np.asarray(points, dtype=float), float(time))

    def child_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return source.transformed_velocity(
            family,
            raw,
            np.asarray(points, dtype=float),
            float(time),
            source.GAIN,
            scale,
        )

    parent = path_engine.measure_material_paths(parent_velocity)
    child = path_engine.measure_material_paths(child_velocity)
    parent_by_radius = by_seed_radius(parent)
    child_by_radius = by_seed_radius(child)
    radius_cmp = {}
    for key in parent_by_radius:
        p = parent_by_radius[key]
        c = child_by_radius[key]
        radius_cmp[key] = {
            "mean_absolute_turns_relative": float(c["mean_absolute_turns"] / p["mean_absolute_turns"] - 1.0),
            "radial_contraction_magnitude_relative": float(
                abs(c["mean_radius_change"]) / abs(p["mean_radius_change"]) - 1.0
            ),
            "mean_abs_z_change_delta": float(c["mean_abs_z_change"] - p["mean_abs_z_change"]),
        }

    report = {
        "task_id": TASK_ID,
        "schema": SCHEMA,
        "base_main_sha": BASE_MAIN_SHA,
        "lineage": {
            "parent_candidate": PARENT_ID,
            "parent_head": PARENT_HEAD,
            "agent7_source_pr": AGENT7_PR,
            "agent7_source_head": AGENT7_HEAD,
            "agent7_task_id": AGENT7_TASK_ID,
            "agent9_path_engine_pr": AGENT9_PATH_PR,
            "agent9_path_engine_head": AGENT9_PATH_HEAD,
        },
        "frozen_transform": {
            "alpha": REDISTRIBUTION_ALPHA,
            "gain": REDISTRIBUTION_GAIN,
            "inner_window": list(INNER_WINDOW),
            "outer_window": list(OUTER_WINDOW),
            "normalization": scale,
            "parent_reference_energy": parent_energy,
            "unscaled_child_reference_energy": unscaled_child_energy,
        },
        "upstream_agent7_receipt": dict(UPSTREAM_AGENT7_RECEIPT),
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
        "st052_parent_material_paths": parent,
        "redistributed_child_material_paths": child,
        "st052_parent_by_seed_radius": parent_by_radius,
        "redistributed_child_by_seed_radius": child_by_radius,
        "child_vs_parent": comparison(child, parent),
        "child_vs_parent_by_seed_radius": radius_cmp,
        "external_method": SCIPY_SOURCE,
        "internal_method": {
            "classification": "direct internal method reuse / independent frozen material-path replay",
            "migrated_scope": "exact PR #528 ST052-M redistribution transform plus exact PR #435 trajectory measurement contract",
            "difference": "Agent 7 #528 assessed Eulerian angular-rate/vorticity/response capacity and explicitly did not integrate material paths; this increment measures the same frozen child on Agent 9's cross-candidate trajectory protocol",
        },
        "truth_boundary": dict(TRUTH_BOUNDARY),
        "direct_contribution_to_final_velocity": (
            "tests whether the already-clean one-dimensional ST052-M winding-control transform improves real-time inward-spiral trajectories before any governed child materialization or pressure/restricted-force rebuild"
        ),
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
