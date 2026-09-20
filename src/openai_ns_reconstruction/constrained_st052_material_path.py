"""Replay frozen ST052-M against its ST051-B parent on Agent-9 material paths.

CR-A9-058 consumes the exact frozen ST052-M candidate from root PR #508 and
compares it with the immutable ST051-B parent under Agent 9's unchanged
48-path material-trajectory protocol.  This answers one narrow routing
question: did the sampled-minimax PDE improvement preserve or damage the
candidate-side inward-spiral trajectory geometry?

This module does not refit coefficients, transfer a parent PDE receipt, define
an OpenAI image target, or promote visualization/PDE truth states.
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

TASK_ID = "CR-A9-058"
SCHEMA = "st052_frozen_material_path_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"

ST052_ID = "ST052-M"
ST052_PR = 508
ST052_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
ST052_PARENT_ID = "ST051-B"
ST052_PARENT_SHA256 = "0071106ef10a5d77b620b942becc00b58c8a6765fd5fb7bc4dcd195ac65c970d"
ST052_REPORTED_RAW_SHA256 = "e078e753fab38ebfa0284d28ba64d26cb8538b700a3afc7705668849c43e12da"

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

TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "production_candidate_selected": False,
    "pressure_or_force_refit_in_this_increment": False,
    "held_out_pde_residual_recomputed_in_this_increment": False,
    "parent_pde_receipt_transferred": False,
    "st052_existing_pde_receipt_promoted_by_this_increment": False,
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


def load_st052(root: str | Path):
    root = Path(root).resolve()
    if _git_head(root) != ST052_HEAD:
        raise ValueError("ST052 checkout is not pinned PR #508 exact head")
    exp = root / "experiments" / "root_st052"
    old_path = list(sys.path)
    names = (
        "replay_st052",
        "minimax_exchange",
        "replay_st051",
        "aligned_continuation",
        "spacetime",
    )
    try:
        sys.path.insert(0, str(exp))
        for name in names:
            sys.modules.pop(name, None)
        mod = importlib.import_module("replay_st052")
    finally:
        sys.path[:] = old_path
    if str(getattr(mod, "PARENT_SHA", "")) != ST052_PARENT_SHA256:
        raise ValueError("ST052 parent identity drift")
    recipe = json.loads((exp / "recipe.json").read_text())
    if recipe.get("parent_id") != ST052_PARENT_ID:
        raise ValueError("ST052 recipe parent id drift")
    if recipe.get("parent_sha256") != ST052_PARENT_SHA256:
        raise ValueError("ST052 recipe parent SHA drift")
    if recipe.get("pde_validated") is not False:
        raise ValueError("ST052 recipe illegally promotes PDE validity")
    if recipe.get("source_correspondence_verified") is not False:
        raise ValueError("ST052 recipe illegally promotes source correspondence")
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


def build_report(st052_root: str | Path, agent9_root: str | Path) -> dict:
    st052 = load_st052(st052_root)
    path_engine = load_path_engine(agent9_root)

    child_family, child_raw = st052.reconstruct()
    parent_family, parent_raw = st052.parent_field()

    def parent_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return parent_family.fields(parent_raw, np.asarray(points, dtype=float), float(time))[0]

    def child_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return child_family.fields(child_raw, np.asarray(points, dtype=float), float(time))[0]

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
            "candidate": ST052_ID,
            "candidate_pr": ST052_PR,
            "candidate_head": ST052_HEAD,
            "parent_candidate": ST052_PARENT_ID,
            "parent_sha256": ST052_PARENT_SHA256,
            "reported_candidate_raw_sha256": ST052_REPORTED_RAW_SHA256,
            "agent9_path_engine_pr": AGENT9_PATH_PR,
            "agent9_path_engine_head": AGENT9_PATH_HEAD,
        },
        "st052_existing_scientific_context": {
            "classification": "upstream receipt only; not rerun or promoted here",
            "reported_holdout_seed_9175291_momentum_max": 0.03681687198184801,
            "reported_holdout_seed_9175291_momentum_volume_L2": 0.05179943624413842,
            "registered_momentum_threshold": 0.001,
            "pde_validated": False,
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
        "parent_material_paths": parent,
        "st052_material_paths": child,
        "parent_by_seed_radius": parent_by_radius,
        "st052_by_seed_radius": child_by_radius,
        "st052_vs_parent": comparison(child, parent),
        "st052_vs_parent_by_seed_radius": radius_cmp,
        "external_method": SCIPY_SOURCE,
        "internal_method": {
            "classification": "direct internal candidate replay / independent frozen-path measurement",
            "migrated_scope": "exact frozen ST052-M recipe and parent loader from PR #508; exact material-path measurement from PR #435",
            "difference": "PR #508 used a separate 36-seed trajectory diagnostic; this increment uses Agent-9's unchanged 48-path / 24-pair contract for cross-candidate routing",
        },
        "truth_boundary": dict(TRUTH_BOUNDARY),
        "direct_contribution_to_final_velocity": (
            "checks whether the residual-improving ST052-M field preserves inward-spiral material-trajectory geometry before visualization/export or further representation work"
        ),
    }
    report["report_sha256"] = _canonical_hash(report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--st052-root", required=True)
    parser.add_argument("--agent9-root", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(args.st052_root, args.agent9_root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
