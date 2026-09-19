"""Replay the frozen ST052-M local-swirl energy child on Agent-9 material paths.

CR-A9-062 consumes Agent 7 PR #559's exact energy-neutral local swirl
compensation on the ST052-M + kappa=.05 redistribution + raw localized taper,
then compares it against the exact untapered redistributed control under Agent
9's unchanged 48-path / 24-pair trajectory contract.

This module performs no coefficient scan, image fitting, pressure/forcing refit,
PDE promotion, production selection, or visual acceptance test.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

import numpy as np

TASK_ID = "CR-A9-062"
SCHEMA = "st052m_local_swirl_energy_material_path_v1"
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
    "clean_local_swirl_energy_capacity": True,
    "energy_root_exists": True,
    "beta": EXPECTED_BETA,
    "post_transform_common_scale": 1.0,
    "central_velocity_identity_abs_max": 0.0,
    "central_angular_rate_retention_r06": 1.0,
    "central_angular_rate_retention_r09": 1.0,
    "tip_radial_rms_relative": [-0.0137325, -0.0149126, -0.0136258],
    "central_radial_rms_relative": [0.0, 0.0, 0.0],
    "global_axial_rms_relative": [0.0106655, 0.0130320, 0.0150734],
    "divergence_fd_max": 2.322e-9,
    "response_rank": 4,
    "response_condition_number": 1.507887,
    "response_max_pairwise_abs_cosine": 0.293159,
    "material_paths_integrated": False,
    "pde_validated": False,
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
    "comparison_is_descriptive_not_acceptance": True,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


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
    if mod.base.SOURCE_ALPHA != REDISTRIBUTION_ALPHA:
        raise ValueError("frozen redistribution alpha drift")
    if mod.base.GAIN != REDISTRIBUTION_GAIN:
        raise ValueError("frozen redistribution gain drift")
    if mod.base.INNER_WINDOW != INNER_WINDOW or mod.base.OUTER_WINDOW != OUTER_WINDOW:
        raise ValueError("frozen redistribution radial-window drift")
    if mod.TRUTH["material_paths_integrated"] is not False:
        raise ValueError("upstream material-path truth boundary drift")
    if mod.TRUTH["pde_validated"] is not False:
        raise ValueError("upstream PDE truth boundary drift")
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


def comparison(child: dict, control: dict) -> dict:
    return {
        "mean_absolute_turns_relative": float(child["mean_absolute_turns"] / control["mean_absolute_turns"] - 1.0),
        "maximum_absolute_turns_relative": float(child["maximum_absolute_turns"] / control["maximum_absolute_turns"] - 1.0),
        "radial_contraction_magnitude_relative": float(
            abs(child["mean_radius_change"]) / abs(control["mean_radius_change"]) - 1.0
        ),
        "mean_pair_axial_separation_change_relative": float(
            child["mean_pair_axial_separation_change"] / control["mean_pair_axial_separation_change"] - 1.0
        ),
        "inward_path_count_delta": int(child["inward_path_count"] - control["inward_path_count"]),
        "pair_growth_count_delta": int(
            child["pair_axial_separation_growth_count"] - control["pair_axial_separation_growth_count"]
        ),
        "pair_shrink_count_delta": int(
            child["pair_axial_separation_shrink_count"] - control["pair_axial_separation_shrink_count"]
        ),
    }


def build_report(agent7_root: str | Path, agent9_root: str | Path) -> dict:
    source = load_agent7_source(agent7_root)
    path_engine = load_path_engine(agent9_root)

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
            family,
            raw,
            np.asarray(points, dtype=float),
            float(time),
            redistribution_scale,
        )

    def child_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return source.compensated_velocity(
            family,
            raw,
            np.asarray(points, dtype=float),
            float(time),
            redistribution_scale=redistribution_scale,
            beta=beta,
        )

    control = path_engine.measure_material_paths(control_velocity)
    child = path_engine.measure_material_paths(child_velocity)
    control_by_radius = by_seed_radius(control)
    child_by_radius = by_seed_radius(child)
    radius_cmp = {}
    for key in control_by_radius:
        p = control_by_radius[key]
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
            "st052_parent_candidate": PARENT_ID,
            "st052_parent_head": PARENT_HEAD,
            "control_candidate": CONTROL_ID,
            "control_head": CONTROL_HEAD,
            "agent7_source_pr": AGENT7_PR,
            "agent7_source_head": AGENT7_HEAD,
            "agent7_task_id": AGENT7_TASK_ID,
            "source_taper_head": SOURCE_TAPER_HEAD,
            "source_failed_compensation_head": SOURCE_FAILED_COMP_HEAD,
            "agent9_path_engine_pr": AGENT9_PATH_PR,
            "agent9_path_engine_head": AGENT9_PATH_HEAD,
        },
        "frozen_transform": {
            "redistribution_alpha": REDISTRIBUTION_ALPHA,
            "redistribution_gain": REDISTRIBUTION_GAIN,
            "redistribution_inner_window": list(INNER_WINDOW),
            "redistribution_outer_window": list(OUTER_WINDOW),
            "redistribution_normalization": redistribution_scale,
            "redistribution_parent_reference_energy": parent_energy,
            "redistribution_unscaled_reference_energy": redistribution_raw_energy,
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
        "redistributed_control_material_paths": control,
        "local_swirl_energy_child_material_paths": child,
        "redistributed_control_by_seed_radius": control_by_radius,
        "local_swirl_energy_child_by_seed_radius": child_by_radius,
        "child_vs_control": comparison(child, control),
        "child_vs_control_by_seed_radius": radius_cmp,
        "external_method": SCIPY_SOURCE,
        "internal_method": {
            "repository": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
            "source_pr": AGENT7_PR,
            "source_commit": AGENT7_HEAD,
            "license": "same-repository internal source; no external code copied",
            "classification": "direct internal method reuse / independent frozen material-path replay",
            "migrated_scope": "exact PR #559 local swirl energy child plus exact PR #435 trajectory measurement contract",
            "difference": "Agent 7 #559 assessed Eulerian morphology, energy closure, support/divergence and response capacity and explicitly did not integrate material paths; this increment measures the same frozen child under Agent 9's cross-candidate trajectory protocol",
        },
        "truth_boundary": dict(TRUTH_BOUNDARY),
        "direct_contribution_to_final_velocity": (
            "tests whether the energy-neutral local swirl compensation removes the localized-taper trajectory penalty before any pressure/restricted-force rebuild or production materialization"
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
