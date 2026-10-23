"""Replay ST051-B frozen redistribution at static gain .05 on material paths.

CR-A9-055 consumes Agent-7 PR #480's preregistered static-gain headroom row and
measures gain=.05 with the unchanged Agent-9 48-material-path contract from
PR #435.  Gain=.025 is replayed in the same run as the already-audited control.
This is visualization-side routing evidence only: no parent pressure/forcing/PDE
receipt is transferred, no hidden OpenAI numerical target is used, and no visual
acceptance threshold is defined.
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

TASK_ID = "CR-A9-055"
SCHEMA = "st051b_static_gain005_material_path_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"
ST051_PARENT_PR = 460
ST051_PARENT_HEAD = "4b784f1b8457af2ead49295631d834d4e882000b"
AGENT7_SOURCE_PR = 480
AGENT7_SOURCE_HEAD = "187157d377b14bae56415648ac100342f2b8bbb0"
AGENT7_TASK_ID = "CR003-ST051B-FROZEN-REDISTRIBUTION-HEADROOM-073"
SOURCE_TRANSFER_PR = 469
SOURCE_PROFILE_PR = 458
PARENT_ID = "ST051-B"
LOWER_GAIN = 0.025
HIGHER_GAIN = 0.05
FROZEN_SOURCE_ALPHA = 2.520520814687742
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
AGENT7_REPORTED_SCALE_0025 = 1.0014791925672812
AGENT7_REPORTED_SCALE_005 = 1.002495582074808
AGENT9_PATH_SOURCE_PR = 435
AGENT9_PATH_SOURCE_HEAD = "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"
SEED_RADII = (0.6, 0.9, 1.2)

SCIPY_SOURCE = {
    "repository": "scipy/scipy",
    "commit": "eff78058ed0d4cb8e1f0b5c5585d62d712e948d2",
    "api": "scipy.integrate.solve_ivp / DOP853",
    "license": "BSD-3-Clause",
    "classification": "direct migration / public API only through frozen Agent-9 path engine",
    "copied_upstream_implementation": False,
}
INTERNAL_SOURCES = [
    {
        "repository": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
        "pr": AGENT7_SOURCE_PR,
        "commit": AGENT7_SOURCE_HEAD,
        "task_id": AGENT7_TASK_ID,
        "classification": "direct internal method reuse",
        "migrated_scope": (
            "ST051-B reconstruction, frozen inner/mid redistribution identity, and the "
            "preregistered static gain=.05 headroom row"
        ),
        "difference": "Agent-7 Eulerian capacity proxy is replaced by cumulative material trajectories",
    },
    {
        "repository": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
        "pr": AGENT9_PATH_SOURCE_PR,
        "commit": AGENT9_PATH_SOURCE_HEAD,
        "classification": "direct internal evidence-method reuse",
        "migrated_scope": "frozen 48-path / 24-pair seed, time, and DOP853 measurement contract",
    },
]
TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "production_candidate_selected": False,
    "production_redistribution_gain_selected": False,
    "redistribution_rebalanced_on_st051b": False,
    "pressure_or_force_changed": False,
    "pressure_or_force_transferred_from_parent": False,
    "held_out_pde_residual_evaluated": False,
    "held_out_pde_residual_transferred_from_parent": False,
    "openai_numeric_target_inferred": False,
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


def _git_head(root: str | Path) -> str:
    done = subprocess.run(
        ["git", "-C", str(Path(root).resolve()), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return done.stdout.strip()


def _hash(payload: dict) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(raw).hexdigest()


def load_path_engine(agent9_root: str | Path):
    root = Path(agent9_root).resolve()
    if _git_head(root) != AGENT9_PATH_SOURCE_HEAD:
        raise ValueError("Agent-9 path checkout is not pinned PR #435 head")
    old_path = list(sys.path)
    try:
        sys.path.insert(0, str(root / "src"))
        module = importlib.import_module(
            "openai_ns_reconstruction.constrained_st048s_piola_swirl_material_path"
        )
    finally:
        sys.path[:] = old_path
    expected = {
        "SEED_RADII": SEED_RADII,
        "SEED_Z": (-0.3, 0.3),
        "SEED_ANGLES": 8,
        "OUTPUT_SAMPLES": 33,
        "SOLVER_METHOD": "DOP853",
        "SOLVER_RTOL": 1e-9,
        "SOLVER_ATOL": 1e-11,
        "SOLVER_MAX_STEP": 0.01,
        "REGISTERED_TIME_INTERVAL": (0.25, 0.75),
        "EVALUATION_BOX": (-2.0, 2.0),
    }
    for key, value in expected.items():
        if getattr(module, key) != value:
            raise ValueError(f"frozen material-path contract drift: {key}")
    return module


def _purge_experiment_modules(root: Path) -> None:
    experiment_root = (root / "experiments").resolve()
    generic = {
        "agent7_st051b_redistribution_headroom",
        "agent7_st051b_frozen_redistribution_transfer",
        "replay_st051",
        "aligned_continuation",
        "replay_recovery",
        "diagnostics",
        "pressure_morph",
        "spacetime",
        "validate",
        "replay_st048",
        "replay_st047",
        "replay_st046",
        "continuation",
        "boundary_shear",
        "localized_model",
        "base",
        "axisreg",
        "bump",
        "support",
        "coupled_fit",
        "structure_and_particles",
    }
    for name, module in list(sys.modules.items()):
        filename = getattr(module, "__file__", None)
        under_experiments = False
        if filename:
            try:
                Path(filename).resolve().relative_to(experiment_root)
                under_experiments = True
            except (OSError, ValueError):
                pass
        if name in generic or under_experiments:
            sys.modules.pop(name, None)


def load_agent7_source_and_parent(candidate_root: str | Path):
    """Load exact #480 and reconstruct its frozen ST051-B parent."""
    root = Path(candidate_root).resolve()
    if _git_head(root) != AGENT7_SOURCE_HEAD:
        raise ValueError("ST051-B/Agent-7 checkout is not pinned PR #480 head")
    experiment = root / "experiments/root_st051"
    old_path = list(sys.path)
    try:
        _purge_experiment_modules(root)
        sys.path.insert(0, str(experiment))
        screen = importlib.import_module("agent7_st051b_redistribution_headroom")
        if screen.TASK_ID != AGENT7_TASK_ID:
            raise ValueError("Agent-7 task drift")
        if screen.PARENT_ID != PARENT_ID or screen.PARENT_HEAD != ST051_PARENT_HEAD:
            raise ValueError("Agent-7 ST051-B lineage drift")
        if screen.SOURCE_PROFILE_PR != SOURCE_PROFILE_PR or screen.TRANSFER_PR != SOURCE_TRANSFER_PR:
            raise ValueError("frozen redistribution lineage drift")
        if float(screen.SOURCE_ALPHA) != FROZEN_SOURCE_ALPHA:
            raise ValueError("frozen redistribution alpha drift")
        if tuple(screen.INNER_WINDOW) != INNER_WINDOW or tuple(screen.OUTER_WINDOW) != OUTER_WINDOW:
            raise ValueError("frozen redistribution window drift")
        if tuple(float(x) for x in screen.GAINS) != (0.0, 0.015, 0.025, 0.035, 0.05):
            raise ValueError("preregistered gain grid drift")
        family, raw = screen.replay_st051.reconstruct(PARENT_ID)
        return screen, family, raw
    finally:
        sys.path[:] = old_path


def _parent_callable(family, raw) -> VelocityCallable:
    def velocity(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("points must be finite shape (n,3)")
        values, _pressure = family.fields(raw, points, np.full(len(points), float(time)))
        values = np.asarray(values, dtype=float)
        if values.shape != points.shape or not np.isfinite(values).all():
            raise ValueError("parent velocity returned invalid values")
        return values
    return velocity


def _summarize_by_seed_radius(measurement: dict) -> dict:
    paths = {f"{r:.1f}": [] for r in SEED_RADII}
    pairs = {f"{r:.1f}": [] for r in SEED_RADII}
    for row in measurement["per_path"]:
        key = f"{float(row['seed']['radius']):.1f}"
        if key not in paths:
            raise ValueError("unexpected seed radius in frozen path receipt")
        paths[key].append(row)
    for row in measurement["pair_rows"]:
        key = f"{float(row['radius']):.1f}"
        if key not in pairs:
            raise ValueError("unexpected pair radius in frozen path receipt")
        pairs[key].append(row)
    result = {}
    for key in paths:
        rows, prows = paths[key], pairs[key]
        if len(rows) != 16 or len(prows) != 8:
            raise ValueError("frozen path population drift")
        result[key] = {
            "path_count": len(rows),
            "pair_count": len(prows),
            "mean_absolute_turns": float(np.mean([r["absolute_turns"] for r in rows])),
            "maximum_absolute_turns": float(np.max([r["absolute_turns"] for r in rows])),
            "mean_radius_change": float(np.mean([r["radius_change"] for r in rows])),
            "mean_pair_axial_separation_change": float(
                np.mean([r["axial_separation_change"] for r in prows])
            ),
        }
    return result


def _relative_change(new: float, old: float) -> float:
    if old == 0.0:
        raise ValueError("zero comparison denominator")
    return float(new / old - 1.0)


def _compare(new: dict, old: dict) -> dict:
    return {
        "radial_contraction_magnitude_relative_change": _relative_change(
            abs(new["mean_radius_change"]), abs(old["mean_radius_change"])
        ),
        "mean_absolute_turns_relative_change": _relative_change(
            new["mean_absolute_turns"], old["mean_absolute_turns"]
        ),
        "maximum_absolute_turns_relative_change": _relative_change(
            new["maximum_absolute_turns"], old["maximum_absolute_turns"]
        ),
        "mean_pair_separation_change_relative_change": _relative_change(
            new["mean_pair_axial_separation_change"], old["mean_pair_axial_separation_change"]
        ),
        "inward_path_count_delta": int(new["inward_path_count"] - old["inward_path_count"]),
        "pair_growth_count_delta": int(
            new["pair_axial_separation_growth_count"] - old["pair_axial_separation_growth_count"]
        ),
        "pair_shrink_count_delta": int(
            new["pair_axial_separation_shrink_count"] - old["pair_axial_separation_shrink_count"]
        ),
    }


def _compare_by_radius(new: dict, old: dict) -> dict:
    out = {}
    for key in old:
        a, b = old[key], new[key]
        out[key] = {
            "mean_absolute_turns_relative_change": _relative_change(
                b["mean_absolute_turns"], a["mean_absolute_turns"]
            ),
            "radial_contraction_magnitude_relative_change": _relative_change(
                abs(b["mean_radius_change"]), abs(a["mean_radius_change"])
            ),
            "mean_pair_separation_change_relative_change": _relative_change(
                b["mean_pair_axial_separation_change"], a["mean_pair_axial_separation_change"]
            ),
        }
    return out


def build_report(candidate_root: str | Path, agent9_root: str | Path) -> dict:
    engine = load_path_engine(agent9_root)
    screen, family, raw = load_agent7_source_and_parent(candidate_root)
    parent = _parent_callable(family, raw)

    balance = screen.base.reference_energy_and_moment(family, raw, order=72)
    parent_energy = float(balance["parent_energy"])
    scales = {}
    raw_energies = {}
    for gain in (LOWER_GAIN, HIGHER_GAIN):
        scale, raw_energy = screen.base.child_scale(
            family, raw, parent_energy, gain=gain, order=72
        )
        scales[gain] = float(scale)
        raw_energies[gain] = float(raw_energy)
    if abs(scales[LOWER_GAIN] - AGENT7_REPORTED_SCALE_0025) > 5e-10:
        raise RuntimeError("Agent-7 gain=.025 normalization cross-check failed")
    if abs(scales[HIGHER_GAIN] - AGENT7_REPORTED_SCALE_005) > 5e-10:
        raise RuntimeError("Agent-7 gain=.05 normalization cross-check failed")

    def child_callable(gain: float) -> VelocityCallable:
        scale = scales[gain]
        def velocity(points: np.ndarray, time: float) -> np.ndarray:
            return np.asarray(
                screen.base.transformed_velocity(
                    family,
                    raw,
                    np.asarray(points, dtype=float),
                    float(time),
                    gain,
                    scale,
                ),
                dtype=float,
            )
        return velocity

    parent_paths = engine.measure_material_paths(parent)
    lower_paths = engine.measure_material_paths(child_callable(LOWER_GAIN))
    higher_paths = engine.measure_material_paths(child_callable(HIGHER_GAIN))
    parent_by_radius = _summarize_by_seed_radius(parent_paths)
    lower_by_radius = _summarize_by_seed_radius(lower_paths)
    higher_by_radius = _summarize_by_seed_radius(higher_paths)

    st006 = engine.ST006_REFERENCE
    higher_vs_st006 = {
        "mean_absolute_turns_relative_change": _relative_change(
            higher_paths["mean_absolute_turns"], st006["mean_absolute_turns"]
        ),
        "maximum_absolute_turns_relative_change": _relative_change(
            higher_paths["maximum_absolute_turns"], st006["maximum_absolute_turns"]
        ),
        "radial_contraction_magnitude_relative_change": _relative_change(
            abs(higher_paths["mean_radius_change"]), abs(st006["mean_radius_change"])
        ),
        "mean_pair_separation_change_relative_change": _relative_change(
            higher_paths["mean_pair_axial_separation_change"], st006["mean_pair_separation_change"]
        ),
        "comparison_role": "descriptive retained-repository baseline only; ST006 is not OpenAI truth",
    }

    report = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_main_sha": BASE_MAIN_SHA,
        "lineage": {
            "parent_candidate": PARENT_ID,
            "parent_pr": ST051_PARENT_PR,
            "parent_head": ST051_PARENT_HEAD,
            "agent7_headroom_pr": AGENT7_SOURCE_PR,
            "agent7_headroom_head": AGENT7_SOURCE_HEAD,
            "frozen_transfer_pr": SOURCE_TRANSFER_PR,
            "frozen_profile_source_pr": SOURCE_PROFILE_PR,
            "agent9_path_engine_pr": AGENT9_PATH_SOURCE_PR,
            "agent9_path_engine_head": AGENT9_PATH_SOURCE_HEAD,
        },
        "frozen_material_path_contract": {
            "time_interval": list(engine.REGISTERED_TIME_INTERVAL),
            "seed_radii": list(engine.SEED_RADII),
            "seed_z": list(engine.SEED_Z),
            "seed_angles": engine.SEED_ANGLES,
            "path_count": 48,
            "paired_material_line_count": 24,
            "output_samples": engine.OUTPUT_SAMPLES,
            "solver": {
                "method": engine.SOLVER_METHOD,
                "rtol": engine.SOLVER_RTOL,
                "atol": engine.SOLVER_ATOL,
                "max_step": engine.SOLVER_MAX_STEP,
            },
        },
        "redistribution": {
            "lower_gain_control": LOWER_GAIN,
            "higher_gain_under_test": HIGHER_GAIN,
            "inner_window": list(INNER_WINDOW),
            "outer_window": list(OUTER_WINDOW),
            "frozen_source_alpha": FROZEN_SOURCE_ALPHA,
            "parent_energy": parent_energy,
            "lower_gain_scale": scales[LOWER_GAIN],
            "higher_gain_scale": scales[HIGHER_GAIN],
            "lower_gain_raw_energy": raw_energies[LOWER_GAIN],
            "higher_gain_raw_energy": raw_energies[HIGHER_GAIN],
            "lower_gain_scale_abs_difference_from_agent7": float(
                abs(scales[LOWER_GAIN] - AGENT7_REPORTED_SCALE_0025)
            ),
            "higher_gain_scale_abs_difference_from_agent7": float(
                abs(scales[HIGHER_GAIN] - AGENT7_REPORTED_SCALE_005)
            ),
            "transform_rebalanced_on_st051b": False,
            "gain_values_outside_preregistered_grid_evaluated": False,
        },
        "parent_material_paths": parent_paths,
        "lower_gain_material_paths": lower_paths,
        "higher_gain_material_paths": higher_paths,
        "parent_by_seed_radius": parent_by_radius,
        "lower_gain_by_seed_radius": lower_by_radius,
        "higher_gain_by_seed_radius": higher_by_radius,
        "higher_vs_lower": _compare(higher_paths, lower_paths),
        "higher_vs_parent": _compare(higher_paths, parent_paths),
        "higher_vs_lower_by_seed_radius": _compare_by_radius(higher_by_radius, lower_by_radius),
        "higher_vs_parent_by_seed_radius": _compare_by_radius(higher_by_radius, parent_by_radius),
        "higher_vs_st006_descriptive": higher_vs_st006,
        "external_method": SCIPY_SOURCE,
        "internal_sources": INTERNAL_SOURCES,
        "truth_boundary": TRUTH_BOUNDARY,
        "direct_contribution": (
            "Tests whether Agent-7's preregistered static gain=.05 headroom improves real cumulative "
            "ST051-B inner/mid material-path winding beyond the already-audited gain=.025 child without "
            "using Eulerian proxy gain as a production or PDE acceptance rule."
        ),
    }
    report["report_sha256"] = _hash(report)
    return report


def _audit_truth(report: dict) -> None:
    if report["task_id"] != TASK_ID or report["base_main_sha"] != BASE_MAIN_SHA:
        raise ValueError("report identity drift")
    if report["frozen_material_path_contract"]["path_count"] != 48:
        raise ValueError("path population drift")
    for key in ("parent_material_paths", "lower_gain_material_paths", "higher_gain_material_paths"):
        if report[key]["path_count"] != 48 or report[key]["paired_material_line_count"] != 24:
            raise ValueError(f"{key} frozen population drift")
    if report["redistribution"]["higher_gain_under_test"] != HIGHER_GAIN:
        raise ValueError("higher-gain identity drift")
    if report["redistribution"]["lower_gain_control"] != LOWER_GAIN:
        raise ValueError("lower-gain control drift")
    if report["redistribution"]["transform_rebalanced_on_st051b"] is not False:
        raise ValueError("transform identity drift")
    if report["redistribution"]["gain_values_outside_preregistered_grid_evaluated"] is not False:
        raise ValueError("gain-grid truth drift")
    for key, value in TRUTH_BOUNDARY.items():
        if report["truth_boundary"].get(key) is not value:
            raise ValueError(f"truth-boundary drift: {key}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", required=True, type=Path)
    parser.add_argument("--agent9-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    report = build_report(args.candidate_root, args.agent9_root)
    _audit_truth(report)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "task_id": TASK_ID,
                "parent_mean_turns": report["parent_material_paths"]["mean_absolute_turns"],
                "lower_gain_mean_turns": report["lower_gain_material_paths"]["mean_absolute_turns"],
                "higher_gain_mean_turns": report["higher_gain_material_paths"]["mean_absolute_turns"],
                "higher_vs_lower": report["higher_vs_lower"],
                "higher_vs_lower_by_seed_radius": report["higher_vs_lower_by_seed_radius"],
                "report_sha256": report["report_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
