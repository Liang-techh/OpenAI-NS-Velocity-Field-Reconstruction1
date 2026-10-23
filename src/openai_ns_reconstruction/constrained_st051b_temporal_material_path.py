"""Replay the ST051-B temporal redistribution on frozen material paths.

CR-A9-056 consumes Agent-7 PR #489's preregistered temporal coefficient
``gamma=.025`` on the already-frozen ST051-B inner/mid redistribution channel.
It compares that one temporal child against static gains .025 and .05 using the
unchanged Agent-9 48-material-path contract from PR #435.

This is visualization-side routing evidence only.  No parent pressure, forcing,
PDE receipt, OpenAI hidden numerical target, or visual acceptance threshold is
transferred or introduced.
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

TASK_ID = "CR-A9-056"
SCHEMA = "st051b_temporal_redistribution_material_path_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"

ST051_PARENT_PR = 460
ST051_PARENT_HEAD = "4b784f1b8457af2ead49295631d834d4e882000b"
PARENT_ID = "ST051-B"

AGENT7_SOURCE_PR = 489
AGENT7_SOURCE_HEAD = "f2828ae259eebc867162ef5842a99f0860053229"
AGENT7_TASK_ID = "CR003-ST051B-TEMPORAL-REDISTRIBUTION-074"
AGENT7_STACK_BASE_PR = 480
AGENT7_STACK_BASE_HEAD = "187157d377b14bae56415648ac100342f2b8bbb0"
SOURCE_PROFILE_PR = 458
SOURCE_TRANSFER_PR = 469

FROZEN_SOURCE_ALPHA = 2.520520814687742
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
STATIC_LOW_GAIN = 0.025
STATIC_HIGH_GAIN = 0.05
TEMPORAL_GAMMA = 0.025
TEMPORAL_GAINS = (0.025, 0.0375, 0.05)
PREREGISTERED_GAMMAS = (0.0, 0.010, 0.015, 0.020, 0.025)

AGENT9_PATH_SOURCE_PR = 435
AGENT9_PATH_SOURCE_HEAD = "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"
SEED_RADII = (0.6, 0.9, 1.2)

# Exact CR-A9-055 static controls.  Replaying these in the same exact-source
# environment prevents a temporal comparison from silently changing the parent
# or path contract.
STATIC_CONTROL_REFERENCE = {
    "0.025": {
        "mean_absolute_turns": 0.02229297191872165,
        "maximum_absolute_turns": 0.027371569482050228,
        "radial_contraction_magnitude": 0.07049762300296604,
        "mean_pair_axial_separation_change": 0.06509150801226499,
    },
    "0.05": {
        "mean_absolute_turns": 0.02249374673546196,
        "maximum_absolute_turns": 0.02789620319169733,
        "radial_contraction_magnitude": 0.07056954974774227,
        "mean_pair_axial_separation_change": 0.06517965283554955,
    },
}

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
        "classification": "direct internal method reuse / independent material-path replay",
        "migrated_scope": (
            "frozen ST051-B redistribution profile and preregistered affine temporal "
            "coefficient gamma=.025"
        ),
        "difference": (
            "Agent-7 Eulerian capacity diagnostics are replaced by cumulative material "
            "trajectories and compared directly with both static controls"
        ),
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
    "production_static_gain_selected": False,
    "production_temporal_coefficient_selected": False,
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
        "agent7_st051b_temporal_redistribution",
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


def load_temporal_source(candidate_root: str | Path):
    root = Path(candidate_root).resolve()
    if _git_head(root) != AGENT7_SOURCE_HEAD:
        raise ValueError("Agent-7 temporal checkout is not pinned PR #489 head")
    experiment = root / "experiments/root_st051"
    old_path = list(sys.path)
    try:
        _purge_experiment_modules(root)
        sys.path.insert(0, str(experiment))
        temporal = importlib.import_module("agent7_st051b_temporal_redistribution")
        if temporal.TASK_ID != AGENT7_TASK_ID:
            raise ValueError("Agent-7 temporal task drift")
        if temporal.PARENT_ID != PARENT_ID or temporal.PARENT_HEAD != ST051_PARENT_HEAD:
            raise ValueError("ST051-B lineage drift")
        if temporal.STACK_BASE_PR != AGENT7_STACK_BASE_PR:
            raise ValueError("Agent-7 static-headroom PR drift")
        if temporal.STACK_BASE_HEAD != AGENT7_STACK_BASE_HEAD:
            raise ValueError("Agent-7 static-headroom head drift")
        if temporal.SOURCE_PROFILE_PR != SOURCE_PROFILE_PR or temporal.TRANSFER_PR != SOURCE_TRANSFER_PR:
            raise ValueError("redistribution source lineage drift")
        if float(temporal.SOURCE_ALPHA) != FROZEN_SOURCE_ALPHA:
            raise ValueError("frozen redistribution alpha drift")
        if tuple(temporal.INNER_WINDOW) != INNER_WINDOW or tuple(temporal.OUTER_WINDOW) != OUTER_WINDOW:
            raise ValueError("frozen redistribution window drift")
        if float(temporal.BASE_GAIN) != STATIC_LOW_GAIN:
            raise ValueError("temporal base-gain drift")
        if tuple(float(x) for x in temporal.GAMMAS) != PREREGISTERED_GAMMAS:
            raise ValueError("preregistered gamma grid drift")
        expected_gains = tuple(float(temporal.gain_at_time(TEMPORAL_GAMMA, t)) for t in temporal.TIMES)
        if not np.allclose(expected_gains, TEMPORAL_GAINS, rtol=0.0, atol=1e-15):
            raise ValueError("temporal gain schedule drift")
        family, raw = temporal.replay_st051.reconstruct(PARENT_ID)
        return temporal, family, raw
    finally:
        sys.path[:] = old_path


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
    result = {}
    for key in old:
        a, b = old[key], new[key]
        result[key] = {
            "mean_absolute_turns_relative_change": _relative_change(
                b["mean_absolute_turns"], a["mean_absolute_turns"]
            ),
            "maximum_absolute_turns_relative_change": _relative_change(
                b["maximum_absolute_turns"], a["maximum_absolute_turns"]
            ),
            "radial_contraction_magnitude_relative_change": _relative_change(
                abs(b["mean_radius_change"]), abs(a["mean_radius_change"])
            ),
            "mean_pair_separation_change_relative_change": _relative_change(
                b["mean_pair_axial_separation_change"], a["mean_pair_axial_separation_change"]
            ),
        }
    return result


def _check_static_control(name: str, measurement: dict) -> None:
    expected = STATIC_CONTROL_REFERENCE[name]
    actual = {
        "mean_absolute_turns": float(measurement["mean_absolute_turns"]),
        "maximum_absolute_turns": float(measurement["maximum_absolute_turns"]),
        "radial_contraction_magnitude": abs(float(measurement["mean_radius_change"])),
        "mean_pair_axial_separation_change": float(measurement["mean_pair_axial_separation_change"]),
    }
    for key, value in expected.items():
        if abs(actual[key] - value) > 5e-12:
            raise RuntimeError(f"static CR-A9-055 control drift: {name} {key}")


def _audit_truth(report: dict) -> None:
    if report["task_id"] != TASK_ID or report["base_main_sha"] != BASE_MAIN_SHA:
        raise ValueError("task/base identity drift")
    redist = report["redistribution"]
    if redist["static_low_gain"] != STATIC_LOW_GAIN or redist["static_high_gain"] != STATIC_HIGH_GAIN:
        raise ValueError("static gain identity drift")
    if redist["temporal_gamma"] != TEMPORAL_GAMMA:
        raise ValueError("temporal gamma identity drift")
    if tuple(redist["temporal_gains_at_registered_times"]) != TEMPORAL_GAINS:
        raise ValueError("temporal schedule identity drift")
    if redist["gamma_values_outside_preregistered_grid_evaluated"]:
        raise ValueError("post-hoc gamma widening")
    if redist["transform_rebalanced_on_st051b"]:
        raise ValueError("transform identity drift")
    for key in ("static_low_material_paths", "static_high_material_paths", "temporal_material_paths"):
        row = report[key]
        if row["path_count"] != 48 or row["paired_material_line_count"] != 24:
            raise ValueError("frozen material-path population drift")
    if report["truth_boundary"] != TRUTH_BOUNDARY:
        raise ValueError("truth-boundary drift")


def build_report(candidate_root: str | Path, agent9_root: str | Path) -> dict:
    engine = load_path_engine(agent9_root)
    temporal, family, raw = load_temporal_source(candidate_root)
    base = temporal.base

    balance = base.reference_energy_and_moment(family, raw, order=72)
    parent_energy = float(balance["parent_energy"])
    low_scale, low_raw_energy = base.child_scale(
        family, raw, parent_energy, gain=STATIC_LOW_GAIN, order=72
    )
    high_scale, high_raw_energy = base.child_scale(
        family, raw, parent_energy, gain=STATIC_HIGH_GAIN, order=72
    )

    def static_callable(gain: float, scale: float) -> VelocityCallable:
        def velocity(points: np.ndarray, time: float) -> np.ndarray:
            return np.asarray(
                base.transformed_velocity(
                    family, raw, np.asarray(points, dtype=float), float(time), float(gain), float(scale)
                ),
                dtype=float,
            )
        return velocity

    def temporal_callable(points: np.ndarray, time: float) -> np.ndarray:
        gain = temporal.gain_at_time(TEMPORAL_GAMMA, float(time))
        return np.asarray(
            base.transformed_velocity(
                family,
                raw,
                np.asarray(points, dtype=float),
                float(time),
                float(gain),
                float(low_scale),
            ),
            dtype=float,
        )

    low_paths = engine.measure_material_paths(static_callable(STATIC_LOW_GAIN, low_scale))
    high_paths = engine.measure_material_paths(static_callable(STATIC_HIGH_GAIN, high_scale))
    temporal_paths = engine.measure_material_paths(temporal_callable)

    _check_static_control("0.025", low_paths)
    _check_static_control("0.05", high_paths)

    low_by_radius = _summarize_by_seed_radius(low_paths)
    high_by_radius = _summarize_by_seed_radius(high_paths)
    temporal_by_radius = _summarize_by_seed_radius(temporal_paths)

    report = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_main_sha": BASE_MAIN_SHA,
        "lineage": {
            "parent_candidate": PARENT_ID,
            "parent_pr": ST051_PARENT_PR,
            "parent_head": ST051_PARENT_HEAD,
            "agent7_temporal_pr": AGENT7_SOURCE_PR,
            "agent7_temporal_head": AGENT7_SOURCE_HEAD,
            "agent7_temporal_task_id": AGENT7_TASK_ID,
            "agent7_static_base_pr": AGENT7_STACK_BASE_PR,
            "agent7_static_base_head": AGENT7_STACK_BASE_HEAD,
            "source_profile_pr": SOURCE_PROFILE_PR,
            "source_transfer_pr": SOURCE_TRANSFER_PR,
            "agent9_path_engine_pr": AGENT9_PATH_SOURCE_PR,
            "agent9_path_engine_head": AGENT9_PATH_SOURCE_HEAD,
        },
        "external_sources": [SCIPY_SOURCE],
        "internal_sources": INTERNAL_SOURCES,
        "frozen_material_path_contract": {
            "time_interval": [0.25, 0.75],
            "seed_radii": list(SEED_RADII),
            "seed_z": [-0.3, 0.3],
            "azimuth_count": 8,
            "path_count": 48,
            "paired_material_line_count": 24,
            "output_samples": 33,
            "solver": "DOP853",
            "rtol": 1e-9,
            "atol": 1e-11,
            "max_step": 0.01,
            "evaluation_box": [-2.0, 2.0],
        },
        "redistribution": {
            "frozen_source_alpha": FROZEN_SOURCE_ALPHA,
            "inner_window": list(INNER_WINDOW),
            "outer_window": list(OUTER_WINDOW),
            "static_low_gain": STATIC_LOW_GAIN,
            "static_high_gain": STATIC_HIGH_GAIN,
            "temporal_gamma": TEMPORAL_GAMMA,
            "temporal_formula": "k(t)=0.025+0.025*2*(t-0.25)",
            "temporal_gains_at_registered_times": list(TEMPORAL_GAINS),
            "preregistered_gamma_grid": list(PREREGISTERED_GAMMAS),
            "gamma_values_outside_preregistered_grid_evaluated": False,
            "transform_rebalanced_on_st051b": False,
            "static_low_reference_energy_scale": float(low_scale),
            "static_high_reference_energy_scale": float(high_scale),
            "temporal_reference_energy_scale": float(low_scale),
            "static_low_raw_reference_energy": float(low_raw_energy),
            "static_high_raw_reference_energy": float(high_raw_energy),
        },
        "static_low_material_paths": low_paths,
        "static_high_material_paths": high_paths,
        "temporal_material_paths": temporal_paths,
        "static_low_by_seed_radius": low_by_radius,
        "static_high_by_seed_radius": high_by_radius,
        "temporal_by_seed_radius": temporal_by_radius,
        "temporal_vs_static_low": _compare(temporal_paths, low_paths),
        "temporal_vs_static_high": _compare(temporal_paths, high_paths),
        "temporal_vs_static_low_by_seed_radius": _compare_by_radius(temporal_by_radius, low_by_radius),
        "temporal_vs_static_high_by_seed_radius": _compare_by_radius(temporal_by_radius, high_by_radius),
        "static_high_vs_static_low": _compare(high_paths, low_paths),
        "comparison_scope": (
            "Descriptive same-contract material-trajectory comparison only. Static .05 and "
            "temporal gamma=.025 are not production-selected, and no public-image numeric target is fitted."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    _audit_truth(report)
    report["report_sha256"] = _hash(report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", required=True, type=Path)
    parser.add_argument("--agent9-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    report = build_report(args.candidate_root, args.agent9_root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
