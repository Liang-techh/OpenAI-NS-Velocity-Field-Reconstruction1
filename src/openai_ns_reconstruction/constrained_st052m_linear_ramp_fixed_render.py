"""Target-free fixed 3-D render of the exact ST052-M linear temporal child.

CR-A9-066 compares Agent-7 PR #587's single frozen linear activation against
its unchanged redistributed control, while retaining the exact static #559
endpoint as a reference. It consumes exact source checkouts and Agent-9's
already-audited CR-A9-064 render machinery. No OpenAI image is fitted, no
visual acceptance threshold is defined, and no pressure/forcing/PDE receipt is
transferred.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType
from typing import Any

import numpy as np

TASK_ID = "CR-A9-066"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"
TEMPORAL_PR = 587
TEMPORAL_HEAD = "0b93819095f6c8576a7571bdc2d2fbef4154944d"
TEMPORAL_CONVERGENCE_PR = 594
TEMPORAL_CONVERGENCE_HEAD = "b8d45d5d04a03e57739b0b3847c0994f4ade8642"
STATIC_CHILD_PR = 559
STATIC_CHILD_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
STATIC_RENDER_PR = 583
STATIC_RENDER_HEAD = "ee80de11775bd45157d114000378b6ee3422ff46"
ST052_PR = 508
ST052_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
REDISTRIBUTION_PR = 528
REDISTRIBUTION_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"
EXPECTED_REDISTRIBUTION_SCALE = 1.0032534663681094
EXPECTED_BETA = 0.08837490297155456
REFERENCE_TIMES = (0.25, 0.375, 0.50, 0.625, 0.75)
GRID_RESOLUTION = 33
ENDPOINT_IDENTITY_ABS_MAX = 1.0e-11

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "visual_pass_threshold_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def expected_activation(time: float) -> float:
    """Frozen autonomous activation used by PR #587."""
    time = float(time)
    if not (REFERENCE_TIMES[0] - 1.0e-12 <= time <= REFERENCE_TIMES[-1] + 1.0e-12):
        raise ValueError("time outside frozen interval")
    return 2.0 * (time - 0.25)


def effect_retention(ramp: float, static: float, *, zero_tol: float = 1.0e-14) -> float | None:
    """Return ramp/static relative-effect retention, or None for zero static effect."""
    if abs(float(static)) <= float(zero_tol):
        return None
    return float(ramp) / float(static)


def _load_file_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module at {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _load_exact_sources(source_root: Path, render_root: Path) -> tuple[ModuleType, ModuleType]:
    source_root = source_root.resolve()
    render_root = render_root.resolve()
    experiment_dir = source_root / "experiments" / "root_st052"
    source_src = source_root / "src"
    helper_path = (
        render_root
        / "src"
        / "openai_ns_reconstruction"
        / "constrained_st052m_local_swirl_fixed_render.py"
    )
    for required in (experiment_dir, source_src, helper_path):
        if not required.exists():
            raise FileNotFoundError(required)

    for value in (str(experiment_dir), str(source_src)):
        if value not in sys.path:
            sys.path.insert(0, value)
    temporal = importlib.import_module("agent7_st052m_linear_temporal_activation")
    helper = _load_file_module(helper_path, "_agent9_cr_a9_064_render_helper")

    if getattr(temporal, "SOURCE_AGENT7_HEAD", None) != STATIC_CHILD_HEAD:
        raise RuntimeError("temporal source static-child identity drift")
    if getattr(temporal, "SOURCE_MORPH_HEAD", None) is None:
        raise RuntimeError("temporal source morphology provenance missing")
    if tuple(getattr(temporal, "TIME_INTERVAL", ())) != (0.25, 0.75):
        raise RuntimeError("temporal source time interval drift")
    if abs(float(getattr(temporal, "TAPER_TAU", np.nan)) - 0.05) > 1.0e-15:
        raise RuntimeError("temporal source taper identity drift")
    if abs(float(getattr(temporal, "EXPECTED_BETA", np.nan)) - EXPECTED_BETA) > 1.0e-15:
        raise RuntimeError("temporal source beta identity drift")

    helper.REFERENCE_TIMES = REFERENCE_TIMES
    helper.GRID_RESOLUTION = GRID_RESOLUTION
    helper.TASK_ID = TASK_ID
    helper.TRUTH = dict(TRUTH)
    return temporal, helper


def _field_digest(axis: np.ndarray, field: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in (axis, field):
        value = np.ascontiguousarray(np.asarray(array, dtype="<f8"))
        digest.update(np.asarray(value.shape, dtype="<i8").tobytes())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _per_time_retention(static_report: dict[str, Any], ramp_report: dict[str, Any]) -> list[dict[str, Any]]:
    static_rows = static_report["per_time"]
    ramp_rows = ramp_report["per_time"]
    if len(static_rows) != len(ramp_rows):
        raise RuntimeError("static/ramp report length mismatch")
    metrics = (
        "mean_absolute_turns",
        "max_absolute_turns",
        "mean_axial_span",
        "mean_radial_span",
        "top_vorticity_axial_rms",
        "top_vorticity_radial_rms",
        "enstrophy_axial_rms",
        "enstrophy_radial_rms",
        "enstrophy_aspect",
        "tip_vorticity_radial_rms",
        "central_vorticity_radial_rms",
        "full_vorticity_rms",
    )
    result = []
    for static_row, ramp_row in zip(static_rows, ramp_rows, strict=True):
        if abs(float(static_row["time"]) - float(ramp_row["time"])) > 1.0e-15:
            raise RuntimeError("static/ramp time mismatch")
        static_rel = static_row["relative"]
        ramp_rel = ramp_row["relative"]
        result.append(
            {
                "time": float(ramp_row["time"]),
                "activation": expected_activation(float(ramp_row["time"])),
                "retention_ramp_over_static": {
                    metric: effect_retention(ramp_rel[metric], static_rel[metric])
                    for metric in metrics
                },
            }
        )
    return result


def run(source_root: Path, render_root: Path, output_dir: Path) -> dict[str, Any]:
    temporal, helper = _load_exact_sources(source_root, render_root)
    field, raw = temporal.morph.prior.replay_st052.reconstruct()
    redistribution_scale, _, _ = temporal.morph.prior.base.child_scale(field, raw)
    energy = temporal.morph.prior.comp.solve_energy_beta(field, raw, redistribution_scale)
    if not energy["root_exists"]:
        raise RuntimeError("frozen #559 energy root disappeared")
    beta = float(energy["beta"])
    if abs(float(redistribution_scale) - EXPECTED_REDISTRIBUTION_SCALE) > 5.0e-10:
        raise RuntimeError("redistribution scale identity drift")
    if abs(beta - EXPECTED_BETA) > 5.0e-10:
        raise RuntimeError("energy-only beta identity drift")

    control = lambda points, time: temporal.control_velocity(
        field, raw, points, time, redistribution_scale
    )
    static = lambda points, time: temporal.static_velocity(
        field, raw, points, time, redistribution_scale, beta
    )
    ramp = lambda points, time: temporal.ramp_velocity(
        field, raw, points, time, redistribution_scale, beta
    )

    axis_c, control_fields = helper.sample_velocity_grid(
        control, times=REFERENCE_TIMES, resolution=GRID_RESOLUTION
    )
    axis_s, static_fields = helper.sample_velocity_grid(
        static, times=REFERENCE_TIMES, resolution=GRID_RESOLUTION
    )
    axis_r, ramp_fields = helper.sample_velocity_grid(
        ramp, times=REFERENCE_TIMES, resolution=GRID_RESOLUTION
    )
    if not (np.array_equal(axis_c, axis_s) and np.array_equal(axis_c, axis_r)):
        raise RuntimeError("grid axis drift")

    start_identity = float(np.max(np.abs(ramp_fields[0] - control_fields[0])))
    end_identity = float(np.max(np.abs(ramp_fields[-1] - static_fields[-1])))
    if start_identity > ENDPOINT_IDENTITY_ABS_MAX:
        raise RuntimeError("ramp no longer equals control at t=.25")
    if end_identity > ENDPOINT_IDENTITY_ABS_MAX:
        raise RuntimeError("ramp no longer equals static child at t=.75")

    static_report = helper.analyze_pair(axis_c, control_fields, static_fields)
    ramp_report = helper.analyze_pair(axis_c, control_fields, ramp_fields)
    retention = _per_time_retention(static_report, ramp_report)

    output_dir.mkdir(parents=True, exist_ok=True)
    render_files = [
        Path(path).name
        for path in helper.render_pair(axis_c, control_fields, ramp_fields, output_dir)
    ]
    np.savez_compressed(output_dir / "control_grid.npz", axis=axis_c, velocity=control_fields)
    np.savez_compressed(output_dir / "static_child_grid.npz", axis=axis_c, velocity=static_fields)
    np.savez_compressed(output_dir / "linear_ramp_grid.npz", axis=axis_c, velocity=ramp_fields)

    report: dict[str, Any] = {
        "task_id": TASK_ID,
        "identities": {
            "base_main_sha": BASE_MAIN_SHA,
            "st052_pr": ST052_PR,
            "st052_head": ST052_HEAD,
            "redistribution_pr": REDISTRIBUTION_PR,
            "redistribution_head": REDISTRIBUTION_HEAD,
            "static_child_pr": STATIC_CHILD_PR,
            "static_child_head": STATIC_CHILD_HEAD,
            "temporal_pr": TEMPORAL_PR,
            "temporal_head": TEMPORAL_HEAD,
            "temporal_convergence_pr": TEMPORAL_CONVERGENCE_PR,
            "temporal_convergence_head": TEMPORAL_CONVERGENCE_HEAD,
            "static_render_pr": STATIC_RENDER_PR,
            "static_render_head": STATIC_RENDER_HEAD,
        },
        "contract": {
            "times": list(REFERENCE_TIMES),
            "grid_resolution": GRID_RESOLUTION,
            "streamline_count": 48,
            "fixed_camera_from_cr_a9_064": True,
            "top_vorticity_quantile_render_only": 0.985,
            "image_or_openai_numeric_target_used": False,
            "visual_acceptance_threshold_used": False,
            "endpoint_identity_abs_max": ENDPOINT_IDENTITY_ABS_MAX,
        },
        "source_replay": {
            "source_checkout_expected_head": TEMPORAL_HEAD,
            "render_helper_expected_head": STATIC_RENDER_HEAD,
            "redistribution_scale": float(redistribution_scale),
            "beta": beta,
            "post_transform_common_scale": 1.0,
            "activation": "g(t)=2*(t-.25)",
            "tau": "0.05*g(t)",
            "beta_time": f"{beta}*g(t)",
        },
        "endpoint_identity": {
            "ramp_vs_control_t025_max_abs": start_identity,
            "ramp_vs_static_t075_max_abs": end_identity,
        },
        "static_vs_control": static_report["per_time"],
        "ramp_vs_control": ramp_report["per_time"],
        "effect_retention": retention,
        "grid_sha256": {
            "control": _field_digest(axis_c, control_fields),
            "static": _field_digest(axis_c, static_fields),
            "ramp": _field_digest(axis_c, ramp_fields),
        },
        "render_files": render_files,
        "truth": dict(TRUTH),
    }
    body = json.dumps(report, indent=2, sort_keys=True) + "\n"
    (output_dir / "report.json").write_text(body)
    summary = {
        "task_id": TASK_ID,
        "report_sha256": hashlib.sha256(body.encode()).hexdigest(),
        "endpoint_identity": report["endpoint_identity"],
        "render_files": render_files,
        "truth": report["truth"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--render-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run(args.source_root, args.render_root, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
