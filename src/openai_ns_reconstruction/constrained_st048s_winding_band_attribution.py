"""Attribute the frozen ST048-S winding deficit by the already-registered seed radii.

This is an evidence-reuse / visualization-governance increment. It consumes the
immutable CR-A9-040 ST006 material-path receipt and CR-A9-048 ST048-S temporal-
Piola/global-swirl receipt, verifies their frozen 48-path contract and truth
boundaries, and stratifies the existing path observables by the preregistered seed
radii 0.6, 0.9 and 1.2.

It does not integrate a new trajectory, select a production radial profile, infer
an OpenAI target, fit pressure/forcing, or evaluate a PDE residual.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

TASK_ID = "CR-A9-049"
SCHEMA = "st048s_winding_band_attribution_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"

ST006_TASK_ID = "CR-A9-040"
ST006_SCHEMA = "st006_material_path_observables_v1"
ST006_HEAD = "3911718884f021f0c0e3ba41ed694b9688764084"
ST006_RUN_ID = 35322722094
ST006_ARTIFACT_ID = 10536949231
ST006_ARTIFACT_DIGEST = "sha256:e73165cbc4ae75e2d9b437ed1b8606a4af9c3e7439ffb50d65c315e33378c88d"
ST006_REPORT_SHA256 = "2595c8dafd0716b3b07cb9dc761b112e2178bf44200b9c50eff3c000c6166772"
ST006_CANDIDATE_SHA256 = "6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3"

SWIRL_TASK_ID = "CR-A9-048"
SWIRL_SCHEMA = "st048s_piola_swirl_material_path_v1"
SWIRL_HEAD = "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"
SWIRL_RUN_ID = 35366730180
SWIRL_ARTIFACT_ID = 10556164980
SWIRL_ARTIFACT_DIGEST = "sha256:15488084eea04a0c4428bc07fbb87dad8dc6d11e719817715ba0c16314f644a0"
SWIRL_REPORT_SHA256 = "9463a67cb556ccd22ff24f776dca0af51cc3b7d3363e8f8b700a8fba434f0f6c"
ST048S_RAW_SHA256 = "6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09"

SEED_RADII = (0.6, 0.9, 1.2)
SEED_Z = (-0.3, 0.3)
SEED_ANGLES = 8
PATHS_PER_RADIUS = 16
TIME_INTERVAL = (0.25, 0.75)
OUTPUT_SAMPLES = 33
SOLVER_METHOD = "DOP853"
SOLVER_RTOL = 1.0e-9
SOLVER_ATOL = 1.0e-11
SOLVER_MAX_STEP = 0.01

EXTERNAL_METHOD = {
    "repository": "scipy/scipy",
    "commit": "eff78058ed0d4cb8e1f0b5c5585d62d712e948d2",
    "api": "scipy.integrate.solve_ivp / DOP853",
    "license": "BSD-3-Clause",
    "classification": "direct migration / public API only (inherited receipts; no new call in CR-A9-049)",
    "copied_upstream_implementation": False,
}

TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "new_trajectory_integration_performed": False,
    "production_candidate_selected": False,
    "production_radial_profile_selected": False,
    "openai_numeric_target_inferred": False,
    "hidden_openai_time_camera_seed_or_velocity_used": False,
    "visual_acceptance_threshold_defined": False,
    "comparison_is_candidate_side_descriptive_routing_only": True,
    "pressure_or_force_fit": False,
    "held_out_pde_residual_evaluated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_hash(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def _load_json(path: str | Path) -> dict[str, Any]:
    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError("receipt must be a JSON object")
    return data


def _same_float_sequence(actual: Any, expected: tuple[float, ...]) -> bool:
    try:
        values = tuple(float(x) for x in actual)
    except (TypeError, ValueError):
        return False
    return values == expected


def _validate_st006(receipt: dict[str, Any]) -> None:
    if receipt.get("task_id") != ST006_TASK_ID or receipt.get("schema") != ST006_SCHEMA:
        raise ValueError("unexpected ST006 receipt identity")
    if receipt.get("report_sha256") != ST006_REPORT_SHA256:
        raise ValueError("ST006 report SHA drift")
    if receipt.get("candidate") != "ST006" or receipt.get("candidate_sha256") != ST006_CANDIDATE_SHA256:
        raise ValueError("ST006 candidate identity drift")
    contract = receipt.get("registered_contract", {})
    if not _same_float_sequence(contract.get("seed_radii"), SEED_RADII):
        raise ValueError("ST006 seed-radii contract drift")
    if not _same_float_sequence(contract.get("seed_z"), SEED_Z):
        raise ValueError("ST006 seed-z contract drift")
    if contract.get("seed_angles") != SEED_ANGLES or contract.get("seed_count") != 48:
        raise ValueError("ST006 seed-count contract drift")
    if not _same_float_sequence(contract.get("time_interval"), TIME_INTERVAL):
        raise ValueError("ST006 time contract drift")
    if contract.get("output_samples") != OUTPUT_SAMPLES:
        raise ValueError("ST006 output-sample contract drift")
    solver = contract.get("solver", {})
    if (
        solver.get("method") != SOLVER_METHOD
        or float(solver.get("rtol")) != SOLVER_RTOL
        or float(solver.get("atol")) != SOLVER_ATOL
        or float(solver.get("max_step")) != SOLVER_MAX_STEP
    ):
        raise ValueError("ST006 solver contract drift")
    boundary = receipt.get("interpretation_boundary", {})
    required_false = (
        "pde_validated", "visualization_ready", "visual_correspondence_verified",
        "openai_field_identified", "paper_exact", "blowup_proved",
        "visual_acceptance_threshold_defined", "image_fit_used",
        "hidden_time_alignment_used", "camera_registration_used", "openai_seed_locations_used",
    )
    if any(boundary.get(k) is not False for k in required_false):
        raise ValueError("ST006 truth boundary was promoted")
    if boundary.get("all_sampled_points_inside_registered_box") is not True:
        raise ValueError("ST006 path contract did not remain inside registered box")
    if len(receipt.get("per_path", [])) != 48:
        raise ValueError("ST006 receipt must contain exactly 48 paths")


def _validate_swirl(receipt: dict[str, Any]) -> None:
    if receipt.get("task_id") != SWIRL_TASK_ID or receipt.get("schema") != SWIRL_SCHEMA:
        raise ValueError("unexpected swirl receipt identity")
    if receipt.get("report_sha256") != SWIRL_REPORT_SHA256:
        raise ValueError("swirl report SHA drift")
    candidate = receipt.get("candidate", {})
    if (
        candidate.get("candidate_id") != "ST048-S"
        or candidate.get("original_raw_candidate_sha256") != ST048S_RAW_SHA256
        or candidate.get("pde_validated") is not False
        or candidate.get("source_correspondence_verified") is not False
    ):
        raise ValueError("ST048-S candidate identity/truth drift")
    contract = receipt.get("frozen_contract", {})
    if not _same_float_sequence(contract.get("seed_radii"), SEED_RADII):
        raise ValueError("swirl seed-radii contract drift")
    if not _same_float_sequence(contract.get("seed_z"), SEED_Z):
        raise ValueError("swirl seed-z contract drift")
    if contract.get("seed_angles") != SEED_ANGLES or contract.get("path_count") != 48:
        raise ValueError("swirl seed-count contract drift")
    if not _same_float_sequence(contract.get("time_interval"), TIME_INTERVAL):
        raise ValueError("swirl time contract drift")
    if (
        contract.get("output_samples") != OUTPUT_SAMPLES
        or contract.get("solver_method") != SOLVER_METHOD
        or float(contract.get("solver_rtol")) != SOLVER_RTOL
        or float(contract.get("solver_atol")) != SOLVER_ATOL
        or float(contract.get("solver_max_step")) != SOLVER_MAX_STEP
    ):
        raise ValueError("swirl numerical contract drift")
    truth = receipt.get("truth_boundary", {})
    required_false = (
        "canonical_velocity_changed", "production_candidate_selected", "production_kappa_selected",
        "held_out_pde_residual_evaluated", "visual_acceptance_threshold_defined",
        "hidden_openai_time_camera_seed_or_velocity_used", "visualization_ready",
        "visual_correspondence_verified", "pde_validated", "source_correspondence_verified",
        "paper_exact", "openai_field_identified", "blowup_proved",
    )
    if any(truth.get(k) is not False for k in required_false):
        raise ValueError("swirl receipt truth boundary was promoted")
    if truth.get("comparison_is_descriptive_not_acceptance") is not True:
        raise ValueError("swirl comparison lost descriptive-only status")
    measurements = receipt.get("measurements", {})
    for key in ("temporal_piola_kappa_0", "temporal_piola_kappa_005"):
        if len(measurements.get(key, {}).get("per_path", [])) != 48:
            raise ValueError(f"{key} must contain exactly 48 paths")


def _normalized_path(row: dict[str, Any], *, st006: bool = False) -> dict[str, float | int]:
    seed = row.get("seed", {})
    radius = float(seed.get("radius"))
    z = float(seed.get("z"))
    angle_index = int(seed.get("angle_index"))
    if radius not in SEED_RADII or z not in SEED_Z or not 0 <= angle_index < SEED_ANGLES:
        raise ValueError("path seed outside frozen contract")
    turns = abs(float(row["signed_turns"])) if st006 else float(row["absolute_turns"])
    if turns < 0.0:
        raise ValueError("absolute turns cannot be negative")
    return {
        "radius": radius,
        "z": z,
        "angle_index": angle_index,
        "absolute_turns": turns,
        "radius_change": float(row["radius_change"]),
        "abs_z_change": float(row["abs_z_change"]),
    }


def _keyed_paths(rows: list[dict[str, Any]], *, st006: bool = False) -> dict[tuple[float, float, int], dict[str, float | int]]:
    result: dict[tuple[float, float, int], dict[str, float | int]] = {}
    for raw in rows:
        row = _normalized_path(raw, st006=st006)
        key = (float(row["radius"]), float(row["z"]), int(row["angle_index"]))
        if key in result:
            raise ValueError("duplicate frozen path seed")
        result[key] = row
    expected = {(r, z, a) for r in SEED_RADII for a in range(SEED_ANGLES) for z in SEED_Z}
    if set(result) != expected:
        raise ValueError("frozen path seed population is incomplete or changed")
    return result


def _mean(values: list[float]) -> float:
    if not values:
        raise ValueError("cannot average empty values")
    return sum(values) / len(values)


def _radius_summary(paths: dict[tuple[float, float, int], dict[str, float | int]], radius: float) -> dict[str, float | int]:
    rows = [row for key, row in paths.items() if key[0] == radius]
    if len(rows) != PATHS_PER_RADIUS:
        raise ValueError("every frozen radius must contain exactly 16 paths")
    return {
        "path_count": len(rows),
        "mean_absolute_turns": _mean([float(x["absolute_turns"]) for x in rows]),
        "mean_radius_change": _mean([float(x["radius_change"]) for x in rows]),
        "mean_abs_z_change": _mean([float(x["abs_z_change"]) for x in rows]),
        "inward_path_count": sum(float(x["radius_change"]) < 0.0 for x in rows),
    }


def build_report(st006_receipt: dict[str, Any], swirl_receipt: dict[str, Any]) -> dict[str, Any]:
    _validate_st006(st006_receipt)
    _validate_swirl(swirl_receipt)

    st006_paths = _keyed_paths(st006_receipt["per_path"], st006=True)
    k0_paths = _keyed_paths(swirl_receipt["measurements"]["temporal_piola_kappa_0"]["per_path"])
    k005_paths = _keyed_paths(swirl_receipt["measurements"]["temporal_piola_kappa_005"]["per_path"])
    if set(st006_paths) != set(k0_paths) or set(k0_paths) != set(k005_paths):
        raise ValueError("upstream receipts do not share the exact same path seeds")

    per_radius: dict[str, Any] = {}
    total_st006_turns = 0.0
    total_child_turns = 0.0
    raw_deficits: dict[float, float] = {}
    for radius in SEED_RADII:
        st006 = _radius_summary(st006_paths, radius)
        k0 = _radius_summary(k0_paths, radius)
        child = _radius_summary(k005_paths, radius)
        if st006["inward_path_count"] != PATHS_PER_RADIUS or k0["inward_path_count"] != PATHS_PER_RADIUS or child["inward_path_count"] != PATHS_PER_RADIUS:
            raise ValueError("expected all frozen paths to remain inward at each radius")
        st006_turn = float(st006["mean_absolute_turns"])
        k0_turn = float(k0["mean_absolute_turns"])
        child_turn = float(child["mean_absolute_turns"])
        deficit = st006_turn - child_turn
        raw_deficits[radius] = deficit * PATHS_PER_RADIUS
        total_st006_turns += st006_turn * PATHS_PER_RADIUS
        total_child_turns += child_turn * PATHS_PER_RADIUS
        per_radius[f"{radius:.1f}"] = {
            "seed_radius": radius,
            "st006": st006,
            "st048s_temporal_piola_kappa_0": k0,
            "st048s_temporal_piola_kappa_005": child,
            "kappa_005_vs_kappa_0": {
                "mean_winding_ratio": child_turn / k0_turn,
                "mean_winding_percent_change": 100.0 * (child_turn / k0_turn - 1.0),
                "radial_contraction_magnitude_ratio": abs(float(child["mean_radius_change"])) / abs(float(k0["mean_radius_change"])),
                "radial_contraction_magnitude_percent_change": 100.0 * (abs(float(child["mean_radius_change"])) / abs(float(k0["mean_radius_change"])) - 1.0),
            },
            "kappa_005_vs_st006": {
                "mean_winding_difference_turns": child_turn - st006_turn,
                "mean_winding_ratio": child_turn / st006_turn,
                "mean_winding_percent_change": 100.0 * (child_turn / st006_turn - 1.0),
                "radial_contraction_magnitude_ratio": abs(float(child["mean_radius_change"])) / abs(float(st006["mean_radius_change"])),
                "radial_contraction_magnitude_percent_change": 100.0 * (abs(float(child["mean_radius_change"])) / abs(float(st006["mean_radius_change"])) - 1.0),
            },
        }

    total_deficit = total_st006_turns - total_child_turns
    if total_deficit <= 0.0:
        raise ValueError("expected frozen child to retain a positive net winding deficit versus ST006")
    for radius in SEED_RADII:
        per_radius[f"{radius:.1f}"]["kappa_005_vs_st006"]["net_winding_deficit_contribution_fraction"] = raw_deficits[radius] / total_deficit

    deficit_radius = max(SEED_RADII, key=lambda r: raw_deficits[r])
    gain_percentages = [per_radius[f"{r:.1f}"]["kappa_005_vs_kappa_0"]["mean_winding_percent_change"] for r in SEED_RADII]
    routing = {
        "largest_candidate_side_winding_deficit_radius": deficit_radius,
        "inner_r06_share_of_net_deficit": per_radius["0.6"]["kappa_005_vs_st006"]["net_winding_deficit_contribution_fraction"],
        "mid_r09_share_of_net_deficit": per_radius["0.9"]["kappa_005_vs_st006"]["net_winding_deficit_contribution_fraction"],
        "outer_r12_share_of_net_deficit": per_radius["1.2"]["kappa_005_vs_st006"]["net_winding_deficit_contribution_fraction"],
        "outer_r12_is_winding_surplus_vs_st006": per_radius["1.2"]["kappa_005_vs_st006"]["mean_winding_difference_turns"] > 0.0,
        "global_kappa_winding_gain_percent_range_across_radii": [min(gain_percentages), max(gain_percentages)],
        "interpretation": (
            "Under the unchanged autonomous 48-path contract and using ST006 only as a candidate-side baseline, "
            "the remaining mean-winding deficit of the kappa=.05 child is concentrated at the inner r=.6 seeds, "
            "secondarily at r=.9, while r=1.2 is already a surplus. The global kappa gain is nearly uniform across "
            "the three seed radii. This routes any later swirl-shape screen toward inner/mid redistribution rather "
            "than adding more outer-ring gain; it is not an OpenAI-derived target or a production coefficient."
        ),
    }

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_main_sha": BASE_MAIN_SHA,
        "upstream_receipts": {
            "st006": {
                "task_id": ST006_TASK_ID,
                "head_sha": ST006_HEAD,
                "workflow_run_id": ST006_RUN_ID,
                "artifact_id": ST006_ARTIFACT_ID,
                "artifact_digest": ST006_ARTIFACT_DIGEST,
                "report_sha256": ST006_REPORT_SHA256,
                "candidate_sha256": ST006_CANDIDATE_SHA256,
            },
            "st048s_swirl": {
                "task_id": SWIRL_TASK_ID,
                "head_sha": SWIRL_HEAD,
                "workflow_run_id": SWIRL_RUN_ID,
                "artifact_id": SWIRL_ARTIFACT_ID,
                "artifact_digest": SWIRL_ARTIFACT_DIGEST,
                "report_sha256": SWIRL_REPORT_SHA256,
                "candidate_raw_sha256": ST048S_RAW_SHA256,
            },
        },
        "external_method": EXTERNAL_METHOD,
        "analysis_contract": {
            "new_integration": False,
            "same_upstream_path_population_required": True,
            "seed_radii": list(SEED_RADII),
            "seed_z": list(SEED_Z),
            "seed_angles": SEED_ANGLES,
            "paths_per_radius": PATHS_PER_RADIUS,
            "time_interval": list(TIME_INTERVAL),
            "output_samples": OUTPUT_SAMPLES,
            "grouping_variable": "autonomous seed radius from frozen receipts",
            "comparison_baseline": "ST006 candidate-side material paths; not OpenAI truth",
        },
        "per_radius": per_radius,
        "routing": routing,
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    report["report_sha256"] = _canonical_hash(report)
    return report


def audit_truth_boundary(report: dict[str, Any]) -> None:
    if report.get("task_id") != TASK_ID or report.get("schema") != SCHEMA:
        raise ValueError("task/schema drift")
    if report.get("truth_boundary") != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drift")
    truth = report["truth_boundary"]
    for key, value in truth.items():
        if key == "comparison_is_candidate_side_descriptive_routing_only":
            if value is not True:
                raise ValueError("comparison must stay descriptive-routing only")
        elif value is not False:
            raise ValueError(f"forbidden truth promotion: {key}")
    if report["analysis_contract"].get("new_integration") is not False:
        raise ValueError("this task may not perform a new trajectory integration")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--st006-report", required=True)
    parser.add_argument("--swirl-report", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    report = build_report(_load_json(args.st006_report), _load_json(args.swirl_report))
    audit_truth_boundary(report)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    compact = {
        "task_id": TASK_ID,
        "report_sha256": report["report_sha256"],
        "largest_deficit_radius": report["routing"]["largest_candidate_side_winding_deficit_radius"],
        "r06_deficit_share": report["routing"]["inner_r06_share_of_net_deficit"],
        "r09_deficit_share": report["routing"]["mid_r09_share_of_net_deficit"],
        "r12_deficit_share": report["routing"]["outer_r12_share_of_net_deficit"],
        "r12_surplus": report["routing"]["outer_r12_is_winding_surplus_vs_st006"],
    }
    print(json.dumps(compact, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
