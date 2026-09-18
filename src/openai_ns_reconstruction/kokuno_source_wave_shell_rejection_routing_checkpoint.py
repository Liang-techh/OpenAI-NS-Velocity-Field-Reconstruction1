"""Kokuno Agent-5 routing for the source-wave-shell Cartesian rejection.

This increment adds no new Kokuno mathematics. It consumes Agent 4 PR #371's
independent source-grid/source-wave-shell Cartesian audit in executable ancestry,
keeps its preregistered local rejection intact, and records the newest Agent-1
and Agent-3 work only as sibling provenance.

The local source-shell audit is not the formal Navier--Stokes gate. Formal
held-out normalized momentum max/L2 remain <=1e-3 and divergence max/L2
<=1e-5. No residual-defined forcing is allowed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_source_wave_shell_cartesian_independent import (
    build_report as build_agent4_report,
)

SCHEMA = "kokuno-agent5-source-wave-shell-rejection-routing-checkpoint-v19"
TASK_ID = "KOKUNO-A5-SOURCE-WAVE-SHELL-REJECTION-ROUTING-019"
BASE_AGENT4_PR = 371
BASE_AGENT4_HEAD = "917db15e5b72a04662cf78cd0345b32a9adeb51c"

FIXED_GATES = {
    "held_out_normalized_full_momentum_max": 1.0e-3,
    "held_out_normalized_full_momentum_l2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_l2": 1.0e-5,
    "changed": False,
}

ST006_REFERENCE = {
    "candidate_sha256": "6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3",
    "validation_seed": 9172801,
    "held_out_cartesian_points": 4096,
    "validation_times": 6,
    "finest_spatial_step": 0.005,
    "momentum_sampled_max": 0.1082289305112118,
    "volume_l2": 0.10758432876230622,
    "pde_validated": False,
}

UPSTREAM = {
    "agent1": {
        "pr": 368,
        "head_sha": "9a8042c89b0d10b94a17df81605338b9aee0f465",
        "standard_run": 35329961966,
        "outer_base_schedule_ready": True,
        "i2_interval_embedded_in_power_law_base_stage": True,
        "cone_modulation_completed": False,
        "all_i1_i2_i3_i4_overlays_completed": False,
        "global_matched_pressure_ready": False,
        "leading_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent2": {
        "pr": 369,
        "head_sha": "6050ea435536d01600d90a9dc02ebbab27293b94",
        "dedicated_run": 35330252022,
        "standard_run": 35330147697,
        "source_wave_shell_wrapper_ready": True,
        "slow_squared_partition_family_contract_ready": True,
        "concrete_source_partition_bumps_reconstructed": False,
        "source_actual_partition_labels_instantiated": False,
        "source_actual_background_path_instantiated": False,
        "public_source_oscillatory_xyz_t_velocity_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent3": {
        "latest_pr": 370,
        "latest_head_sha": "251ad3f3084486e6a0c2fb0e4ef9b617144259a9",
        "last_verified_pr": 357,
        "last_verified_head_sha": "169558a6994a255119d6af97430e717a221c034c",
        "radial_force_side_effect_required": True,
        "finite_correction_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "pr": BASE_AGENT4_PR,
        "head_sha": BASE_AGENT4_HEAD,
        "dedicated_run": 35331769250,
        "standard_run": 35331769237,
        "artifact_id": 10541526034,
        "artifact_digest": "sha256:5f24d5b153b8d501b4c8f2cb312604392d1a922dd42723c15cf1644fcf1c8a92",
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5": {
        "pr": 362,
        "head_sha": "289e49a119d96266c87a68429574c274703442d2",
        "consumed_in_executable_ancestry": False,
    },
}

FUTURE_CANDIDATE_ARTIFACT_CONTRACT = {
    "schema_name": "kokuno-candidate-artifact-v1-reserved",
    "status": "reserved_not_instantiated",
    "required_provenance": [
        "source_repository",
        "source_commit_or_version",
        "source_document",
        "autonomous_choices",
        "stage_artifact_hashes",
    ],
    "required_public_api": {
        "velocity": "velocity(points, t) -> (..., 3)",
        "pressure": "pressure(points, t) -> (...)",
        "forcing": "forcing(points, t) -> (..., 3) under the fixed/restricted contract only",
        "save_load": "deterministic metadata plus hashed numeric payloads",
    },
    "required_stage_receipts": ["leading", "oscillatory", "correction", "independent_validation"],
}

STATES = {
    "outer_base_backbone_ready": True,
    "slow_squared_partition_family_contract_ready": True,
    "source_wave_shell_independent_cartesian_preflight_assessed": True,
    "source_wave_shell_independent_cartesian_preflight_passed": False,
    "leading_ready": False,
    "oscillatory_machinery_ready": True,
    "public_source_oscillatory_xyz_t_velocity_ready": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "correction_ready": False,
    "candidate_artifact_instantiated": False,
    "complete_kokuno_composite_velocity_ready": False,
    "velocity_export_ready": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _validate_agent4_report(report: dict[str, Any]) -> None:
    if report.get("task_id") != "KOKUNO-A4-SOURCE-WAVE-SHELL-CARTESIAN-AUDIT-019":
        raise ValueError("unexpected Agent-4 source-wave-shell task")
    dependency = report.get("dependency", {})
    if dependency.get("agent2_parent_pr") != 369:
        raise ValueError("Agent-4 audit is not bound to Agent-2 PR #369")
    formal = report.get("formal_project_gates", {})
    if formal.get("normalized_momentum_max_l2") != 1.0e-3:
        raise ValueError("formal momentum gate changed")
    if formal.get("divergence_max_l2") != 1.0e-5:
        raise ValueError("formal divergence gate changed")
    if formal.get("formal_full_domain_pde_gate_assessed") is not False:
        raise ValueError("local source-shell audit cannot be formal PDE validation")
    truth = report.get("truth_boundary", {})
    if truth.get("pde_validated") is not False:
        raise ValueError("local source-shell audit cannot promote pde_validated")

    guards = report.get("local_guards", {})
    thresholds = guards.get("thresholds", {})
    if thresholds.get("finest_curl_relative_rms_max") != 2.0e-4:
        raise ValueError("Agent-4 local curl guard changed")
    if thresholds.get("finest_divergence_max_abs_max") != 2.0e-4:
        raise ValueError("Agent-4 local divergence guard changed")
    if thresholds.get("minimum_refinement_ratio") != 2.5:
        raise ValueError("Agent-4 refinement guard changed")
    if guards.get("same_thresholds_as_agent4_prs_353_360") is not True:
        raise ValueError("historical local guards were not preserved")
    if guards.get("source_wave_shell_independent_cartesian_preflight_passed") is not False:
        raise ValueError("expected source-wave-shell local rejection is absent")
    checks = guards.get("checks", {})
    if checks.get("divergence_refinement_ratio") is not False:
        raise ValueError("expected divergence-refinement rejection is absent")
    other_checks = {key: value for key, value in checks.items() if key != "divergence_refinement_ratio"}
    if not other_checks or not all(value is True for value in other_checks.values()):
        raise ValueError("source-shell rejection is not isolated to divergence refinement")


def _finest(report: dict[str, Any]) -> dict[str, Any]:
    ladder = report.get("resolution_ladder")
    if not isinstance(ladder, list) or len(ladder) != 3:
        raise ValueError("unexpected Agent-4 source-shell resolution ladder")
    row = ladder[-1]
    return {
        "step": float(row["step"]),
        "sample_count": int(row["sample_count"]),
        "curl_relative_rms": float(row["curl_relative_rms"]),
        "divergence_max_abs": float(row["divergence_max_abs"]),
        "divergence_rms_abs": float(row["divergence_rms_abs"]),
        "regions": row["regions"],
    }


def build_checkpoint(report: dict[str, Any] | None = None) -> dict[str, Any]:
    if report is None:
        report = build_agent4_report()
    _validate_agent4_report(report)
    finest = _finest(report)
    divergence_ratios = [float(x) for x in report["refinement"]["divergence_rms_coarse_to_fine_ratios"]]
    curl_ratios = [float(x) for x in report["refinement"]["curl_relative_rms_coarse_to_fine_ratios"]]
    failed_checks = [key for key, value in report["local_guards"]["checks"].items() if not value]

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base": {"agent4_pr": BASE_AGENT4_PR, "agent4_head": BASE_AGENT4_HEAD},
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "upstream": UPSTREAM,
        "future_candidate_artifact_contract": FUTURE_CANDIDATE_ARTIFACT_CONTRACT,
        "agent4_source_wave_shell_audit": {
            "seed": int(report["seed"]),
            "operator": report["operator"],
            "sampling": report["sampling"],
            "finest": finest,
            "curl_refinement_ratios": curl_ratios,
            "divergence_refinement_ratios": divergence_ratios,
            "mutation": report["mutation"],
            "local_guards": report["local_guards"],
            "failed_checks": failed_checks,
            "source_wave_shell_independent_cartesian_preflight_assessed": True,
            "source_wave_shell_independent_cartesian_preflight_passed": False,
            "formal_full_domain_verdict": None,
            "directly_comparable_to_ST006": False,
        },
        "pipeline_frontier": [
            {
                "stage": "source_profile_ingest",
                "ready": False,
                "status": "partial",
                "blockers": ["concrete source partition bumps/labels", "actual positive-order background", "remaining leading cone/overlays"],
            },
            {
                "stage": "leading_candidate",
                "ready": False,
                "status": "blocked",
                "blockers": ["cone modulation", "I1/I2/I3/I4 composition", "global matched pressure"],
            },
            {
                "stage": "oscillatory_augmentation",
                "ready": False,
                "status": "blocked_by_local_preflight",
                "blockers": ["Agent-4 #371 divergence refinement ratio < 2.5", "concrete source partition/background", "public actual-source Q-scaled velocity"],
            },
            {
                "stage": "mean_radial_corrections",
                "ready": False,
                "status": "blocked",
                "blockers": ["genuinely independent second public covariance column"],
            },
            {"stage": "candidate_artifact", "ready": False, "status": "blocked"},
            {"stage": "independent_full_pde_validation", "ready": False, "status": "waiting_for_candidate"},
            {"stage": "report_and_export", "ready": False, "status": "blocked"},
        ],
        "routing": {
            "new_fact": (
                "Agent-4 #371 assessed the actual source-grid/source-wave-shell wrapper and retained a local rejection: every frozen local check passes except the divergence refinement-ratio guard."
            ),
            "numerical_signature": (
                "Absolute finest divergence is already about 2.19e-10 and curl relative RMS about 1.32e-12, while divergence RMS refines about 11.79x then only 1.19x; this is consistent with reaching a finite-difference/roundoff floor but is not promoted without a passing independent guard."
            ),
            "do_not_do": [
                "do not weaken the 2e-4 absolute guards or the 2.5 refinement-ratio guard after observing #371",
                "do not call the tiny absolute divergence a pass while the frozen refinement check fails",
                "do not bind concrete source partition/background into a promoted public oscillatory field until the local source-shell preflight passes",
                "do not rerun Agent-3 finite correction before a genuinely independent public second column exists",
                "do not instantiate the final candidate artifact or export state before a global velocity/pressure/restricted-forcing composite exists",
            ],
            "shortest_next_closure": [
                "Agent 4: isolate #371's divergence-refinement failure with a preregistered roundoff-aware independent operator or step ladder while keeping the same 2e-4 absolute and 2.5 refinement guards",
                "Agent 2: only after that independent pass, bind concrete source partition values/derivatives plus actual positive-order background and expose Q-scaled velocity(x,y,z,t) and a genuine second covariance column",
                "Agent 3: screen that exact second column, retain radial d_z sigma_1, then run the finite correction cycle only if frozen rank/bounded-inverse guards pass",
                "Agent 1: compose #368 outer base with cone modulation and I1/I2/I3/I4 overlays and provide global matched pressure/admissible fixed forcing",
                "Agent 5: instantiate the reserved candidate artifact only after all stage inputs are executable together",
                "Agent 4: run the formal held-out 1e-3 momentum and 1e-5 divergence gate on the frozen complete candidate",
            ],
        },
        "baseline_vs_kokuno": [
            {
                "name": "ST006",
                "scope": "held-out full-domain normalized NS momentum",
                "momentum_sampled_max": ST006_REFERENCE["momentum_sampled_max"],
                "volume_l2": ST006_REFERENCE["volume_l2"],
                "pde_validated": False,
            },
            {
                "name": "Kokuno source-wave-shell seam #371",
                "scope": "local source-grid/source-shell Cartesian curl/divergence implementation audit",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "directly_comparable_to_ST006": False,
                "local_preflight_passed": False,
                "failed_checks": failed_checks,
            },
        ],
        "states": STATES,
        "truth_boundary": {
            "threshold_relaxed": False,
            "surrogate_data_promoted": False,
            "free_residual_defined_forcing_used": False,
            "kokuno_replay_used_as_independent_pde_validation": False,
            "local_absolute_error_used_to_override_failed_refinement_guard": False,
            "open_pr_sibling_provenance_laundered_into_executable_ancestry": False,
            "candidate_artifact_fabricated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }
    payload["checkpoint_sha256"] = _sha(payload)
    return payload


def validate_checkpoint(payload: dict[str, Any]) -> dict[str, Any]:
    if payload.get("schema") != SCHEMA or payload.get("task_id") != TASK_ID:
        raise ValueError("unexpected Agent-5 source-shell routing schema/task")
    signature = payload.get("checkpoint_sha256")
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    if signature != _sha(unsigned):
        raise ValueError("checkpoint SHA mismatch")
    if payload.get("fixed_gates") != FIXED_GATES:
        raise ValueError("fixed gates changed")
    if payload.get("st006_reference") != ST006_REFERENCE:
        raise ValueError("ST006 reference changed")
    if payload.get("states") != STATES:
        raise ValueError("routing states changed")
    if payload.get("future_candidate_artifact_contract") != FUTURE_CANDIDATE_ARTIFACT_CONTRACT:
        raise ValueError("future candidate artifact contract changed")
    audit = payload.get("agent4_source_wave_shell_audit", {})
    if audit.get("source_wave_shell_independent_cartesian_preflight_assessed") is not True:
        raise ValueError("source-wave-shell independent audit disappeared")
    if audit.get("source_wave_shell_independent_cartesian_preflight_passed") is not False:
        raise ValueError("source-wave-shell rejection was falsely promoted")
    if audit.get("failed_checks") != ["divergence_refinement_ratio"]:
        raise ValueError("source-wave-shell rejection reason changed")
    if audit.get("directly_comparable_to_ST006") is not False:
        raise ValueError("local source-shell audit cannot be ST006-comparable")
    stages = payload.get("pipeline_frontier", [])
    if not stages or any(stage.get("ready") is not False for stage in stages):
        raise ValueError("an incomplete pipeline stage was promoted")
    truth = payload.get("truth_boundary", {})
    for key, value in truth.items():
        if value is not False:
            raise ValueError(f"truth boundary promoted: {key}")
    return payload


def write_bundle(output_dir: str | Path) -> dict[str, Path]:
    root = Path(output_dir)
    if root.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {root}")
    root.mkdir(parents=True)
    report = build_agent4_report()
    checkpoint = build_checkpoint(report)
    validate_checkpoint(checkpoint)
    report_path = root / "agent4_source_wave_shell_cartesian_audit.json"
    checkpoint_path = root / "source_wave_shell_rejection_routing_checkpoint.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    checkpoint_path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema": SCHEMA,
        "files": {
            report_path.name: hashlib.sha256(report_path.read_bytes()).hexdigest(),
            checkpoint_path.name: hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
        },
        "checkpoint_payload_sha256": checkpoint["checkpoint_sha256"],
        "pde_validated": False,
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"report": report_path, "checkpoint": checkpoint_path, "manifest": manifest_path}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    paths = write_bundle(args.output_dir)
    print(json.dumps({key: str(value) for key, value in paths.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
