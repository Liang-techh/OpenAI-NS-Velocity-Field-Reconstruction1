"""Kokuno Agent-5 cross-lane integration frontier after A1 #368 / A2 #369.

This module adds no new Kokuno mathematics.  It turns the latest independently
usable handoffs into a fail-closed pipeline manifest while preserving Agent-4's
unresolved source-wave-shell validation dependency.

The manifest is deliberately not a candidate artifact: there is still no global
matched leading pressure, no actual-source public oscillatory xyz,t velocity,
no genuinely independent second covariance column, and no promoted finite
correction cycle.  The formal held-out gates remain normalized momentum max/L2
<=1e-3 and divergence max/L2 <=1e-5.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_support_localization_resolution_routing_checkpoint import (
    FIXED_GATES,
    ST006_REFERENCE,
    build_checkpoint as build_v18_checkpoint,
    validate_checkpoint as validate_v18_checkpoint,
)

SCHEMA = "kokuno-agent5-multilane-frontier-checkpoint-v19"
TASK_ID = "KOKUNO-A5-MULTILANE-FRONTIER-019"
BASE_AGENT5_PR = 362
BASE_AGENT5_HEAD = "289e49a119d96266c87a68429574c274703442d2"

LATEST = {
    "agent1": {
        "pr": 368,
        "head_sha": "9a8042c89b0d10b94a17df81605338b9aee0f465",
        "verification": {"workflow_run": 35329961966, "conclusion": "success"},
        "outer_base_schedule_ready": True,
        "i2_interval_embedded_in_power_law_base_stage": True,
        "cone_modulation_completed": False,
        "all_i1_i2_i3_i4_overlays_completed": False,
        "global_matched_pressure_ready": False,
        "global_leading_profile_reconstructed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent2": {
        "pr": 369,
        "head_sha": "6050ea435536d01600d90a9dc02ebbab27293b94",
        "verification": {
            "dedicated_run": 35330252022,
            "standard_run": 35330147697,
            "conclusion": "success",
        },
        "source_wave_shell_wrapper_ready": True,
        "slow_squared_partition_family_contract_ready": True,
        "concrete_source_partition_bumps_reconstructed": False,
        "source_actual_partition_labels_instantiated": False,
        "source_actual_background_path_instantiated": False,
        "public_source_oscillatory_xyz_t_velocity_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "source_wave_shell_independent_cartesian_preflight_passed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent3": {
        "latest_pr": 370,
        "latest_head_sha": "251ad3f3084486e6a0c2fb0e4ef9b617144259a9",
        "base_verified_pr": 357,
        "base_verified_head_sha": "169558a6994a255119d6af97430e717a221c034c",
        "latest_spacetime_envelope_ci_complete_at_cut": False,
        "radial_force_side_effect_required": True,
        "finite_correction_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "latest_verified_pr": 360,
        "head_sha": "948200cfac9c95e5bed2cf08371c5d113f00188f",
        "dedicated_run": 35326140220,
        "standard_run": 35326140228,
        "manufactured_support_localization_axis_refinement_passed": True,
        "source_wave_shell_independent_cartesian_preflight_passed": False,
        "next_required_dependency_pr": 369,
        "consumed_via_previous_agent5_checkpoint": True,
    },
    "previous_agent5": {
        "pr": BASE_AGENT5_PR,
        "head_sha": BASE_AGENT5_HEAD,
        "consumed_in_executable_ancestry": True,
    },
}

PIPELINE_STAGES = [
    {
        "stage": "source_profile_ingest",
        "status": "partial",
        "ready": False,
        "evidence": ["A1 #368 outer/base schedule", "A2 #369 partition-family contract"],
        "blockers": [
            "concrete source partition bumps/labels",
            "actual positive-order background path",
            "remaining leading cone/overlay data",
        ],
    },
    {
        "stage": "leading_candidate",
        "status": "blocked",
        "ready": False,
        "evidence": ["A1 #368 executable outer base backbone", "independently preflighted heat repair from prior lane"],
        "blockers": ["cone modulation", "I1/I2/I3/I4 overlay composition", "global matched pressure"],
    },
    {
        "stage": "oscillatory_augmentation",
        "status": "blocked",
        "ready": False,
        "evidence": ["complete-curl machinery", "A2 #369 squared-partition family contract"],
        "blockers": [
            "A4 independent Cartesian audit on source-grid/source-wave-shell wrapper",
            "concrete source partition/background binding",
            "public actual-source Q-scaled velocity(x,y,z,t)",
        ],
    },
    {
        "stage": "mean_radial_corrections",
        "status": "blocked",
        "ready": False,
        "evidence": ["A3 #357 retained radial d_z sigma_1 side effect"],
        "blockers": ["genuinely independent second public covariance column", "A3 #370 exact-head completion at this cut"],
    },
    {
        "stage": "candidate_artifact",
        "status": "blocked",
        "ready": False,
        "blockers": ["leading_candidate", "oscillatory_augmentation", "mean_radial_corrections"],
    },
    {
        "stage": "independent_validation",
        "status": "waiting_for_candidate",
        "ready": False,
        "validator_owner": "Kokuno Agent 4",
        "fixed_momentum_gate": 1.0e-3,
        "fixed_divergence_gate": 1.0e-5,
    },
    {
        "stage": "report_and_export",
        "status": "blocked",
        "ready": False,
        "blockers": ["candidate_artifact", "independent_validation"],
    },
]

CANDIDATE_ARTIFACT_CONTRACT = {
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
        "pressure": "pressure(points, t) -> (...) or explicit unavailable=false/true state",
        "forcing": "forcing(points, t) -> (..., 3) under the fixed/restricted contract only",
        "save_load": "deterministic JSON metadata plus hashed numeric payloads",
    },
    "required_stage_receipts": ["leading", "oscillatory", "correction", "independent_validation"],
    "required_truth_flags": [
        "leading_ready",
        "public_source_oscillatory_xyz_t_velocity_ready",
        "genuinely_independent_second_covariance_column_ready",
        "correction_ready",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
        "paper_exact",
    ],
}

STATES = {
    "outer_base_backbone_ready": True,
    "slow_squared_partition_family_contract_ready": True,
    "manufactured_support_localization_axis_refinement_passed": True,
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


def build_checkpoint(previous: dict[str, Any] | None = None) -> dict[str, Any]:
    if previous is None:
        previous = build_v18_checkpoint()
    validate_v18_checkpoint(previous)
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base": {"agent5_pr": BASE_AGENT5_PR, "agent5_head": BASE_AGENT5_HEAD},
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "previous_checkpoint_sha256": previous["checkpoint_sha256"],
        "latest": LATEST,
        "pipeline_stages": PIPELINE_STAGES,
        "candidate_artifact_contract": CANDIDATE_ARTIFACT_CONTRACT,
        "states": STATES,
        "shortest_closure": [
            "Agent 4: independently audit A2 #369 ancestry on a source-grid phase inside the source wave shell using the unchanged Cartesian curl/divergence guard",
            "Agent 2: after that pass, instantiate concrete source partition/background and expose actual-source Q-scaled velocity(x,y,z,t) plus a genuinely independent covariance column",
            "Agent 3: admit that exact second column, retain radial d_z sigma_1, then rerun the finite correction cycle only if the frozen rank/bounded-inverse guards pass",
            "Agent 1: compose #368 base backbone with cone modulation and I1/I2/I3/I4 overlays and provide global matched pressure/admissible fixed forcing",
            "Agent 5: instantiate the reserved candidate artifact only after those four inputs are executable in one ancestry",
            "Agent 4: run the formal held-out 1e-3 momentum and 1e-5 divergence gate on the frozen composite artifact",
        ],
        "baseline_vs_kokuno": [
            {
                "name": "ST006",
                "scope": "held-out full-domain normalized NS momentum",
                "momentum_sampled_max": ST006_REFERENCE["momentum_sampled_max"],
                "volume_l2": ST006_REFERENCE["volume_l2"],
                "pde_validated": False,
            },
            {
                "name": "Kokuno multilane frontier v19",
                "scope": "integration/readiness manifest; no complete velocity/pressure/forcing composite",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "directly_comparable_to_ST006": False,
                "formal_full_domain_pde_gate_assessed": False,
            },
        ],
        "truth_boundary": {
            "threshold_relaxed": False,
            "surrogate_data_promoted": False,
            "free_residual_defined_forcing_used": False,
            "kokuno_replay_used_as_independent_validation": False,
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
        raise ValueError("unexpected Agent-5 frontier schema/task")
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
        raise ValueError("frontier states changed")
    if payload.get("candidate_artifact_contract") != CANDIDATE_ARTIFACT_CONTRACT:
        raise ValueError("candidate artifact contract changed")
    stages = payload.get("pipeline_stages")
    if stages != PIPELINE_STAGES or any(stage.get("ready") is not False for stage in stages):
        raise ValueError("an incomplete pipeline stage was promoted")
    if payload["latest"]["agent2"]["source_wave_shell_independent_cartesian_preflight_passed"] is not False:
        raise ValueError("source-wave-shell independent preflight was falsely promoted")
    if payload["latest"]["agent3"]["finite_correction_cycle_rerun_allowed"] is not False:
        raise ValueError("finite correction cycle was falsely authorized")
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
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)
    checkpoint_path = root / "multilane_frontier_checkpoint.json"
    checkpoint_path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema": SCHEMA,
        "checkpoint_file": checkpoint_path.name,
        "checkpoint_file_sha256": hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
        "checkpoint_payload_sha256": checkpoint["checkpoint_sha256"],
        "pde_validated": False,
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"checkpoint": checkpoint_path, "manifest": manifest_path}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    paths = write_bundle(args.output_dir)
    print(json.dumps({key: str(value) for key, value in paths.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
