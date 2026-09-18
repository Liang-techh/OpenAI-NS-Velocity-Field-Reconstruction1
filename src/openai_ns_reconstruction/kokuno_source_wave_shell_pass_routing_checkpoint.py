"""Kokuno Agent-5 routing after the independent source-wave-shell roundoff ladder.

This module adds no Kokuno mathematics and does not reinterpret the historical
Agent-4 #371 rejection.  It binds the separately preregistered Agent-4 #383
receipt, where a coarser truncation-dominated FD4 ladder passes the *same*
local guards while #371 remains recorded as a failed finer-ladder experiment.

The result clears only the source-wave-shell implementation preflight.  It does
not create concrete source partitions/background data, a public actual-source
oscillatory velocity, a second covariance column, a global leading pressure,
or a full Kokuno candidate.  Formal held-out normalized NS momentum max/L2
remain <=1e-3 and divergence max/L2 <=1e-5; no residual-defined forcing is
allowed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-source-wave-shell-pass-routing-checkpoint-v20"
TASK_ID = "KOKUNO-A5-SOURCE-WAVE-SHELL-PASS-ROUTING-020"
BASE_AGENT4_PR = 383
BASE_AGENT4_HEAD = "71c6525ef710cc7bc2bf3599777da3b59af4b900"

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

HISTORICAL_AGENT4_REJECTION = {
    "pr": 371,
    "head_sha": "917db15e5b72a04662cf78cd0345b32a9adeb51c",
    "seed": 9173101,
    "fd4_steps": [0.001, 0.0005, 0.00025],
    "finest_curl_relative_rms": 1.3245337639271412e-12,
    "finest_divergence_max_abs": 2.1924837789035134e-10,
    "divergence_rms_refinement_ratios": [11.789475096192975, 1.193958216011321],
    "minimum_refinement_ratio": 2.5,
    "source_wave_shell_independent_cartesian_preflight_passed": False,
    "failed_checks": ["divergence_refinement_ratio"],
    "overwritten": False,
}

AGENT4_PASS_RECEIPT = {
    "task_id": "KOKUNO-A4-SOURCE-WAVE-SHELL-ROUNDOFF-LADDER-AUDIT-020",
    "source_pr": 383,
    "source_head": BASE_AGENT4_HEAD,
    "dedicated_run": 35336259814,
    "standard_run": 35336259839,
    "artifact_id": 10542893497,
    "artifact_digest": "sha256:4350fbeba26ca8a54ee08386f48c33709af034fd96df443eb4a9ec904e1e91a1",
    "seed": 9173111,
    "sample_count": 36,
    "fd4_steps": [0.004, 0.002, 0.001],
    "finest": {
        "step": 0.001,
        "curl_relative_rms": 2.755788271978995e-10,
        "divergence_max_abs": 3.580782588711177e-09,
        "divergence_rms_abs": 1.226740480615126e-09,
        "source_shell_interior_divergence_max_abs": 2.176313019047129e-09,
        "near_inner_shell_edge_divergence_max_abs": 3.580782588711177e-09,
    },
    "curl_refinement_ratios": [15.995104797103064, 15.998167543805044],
    "divergence_refinement_ratios": [15.995372546512321, 16.024253341545737],
    "mutation": {
        "curl_point_relative_min": 0.010739225368083099,
        "curl_point_relative_rms": 0.024512254413883293,
        "mutated_divergence_max_abs": 0.05701003361962002,
    },
    "local_thresholds": {
        "finest_curl_relative_rms_max": 2.0e-4,
        "finest_divergence_max_abs_max": 2.0e-4,
        "minimum_refinement_ratio": 2.5,
        "missing_gradient_relative_error_min": 1.0e-3,
        "missing_gradient_divergence_max_abs_min": 1.0e-4,
    },
    "all_local_checks_passed": True,
    "source_wave_shell_roundoff_ladder_preflight_passed": True,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}

UPSTREAM = {
    "agent1": {
        "latest_pr": 381,
        "head_sha": "0ac934c640a0520f70bdd66d69bd42fdc80ccf5d",
        "latest_increment": "PA.17 I1 shear five-moment repair API",
        "cone_modulation_completed": False,
        "i1_repair_applied_to_actual_modulation": False,
        "global_matched_pressure_ready": False,
        "leading_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent2": {
        "latest_pr": 380,
        "head_sha": "43f021e904dbb09c9bc2aba1d9ea2954083937e1",
        "dedicated_run": 35335492487,
        "standard_run": 35335492324,
        "localized_real_pair_family_ready": True,
        "per_label_q_scaling_ready": True,
        "concrete_source_partition_bumps_reconstructed": False,
        "source_actual_partition_labels_instantiated": False,
        "source_actual_background_path_instantiated": False,
        "public_source_oscillatory_xyz_t_velocity_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent3": {
        "latest_pr": 382,
        "head_sha": "c80d3d0993c0cc62e2f24cab5d427c5fcb12138d",
        "dedicated_run": 35336161949,
        "artifact_id": 10542888440,
        "artifact_digest": "sha256:359336dc7e83da01052f06d72513d18d4dd65a6bad73f569f36ad78baed688a1",
        "family_cross_term_contract_executable": True,
        "actual_source_multilabel_family_consumed": False,
        "source_coefficient_unit_mapping_available": False,
        "bounded_inverse_ready": False,
        "finite_correction_cycle_rerun_allowed": False,
        "radial_force_side_effect_required": True,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "pr": BASE_AGENT4_PR,
        "head_sha": BASE_AGENT4_HEAD,
        "dedicated_run": 35336259814,
        "standard_run": 35336259839,
        "artifact_id": 10542893497,
        "artifact_digest": "sha256:4350fbeba26ca8a54ee08386f48c33709af034fd96df443eb4a9ec904e1e91a1",
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5": {
        "pr": 374,
        "head_sha": "1f1a3d3e93245ba7bec8168259496ac2ebffc3d3",
        "historical_rejection_routed": True,
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
    "source_wave_shell_independent_cartesian_preflight_assessed": True,
    "source_wave_shell_independent_cartesian_preflight_passed": True,
    "historical_source_wave_shell_rejection_preserved": True,
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

PIPELINE_FRONTIER = [
    {
        "stage": "source_profile_ingest",
        "ready": False,
        "status": "partial",
        "blockers": [
            "concrete source partition bumps/labels",
            "actual positive-order background path",
            "remaining leading cone/overlay inputs",
        ],
    },
    {
        "stage": "leading_candidate",
        "ready": False,
        "status": "blocked",
        "blockers": [
            "actual cone/radial modulation and its five-moment discrepancy",
            "eta-smooth I1 repair coefficient family",
            "remaining I1/I2/I3/I4 composition",
            "global matched pressure and admissible fixed forcing",
        ],
    },
    {
        "stage": "oscillatory_augmentation",
        "ready": False,
        "status": "source_shell_preflight_cleared_waiting_for_actual_source_inputs",
        "evidence": ["Agent-4 #383 passes the separately preregistered roundoff ladder with the original local guards"],
        "blockers": [
            "concrete source partition values/derivatives and labels",
            "actual positive-order background",
            "public Q-scaled velocity(x,y,z,t)",
            "genuinely independent second covariance column",
        ],
    },
    {
        "stage": "mean_radial_corrections",
        "ready": False,
        "status": "blocked",
        "evidence": [
            "Agent-3 #382 cross-label covariance contract is executable",
            "radial d_z sigma_1 side effect remains required",
        ],
        "blockers": [
            "actual source multi-label family",
            "source coefficient-unit mapping",
            "second public covariance column passing the frozen spacetime bounded-inverse screen",
        ],
    },
    {"stage": "candidate_artifact", "ready": False, "status": "blocked"},
    {"stage": "independent_full_pde_validation", "ready": False, "status": "waiting_for_candidate"},
    {"stage": "report_and_export", "ready": False, "status": "blocked"},
]


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def build_checkpoint() -> dict[str, Any]:
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base": {"agent4_pr": BASE_AGENT4_PR, "agent4_head": BASE_AGENT4_HEAD},
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "historical_agent4_rejection": HISTORICAL_AGENT4_REJECTION,
        "agent4_pass_receipt": AGENT4_PASS_RECEIPT,
        "upstream": UPSTREAM,
        "future_candidate_artifact_contract": FUTURE_CANDIDATE_ARTIFACT_CONTRACT,
        "pipeline_frontier": PIPELINE_FRONTIER,
        "routing": {
            "new_fact": (
                "Agent-4 #383 independently clears the source-grid/source-wave-shell complete-curl preflight on a separately preregistered coarser FD4 ladder while preserving #371 as a historical failed finer-ladder experiment."
            ),
            "numerical_signature": (
                "At h=1e-3, curl relative RMS is 2.76e-10 and divergence max is 3.58e-9; curl refinement is 15.995x/15.998x and divergence refinement is 15.995x/16.024x under the unchanged 2e-4 absolute and 2.5 refinement guards."
            ),
            "interpretation": (
                "The source-shell implementation seam is now independently cleared for routing. This does not retroactively turn #371 into a pass and does not create actual source partitions/background or a public oscillatory field."
            ),
            "do_not_do": [
                "do not overwrite or reinterpret Agent-4 #371's historical rejection",
                "do not promote supplied/manufactured partition/background data as actual source inputs",
                "do not call the supplied-data real-pair family a public actual-source velocity(x,y,z,t)",
                "do not rerun Agent-3 finite correction until a genuinely independent public second column and source coefficient units exist",
                "do not instantiate/export the final candidate before global velocity/pressure/restricted-forcing composition exists",
                "do not weaken the formal 1e-3 momentum or 1e-5 divergence gates",
            ],
            "shortest_next_closure": [
                "Agent 2: consume the cleared source-shell seam, bind concrete source partition values/derivatives and actual positive-order background into #380's localized real-pair family, then expose Q-scaled public velocity(x,y,z,t) plus a genuinely independent second covariance column",
                "Agent 3: route that exact multi-label family through the #382 cross-term-aware covariance contract and frozen spacetime rank/bounded-inverse screen, retain radial d_z sigma_1, and only then run the finite correction cycle",
                "Agent 1: complete actual cone/radial modulation, its PA.17 I1 discrepancy and eta-smooth repair family, remaining overlays, global matched pressure and admissible fixed forcing",
                "Agent 5: instantiate the reserved candidate artifact only once leading, oscillatory and correction stages coexist in one executable ancestry",
                "Agent 4: run the formal held-out normalized momentum <=1e-3 and divergence <=1e-5 gate on that frozen complete candidate, with same-protocol ST006 comparison",
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
                "name": "Kokuno source-wave-shell seam #383",
                "scope": "local source-grid/source-shell Cartesian curl/divergence implementation audit",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "directly_comparable_to_ST006": False,
                "local_preflight_passed": True,
                "historical_rejection_preserved": True,
            },
        ],
        "states": STATES,
        "truth_boundary": {
            "threshold_relaxed": False,
            "surrogate_data_promoted": False,
            "free_residual_defined_forcing_used": False,
            "kokuno_replay_used_as_independent_pde_validation": False,
            "historical_rejection_overwritten": False,
            "supplied_data_family_promoted_to_actual_source": False,
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
        raise ValueError("unexpected Agent-5 source-shell pass routing schema/task")
    signature = payload.get("checkpoint_sha256")
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    if signature != _sha(unsigned):
        raise ValueError("checkpoint SHA mismatch")
    if payload.get("fixed_gates") != FIXED_GATES:
        raise ValueError("fixed gates changed")
    if payload.get("st006_reference") != ST006_REFERENCE:
        raise ValueError("ST006 reference changed")
    if payload.get("historical_agent4_rejection") != HISTORICAL_AGENT4_REJECTION:
        raise ValueError("historical Agent-4 rejection changed")
    if payload.get("agent4_pass_receipt") != AGENT4_PASS_RECEIPT:
        raise ValueError("Agent-4 pass receipt changed")
    if payload.get("upstream") != UPSTREAM:
        raise ValueError("upstream routing receipt changed")
    if payload.get("states") != STATES:
        raise ValueError("routing states changed")
    if payload.get("pipeline_frontier") != PIPELINE_FRONTIER:
        raise ValueError("pipeline frontier changed")
    if payload.get("future_candidate_artifact_contract") != FUTURE_CANDIDATE_ARTIFACT_CONTRACT:
        raise ValueError("future candidate artifact contract changed")
    historical = payload["historical_agent4_rejection"]
    if historical.get("source_wave_shell_independent_cartesian_preflight_passed") is not False:
        raise ValueError("historical source-shell rejection was overwritten")
    if historical.get("overwritten") is not False:
        raise ValueError("historical source-shell rejection was marked overwritten")
    receipt = payload["agent4_pass_receipt"]
    if receipt.get("source_wave_shell_roundoff_ladder_preflight_passed") is not True:
        raise ValueError("Agent-4 roundoff-ladder pass disappeared")
    if receipt.get("all_local_checks_passed") is not True:
        raise ValueError("Agent-4 local pass receipt incomplete")
    if payload["states"]["source_wave_shell_independent_cartesian_preflight_passed"] is not True:
        raise ValueError("cleared source-shell routing state disappeared")
    for key in (
        "leading_ready",
        "public_source_oscillatory_xyz_t_velocity_ready",
        "genuinely_independent_second_covariance_column_ready",
        "correction_ready",
        "candidate_artifact_instantiated",
        "complete_kokuno_composite_velocity_ready",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        if payload["states"].get(key) is not False:
            raise ValueError(f"global readiness was falsely promoted: {key}")
    if payload["upstream"]["agent3"]["bounded_inverse_ready"] is not False:
        raise ValueError("Agent-3 bounded inverse was falsely promoted")
    if payload["upstream"]["agent3"]["finite_correction_cycle_rerun_allowed"] is not False:
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
    checkpoint_path = root / "source_wave_shell_pass_routing_checkpoint.json"
    checkpoint_path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema": SCHEMA,
        "checkpoint_file": checkpoint_path.name,
        "checkpoint_file_sha256": hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
        "checkpoint_payload_sha256": checkpoint["checkpoint_sha256"],
        "bound_agent4_artifact_id": AGENT4_PASS_RECEIPT["artifact_id"],
        "bound_agent4_artifact_digest": AGENT4_PASS_RECEIPT["artifact_digest"],
        "bound_agent3_sibling_artifact_id": UPSTREAM["agent3"]["artifact_id"],
        "bound_agent3_sibling_artifact_digest": UPSTREAM["agent3"]["artifact_digest"],
        "historical_agent4_rejection_preserved": True,
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
