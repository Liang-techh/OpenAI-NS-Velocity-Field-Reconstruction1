"""Agent-5 routing checkpoint after the concrete slow-partition independent audit.

This module adds no Kokuno mathematics.  It binds Agent 2 PR #392's concrete,
repository-autonomous source-compatible slow squared partition to Agent 4 PR
#394's independent black-box PASS and records the current Agent 1/3 sibling
frontiers.  The key truth boundary is explicit: a source-compatible autonomous
realization may be used in a *Kokuno-derived* candidate when its provenance is
retained, but it is not recovered Kokuno cutoff data and cannot support a
paper-exact claim.

The formal held-out project gates remain normalized momentum max/L2 <= 1e-3
and divergence max/L2 <= 1e-5.  This local partition audit is not a PDE gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-source-compatible-partition-routing-checkpoint-v21"
TASK_ID = "KOKUNO-A5-SOURCE-COMPATIBLE-PARTITION-ROUTING-021"
BASE_AGENT4_PR = 394
BASE_AGENT4_HEAD = "69a7ca724db377510815eddcc2b152b05d7193ac"

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

AGENT4_PARTITION_AUDIT_RECEIPT = {
    "task_id": "KOKUNO-A4-SOURCE-COMPATIBLE-PARTITION-INDEPENDENT-AUDIT-021",
    "source_pr": 394,
    "source_head": BASE_AGENT4_HEAD,
    "dedicated_run": 35341494077,
    "standard_run": 35341494053,
    "artifact_id": 10545247621,
    "artifact_digest": "sha256:6b2cc78a2f0ba6977b938fdafb6aa62a4059964727e14ad9fa681c7246c09765",
    "seed": 9173121,
    "held_out_points": 24,
    "grid_origins": [[0.0, 0.0, 0.0], [1.3e-5, -0.7e-5, 2.1e-5]],
    "fd6_steps": [0.04, 0.02, 0.01],
    "worst_partition_max_abs": 4.440892098500626e-16,
    "worst_differentiated_closure_max_abs": 1.3877787807814457e-16,
    "maximum_active_ell_span": 1,
    "source_active_ell_span_bound": 2,
    "finest_radial_relative_rms_by_origin": [8.118004557177425e-08, 8.934172618602344e-08],
    "finest_axial_relative_rms_by_origin": [1.9270369065891428e-07, 1.7950991176237213e-07],
    "radial_refinement_ratios_by_origin": [
        [23.956051695962998, 29.852541274973195],
        [9.909422062012535, 39.66177474924965],
    ],
    "axial_refinement_ratios_by_origin": [
        [35.19219990902432, 20.074493295053404],
        [17.995403496790125, 33.36958924722895],
    ],
    "dropped_label_partition_error_by_origin": [0.9414901976120158, 0.9714764715732194],
    "dropped_radial_derivative_relative_error_by_origin": [0.3156638879506775, 0.28572745941576727],
    "local_guards": {
        "partition_max_abs": 5.0e-12,
        "differentiated_closure_max_abs": 5.0e-10,
        "finest_fd6_relative_rms": 1.0e-6,
        "minimum_refinement_ratio": 8.0,
        "dropped_label_partition_error_min": 1.0e-4,
        "dropped_derivative_relative_error_min": 1.0e-3,
    },
    "all_local_checks_passed": True,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}

UPSTREAM = {
    "agent1": {
        "latest_pr": 391,
        "head_sha": "c52002793547f18cb5e433050e0478a06a57ffe3",
        "standard_run": 35339897269,
        "finite_frequency_five_moment_discrepancy_executable": True,
        "autonomous_diagnostic_loop_available": True,
        "source_admissible_loop_reconstructed": False,
        "eta_smooth_repair_coefficient_family_reconstructed": False,
        "global_leading_profile_reconstructed": False,
        "global_matched_pressure_ready": False,
        "leading_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent2": {
        "latest_pr": 392,
        "head_sha": "e2f73ccdc4eccaebf03fd20e3b3d9c3d9df53977",
        "dedicated_run": 35340907580,
        "standard_run": 35340907632,
        "source_compatible_partition_realization_ready": True,
        "partition_values_and_derivatives_materialized": True,
        "partition_independently_audited": True,
        "partition_is_repository_autonomous": True,
        "source_actual_partition_recovered": False,
        "source_actual_background_path_instantiated": False,
        "actual_positive_order_mode_binding_ready": False,
        "public_xyz_t_by_beta_velocity_ready": False,
        "source_coefficient_unit_mapping_available": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent3": {
        "latest_pr": 393,
        "head_sha": "bef9206d3aac61f2de8f0c479dfa74cc40938746",
        "dedicated_run": 35341173639,
        "standard_run": 35341173484,
        "artifact_id": 10545117446,
        "artifact_digest": "sha256:4781b7bc0c5bb7fc493c2caae6f8bdda06781f98805330752f4a4da67dfd8feb",
        "family_bounded_inverse_machinery_ready": True,
        "real_active_target_nodes": 27,
        "nodes_requiring_second_direction": 25,
        "duplicate_control_rank_two_coverage": "0/25",
        "missing_relative_stress_response": 0.210194502136004,
        "signed_coefficient_budget": 0.012513055889617838,
        "actual_source_multilabel_family_consumed": False,
        "source_coefficient_unit_mapping_available": False,
        "actual_source_family_bounded_inverse_assessed": False,
        "finite_correction_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "latest_pr": BASE_AGENT4_PR,
        "head_sha": BASE_AGENT4_HEAD,
        "dedicated_run": 35341494077,
        "standard_run": 35341494053,
        "artifact_id": 10545247621,
        "artifact_digest": "sha256:6b2cc78a2f0ba6977b938fdafb6aa62a4059964727e14ad9fa681c7246c09765",
        "source_compatible_partition_independent_preflight_passed": True,
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5": {
        "latest_pr": 385,
        "head_sha": "090514b1cf87b5a16e710e1867004d83eb77014e",
        "source_wave_shell_preflight_routed": True,
        "consumed_in_executable_ancestry": False,
    },
}

AUTONOMOUS_ROUTE_POLICY = {
    "kokuno_derived_autonomous_choices_allowed": True,
    "condition": (
        "Every autonomous choice must remain explicitly labeled, versioned and hashed; "
        "it may support a Kokuno-derived candidate but never a recovered-source or paper-exact claim."
    ),
    "source_compatible_partition_may_be_used_as_scaffold": True,
    "source_compatible_partition_may_be_called_recovered_source": False,
    "autonomous_agent1_loop_may_be_used_as_kokuno_derived_scaffold": True,
    "autonomous_agent1_loop_may_be_called_source_hidden_loop": False,
    "paper_exact": False,
}

FUTURE_CANDIDATE_ARTIFACT_CONTRACT = {
    "schema_name": "kokuno-candidate-artifact-v1-reserved",
    "status": "reserved_not_instantiated",
    "required_provenance": [
        "source_repository",
        "source_commit_or_version",
        "source_document",
        "source_derived_formulas",
        "autonomous_choices",
        "stage_artifact_hashes",
    ],
    "required_public_api": {
        "velocity": "velocity(points, t) -> (..., 3)",
        "pressure": "pressure(points, t) -> (...)",
        "forcing": "forcing(points, t) -> (..., 3) under the fixed/restricted contract only",
        "save_load": "deterministic metadata plus hashed numeric payloads",
    },
    "required_stage_receipts": [
        "source_profile_ingest",
        "leading",
        "oscillatory",
        "mean_radial_correction",
        "independent_validation",
    ],
    "required_stage_metrics": [
        "normalized_momentum_max_or_not_assessed",
        "normalized_momentum_l2_or_not_assessed",
        "divergence_max_or_not_assessed",
        "divergence_l2_or_not_assessed",
        "nontriviality_metric",
    ],
}

STATES = {
    "source_compatible_partition_realization_ready": True,
    "source_compatible_partition_independently_audited": True,
    "source_actual_partition_recovered": False,
    "leading_ready": False,
    "oscillatory_machinery_ready": True,
    "oscillatory_ready": False,
    "public_kokuno_derived_oscillatory_xyz_t_velocity_ready": False,
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
        "status": "source_compatible_autonomous_partition_independently_cleared",
        "evidence": [
            "Agent-2 #392 materializes eta_beta, D_r eta_beta and D_z eta_beta",
            "Agent-4 #394 independently passes the value/derivative seam under frozen guards",
        ],
        "blockers": [
            "actual positive-order/background path",
            "actual source mode/label binding or an explicitly autonomous Kokuno-derived substitute",
            "explicit coefficient-unit mapping for the public by-beta family",
        ],
    },
    {
        "stage": "leading_candidate",
        "ready": False,
        "status": "autonomous_diagnostic_available_global_assembly_missing",
        "evidence": ["Agent-1 #391 supplies an executable finite-frequency five-moment diagnostic target"],
        "blockers": [
            "eta-smooth I1 repair coefficient family for the selected modulation",
            "complete radial overlay composition",
            "global matched pressure",
            "admissible fixed/restricted forcing contract on the assembled leading field",
        ],
    },
    {
        "stage": "oscillatory_augmentation",
        "ready": False,
        "status": "partition_preflight_cleared_waiting_for_background_mode_and_public_family",
        "blockers": [
            "positive-order/background and mode binding",
            "public Q-scaled xyz,t velocity by beta and total",
            "explicit source/autonomous coefficient units",
            "genuinely independent covariance direction",
        ],
    },
    {
        "stage": "mean_radial_corrections",
        "ready": False,
        "status": "bounded_inverse_machinery_ready_waiting_for_actual_family",
        "evidence": [
            "Agent-3 #393 enforces rank, algebraic residual, unit compatibility and aggregate family L1 budget",
            "the current duplicate-family negative control remains 0/25 rank-two with 0.2101945 missing relative stress",
            "the radial d_z sigma_1 side effect remains part of the future materialized correction",
        ],
        "blockers": [
            "public multi-label oscillatory family",
            "coefficient-unit mapping",
            "frozen spacetime/radial bounded-inverse PASS",
        ],
    },
    {"stage": "candidate_artifact", "ready": False, "status": "waiting_for_leading_oscillatory_correction"},
    {"stage": "independent_full_pde_validation", "ready": False, "status": "waiting_for_frozen_complete_candidate"},
    {"stage": "report_and_export", "ready": False, "status": "waiting_for_validation_receipt"},
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
        "agent4_partition_audit_receipt": AGENT4_PARTITION_AUDIT_RECEIPT,
        "upstream": UPSTREAM,
        "autonomous_route_policy": AUTONOMOUS_ROUTE_POLICY,
        "future_candidate_artifact_contract": FUTURE_CANDIDATE_ARTIFACT_CONTRACT,
        "pipeline_frontier": PIPELINE_FRONTIER,
        "states": STATES,
        "routing": {
            "new_fact": (
                "Agent-4 #394 independently clears Agent-2 #392's concrete source-compatible slow-partition value/derivative seam; no remaining numerical blocker is retained at this partition interface."
            ),
            "truth_boundary_unlock": (
                "The cleared #392 partition is repository-autonomous. It can now be routed as a provenance-labeled Kokuno-derived scaffold without waiting for unrecoverable hidden cutoff data, but it must never be labeled recovered source or paper-exact. The same distinction applies to Agent-1 #391's autonomous modulation diagnostic."
            ),
            "do_not_do": [
                "do not repeat generic source-partition derivative validation after the independent #394 PASS unless the implementation changes",
                "do not rename #392's autonomous cutoff/grid choices as recovered Kokuno data",
                "do not fabricate a second covariance direction from duplicated columns",
                "do not rerun the finite correction cycle before #393's actual-family bounded-inverse screen passes",
                "do not instantiate/export the final Kokuno candidate before leading, oscillatory and correction stages coexist",
                "do not weaken the fixed 1e-3 momentum or 1e-5 divergence gates",
                "do not use residual-defined free forcing",
            ],
            "shortest_next_closure": [
                "Agent 2: stop revisiting the cleared partition seam; bind the positive-order/background and mode data into #392+#380, expose public Q-scaled xyz,t velocity by beta and total, and publish an explicit coefficient-unit mapping",
                "Agent 3: consume that exact family through #382 cross-label covariance and #393 aggregate-budget bounded inverse; only a frozen spacetime/radial PASS may authorize correction materialization and the finite cycle, retaining radial d_z sigma_1",
                "Agent 1: use a provenance-labeled selected modulation (source-recovered if ever available, otherwise explicitly autonomous Kokuno-derived) to build an eta-smooth PA.17 repair family, then assemble a global leading velocity/pressure with admissible fixed/restricted forcing",
                "Agent 5: once those three stages coexist in one executable ancestry, instantiate the reserved deterministic candidate artifact and add Python/MATLAB smoke export",
                "Agent 4: run the formal held-out normalized momentum <=1e-3 and divergence <=1e-5 gate on the frozen complete candidate with a same-protocol ST006 comparison",
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
                "name": "Kokuno source-compatible partition #392/#394",
                "scope": "local slow-partition value/derivative implementation audit",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "directly_comparable_to_ST006": False,
                "local_preflight_passed": True,
                "source_recovered": False,
            },
        ],
        "truth_boundary": {
            "threshold_relaxed": False,
            "surrogate_data_promoted": False,
            "free_residual_defined_forcing_used": False,
            "kokuno_replay_used_as_independent_pde_validation": False,
            "autonomous_partition_claimed_as_recovered_source": False,
            "autonomous_modulation_claimed_as_source_hidden_loop": False,
            "duplicate_covariance_columns_promoted_as_independent": False,
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
        raise ValueError("unexpected Agent-5 source-compatible-partition routing schema/task")
    signature = payload.get("checkpoint_sha256")
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    if signature != _sha(unsigned):
        raise ValueError("checkpoint SHA mismatch")
    if payload.get("fixed_gates") != FIXED_GATES:
        raise ValueError("fixed gates changed")
    if payload.get("st006_reference") != ST006_REFERENCE:
        raise ValueError("ST006 reference changed")
    if payload.get("agent4_partition_audit_receipt") != AGENT4_PARTITION_AUDIT_RECEIPT:
        raise ValueError("Agent-4 partition audit receipt changed")
    if payload.get("upstream") != UPSTREAM:
        raise ValueError("upstream routing receipt changed")
    if payload.get("autonomous_route_policy") != AUTONOMOUS_ROUTE_POLICY:
        raise ValueError("autonomous route policy changed")
    if payload.get("future_candidate_artifact_contract") != FUTURE_CANDIDATE_ARTIFACT_CONTRACT:
        raise ValueError("future candidate artifact contract changed")
    if payload.get("pipeline_frontier") != PIPELINE_FRONTIER:
        raise ValueError("pipeline frontier changed")
    if payload.get("states") != STATES:
        raise ValueError("routing states changed")

    receipt = payload["agent4_partition_audit_receipt"]
    if receipt.get("all_local_checks_passed") is not True:
        raise ValueError("independent partition preflight lost its PASS")
    if receipt.get("formal_full_domain_pde_gate_assessed") is not False:
        raise ValueError("local partition audit was promoted to a PDE gate")
    a2 = payload["upstream"]["agent2"]
    if a2.get("source_compatible_partition_realization_ready") is not True:
        raise ValueError("source-compatible partition realization disappeared")
    if a2.get("source_actual_partition_recovered") is not False:
        raise ValueError("autonomous partition was promoted to recovered source")
    if a2.get("public_xyz_t_by_beta_velocity_ready") is not False:
        raise ValueError("public oscillatory family was falsely promoted")
    a3 = payload["upstream"]["agent3"]
    if a3.get("finite_correction_cycle_rerun_allowed") is not False:
        raise ValueError("finite correction cycle was falsely authorized")
    policy = payload["autonomous_route_policy"]
    if policy.get("kokuno_derived_autonomous_choices_allowed") is not True:
        raise ValueError("truthful Kokuno-derived autonomous route was disabled")
    if policy.get("source_compatible_partition_may_be_called_recovered_source") is not False:
        raise ValueError("autonomous partition may not be called recovered source")
    if policy.get("paper_exact") is not False:
        raise ValueError("paper-exact was falsely promoted")
    for key in (
        "source_actual_partition_recovered",
        "leading_ready",
        "oscillatory_ready",
        "public_kokuno_derived_oscillatory_xyz_t_velocity_ready",
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
    checkpoint_path = root / "source_compatible_partition_routing_checkpoint.json"
    checkpoint_path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema": SCHEMA,
        "checkpoint_file": checkpoint_path.name,
        "checkpoint_file_sha256": hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
        "checkpoint_payload_sha256": checkpoint["checkpoint_sha256"],
        "bound_agent4_artifact_id": AGENT4_PARTITION_AUDIT_RECEIPT["artifact_id"],
        "bound_agent4_artifact_digest": AGENT4_PARTITION_AUDIT_RECEIPT["artifact_digest"],
        "bound_agent3_sibling_artifact_id": UPSTREAM["agent3"]["artifact_id"],
        "bound_agent3_sibling_artifact_digest": UPSTREAM["agent3"]["artifact_digest"],
        "source_compatible_partition_independently_audited": True,
        "source_actual_partition_recovered": False,
        "candidate_artifact_instantiated": False,
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
