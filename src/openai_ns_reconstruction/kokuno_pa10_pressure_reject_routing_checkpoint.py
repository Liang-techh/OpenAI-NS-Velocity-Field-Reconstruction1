"""Kokuno Agent-5 routing after the independent PA.10 pressure audit.

This checkpoint performs one narrow integration step.  It reruns Agent 4 PR #608's
public-only audit of Agent 1 PR #606's selected PA.10 pressure primitive and routes
its scientific REJECT without weakening any gate.  Existing admitted oscillatory
and radial-stress seams remain usable, but the rejected pressure value path is not
allowed to enter the full same-cycle defect/correction pipeline.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_agent4_pa10_pressure_primitive_audit import run_audit

TASK = "KOKUNO-A5-PA10-PRESSURE-REJECT-ROUTING-045"
SCHEMA = "kokuno-agent5-pa10-pressure-reject-routing-checkpoint-v45"

ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622

AGENT1_PR = 606
AGENT1_HEAD = "56d7e47e4033d869f802ee51ed380b1daa4a34a1"
AGENT1_CANDIDATE_SHA256 = "8b5cdfb749d395f25d7c5f09f9200a3b281dd346fcb4ead6c8c5b1583d9ff675"
AGENT1_STANDARD_RUN = 35434957295

AGENT2_VELOCITY_PR = 561
AGENT2_VELOCITY_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
AGENT2_TIME_DERIVATIVE_PR = 579
AGENT2_TIME_DERIVATIVE_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"

AGENT3_SPACETIME_ADMISSION_PR = 607
AGENT3_SPACETIME_ADMISSION_HEAD = "a2fc36a44771cffb173e391d53fb04df220ab638"
AGENT3_SPACETIME_DEDICATED_RUN = 35435350975
AGENT3_SPACETIME_STANDARD_RUN = 35435350972
AGENT3_SPACETIME_ARTIFACT_ID = 10581334240
AGENT3_SPACETIME_ARTIFACT_DIGEST = (
    "sha256:68dbb652b633092b766d7ec7ea089d937ed26c4d05d2b45ba4d29fc3077b5d3a"
)

AGENT4_PRESSURE_AUDIT_PR = 608
AGENT4_PRESSURE_AUDIT_HEAD = "9a36e257691b4757c37194aa09c11205f21a1ddd"
AGENT4_PRESSURE_DEDICATED_RUN = 35435721364
AGENT4_PRESSURE_STANDARD_RUN = 35435750954
AGENT4_PRESSURE_ARTIFACT_ID = 10582059362
AGENT4_PRESSURE_ARTIFACT_DIGEST = (
    "sha256:408892549b47495142aab4c1055a7a4a40cca035a75e189f87757b68398b614e"
)

EXPECTED_FAILED_GUARDS = {
    "p_y_refinement",
    "Y_Phi_Y_finest_relative_rms",
    "Y_Phi_Y_finest_relative_max",
    "Y_Phi_Y_refinement",
    "p_eta_finest_relative_rms",
    "p_eta_finest_relative_max",
    "p_eta_refinement",
}

FORMAL_GATES = {
    "held_out_normalized_momentum_max": 1.0e-3,
    "held_out_normalized_momentum_l2": 1.0e-3,
    "held_out_divergence_max": 1.0e-5,
    "held_out_divergence_l2": 1.0e-5,
}


def _canonical_sha256(payload: dict[str, Any]) -> str:
    unsigned = copy.deepcopy(payload)
    unsigned.pop("checkpoint_sha256", None)
    encoded = json.dumps(
        unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _audit_selected_pressure() -> dict[str, Any]:
    audit = json.loads(json.dumps(run_audit(), sort_keys=True, allow_nan=False))
    provenance = audit["provenance"]
    if provenance.get("agent1_pr") != AGENT1_PR:
        raise ValueError("wrong Agent-1 PR in pressure audit")
    if provenance.get("agent1_exact_head") != AGENT1_HEAD:
        raise ValueError("wrong Agent-1 head in pressure audit")
    if provenance.get("agent1_candidate_sha256") != AGENT1_CANDIDATE_SHA256:
        raise ValueError("wrong Agent-1 candidate identity in pressure audit")

    project_gates = audit["frozen_protocol"]["final_project_gates_unchanged"]
    expected_project_gates = {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_L2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
    }
    if project_gates != expected_project_gates:
        raise ValueError("Agent-4 pressure audit changed final project gates")

    if audit.get("selected_pressure_primitive_independent_preflight_passed") is not False:
        raise ValueError("expected Agent-4 pressure audit to remain scientific REJECT")
    if set(audit.get("failed_guards", [])) != EXPECTED_FAILED_GUARDS:
        raise ValueError("Agent-4 pressure failed-guard set changed")

    truth = audit["truth_boundary"]
    for key in (
        "source_pressure_radius_one_ball_bound_assessed",
        "global_pressure_matched",
        "global_leading_profile_reconstructed",
        "full_same_cycle_requested_stress_materialized",
        "public_velocity_correction_materialized",
        "heldout_ns_momentum_residual_assessed",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"pressure audit over-promoted downstream truth: {key}")

    p_eta = audit["results"]["p_eta"]
    levels = p_eta["levels"]
    if len(levels) != 3:
        raise ValueError("unexpected p_eta refinement ladder")
    if levels[-1]["relative_rms"] <= audit["frozen_protocol"]["local_gates"]["p_eta_relative_rms"]:
        raise ValueError("p_eta blocker unexpectedly disappeared")
    if audit["results"]["nontrivial_reference_rms"]["p_eta"] <= 1.0e8:
        raise ValueError("p_eta reference is no longer demonstrably nontrivial")
    return audit


def build_checkpoint() -> dict[str, Any]:
    pressure = _audit_selected_pressure()
    p_y = pressure["results"]["p_Y"]
    phi_eta = pressure["results"]["Phi_eta"]
    p_eta = pressure["results"]["p_eta"]
    fto = pressure["results"]["fundamental_theorem_Y"]

    checkpoint: dict[str, Any] = {
        "task": TASK,
        "schema": SCHEMA,
        "round": 45,
        "purpose": (
            "route the independent scientific rejection of the selected PA.10 pressure value path, "
            "while preserving already-admitted oscillatory/radial-stress seams and keeping full "
            "leading, correction, candidate-artifact and PDE gates fail closed"
        ),
        "upstream": {
            "agent1_selected_pressure": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "candidate_sha256": AGENT1_CANDIDATE_SHA256,
                "standard_workflow_run": AGENT1_STANDARD_RUN,
                "public_pressure_interface_exposed": True,
                "source_pressure_radius_one_ball_bound_machine_certified": False,
                "source_R1_R2_machine_bound": False,
                "source_M_K_machine_bound": False,
                "source_B0_T_sh_machine_bound": False,
                "global_matched_pressure_ready": False,
                "global_leading_velocity_ready": False,
            },
            "agent2_admitted_oscillatory": {
                "velocity_pr": AGENT2_VELOCITY_PR,
                "velocity_head": AGENT2_VELOCITY_HEAD,
                "time_derivative_pr": AGENT2_TIME_DERIVATIVE_PR,
                "time_derivative_head": AGENT2_TIME_DERIVATIVE_HEAD,
                "velocity_candidate_changed_in_this_round": False,
            },
            "agent3_spacetime_radial_stress_admission": {
                "pr": AGENT3_SPACETIME_ADMISSION_PR,
                "head": AGENT3_SPACETIME_ADMISSION_HEAD,
                "dedicated_workflow_run": AGENT3_SPACETIME_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT3_SPACETIME_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "artifact_id": AGENT3_SPACETIME_ARTIFACT_ID,
                "artifact_zip_digest": AGENT3_SPACETIME_ARTIFACT_DIGEST,
                "independent_spacetime_radial_stress_generalization_admitted": True,
                "finite_correction_cycle_radial_operator_reuse_allowed_without_retuning": True,
                "actual_same_cycle_full_defect_still_required": True,
            },
            "agent4_selected_pressure_audit": {
                "pr": AGENT4_PRESSURE_AUDIT_PR,
                "head": AGENT4_PRESSURE_AUDIT_HEAD,
                "dedicated_workflow_run": AGENT4_PRESSURE_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT4_PRESSURE_STANDARD_RUN,
                "standard_workflow_status_at_checkpoint_freeze": "in_progress",
                "artifact_id": AGENT4_PRESSURE_ARTIFACT_ID,
                "artifact_zip_digest": AGENT4_PRESSURE_ARTIFACT_DIGEST,
                "selected_pressure_primitive_independent_preflight_passed": False,
                "failed_guards": pressure["failed_guards"],
                "p_Y": p_y,
                "Phi_eta": phi_eta,
                "p_eta": p_eta,
                "fundamental_theorem_Y": fto,
                "nontrivial_reference_rms": pressure["results"]["nontrivial_reference_rms"],
                "mutation_relative_rms": pressure["results"]["mutation_relative_rms"],
            },
        },
        "states": {
            "leading_ready": False,
            "selected_pressure_primitive_interface_exposed": True,
            "selected_pressure_primitive_independent_preflight_passed": False,
            "selected_pressure_public_value_path_repair_required": True,
            "oscillatory_ready": True,
            "public_oscillatory_time_derivative_independently_validated": True,
            "correction_ingest_allowed": True,
            "oscillatory_requested_stress_component_materialized": True,
            "oscillatory_component_finite_head_mean_debt_materialized": True,
            "compact_radial_stress_operator_independently_admitted": True,
            "independent_spacetime_radial_stress_generalization_admitted": True,
            "finite_correction_cycle_radial_operator_reuse_allowed_without_retuning": True,
            "same_cycle_requested_stress_materialized": False,
            "candidate_numeric_finite_head_mean_debt_materialized": False,
            "real_full_candidate_defect_consumed": False,
            "signed_mean_inverse_input_ready": False,
            "public_velocity_correction_materialized": False,
            "correction_ready": False,
            "finite_correction_cycle_run": False,
            "candidate_artifact_instantiated": False,
            "velocity_export_ready": False,
            "formal_full_domain_pde_gate_assessed": False,
            "heldout_ns_residual_assessed": False,
            "pde_validated": False,
        },
        "typed_handoffs": {
            "leading": {
                "producer": "agent1",
                "validator": "agent4",
                "ready": False,
                "public_selected_pressure_interface_available": True,
                "independent_selected_pressure_preflight_passed": False,
                "next_required": (
                    "repair the public g(eta)/pressure value path on the Lambda^{-1/2} scale using a local "
                    "phase difference integral or equivalent precision-safe stationary-point representation; "
                    "rerun the unchanged Agent-4 public-value audit; only after that PASS derive coefficient-space "
                    "radius-one-ball pressure/mixed bounds and Y*u_Y, then close R1/R2 -> M,K -> B0/T_sh -> "
                    "PA.16 -> global matched pressure/leading"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "validator": "agent4",
                "ready": True,
                "velocity_api": "velocity_osc(x,y,z,t)->[...,3]",
                "time_derivative_api": "velocity_osc_dt(x,y,z,t)->[...,3]",
                "retune_allowed_by_this_checkpoint": False,
            },
            "correction": {
                "producer": "agent3",
                "ingest_allowed": True,
                "compact_radial_stress_operator_independently_admitted": True,
                "spacetime_reuse_without_retuning_admitted": True,
                "ready": False,
                "next_required": (
                    "do not consume Agent-1 #606 pressure as a full-defect input; wait for the repaired, "
                    "independently accepted global leading/cross + matched-pressure + restricted-forcing handoff, "
                    "then form the actual same-cycle theta/axial defect and continue the already-guarded correction cycle"
                ),
            },
            "final_candidate": {
                "ready": False,
                "artifact_schema_ready": False,
                "save_load_ready": False,
                "python_visualization_smoke_ready": False,
                "matlab_visualization_smoke_ready": False,
                "independent_pde_validation_ready": False,
            },
        },
        "pressure_audit_interpretation": {
            "p_Y_value_derivative_identity_numerically_accurate": True,
            "p_Y_refinement_guard_failed_at_machine_floor": True,
            "Phi_eta_independent_fd_passed": True,
            "Y_Phi_Y_stationary_center_cloud_not_binary64_identifiable": True,
            "p_eta_is_substantive_blocker": True,
            "p_eta_failure_attributed_to_formula_error": False,
            "likely_integration_issue": (
                "public log_g forms Lambda*(phase(eta)-phase_max) from separately materialized binary64 global "
                "phases; at O(1/Lambda) differences this can flatten g while declared analytic p_eta stays nonzero"
            ),
        },
        "baseline_vs_kokuno": {
            "st006": {
                "same_protocol_full_domain_baseline": True,
                "momentum_sampled_max": ST006_MOMENTUM_MAX,
                "momentum_volume_l2": ST006_MOMENTUM_L2,
                "pde_validated": False,
            },
            "kokuno": {
                "comparable_full_domain_momentum_receipt_exists": False,
                "pressure_interface_metrics_are_not_ns_residuals": True,
                "p_eta_finest_relative_rms": p_eta["levels"][-1]["relative_rms"],
                "p_eta_finest_relative_max": p_eta["levels"][-1]["relative_max"],
            },
        },
        "formal_gates_unchanged": copy.deepcopy(FORMAL_GATES),
        "truth_boundary": {
            "kokuno_reconstruction_is_paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "free_residual_defined_forcing_used": False,
            "surrogate_defect_used_for_promotion": False,
            "validation_thresholds_relaxed": False,
            "green_ci_laundered_as_scientific_pressure_pass": False,
            "selected_center_interface_laundered_as_source_ball_certificate": False,
            "rejected_pressure_laundered_into_full_same_cycle_defect": False,
            "component_operator_metric_laundered_as_pde_residual": False,
            "agent4_replay_laundered_as_final_independent_pde_validation": False,
        },
        "blockers": [
            "Agent-1 selected pressure public value path fails the independent p_eta audit on the Lambda^{-1/2} scale",
            "source radius-one-ball pressure/mixed and Y*u_Y bounds are not machine certified",
            "source R1/R2 -> M,K -> B0/T_sh -> PA.16 -> global matched pressure/leading are not closed",
            "full same-cycle composite theta/axial defect has not been materialized",
            "signed mean inverse and public correction velocity remain closed",
            "no composite velocity/pressure/restricted-forcing artifact exists",
            "no formal held-out full-domain NS residual has been assessed",
        ],
    }
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    return checkpoint


def validate_checkpoint(checkpoint: dict[str, Any]) -> None:
    if checkpoint.get("schema") != SCHEMA or checkpoint.get("task") != TASK:
        raise ValueError("wrong Agent-5 v45 schema/task")
    if checkpoint.get("checkpoint_sha256") != _canonical_sha256(checkpoint):
        raise ValueError("checkpoint SHA256 mismatch")

    states = checkpoint["states"]
    required_true = (
        "selected_pressure_primitive_interface_exposed",
        "selected_pressure_public_value_path_repair_required",
        "oscillatory_ready",
        "public_oscillatory_time_derivative_independently_validated",
        "correction_ingest_allowed",
        "compact_radial_stress_operator_independently_admitted",
        "independent_spacetime_radial_stress_generalization_admitted",
        "finite_correction_cycle_radial_operator_reuse_allowed_without_retuning",
    )
    for key in required_true:
        if states.get(key) is not True:
            raise ValueError(f"required retained/admitted state is false: {key}")

    required_false = (
        "leading_ready",
        "selected_pressure_primitive_independent_preflight_passed",
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "real_full_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "public_velocity_correction_materialized",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    )
    for key in required_false:
        if states.get(key) is not False:
            raise ValueError(f"downstream/rejected state promoted too early: {key}")

    if checkpoint["formal_gates_unchanged"] != FORMAL_GATES:
        raise ValueError("formal final gates changed")

    audit = checkpoint["upstream"]["agent4_selected_pressure_audit"]
    if audit.get("head") != AGENT4_PRESSURE_AUDIT_HEAD:
        raise ValueError("wrong Agent-4 pressure audit head")
    if audit.get("artifact_id") != AGENT4_PRESSURE_ARTIFACT_ID:
        raise ValueError("wrong Agent-4 pressure audit artifact")
    if audit.get("artifact_zip_digest") != AGENT4_PRESSURE_ARTIFACT_DIGEST:
        raise ValueError("wrong Agent-4 pressure audit artifact digest")
    if audit.get("selected_pressure_primitive_independent_preflight_passed") is not False:
        raise ValueError("scientific pressure rejection was laundered into PASS")
    if set(audit.get("failed_guards", [])) != EXPECTED_FAILED_GUARDS:
        raise ValueError("pressure failed-guard set changed")

    leading = checkpoint["upstream"]["agent1_selected_pressure"]
    for key in (
        "source_pressure_radius_one_ball_bound_machine_certified",
        "source_R1_R2_machine_bound",
        "source_M_K_machine_bound",
        "source_B0_T_sh_machine_bound",
        "global_matched_pressure_ready",
        "global_leading_velocity_ready",
    ):
        if leading.get(key) is not False:
            raise ValueError(f"selected-center pressure seam laundered into source/global leading: {key}")

    truth = checkpoint["truth_boundary"]
    for key in (
        "free_residual_defined_forcing_used",
        "surrogate_defect_used_for_promotion",
        "validation_thresholds_relaxed",
        "green_ci_laundered_as_scientific_pressure_pass",
        "selected_center_interface_laundered_as_source_ball_certificate",
        "rejected_pressure_laundered_into_full_same_cycle_defect",
        "component_operator_metric_laundered_as_pde_residual",
        "agent4_replay_laundered_as_final_independent_pde_validation",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"truth-boundary violation: {key}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(checkpoint, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
