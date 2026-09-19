"""Agent-5 fail-closed routing after the repaired PA.10 pressure re-audit.

Agent 4 PR #617 independently confirms that Agent 1 PR #615 repaired the old
stationary-local public pressure value/derivative mismatch: the unchanged local
FD gate now passes with strong refinement.  The overall pressure preflight still
rejects because the eta fundamental-theorem check hits a floor comparable with
the binary64 eta-coordinate resolution on the Lambda**(-1/2) layer.

This checkpoint promotes only the local repair evidence and routes the next
interface requirement: expose the narrow layer in a typed O(1) local coordinate
(or an equivalent compensated/high-precision coordinate contract), then rerun
an independent audit.  It does not weaken the failed eta-only guard and does not
promote leading/correction/PDE readiness.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_agent4_pa10_local_phase_pressure_reaudit import run_audit

TASK = "KOKUNO-A5-PA10-LOCAL-COORDINATE-ROUTING-046"
SCHEMA = "kokuno-agent5-pa10-local-coordinate-routing-checkpoint-v46"
ROUND = 46

ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622

AGENT1_PR = 615
AGENT1_HEAD = "7bd597398aef8b4ef1671f79e39eb352adacdd0a"
AGENT1_CANDIDATE_SHA256 = "57ad627b55d9efd4c3669796dc80b9f67bc2af5dc42ce65d5a019d19392930e6"

AGENT2_VELOCITY_PR = 561
AGENT2_VELOCITY_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
AGENT2_TIME_DERIVATIVE_PR = 579
AGENT2_TIME_DERIVATIVE_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"
AGENT2_VECTOR_POTENTIAL_PR = 616
AGENT2_VECTOR_POTENTIAL_HEAD = "aa79cb5eb601ad179c79fb8a4f3e21a093d94c31"

AGENT3_SPACETIME_ADMISSION_PR = 607
AGENT3_SPACETIME_ADMISSION_HEAD = "a2fc36a44771cffb173e391d53fb04df220ab638"
AGENT3_DEDICATED_RUN = 35435350975
AGENT3_STANDARD_RUN = 35435350972
AGENT3_ARTIFACT_ID = 10581334240
AGENT3_ARTIFACT_DIGEST = (
    "sha256:68dbb652b633092b766d7ec7ea089d937ed26c4d05d2b45ba4d29fc3077b5d3a"
)

AGENT4_REAUDIT_PR = 617
AGENT4_REAUDIT_HEAD = "31f8545d8f26bac695f1eb39ee379388c6d283a2"
AGENT4_DEDICATED_RUN = 35438452570
AGENT4_STANDARD_RUN = 35438466180
AGENT4_ARTIFACT_ID = 10582813035
AGENT4_ARTIFACT_DIGEST = (
    "sha256:0b290d6108e8e07e45ab67c409360c2b775447cbf79e2a6967c18bf1108e2b0f"
)

EXPECTED_FAILED_GUARDS = {
    "eta_fundamental_theorem_finest_relative_rms",
    "eta_fundamental_theorem_finest_relative_max",
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


def _checked_reaudit() -> dict[str, Any]:
    audit = json.loads(json.dumps(run_audit(), sort_keys=True, allow_nan=False))
    if audit.get("agent1_pr") != AGENT1_PR:
        raise ValueError("wrong Agent-1 PR in repaired pressure re-audit")
    if audit.get("agent1_head") != AGENT1_HEAD:
        raise ValueError("wrong Agent-1 head in repaired pressure re-audit")
    if audit.get("candidate_sha256") != AGENT1_CANDIDATE_SHA256:
        raise ValueError("wrong Agent-1 candidate identity in repaired pressure re-audit")

    expected_project_gates = {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_L2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
    }
    if audit["frozen_protocol"]["final_project_gates_unchanged"] != expected_project_gates:
        raise ValueError("Agent-4 re-audit changed final project gates")
    if audit.get("repaired_public_p_eta_independent_preflight_passed") is not False:
        raise ValueError("expected current eta-only pressure preflight to remain REJECT")
    if set(audit.get("failed_guards", [])) != EXPECTED_FAILED_GUARDS:
        raise ValueError("Agent-4 repaired-pressure failed-guard set changed")

    fd_levels = audit["p_eta_fd6_levels"]
    ratios = audit["p_eta_fd6_refinement_ratios"]
    local_gates = audit["frozen_protocol"]["local_gates"]
    if fd_levels[-1]["relative_rms"] > local_gates["p_eta_relative_rms"]:
        raise ValueError("repaired local p_eta RMS gate regressed")
    if fd_levels[-1]["relative_max"] > local_gates["p_eta_relative_max"]:
        raise ValueError("repaired local p_eta max gate regressed")
    if min(ratios) < local_gates["p_eta_min_refinement_ratio"]:
        raise ValueError("repaired local p_eta refinement regressed")

    fto = audit["eta_fundamental_theorem_levels"][-1]
    if fto["relative_rms"] <= local_gates["eta_fundamental_theorem_relative_rms"]:
        raise ValueError("eta fundamental-theorem RMS blocker unexpectedly disappeared")
    if fto["relative_max"] <= local_gates["eta_fundamental_theorem_relative_max"]:
        raise ValueError("eta fundamental-theorem max blocker unexpectedly disappeared")

    truth = audit["truth_boundary"]
    for key in (
        "source_pressure_radius_one_ball_bound_assessed",
        "source_R1_R2_bounds_assessed",
        "source_MK_assessed",
        "global_pressure_matched",
        "global_leading_profile_reconstructed",
        "full_same_cycle_requested_stress_materialized",
        "public_velocity_correction_materialized",
        "heldout_ns_momentum_residual_assessed",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"Agent-4 re-audit over-promoted downstream truth: {key}")
    return audit


def build_checkpoint() -> dict[str, Any]:
    audit = _checked_reaudit()
    protocol = audit["frozen_protocol"]
    width = float(protocol["lambda_minus_half_width"])
    eta_star = float(protocol["phase_stationary_eta"])
    eta_ulp = float(abs(math.nextafter(eta_star, math.inf) - eta_star))
    ulps_per_width = width / max(eta_ulp, 1.0e-300)
    ulp_over_width = eta_ulp / max(width, 1.0e-300)

    fd_levels = audit["p_eta_fd6_levels"]
    fto_levels = audit["eta_fundamental_theorem_levels"]
    local_gates = protocol["local_gates"]

    checkpoint: dict[str, Any] = {
        "task": TASK,
        "schema": SCHEMA,
        "round": ROUND,
        "purpose": (
            "admit the independently confirmed stationary-local pressure repair, retain the current eta-only "
            "global identity rejection, and route Agent 1 to a typed O(1) local-coordinate pressure interface "
            "without changing any scientific or final PDE gate"
        ),
        "upstream": {
            "agent1_repaired_selected_pressure": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "candidate_sha256": AGENT1_CANDIDATE_SHA256,
                "stationary_local_phase_repair_present": True,
                "source_radius_one_ball_pressure_bound_machine_certified": False,
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
                "new_vector_potential_pr": AGENT2_VECTOR_POTENTIAL_PR,
                "new_vector_potential_head": AGENT2_VECTOR_POTENTIAL_HEAD,
                "vector_potential_self_check_passed": True,
                "vector_potential_independent_agent4_audit_required": True,
                "admitted_velocity_candidate_changed": False,
            },
            "agent3_spacetime_radial_stress_admission": {
                "pr": AGENT3_SPACETIME_ADMISSION_PR,
                "head": AGENT3_SPACETIME_ADMISSION_HEAD,
                "dedicated_workflow_run": AGENT3_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT3_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "artifact_id": AGENT3_ARTIFACT_ID,
                "artifact_zip_digest": AGENT3_ARTIFACT_DIGEST,
                "compact_radial_stress_operator_independently_admitted": True,
                "finite_correction_cycle_radial_operator_reuse_allowed_without_retuning": True,
                "actual_same_cycle_full_defect_still_required": True,
            },
            "agent4_repaired_pressure_reaudit": {
                "pr": AGENT4_REAUDIT_PR,
                "head": AGENT4_REAUDIT_HEAD,
                "dedicated_workflow_run": AGENT4_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT4_STANDARD_RUN,
                "standard_workflow_status_at_checkpoint_freeze": "in_progress",
                "artifact_id": AGENT4_ARTIFACT_ID,
                "artifact_zip_digest": AGENT4_ARTIFACT_DIGEST,
                "repaired_public_p_eta_independent_preflight_passed": False,
                "failed_guards": audit["failed_guards"],
                "p_eta_fd6_levels": fd_levels,
                "p_eta_fd6_refinement_ratios": audit["p_eta_fd6_refinement_ratios"],
                "eta_fundamental_theorem_levels": fto_levels,
                "p_eta_public_reference": audit["p_eta_public_reference"],
                "p_eta_sign_flip_relative_rms": audit["p_eta_sign_flip_relative_rms"],
            },
        },
        "states": {
            "leading_ready": False,
            "selected_pressure_primitive_interface_exposed": True,
            "selected_pressure_stationary_local_value_path_repair_confirmed": True,
            "selected_pressure_local_p_eta_fd_guard_passed": True,
            "selected_pressure_eta_fundamental_theorem_guard_passed": False,
            "selected_pressure_current_eta_api_independent_preflight_passed": False,
            "selected_pressure_old_global_phase_subtraction_blocker_closed": True,
            "selected_pressure_local_coordinate_contract_required": True,
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
                "current_selected_pressure_api": "evaluate(Y, eta)",
                "local_value_derivative_repair_confirmed": True,
                "current_eta_only_global_identity_passed": False,
                "next_required_api": "evaluate_local(Y, s) or equivalent typed compensated coordinate",
                "local_coordinate_definition": "s=sqrt(Lambda)*(eta-eta_*)",
                "physical_coordinate_reconstruction": "eta=eta_*+s/sqrt(Lambda)",
                "required_chain_rule": "p_eta=sqrt(Lambda)*p_s",
                "required_validation": (
                    "preregister and independently audit local-coordinate value/derivative and fundamental-theorem "
                    "identities, plus the physical chain-rule mapping; retain #617 as a REJECT for the current "
                    "binary64 eta-only API rather than deleting or weakening its failed guards"
                ),
                "after_local_coordinate_pass": (
                    "derive coefficient-space radius-one-ball pressure/mixed and Y*u_Y bounds, then close "
                    "R1/R2 -> M,K -> B0/T_sh -> PA.16 -> global matched pressure/leading"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "validator": "agent4",
                "ready": True,
                "velocity_api": "velocity_osc(x,y,z,t)->[...,3]",
                "time_derivative_api": "velocity_osc_dt(x,y,z,t)->[...,3]",
                "candidate_retune_allowed_by_this_checkpoint": False,
                "vector_potential_new_self_check_only": True,
            },
            "correction": {
                "producer": "agent3",
                "ready": False,
                "ingest_allowed": True,
                "radial_operator_reuse_without_retuning": True,
                "next_required": (
                    "wait for the independently accepted global leading/cross + matched-pressure + restricted-forcing "
                    "handoff, then form the actual same-cycle full theta/axial defect before signed inverse/correction"
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
        "pressure_coordinate_diagnosis": {
            "lambda_minus_half_width": width,
            "eta_star": eta_star,
            "eta_binary64_ulp": eta_ulp,
            "representable_ulps_per_layer_width": ulps_per_width,
            "ulp_over_layer_width": ulp_over_width,
            "finest_local_fd_relative_rms": fd_levels[-1]["relative_rms"],
            "finest_local_fd_relative_max": fd_levels[-1]["relative_max"],
            "local_fd_refinement_ratios": audit["p_eta_fd6_refinement_ratios"],
            "finest_eta_fundamental_theorem_relative_rms": fto_levels[-1]["relative_rms"],
            "finest_eta_fundamental_theorem_relative_max": fto_levels[-1]["relative_max"],
            "eta_fundamental_theorem_rms_gate": local_gates[
                "eta_fundamental_theorem_relative_rms"
            ],
            "eta_fundamental_theorem_max_gate": local_gates[
                "eta_fundamental_theorem_relative_max"
            ],
            "coordinate_resolution_hypothesis_only": True,
            "observed_fto_floor_is_same_order_as_ulp_over_width": (
                0.2 <= fto_levels[-1]["relative_rms"] / max(ulp_over_width, 1.0e-300) <= 5.0
            ),
            "formula_error_claimed": False,
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
                "selected_pressure_local_fd_gate_passed": True,
                "selected_pressure_eta_global_identity_passed": False,
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
            "local_fd_pass_laundered_as_full_pressure_preflight_pass": False,
            "eta_fundamental_theorem_failure_deleted_or_weakened": False,
            "coordinate_resolution_hypothesis_laundered_as_proof": False,
            "selected_center_interface_laundered_as_source_ball_certificate": False,
            "component_operator_metric_laundered_as_pde_residual": False,
            "agent4_replay_laundered_as_final_independent_pde_validation": False,
        },
        "blockers": [
            "current binary64 eta-only selected pressure API fails the preregistered global eta fundamental-theorem guards",
            "typed cancellation-safe O(1) local-coordinate pressure API has not yet been independently validated",
            "source radius-one-ball pressure/mixed and Y*u_Y bounds are not machine certified",
            "source R1/R2 -> M,K -> B0/T_sh -> PA.16 -> global matched pressure/leading are not closed",
            "full same-cycle composite theta/axial defect has not been materialized",
            "signed mean inverse and public correction velocity remain closed",
            "no Kokuno composite candidate artifact exists for formal held-out NS validation",
        ],
        "shortest_closure": [
            "Agent 1: expose a versioned local-coordinate or compensated-coordinate selected pressure API without retuning source formulas",
            "Agent 4: independently audit that new API with preregistered local-coordinate FTO + chain-rule checks while retaining #617 as current eta-only REJECT evidence",
            "Agent 1: after independent pressure-interface PASS, close pressure/mixed + Y*u_Y ball bounds through R1/R2 -> M,K -> B0/T_sh -> PA.16 -> global matched leading/pressure",
            "Agent 3: consume the actual same-cycle full defect and reuse the already admitted radial-stress operator; do not rerun radial-operator research",
            "Agent 5: instantiate/save/load/export the composite only after leading and guarded correction coexist, then hand the frozen artifact to Agent 4 for the unchanged formal PDE gate",
        ],
    }
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    return checkpoint


def validate_checkpoint(checkpoint: dict[str, Any]) -> None:
    if checkpoint.get("task") != TASK or checkpoint.get("schema") != SCHEMA:
        raise ValueError("wrong Agent-5 local-coordinate routing identity")
    if checkpoint.get("round") != ROUND:
        raise ValueError("wrong Agent-5 routing round")
    if checkpoint.get("checkpoint_sha256") != _canonical_sha256(checkpoint):
        raise ValueError("checkpoint SHA mismatch")
    if checkpoint.get("formal_gates_unchanged") != FORMAL_GATES:
        raise ValueError("formal final gates changed")

    states = checkpoint["states"]
    required_true = (
        "selected_pressure_primitive_interface_exposed",
        "selected_pressure_stationary_local_value_path_repair_confirmed",
        "selected_pressure_local_p_eta_fd_guard_passed",
        "selected_pressure_old_global_phase_subtraction_blocker_closed",
        "selected_pressure_local_coordinate_contract_required",
        "oscillatory_ready",
        "public_oscillatory_time_derivative_independently_validated",
        "correction_ingest_allowed",
        "compact_radial_stress_operator_independently_admitted",
        "independent_spacetime_radial_stress_generalization_admitted",
        "finite_correction_cycle_radial_operator_reuse_allowed_without_retuning",
    )
    for key in required_true:
        if states.get(key) is not True:
            raise ValueError(f"required admitted state changed: {key}")

    required_false = (
        "leading_ready",
        "selected_pressure_eta_fundamental_theorem_guard_passed",
        "selected_pressure_current_eta_api_independent_preflight_passed",
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
            raise ValueError(f"downstream state over-promoted: {key}")

    a4 = checkpoint["upstream"]["agent4_repaired_pressure_reaudit"]
    if a4.get("head") != AGENT4_REAUDIT_HEAD:
        raise ValueError("wrong Agent-4 repaired-pressure head")
    if a4.get("artifact_id") != AGENT4_ARTIFACT_ID:
        raise ValueError("wrong Agent-4 repaired-pressure artifact")
    if a4.get("artifact_zip_digest") != AGENT4_ARTIFACT_DIGEST:
        raise ValueError("wrong Agent-4 repaired-pressure artifact digest")
    if set(a4.get("failed_guards", [])) != EXPECTED_FAILED_GUARDS:
        raise ValueError("repaired-pressure failed-guard set changed")

    diagnosis = checkpoint["pressure_coordinate_diagnosis"]
    if diagnosis.get("coordinate_resolution_hypothesis_only") is not True:
        raise ValueError("coordinate-resolution diagnosis truth boundary changed")
    if diagnosis.get("observed_fto_floor_is_same_order_as_ulp_over_width") is not True:
        raise ValueError("eta-coordinate resolution evidence moved outside routed order-of-magnitude relation")
    if diagnosis["finest_local_fd_relative_rms"] > 2.0e-3:
        raise ValueError("local p_eta repair regressed")
    if diagnosis["finest_eta_fundamental_theorem_relative_rms"] <= 1.0e-8:
        raise ValueError("expected current eta global-identity blocker disappeared")

    truth = checkpoint["truth_boundary"]
    forbidden_true = (
        "kokuno_reconstruction_is_paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "free_residual_defined_forcing_used",
        "surrogate_defect_used_for_promotion",
        "validation_thresholds_relaxed",
        "green_ci_laundered_as_scientific_pressure_pass",
        "local_fd_pass_laundered_as_full_pressure_preflight_pass",
        "eta_fundamental_theorem_failure_deleted_or_weakened",
        "coordinate_resolution_hypothesis_laundered_as_proof",
        "selected_center_interface_laundered_as_source_ball_certificate",
        "component_operator_metric_laundered_as_pde_residual",
        "agent4_replay_laundered_as_final_independent_pde_validation",
    )
    for key in forbidden_true:
        if truth.get(key) is not False:
            raise ValueError(f"truth-boundary violation: {key}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n")
    print(json.dumps(checkpoint, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
