"""Kokuno Agent-5 routing checkpoint after real oscillatory mean-debt materialization.

This integration increment admits only the Agent-3 oscillatory-component finite-head
mean debt produced from the already admitted Agent-2 public field.  It does not
promote a full same-cycle composite defect, signed mean inverse, correction velocity,
finite correction cycle, full candidate artifact, or PDE validation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

TASK = "KOKUNO-A5-OSCILLATORY-MEAN-DEBT-ROUTING-043"
SCHEMA = "kokuno-agent5-oscillatory-mean-debt-routing-checkpoint-v43"

ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622

# Previously admitted public oscillatory field / time derivative.
AGENT2_VELOCITY_PR = 561
AGENT2_VELOCITY_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
AGENT2_TIME_DERIVATIVE_PR = 579
AGENT2_TIME_DERIVATIVE_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"
AGENT4_TIME_DERIVATIVE_PR = 580
AGENT4_TIME_DERIVATIVE_HEAD = "6f7edb4d66dcb162a8510a40e4227f48f3c57e28"
AGENT4_TIME_DERIVATIVE_DEDICATED_RUN = 35427466679
AGENT4_TIME_DERIVATIVE_STANDARD_RUN = 35427466738

# Exact Agent-3 receipt admitted in this round.
AGENT3_MEAN_DEBT_PR = 586
AGENT3_MEAN_DEBT_HEAD = "ca4da20c90ccd7a5046f79fb2e7965ec05b6086f"
AGENT3_MEAN_DEBT_DEDICATED_RUN = 35428770198
AGENT3_MEAN_DEBT_STANDARD_RUN = 35428770163
AGENT3_MEAN_DEBT_ARTIFACT_ID = 10580380247
AGENT3_MEAN_DEBT_ARTIFACT_DIGEST = (
    "sha256:1cdd166f0deea4f90f57c502275abcf54eb55aa9c731a2d3bf58c9a4a57258ce"
)
AGENT3_MEAN_DEBT_RECEIPT_RAW_SHA256 = (
    "de9b4361baef6c9a524bacd8146f5f608fc9b7af98fa10bcc1d7225049905a62"
)
AGENT3_MEAN_DEBT_RECEIPT_CANONICAL_SHA256 = (
    "b7f3dcc9564de1d3eb7f9130686d204353ecf5d37a7b13b9e1554601d216048e"
)
AUTONOMOUS_MISSING_WEIGHT = 0.28410624176179694
OSC_REQUESTED_STRESS_RMS = 4193262.2244209745
OSC_REQUESTED_STRESS_MAX = 9544356.863937614
OSC_COMPONENT_DEBT_RMS = 1191331.9713019556
OSC_COMPONENT_DEBT_MAX = 2711611.358646726
OSC_COMPONENT_IDENTITY_CLOSURE_MAX = 0.0

# Latest leading sibling.  Its selected-data chi certificate is useful planning
# evidence only; the source rho/M_chi and downstream M,K remain unbound.
AGENT1_PR = 588
AGENT1_HEAD = "f9aac44bbef00ec22592d29010a6250a8e674336"
AGENT1_DEDICATED_RUN = 35429769039
AGENT1_SELECTED_RHO = 5.0e-4
AGENT1_SELECTED_M_CHI = 183.72218879365312
AGENT1_SOURCE_RHO_MACHINE_BOUND = False
AGENT1_SOURCE_M_CHI_MACHINE_BOUND = False
AGENT1_SOURCE_M_K_MACHINE_BOUND = False

# A3 #589 is a later component-level radial-closure audit, but its exact-head
# Actions were still in progress when this checkpoint was authored.  It is
# recorded fail-closed and is not used for any promotion in this round.
AGENT3_RADIAL_AUDIT_PR = 589
AGENT3_RADIAL_AUDIT_HEAD = "2271c52412922c8255cede8aee9a01e68ba55a55"
AGENT3_RADIAL_AUDIT_DEDICATED_RUN = 35430099230
AGENT3_RADIAL_AUDIT_STANDARD_RUN = 35430099232
AGENT3_RADIAL_AUDIT_STATUS_AT_FREEZE = "in_progress"

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


def build_checkpoint() -> dict[str, Any]:
    if OSC_COMPONENT_DEBT_RMS <= 0.0 or OSC_COMPONENT_DEBT_MAX <= 0.0:
        raise ValueError("oscillatory component mean debt must be nontrivial")
    if OSC_COMPONENT_IDENTITY_CLOSURE_MAX > 1.0e-12:
        raise ValueError("oscillatory component mean-debt identity no longer closes")
    if abs(AUTONOMOUS_MISSING_WEIGHT - 0.28410624176179694) > 5.0e-15:
        raise ValueError("autonomous finite-head factor drifted")

    checkpoint: dict[str, Any] = {
        "task": TASK,
        "schema": SCHEMA,
        "round": 43,
        "purpose": (
            "pin the first successful candidate-generated oscillatory-component "
            "finite-head mean debt while keeping the full same-cycle correction and PDE gates closed"
        ),
        "upstream": {
            "agent1_latest_leading": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "dedicated_workflow_run": AGENT1_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "selected_rho": AGENT1_SELECTED_RHO,
                "selected_M_chi_upper_bound": AGENT1_SELECTED_M_CHI,
                "selected_data_certificate_only": True,
                "source_rho_machine_bound": AGENT1_SOURCE_RHO_MACHINE_BOUND,
                "source_chi_multiplier_norm_machine_bound": AGENT1_SOURCE_M_CHI_MACHINE_BOUND,
                "source_operator_M_K_machine_bound": AGENT1_SOURCE_M_K_MACHINE_BOUND,
                "global_leading_velocity_pressure_ready": False,
            },
            "agent2_admitted_oscillatory": {
                "velocity_pr": AGENT2_VELOCITY_PR,
                "velocity_head": AGENT2_VELOCITY_HEAD,
                "time_derivative_pr": AGENT2_TIME_DERIVATIVE_PR,
                "time_derivative_head": AGENT2_TIME_DERIVATIVE_HEAD,
                "public_velocity_api": "velocity_osc(x,y,z,t)->[...,3]",
                "public_time_derivative_api": "velocity_osc_dt(x,y,z,t)->[...,3]",
                "velocity_candidate_changed_in_this_round": False,
            },
            "agent4_time_derivative_audit": {
                "pr": AGENT4_TIME_DERIVATIVE_PR,
                "head": AGENT4_TIME_DERIVATIVE_HEAD,
                "dedicated_workflow_run": AGENT4_TIME_DERIVATIVE_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT4_TIME_DERIVATIVE_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "public_time_derivative_preflight_passed": True,
                "full_ns_momentum_assessed": False,
            },
            "agent3_oscillatory_component_mean_debt": {
                "pr": AGENT3_MEAN_DEBT_PR,
                "head": AGENT3_MEAN_DEBT_HEAD,
                "dedicated_workflow_run": AGENT3_MEAN_DEBT_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT3_MEAN_DEBT_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "artifact_id": AGENT3_MEAN_DEBT_ARTIFACT_ID,
                "artifact_zip_digest": AGENT3_MEAN_DEBT_ARTIFACT_DIGEST,
                "receipt_raw_sha256": AGENT3_MEAN_DEBT_RECEIPT_RAW_SHA256,
                "receipt_canonical_sha256": AGENT3_MEAN_DEBT_RECEIPT_CANONICAL_SHA256,
                "autonomous_missing_weight": AUTONOMOUS_MISSING_WEIGHT,
                "requested_stress_rms": OSC_REQUESTED_STRESS_RMS,
                "requested_stress_max_abs": OSC_REQUESTED_STRESS_MAX,
                "oscillatory_component_debt_rms": OSC_COMPONENT_DEBT_RMS,
                "oscillatory_component_debt_max_abs": OSC_COMPONENT_DEBT_MAX,
                "identity_closure_max_abs": OSC_COMPONENT_IDENTITY_CLOSURE_MAX,
                "real_oscillatory_self_defect_component_consumed": True,
                "oscillatory_requested_stress_component_materialized": True,
                "oscillatory_component_finite_head_mean_debt_materialized": True,
                "autonomous_factor_is_formal_theorem_missing_weight": False,
                "full_same_cycle_composite_requested_stress_materialized": False,
                "candidate_full_finite_head_mean_debt_materialized": False,
                "signed_mean_inverse_input_ready": False,
            },
            "agent3_later_radial_closure_audit": {
                "pr": AGENT3_RADIAL_AUDIT_PR,
                "head": AGENT3_RADIAL_AUDIT_HEAD,
                "dedicated_workflow_run": AGENT3_RADIAL_AUDIT_DEDICATED_RUN,
                "standard_workflow_run": AGENT3_RADIAL_AUDIT_STANDARD_RUN,
                "status_at_checkpoint_freeze": AGENT3_RADIAL_AUDIT_STATUS_AT_FREEZE,
                "scientific_promotion_allowed": False,
            },
        },
        "states": {
            "leading_ready": False,
            "public_provenance_labelled_xyz_t_oscillatory_velocity_ready": True,
            "public_oscillatory_preflight_passed": True,
            "oscillatory_ready": True,
            "public_oscillatory_time_derivative_handoff_ready": True,
            "public_oscillatory_time_derivative_independently_validated": True,
            "correction_receipt_identity_pinned_to_current_pass": True,
            "correction_ingest_allowed": True,
            "oscillatory_requested_stress_component_materialized": True,
            "oscillatory_component_finite_head_mean_debt_materialized": True,
            "same_cycle_requested_stress_materialized": False,
            "candidate_numeric_finite_head_mean_debt_materialized": False,
            "real_full_candidate_defect_consumed": False,
            "signed_mean_inverse_input_ready": False,
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
                "ready": False,
                "next_required": (
                    "machine-bind exact displayed R1/R2 radius-one-ball norm/Lipschitz bounds, "
                    "then source M,K -> B0/T_sh -> PA.16 -> global matched pressure"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "validator": "agent4",
                "ready": True,
                "velocity_api": "velocity_osc(x,y,z,t)->[...,3]",
                "time_derivative_api": "velocity_osc_dt(x,y,z,t)->[...,3]",
            },
            "correction": {
                "producer": "agent3",
                "ingest_allowed": True,
                "oscillatory_component_requested_stress_materialized": True,
                "oscillatory_component_mean_debt_materialized": True,
                "full_same_cycle_requested_stress_materialized": False,
                "full_candidate_mean_debt_materialized": False,
                "signed_mean_inverse_input_ready": False,
                "ready": False,
                "next_required": (
                    "wait for Agent-1 leading/cross + matched-pressure handoff, form the full same-cycle "
                    "theta/axial defect, then run compact radial stress and only then enter DeltaC/epsilon "
                    "signed inverse, bounded-inverse, budget, spacetime, radial, quadratic and joint-gain guards"
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
        "baseline_vs_kokuno": {
            "st006": {
                "same_protocol_full_domain_baseline": True,
                "momentum_sampled_max": ST006_MOMENTUM_MAX,
                "momentum_volume_l2": ST006_MOMENTUM_L2,
                "pde_validated": False,
            },
            "kokuno": {
                "comparable_full_domain_momentum_receipt_exists": False,
                "component_requested_stress_rms": OSC_REQUESTED_STRESS_RMS,
                "component_mean_debt_rms": OSC_COMPONENT_DEBT_RMS,
                "component_metrics_are_not_ns_residuals": True,
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
            "component_debt_laundered_as_full_same_cycle_debt": False,
            "component_metric_laundered_as_pde_residual": False,
        },
        "blockers": [
            "Agent-1 global leading/cross terms and matched pressure are not ready",
            "full same-cycle composite requestedStress has not been materialized",
            "candidate full finite-head mean debt has not been materialized",
            "signed mean inverse and public correction velocity remain closed",
            "no composite velocity/pressure/restricted-forcing artifact exists",
            "no formal held-out full-domain NS residual has been assessed",
        ],
    }
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    return checkpoint


def validate_checkpoint(checkpoint: dict[str, Any]) -> None:
    if checkpoint.get("schema") != SCHEMA or checkpoint.get("task") != TASK:
        raise ValueError("wrong Agent-5 v43 schema/task")
    expected = _canonical_sha256(checkpoint)
    if checkpoint.get("checkpoint_sha256") != expected:
        raise ValueError("checkpoint SHA256 mismatch")
    s = checkpoint["states"]
    required_true = (
        "oscillatory_ready",
        "public_oscillatory_time_derivative_independently_validated",
        "correction_ingest_allowed",
        "oscillatory_requested_stress_component_materialized",
        "oscillatory_component_finite_head_mean_debt_materialized",
    )
    for key in required_true:
        if s.get(key) is not True:
            raise ValueError(f"required admitted state is false: {key}")
    required_false = (
        "leading_ready",
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "real_full_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    )
    for key in required_false:
        if s.get(key) is not False:
            raise ValueError(f"downstream state promoted too early: {key}")
    if checkpoint["formal_gates_unchanged"] != FORMAL_GATES:
        raise ValueError("formal final gates changed")
    a3 = checkpoint["upstream"]["agent3_oscillatory_component_mean_debt"]
    if a3["head"] != AGENT3_MEAN_DEBT_HEAD:
        raise ValueError("wrong Agent-3 mean-debt head")
    if a3["artifact_zip_digest"] != AGENT3_MEAN_DEBT_ARTIFACT_DIGEST:
        raise ValueError("wrong Agent-3 mean-debt artifact digest")
    if a3["receipt_raw_sha256"] != AGENT3_MEAN_DEBT_RECEIPT_RAW_SHA256:
        raise ValueError("wrong Agent-3 mean-debt receipt hash")
    if a3["full_same_cycle_composite_requested_stress_materialized"] is not False:
        raise ValueError("component requestedStress laundered into full composite")
    later = checkpoint["upstream"]["agent3_later_radial_closure_audit"]
    if later["status_at_checkpoint_freeze"] != "in_progress":
        raise ValueError("later A3 audit status unexpectedly changed in frozen v43")
    if later["scientific_promotion_allowed"] is not False:
        raise ValueError("in-progress sibling audit cannot promote science")


def _write(path: Path) -> dict[str, Any]:
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n")
    return checkpoint


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    if bool(args.output) == bool(args.check):
        parser.error("provide exactly one of --output or --check")
    if args.output:
        checkpoint = _write(args.output)
        print(json.dumps(checkpoint, indent=2, sort_keys=True))
        return
    checkpoint = json.loads(args.check.read_text())
    validate_checkpoint(checkpoint)
    print(checkpoint["checkpoint_sha256"])


if __name__ == "__main__":
    main()
