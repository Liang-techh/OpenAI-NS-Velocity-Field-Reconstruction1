"""Kokuno Agent-5 routing for the independently rejected PA.10 pressure seam.

This is one narrow integration increment.  It reruns Agent 4's public-only audit
of Agent 1's selected-center PA.10 pressure primitive, pins the exact execution
identity, and routes the scientific REJECT into the leading-lane contract.

A green workflow here means the rejection was reproduced and the truth boundary
was preserved.  It does not make the selected pressure primitive independently
admitted, does not promote ``leading_ready``, and does not assess a full-domain
Navier--Stokes residual.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_agent4_pa10_pressure_primitive_audit import run_audit

TASK = "KOKUNO-A5-PA10-PRESSURE-REJECTION-ROUTING-045"
SCHEMA = "kokuno-agent5-pa10-pressure-rejection-routing-v45"

ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622

AGENT1_PR = 606
AGENT1_HEAD = "56d7e47e4033d869f802ee51ed380b1daa4a34a1"
AGENT1_CANDIDATE_SHA256 = "8b5cdfb749d395f25d7c5f09f9200a3b281dd346fcb4ead6c8c5b1583d9ff675"
AGENT1_DEDICATED_RUN = 35434955050
AGENT1_ARTIFACT_ID = 10582048255
AGENT1_ARTIFACT_DIGEST = (
    "sha256:524b0f2b9bb3f59f848fe365ece55579716e63f290ab5f8cc4858b67ad7a7e42"
)

AGENT4_PR = 608
AGENT4_HEAD = "9a36e257691b4757c37194aa09c11205f21a1ddd"
AGENT4_DEDICATED_RUN = 35435721364
AGENT4_ARTIFACT_ID = 10582059362
AGENT4_ARTIFACT_DIGEST = (
    "sha256:408892549b47495142aab4c1055a7a4a40cca035a75e189f87757b68398b614e"
)

# Accepted/frozen interfaces from earlier routing rounds.  Nothing in this
# checkpoint retunes or revalidates them.
AGENT2_VELOCITY_PR = 561
AGENT2_VELOCITY_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
AGENT2_TIME_DERIVATIVE_PR = 579
AGENT2_TIME_DERIVATIVE_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"

EXPECTED_FAILED_GUARDS = [
    "p_y_refinement",
    "Y_Phi_Y_finest_relative_rms",
    "Y_Phi_Y_finest_relative_max",
    "Y_Phi_Y_refinement",
    "p_eta_finest_relative_rms",
    "p_eta_finest_relative_max",
    "p_eta_refinement",
]

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


def _validate_agent4_receipt(audit: dict[str, Any]) -> None:
    if audit.get("schema") != "kokuno-agent4-pa10-pressure-primitive-independent-audit-v1":
        raise ValueError("unexpected Agent-4 pressure audit schema")
    provenance = audit.get("provenance", {})
    expected_provenance = {
        "agent1_pr": AGENT1_PR,
        "agent1_exact_head": AGENT1_HEAD,
        "agent1_candidate_sha256": AGENT1_CANDIDATE_SHA256,
    }
    if provenance != expected_provenance:
        raise ValueError("Agent-4 receipt is not bound to the admitted Agent-1 candidate")

    independence = audit.get("independence_contract", {})
    if independence.get("agent1_private_integral_helper_called") is not False:
        raise ValueError("Agent-4 audit reused Agent-1 private integral helper")
    if independence.get("agent1_axis_state_called_by_validator") is not False:
        raise ValueError("Agent-4 audit reused Agent-1 axis-state oracle")
    if independence.get("agent1_f0_prime_called_by_validator") is not False:
        raise ValueError("Agent-4 audit reused Agent-1 f0-prime oracle")
    if independence.get("agent1_quadrature_used_as_validation_operator") is not False:
        raise ValueError("Agent-4 audit reused Agent-1 quadrature as validator")

    project_gates = audit["frozen_protocol"]["final_project_gates_unchanged"]
    if project_gates != {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_L2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
    }:
        raise ValueError("final project PDE gates changed")

    # This checkpoint intentionally admits the rejection, not the primitive.
    if audit.get("selected_pressure_primitive_independent_preflight_passed") is not False:
        raise ValueError("expected frozen Agent-4 scientific REJECT did not reproduce")
    if audit.get("failed_guards") != EXPECTED_FAILED_GUARDS:
        raise ValueError("Agent-4 pressure failure signature changed")

    truth = audit.get("truth_boundary", {})
    required_false = (
        "selected_center_pressure_primitive_independently_audited",
        "source_pressure_radius_one_ball_bound_assessed",
        "global_pressure_matched",
        "global_leading_profile_reconstructed",
        "full_same_cycle_requested_stress_materialized",
        "public_velocity_correction_materialized",
        "leading_only_ns_residual_assessed",
        "leading_plus_oscillatory_ns_residual_assessed",
        "after_correction_ns_residual_assessed",
        "heldout_ns_momentum_residual_assessed",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    )
    for key in required_false:
        if truth.get(key) is not False:
            raise ValueError(f"Agent-4 receipt over-promoted truth-boundary field: {key}")


def build_checkpoint(audit: dict[str, Any] | None = None) -> dict[str, Any]:
    if audit is None:
        audit = run_audit()
    audit = json.loads(json.dumps(audit, sort_keys=True, allow_nan=False))
    _validate_agent4_receipt(audit)

    results = audit["results"]
    p_y = results["p_Y"]
    y_phi_y = results["Y_Phi_Y"]
    phi_eta = results["Phi_eta"]
    p_eta = results["p_eta"]
    fto = results["fundamental_theorem_Y"]

    checkpoint: dict[str, Any] = {
        "task": TASK,
        "schema": SCHEMA,
        "round": 45,
        "purpose": (
            "route Agent-4's exact public-only rejection of the selected PA.10 pressure primitive "
            "into the typed leading-lane contract without weakening any guard or promoting downstream state"
        ),
        "upstream": {
            "agent1_selected_pressure_primitive": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "candidate_sha256": AGENT1_CANDIDATE_SHA256,
                "dedicated_workflow_run": AGENT1_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "artifact_id": AGENT1_ARTIFACT_ID,
                "artifact_zip_digest": AGENT1_ARTIFACT_DIGEST,
                "public_api": "KokunoPA10SelectedPressurePrimitive.evaluate(Y,eta)",
                "public_api_materialized": True,
                "source_radius_one_ball_pressure_bound_machine_bound": False,
                "source_operator_M_K_machine_bound": False,
                "source_B0_T_sh_machine_bound": False,
                "global_matched_pressure_ready": False,
                "global_leading_velocity_pressure_ready": False,
            },
            "agent4_independent_pressure_audit": {
                "pr": AGENT4_PR,
                "head": AGENT4_HEAD,
                "dedicated_workflow_run": AGENT4_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "artifact_id": AGENT4_ARTIFACT_ID,
                "artifact_zip_digest": AGENT4_ARTIFACT_DIGEST,
                "scientific_verdict": "REJECT",
                "failed_guards": list(audit["failed_guards"]),
                "p_Y": {
                    "relative_rms_ladder": [level["relative_rms"] for level in p_y["levels"]],
                    "finest_relative_max": p_y["levels"][-1]["relative_max"],
                    "refinement_ratios": list(p_y["refinement_ratios"]),
                },
                "Y_Phi_Y": {
                    "relative_rms_ladder": [level["relative_rms"] for level in y_phi_y["levels"]],
                    "finest_absolute_rms": y_phi_y["levels"][-1]["absolute_rms"],
                    "reference_rms": y_phi_y["levels"][-1]["reference_rms"],
                    "refinement_ratios": list(y_phi_y["refinement_ratios"]),
                    "binary64_identifiable_under_frozen_stationary_cloud": False,
                },
                "Phi_eta": {
                    "relative_rms_ladder": [level["relative_rms"] for level in phi_eta["levels"]],
                    "finest_relative_max": phi_eta["levels"][-1]["relative_max"],
                    "refinement_ratios": list(phi_eta["refinement_ratios"]),
                    "independent_channel_passed": True,
                },
                "p_eta": {
                    "relative_rms_ladder": [level["relative_rms"] for level in p_eta["levels"]],
                    "relative_max_ladder": [level["relative_max"] for level in p_eta["levels"]],
                    "refinement_ratios": list(p_eta["refinement_ratios"]),
                    "reference_rms": results["nontrivial_reference_rms"]["p_eta"],
                    "substantive_blocker": True,
                },
                "fundamental_theorem_Y": {
                    "relative_rms_ladder": [level["relative_rms"] for level in fto["levels"]],
                    "finest_relative_max": fto["levels"][-1]["relative_max"],
                },
                "origin_zero_abs_max": results["origin_zero_abs_max"],
                "mutation_relative_rms": copy.deepcopy(results["mutation_relative_rms"]),
            },
            "agent2_frozen_oscillatory": {
                "velocity_pr": AGENT2_VELOCITY_PR,
                "velocity_head": AGENT2_VELOCITY_HEAD,
                "time_derivative_pr": AGENT2_TIME_DERIVATIVE_PR,
                "time_derivative_head": AGENT2_TIME_DERIVATIVE_HEAD,
                "retuned_in_this_round": False,
            },
        },
        "states": {
            "selected_pressure_primitive_public_api_ready": True,
            "selected_pressure_primitive_independently_audited": False,
            "selected_pressure_primitive_independent_preflight_passed": False,
            "agent1_public_pressure_value_path_repair_required": True,
            "leading_ready": False,
            "oscillatory_ready": True,
            "correction_ingest_allowed": True,
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
                "public_pressure_api_exists": True,
                "independent_pressure_preflight_passed": False,
                "ready": False,
                "next_required": (
                    "repair the public g(eta)/pressure value path on the Lambda^{-1/2} scale by evaluating "
                    "the local stationary-point phase difference directly (or an equivalent high-precision/scaled "
                    "representation), rerun the unchanged Agent-4 p_eta value-vs-derivative audit, then separately "
                    "close source radius-one-ball pressure/mixed bounds and Y u_Y before R1/R2 -> M,K -> "
                    "B0/T_sh -> PA.16 -> global matched pressure"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "ready": True,
                "velocity_api": "velocity_osc(x,y,z,t)->[...,3]",
                "time_derivative_api": "velocity_osc_dt(x,y,z,t)->[...,3]",
                "retune_allowed_by_this_checkpoint": False,
            },
            "correction": {
                "producer": "agent3",
                "ingest_allowed": True,
                "ready": False,
                "blocked_by_missing_full_same_cycle_leading_cross_pressure_forcing_defect": True,
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
        "diagnosis": {
            "p_y_formula_identity_substantively_failed": False,
            "p_y_refinement_guard_failed_at_machine_floor": True,
            "Y_Phi_Y_stationary_cloud_is_inconclusive_due_to_binary64_identifiability": True,
            "p_eta_is_substantive_numeric_blocker": True,
            "likely_public_value_path_issue": (
                "global binary64 phase subtraction followed by multiplication by Lambda loses the O(1/Lambda) "
                "stationary-point phase difference and can create a machine-flat clipped g=1 plateau"
            ),
            "gate_lowering_allowed": False,
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
                "pressure_preflight_metrics_are_not_ns_momentum_residuals": True,
            },
        },
        "formal_gates_unchanged": copy.deepcopy(FORMAL_GATES),
        "truth_boundary": {
            "kokuno_reconstruction_is_paper_exact": False,
            "selected_pressure_reject_laundered_as_pass": False,
            "source_pressure_radius_one_ball_bound_claimed": False,
            "global_matched_pressure_claimed": False,
            "global_leading_field_claimed": False,
            "free_residual_defined_forcing_used": False,
            "validation_thresholds_relaxed": False,
            "surrogate_data_promoted": False,
            "agent4_replay_laundered_as_final_independent_pde_validation": False,
        },
        "blockers": [
            "Agent-1 public p_eta value/derivative seam is independently rejected on the frozen Lambda^-1/2-scale audit",
            "source radius-one-ball pressure/mixed bounds and Y u_Y are not machine-bound",
            "source R1/R2 -> M,K -> B0/T_sh -> PA.16 -> global matched pressure is still open",
            "full same-cycle leading/cross + matched-pressure + restricted-forcing defect does not exist",
            "public correction velocity and composite candidate artifact do not exist",
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
    if states.get("selected_pressure_primitive_public_api_ready") is not True:
        raise ValueError("public selected pressure API unexpectedly unavailable")
    if states.get("selected_pressure_primitive_independently_audited") is not False:
        raise ValueError("rejected pressure primitive was laundered as independently admitted")
    if states.get("selected_pressure_primitive_independent_preflight_passed") is not False:
        raise ValueError("rejected pressure primitive was laundered as PASS")
    if states.get("agent1_public_pressure_value_path_repair_required") is not True:
        raise ValueError("pressure value-path repair requirement was dropped")
    if states.get("oscillatory_ready") is not True:
        raise ValueError("previously admitted oscillatory lane regressed")

    required_false = (
        "leading_ready",
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

    if checkpoint.get("formal_gates_unchanged") != FORMAL_GATES:
        raise ValueError("formal project gates changed")
    if checkpoint["upstream"]["agent4_independent_pressure_audit"]["scientific_verdict"] != "REJECT":
        raise ValueError("Agent-4 scientific rejection was not preserved")
    if checkpoint["upstream"]["agent4_independent_pressure_audit"]["failed_guards"] != EXPECTED_FAILED_GUARDS:
        raise ValueError("pressure failure signature changed")

    truth = checkpoint["truth_boundary"]
    for key, value in truth.items():
        if value is not False:
            raise ValueError(f"truth-boundary prohibition weakened: {key}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)
    text = json.dumps(checkpoint, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
