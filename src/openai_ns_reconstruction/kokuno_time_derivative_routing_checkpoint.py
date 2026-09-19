"""Kokuno Agent-5 routing checkpoint after independent time-derivative admission.

This is a narrow integration increment. Agent 2 PR #579 exposes an exact public
``velocity_osc_dt(x,y,z,t)`` for the already-admitted oscillatory velocity, and
Agent 4 PR #580 independently validates that interface using only public values.

The latest Agent-3 PR #581 attempts the next real candidate-generated
oscillatory mean-stress materialization, but its exact-head dedicated workflow
failed while generating the numerical receipt. This checkpoint therefore
admits only the independently validated time-derivative handoff and keeps every
correction-materialization / full-PDE state fail closed.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

TASK = "KOKUNO-A5-TIME-DERIVATIVE-ROUTING-042"
SCHEMA = "kokuno-agent5-time-derivative-routing-checkpoint-v42"

ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622

# Preserve the previously green Agent-5 correction-ingest checkpoint as the
# routing parent. This round does not revoke its permission to ingest real
# same-cycle correction inputs; it records that materialization still has not
# completed successfully.
PARENT_AGENT5_PR = 573
PARENT_AGENT5_HEAD = "b63d7d534e6453576d8f4d4da8c607f102e9dfe0"
PARENT_AGENT5_DEDICATED_RUN = 35425277883
PARENT_AGENT5_STANDARD_RUN = 35425277852

# Leading lane remains independently blocked on source-dependent operator
# constants and the B0/T_sh -> PA.16 -> matched-pressure closure.
AGENT1_PR = 569
AGENT1_HEAD = "8836af19f77518f4e160fd1e8d29b2eb57dc61e5"
AGENT1_DEDICATED_RUN = 35424195154
AGENT1_ARTIFACT_ID = 10579052096
AGENT1_ARTIFACT_DIGEST = (
    "sha256:3f8a12512b49121a7ed788bf9b0e0bc165baec54d57c88f0543e789101f93e45"
)

# The admitted oscillatory velocity is unchanged from Agent 2 #561. Agent 2
# #579 adds only its exact autonomous time derivative.
ADMITTED_VELOCITY_AGENT2_PR = 561
ADMITTED_VELOCITY_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
AGENT2_TIME_DERIVATIVE_PR = 579
AGENT2_TIME_DERIVATIVE_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"
AGENT2_TIME_DERIVATIVE_DEDICATED_RUN = 35426870214
AGENT2_TIME_DERIVATIVE_STANDARD_RUN = 35426870184
AGENT2_TIME_DERIVATIVE_ARTIFACT_ID = 10579456478
AGENT2_TIME_DERIVATIVE_ARTIFACT_DIGEST = (
    "sha256:460269b1e07ba31ce6b1759a51214e60d70e86be0e9603a12eaaf0e7feaa3b2f"
)

# Agent 4 #580 independently audits only the public velocity/time-derivative
# interfaces. The exact receipt bytes were retrieved from this workflow artifact.
AGENT4_TIME_DERIVATIVE_PR = 580
AGENT4_TIME_DERIVATIVE_HEAD = "6f7edb4d66dcb162a8510a40e4227f48f3c57e28"
AGENT4_TIME_DERIVATIVE_DEDICATED_RUN = 35427466679
AGENT4_TIME_DERIVATIVE_ARTIFACT_ID = 10579707150
AGENT4_TIME_DERIVATIVE_ARTIFACT_DIGEST = (
    "sha256:0fbf8b0666b130f22a3a9c9c8afed9b6665f0b934fffecb89b51f67d40d9395f"
)
AGENT4_RECEIPT_RAW_SHA256 = (
    "4dc74aee7469f14b97973cb2a9bc69f085a5caa5c1d13db6e0ef8fef3213f437"
)
AGENT4_RECEIPT_CANONICAL_SHA256 = (
    "e3f2c9cba6b476070acbca5d5fd2ef7a585e2cba30dcbb6456c9a1ff0f294f1d"
)

TIME_DERIVATIVE_GUARDS = {
    "finest_fd4_relative_rms": 5.0e-6,
    "finest_fd4_relative_max": 1.0e-5,
    "minimum_fd4_refinement_ratio": 8.0,
    "finest_integral_relative_rms": 1.0e-9,
    "finest_integral_relative_max": 1.0e-9,
    "support_exterior_absolute_max": 1.0e-12,
    "minimum_derivative_rms": 1.0e-8,
    "minimum_sign_flip_integral_mismatch": 0.5,
}

TIME_DERIVATIVE_METRICS = {
    "seed": 9173261,
    "local_point_count": 36,
    "integral_point_count": 12,
    "fd4_steps": [0.012, 0.006, 0.003],
    "fd4_relative_rms": [
        2.0047261979050232e-5,
        1.257899731096725e-6,
        7.869644791783462e-8,
    ],
    "fd4_relative_max": [
        2.9449145250738712e-5,
        1.8506934621223007e-6,
        1.158273501818915e-7,
    ],
    "fd4_refinement_ratios": [15.937090599082667, 15.984199597039918],
    "quadrature_orders": [16, 32, 64],
    "integral_relative_rms": [
        3.549257807868084e-14,
        3.5451549812014814e-14,
        3.100674448199227e-14,
    ],
    "integral_relative_max": [
        2.7354502843124008e-14,
        3.475417060545836e-14,
        3.4001632157038114e-14,
    ],
    "support_exterior_absolute_max": 0.0,
    "public_derivative_rms": 1193.7060017755762,
    "public_derivative_sampled_max": 3235.6833663573893,
    "sign_flip_integral_mismatch": 1.9999999999999967,
}

# Latest Agent-3 exact-head attempt. Focused tests compiled and passed, but the
# actual numerical receipt generation step failed and no artifact was uploaded.
# This is a plumbing/runtime blocker until Agent 3 produces a successful exact
# receipt; it is not treated as a scientific PASS or as a materialized defect.
AGENT3_MEAN_STRESS_PR = 581
AGENT3_MEAN_STRESS_HEAD = "c9b9adb5d81960bd39f41e1a242ff122b8fb30b5"
AGENT3_MEAN_STRESS_DEDICATED_RUN = 35427565520
AGENT3_MEAN_STRESS_DEDICATED_CONCLUSION = "failure"
AGENT3_FAILED_STEP = "Generate actual oscillatory mean-stress receipt"
AGENT3_FOCUSED_TESTS_PASSED = True
AGENT3_RECEIPT_GENERATED = False
AGENT3_ARTIFACT_UPLOADED = False

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


def _time_derivative_passes(
    metrics: dict[str, Any] | None = None,
) -> dict[str, bool]:
    m = TIME_DERIVATIVE_METRICS if metrics is None else metrics
    return {
        "fd4_relative_rms": (
            m["fd4_relative_rms"][-1]
            <= TIME_DERIVATIVE_GUARDS["finest_fd4_relative_rms"]
        ),
        "fd4_relative_max": (
            m["fd4_relative_max"][-1]
            <= TIME_DERIVATIVE_GUARDS["finest_fd4_relative_max"]
        ),
        "fd4_refinement": (
            min(m["fd4_refinement_ratios"])
            >= TIME_DERIVATIVE_GUARDS["minimum_fd4_refinement_ratio"]
        ),
        "integral_relative_rms": (
            m["integral_relative_rms"][-1]
            <= TIME_DERIVATIVE_GUARDS["finest_integral_relative_rms"]
        ),
        "integral_relative_max": (
            m["integral_relative_max"][-1]
            <= TIME_DERIVATIVE_GUARDS["finest_integral_relative_max"]
        ),
        "support": (
            m["support_exterior_absolute_max"]
            <= TIME_DERIVATIVE_GUARDS["support_exterior_absolute_max"]
        ),
        "nontriviality": (
            m["public_derivative_rms"]
            >= TIME_DERIVATIVE_GUARDS["minimum_derivative_rms"]
        ),
        "mutation_detection": (
            m["sign_flip_integral_mismatch"]
            >= TIME_DERIVATIVE_GUARDS["minimum_sign_flip_integral_mismatch"]
        ),
    }


def build_checkpoint(
    metrics: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build v42 and fail closed if the frozen A4 derivative receipt no longer passes."""
    used_metrics = copy.deepcopy(TIME_DERIVATIVE_METRICS if metrics is None else metrics)
    derivative_passes = _time_derivative_passes(used_metrics)
    if not all(derivative_passes.values()):
        failed = sorted(k for k, value in derivative_passes.items() if not value)
        raise ValueError(f"independent public time-derivative handoff rejected: {failed}")

    checkpoint: dict[str, Any] = {
        "task": TASK,
        "schema": SCHEMA,
        "round": 42,
        "purpose": (
            "admit the independently audited public oscillatory time-derivative handoff "
            "while recording the failed exact-head Agent-3 mean-stress receipt generation "
            "and keeping correction/composite/PDE states fail closed"
        ),
        "upstream": {
            "parent_agent5_correction_ingest": {
                "pr": PARENT_AGENT5_PR,
                "head": PARENT_AGENT5_HEAD,
                "dedicated_workflow_run": PARENT_AGENT5_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": PARENT_AGENT5_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "correction_ingest_allowed": True,
            },
            "agent1_leading": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "dedicated_workflow_run": AGENT1_DEDICATED_RUN,
                "artifact_id": AGENT1_ARTIFACT_ID,
                "artifact_zip_digest": AGENT1_ARTIFACT_DIGEST,
                "source_operator_M_K_machine_bound": False,
                "source_B0_T_sh_verified": False,
                "selected_pa16_handoff_allowed": False,
                "global_leading_velocity_pressure_ready": False,
            },
            "agent2_public_time_derivative": {
                "pr": AGENT2_TIME_DERIVATIVE_PR,
                "head": AGENT2_TIME_DERIVATIVE_HEAD,
                "admitted_velocity_pr": ADMITTED_VELOCITY_AGENT2_PR,
                "admitted_velocity_head": ADMITTED_VELOCITY_AGENT2_HEAD,
                "public_velocity_api": (
                    "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc"
                ),
                "public_time_derivative_api": (
                    "openai_ns_reconstruction.kokuno_public_oscillatory_time_derivative:velocity_osc_dt"
                ),
                "dedicated_workflow_run": AGENT2_TIME_DERIVATIVE_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT2_TIME_DERIVATIVE_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "artifact_id": AGENT2_TIME_DERIVATIVE_ARTIFACT_ID,
                "artifact_zip_digest": AGENT2_TIME_DERIVATIVE_ARTIFACT_DIGEST,
                "velocity_candidate_changed": False,
                "repository_autonomous_time_law": True,
                "paper_exact": False,
            },
            "agent4_independent_time_derivative_audit": {
                "pr": AGENT4_TIME_DERIVATIVE_PR,
                "head": AGENT4_TIME_DERIVATIVE_HEAD,
                "audited_agent2_head": AGENT2_TIME_DERIVATIVE_HEAD,
                "dedicated_workflow_run": AGENT4_TIME_DERIVATIVE_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "artifact_id": AGENT4_TIME_DERIVATIVE_ARTIFACT_ID,
                "artifact_zip_digest": AGENT4_TIME_DERIVATIVE_ARTIFACT_DIGEST,
                "receipt_raw_sha256": AGENT4_RECEIPT_RAW_SHA256,
                "receipt_canonical_sha256": AGENT4_RECEIPT_CANONICAL_SHA256,
                "guards": copy.deepcopy(TIME_DERIVATIVE_GUARDS),
                "metrics": used_metrics,
                "derived_guard_passes": derivative_passes,
                "scientific_verdict": "PASS_public_time_derivative_interface_only",
                "pressure_assessed": False,
                "forcing_assessed": False,
                "heldout_ns_momentum_residual_assessed": False,
            },
            "agent3_latest_mean_stress_attempt": {
                "pr": AGENT3_MEAN_STRESS_PR,
                "head": AGENT3_MEAN_STRESS_HEAD,
                "dedicated_workflow_run": AGENT3_MEAN_STRESS_DEDICATED_RUN,
                "dedicated_workflow_conclusion": AGENT3_MEAN_STRESS_DEDICATED_CONCLUSION,
                "focused_tests_passed": AGENT3_FOCUSED_TESTS_PASSED,
                "failed_step": AGENT3_FAILED_STEP,
                "receipt_generated": AGENT3_RECEIPT_GENERATED,
                "artifact_uploaded": AGENT3_ARTIFACT_UPLOADED,
                "scientific_promotion_allowed": False,
                "blocker": "exact_head_numeric_receipt_generation_failed_no_artifact",
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
            "oscillatory_requested_stress_component_materialized": False,
            "same_cycle_requested_stress_materialized": False,
            "candidate_numeric_finite_head_mean_debt_materialized": False,
            "real_candidate_defect_consumed": False,
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
                    "machine-bind source-dependent PA.10 inputs and complete "
                    "M,K -> B0/T_sh -> PA.16 -> global matched pressure"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "validator": "agent4",
                "velocity_candidate_head": ADMITTED_VELOCITY_AGENT2_HEAD,
                "velocity_api": "velocity_osc(x,y,z,t)->[...,3]",
                "time_derivative_head": AGENT2_TIME_DERIVATIVE_HEAD,
                "time_derivative_api": "velocity_osc_dt(x,y,z,t)->[...,3]",
                "velocity_admitted": True,
                "time_derivative_independently_admitted": True,
                "candidate_changed_by_derivative_api": False,
                "frozen_do_not_retune": True,
            },
            "correction": {
                "producer": "agent3",
                "ingest_allowed": True,
                "materialization_successful": False,
                "current_blocker": (
                    "Agent-3 #581 exact-head receipt generation failed after focused tests; "
                    "produce a deterministic candidate-generated mean-stress receipt/artifact "
                    "without changing the admitted oscillatory field or scientific gates"
                ),
                "required_next_outputs": [
                    "successful_oscillatory_requested_stress_component_receipt",
                    "candidate_specific_finite_head_mean_debt",
                    "then_DeltaC_over_epsilon_signed_inverse_and_existing_guards",
                ],
                "full_same_cycle_composite_target_blocked_by_leading_cross_terms": True,
            },
            "final_candidate": {
                "ready": False,
                "required_before_instantiation": [
                    "global_leading_velocity_and_matched_pressure",
                    "admitted_oscillatory_velocity_and_time_derivative",
                    "guarded_materialized_mean_radial_correction",
                ],
                "required_artifact_api": [
                    "velocity(x,y,z,t)",
                    "pressure(x,y,z,t)",
                    "restricted_forcing(x,y,z,t)",
                    "save_load",
                    "python_matlab_smoke",
                ],
            },
        },
        "baseline_comparison": {
            "st006": {
                "full_domain_comparable": True,
                "momentum_sampled_max": ST006_MOMENTUM_MAX,
                "momentum_volume_l2": ST006_MOMENTUM_L2,
                "pde_validated": False,
            },
            "kokuno_current": {
                "full_domain_comparable": False,
                "full_ns_momentum_residual_available": False,
                "component_time_derivative_interface_validated": True,
                "do_not_compare_component_derivative_error_to_st006_momentum": True,
            },
        },
        "formal_gates_unchanged": copy.deepcopy(FORMAL_GATES),
        "truth_boundary": {
            "source_or_paper_exact_claimed": False,
            "openai_field_identified": False,
            "free_residual_forcing_allowed": False,
            "surrogate_mean_defect_accepted": False,
            "thresholds_relaxed": False,
            "component_validation_promoted_to_full_pde": False,
            "pde_validated": False,
        },
    }
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    validate_checkpoint(checkpoint)
    return checkpoint


def validate_checkpoint(payload: dict[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or payload.get("task") != TASK:
        raise ValueError("wrong Agent-5 v42 checkpoint identity")
    expected = _canonical_sha256(payload)
    if payload.get("checkpoint_sha256") != expected:
        raise ValueError("checkpoint SHA mismatch")

    states = payload.get("states", {})
    required_true = (
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready",
        "public_oscillatory_preflight_passed",
        "oscillatory_ready",
        "public_oscillatory_time_derivative_handoff_ready",
        "public_oscillatory_time_derivative_independently_validated",
        "correction_receipt_identity_pinned_to_current_pass",
        "correction_ingest_allowed",
    )
    required_false = (
        "leading_ready",
        "oscillatory_requested_stress_component_materialized",
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "real_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    )
    for key in required_true:
        if states.get(key) is not True:
            raise ValueError(f"required admitted state lost: {key}")
    for key in required_false:
        if states.get(key) is not False:
            raise ValueError(f"downstream state promoted prematurely: {key}")

    a4 = payload["upstream"]["agent4_independent_time_derivative_audit"]
    if a4["head"] != AGENT4_TIME_DERIVATIVE_HEAD:
        raise ValueError("Agent-4 audited head changed")
    if a4["artifact_id"] != AGENT4_TIME_DERIVATIVE_ARTIFACT_ID:
        raise ValueError("Agent-4 artifact identity changed")
    if not all(a4["derived_guard_passes"].values()):
        raise ValueError("Agent-4 time derivative guard no longer passes")

    a3 = payload["upstream"]["agent3_latest_mean_stress_attempt"]
    if a3["dedicated_workflow_conclusion"] != "failure":
        raise ValueError("A3 v42 snapshot must preserve observed exact-head failure")
    if a3["receipt_generated"] is not False or a3["artifact_uploaded"] is not False:
        raise ValueError("failed A3 attempt cannot be laundered into materialized evidence")
    if a3["scientific_promotion_allowed"] is not False:
        raise ValueError("failed A3 attempt cannot authorize promotion")

    gates = payload.get("formal_gates_unchanged", {})
    if gates != FORMAL_GATES:
        raise ValueError("formal PDE gates changed")
    truth = payload.get("truth_boundary", {})
    for key in (
        "source_or_paper_exact_claimed",
        "openai_field_identified",
        "free_residual_forcing_allowed",
        "surrogate_mean_defect_accepted",
        "thresholds_relaxed",
        "component_validation_promoted_to_full_pde",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"truth boundary promoted: {key}")


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args(argv)
    if args.check is not None:
        payload = json.loads(args.check.read_text(encoding="utf-8"))
        validate_checkpoint(payload)
        print(payload["checkpoint_sha256"])
        return 0
    payload = build_checkpoint()
    if args.output is not None:
        _write_json(args.output, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
