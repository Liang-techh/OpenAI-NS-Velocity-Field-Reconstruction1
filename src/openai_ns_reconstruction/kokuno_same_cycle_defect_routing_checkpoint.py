"""Agent-5 fail-closed routing for the typed same-cycle momentum-defect seam.

Agent 3 PR #627 exposes a narrow engineering contract that computes the raw
Navier--Stokes momentum defect directly from same-cycle velocity, velocity-time-
derivative, scalar-pressure and restricted-forcing providers.  The public API
does not accept a caller-supplied residual, stress, target, pressure gradient or
gain, and it enforces disjoint held-in / held-out point sets.

This checkpoint admits only that typed integration seam.  It does not claim that
Agent 1 has supplied the missing global leading/cross + matched-pressure +
restricted-forcing providers, does not consume a real full candidate defect,
and does not promote correction or PDE validation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_same_cycle_defect_contract import deterministic_receipt, truth_boundary

TASK = "KOKUNO-A5-SAME-CYCLE-DEFECT-ROUTING-047"
SCHEMA = "kokuno-agent5-same-cycle-defect-routing-checkpoint-v47"
ROUND = 47

ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622

AGENT1_PR = 624
AGENT1_HEAD = "e036b057d9c9fd320e2bb87aeb0fcdb1395f2387"
AGENT1_CANDIDATE_SHA256 = "2c41504eeb448b840280dbb80f40f1209f0efff9d0db2ca0e8826820e73944f6"
AGENT1_DEDICATED_RUN = 35440168657
AGENT1_ARTIFACT_ID = 10583965457
AGENT1_ARTIFACT_DIGEST = (
    "sha256:759fcfc9e6f346969aa9ec3eb5de8f5fdb7c2c2d4a2560cc26f2c915099e8c05"
)

AGENT2_ADMITTED_VELOCITY_PR = 561
AGENT2_ADMITTED_VELOCITY_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
AGENT2_ADMITTED_DT_PR = 579
AGENT2_ADMITTED_DT_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"
AGENT2_VECTOR_POTENTIAL_DT_PR = 625
AGENT2_VECTOR_POTENTIAL_DT_HEAD = "e12bfb0a65031e44b1556443d88d11f2e9e92904"

AGENT3_PR = 627
AGENT3_HEAD = "0504056a6a5c9fb61193c68c1227b83dbd0a63ea"
AGENT3_DEDICATED_RUN = 35441005166
AGENT3_ARTIFACT_ID = 10583039448
AGENT3_ARTIFACT_DIGEST = (
    "sha256:b5fb108742501d80fe1eecccce40e587356ccccc3ee17b8c164f7ebbfa5b346e"
)

AGENT4_LOCAL_PRESSURE_AUDIT_PR = 628
AGENT4_LOCAL_PRESSURE_AUDIT_HEAD = "21ec7cb03fe91dac008160419207d9cbe654eadf"

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


def _checked_agent3_receipt() -> dict[str, Any]:
    receipt = json.loads(json.dumps(deterministic_receipt(), sort_keys=True, allow_nan=False))
    if receipt.get("schema") != "kokuno-a3-same-cycle-defect-contract-v1":
        raise ValueError("unexpected Agent-3 same-cycle-defect receipt schema")
    if receipt.get("task") != "KOKUNO-A3-SAME-CYCLE-DEFECT-CONTRACT-046":
        raise ValueError("unexpected Agent-3 same-cycle-defect task")
    provenance = receipt.get("provenance", {})
    if provenance.get("parent_agent3_pr") != 607:
        raise ValueError("Agent-3 contract parent PR changed")
    if provenance.get("parent_agent3_head") != "a2fc36a44771cffb173e391d53fb04df220ab638":
        raise ValueError("Agent-3 contract parent head changed")

    regression = receipt.get("formula_regression", {})
    if regression.get("expected_residual") != [2.0, 0.0, 3.0]:
        raise ValueError("Agent-3 analytic residual regression changed")
    if not math.isclose(
        float(regression.get("expected_vector_rms", math.nan)),
        math.sqrt(13.0),
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise ValueError("Agent-3 analytic vector RMS regression changed")
    if float(regression.get("maximum_absolute_error", math.inf)) > 2.0e-11:
        raise ValueError("Agent-3 FD4 manufactured residual regression exceeded guard")

    held_in = regression.get("held_in", {})
    held_out = regression.get("held_out", {})
    if held_in.get("point_count") != 2 or held_out.get("point_count") != 2:
        raise ValueError("Agent-3 deterministic split size changed")
    if held_in.get("identity") != held_out.get("identity"):
        raise ValueError("Agent-3 held-in/out cycle identities differ")
    if held_in.get("identity") != {
        "cycle_id": "analytic-contract-regression-v1",
        "cycle_index": 0,
        "state_token": "fixed",
    }:
        raise ValueError("Agent-3 deterministic cycle identity changed")

    expected_truth = truth_boundary()
    if receipt.get("truth_boundary") != expected_truth:
        raise ValueError("Agent-3 truth boundary changed")
    if expected_truth["actual_momentum_defect_formula_executable"] is not True:
        raise ValueError("Agent-3 momentum defect formula is not executable")
    for key in (
        "caller_supplied_residual_allowed",
        "caller_supplied_stress_allowed",
        "caller_supplied_target_allowed",
        "caller_supplied_pressure_gradient_allowed",
        "caller_supplied_gain_allowed",
        "restricted_forcing_semantics_independently_validated_here",
        "full_same_cycle_composite_requested_stress_materialized",
        "candidate_finite_head_mean_debt_materialized",
        "real_full_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "public_velocity_correction_materialized",
        "finite_correction_cycle_run",
        "heldout_ns_momentum_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
    ):
        if expected_truth.get(key) is not False:
            raise ValueError(f"Agent-3 truth boundary over-promoted: {key}")
    for key in (
        "same_cycle_identity_enforced",
        "restricted_forcing_provider_required",
        "heldin_heldout_disjointness_enforced",
    ):
        if expected_truth.get(key) is not True:
            raise ValueError(f"Agent-3 contract prerequisite regressed: {key}")
    if expected_truth.get("final_normalized_momentum_gate") != 1.0e-3:
        raise ValueError("Agent-3 momentum gate changed")
    if expected_truth.get("final_normalized_divergence_gate") != 1.0e-5:
        raise ValueError("Agent-3 divergence gate changed")
    return receipt


def build_checkpoint() -> dict[str, Any]:
    receipt = _checked_agent3_receipt()
    regression = receipt["formula_regression"]
    checkpoint: dict[str, Any] = {
        "task": TASK,
        "schema": SCHEMA,
        "round": ROUND,
        "purpose": (
            "admit Agent 3's typed same-cycle raw momentum-defect contract as reusable integration glue while "
            "keeping the missing full candidate providers, full correction and final PDE gate fail-closed"
        ),
        "upstream": {
            "agent1_local_pressure_coordinate": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "candidate_sha256": AGENT1_CANDIDATE_SHA256,
                "dedicated_workflow_run": AGENT1_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "artifact_id": AGENT1_ARTIFACT_ID,
                "artifact_zip_digest": AGENT1_ARTIFACT_DIGEST,
                "local_coordinate_self_check_passed": True,
                "independent_agent4_local_coordinate_audit_pr": AGENT4_LOCAL_PRESSURE_AUDIT_PR,
                "independent_agent4_local_coordinate_audit_head": AGENT4_LOCAL_PRESSURE_AUDIT_HEAD,
                "independent_agent4_local_coordinate_audit_resolved_at_checkpoint_freeze": False,
                "source_radius_one_ball_pressure_bounds_ready": False,
                "global_matched_pressure_ready": False,
                "global_leading_velocity_ready": False,
            },
            "agent2_admitted_oscillatory": {
                "velocity_pr": AGENT2_ADMITTED_VELOCITY_PR,
                "velocity_head": AGENT2_ADMITTED_VELOCITY_HEAD,
                "velocity_time_derivative_pr": AGENT2_ADMITTED_DT_PR,
                "velocity_time_derivative_head": AGENT2_ADMITTED_DT_HEAD,
                "new_vector_potential_time_derivative_pr": AGENT2_VECTOR_POTENTIAL_DT_PR,
                "new_vector_potential_time_derivative_head": AGENT2_VECTOR_POTENTIAL_DT_HEAD,
                "new_vector_potential_time_derivative_is_self_check_only": True,
                "admitted_velocity_candidate_changed": False,
            },
            "agent3_same_cycle_defect_contract": {
                "pr": AGENT3_PR,
                "head": AGENT3_HEAD,
                "dedicated_workflow_run": AGENT3_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_status_at_checkpoint_freeze": "in_progress",
                "artifact_id": AGENT3_ARTIFACT_ID,
                "artifact_zip_digest": AGENT3_ARTIFACT_DIGEST,
                "receipt_schema": receipt["schema"],
                "formula": receipt["truth_boundary"]["momentum_defect_formula"],
                "spatial_operator": receipt["truth_boundary"]["spatial_derivative_operator"],
                "analytic_expected_residual": regression["expected_residual"],
                "analytic_expected_vector_rms": regression["expected_vector_rms"],
                "analytic_maximum_absolute_error": regression["maximum_absolute_error"],
            },
        },
        "states": {
            "leading_ready": False,
            "selected_pressure_local_coordinate_interface_exposed": True,
            "selected_pressure_local_coordinate_independent_preflight_passed": False,
            "oscillatory_ready": True,
            "public_oscillatory_time_derivative_independently_validated": True,
            "correction_ingest_allowed": True,
            "oscillatory_requested_stress_component_materialized": True,
            "oscillatory_component_finite_head_mean_debt_materialized": True,
            "compact_radial_stress_operator_independently_admitted": True,
            "independent_spacetime_radial_stress_generalization_admitted": True,
            "finite_correction_cycle_radial_operator_reuse_allowed_without_retuning": True,
            "same_cycle_momentum_defect_contract_ready": True,
            "actual_momentum_defect_formula_executable": True,
            "same_cycle_identity_enforced": True,
            "restricted_forcing_provider_required": True,
            "heldin_heldout_disjointness_enforced": True,
            "restricted_forcing_semantics_independently_validated": False,
            "real_full_candidate_defect_consumed": False,
            "same_cycle_requested_stress_materialized": False,
            "candidate_numeric_finite_head_mean_debt_materialized": False,
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
                "current_api": "evaluate_local(Y,s) selected-center pressure coordinate",
                "independent_local_coordinate_audit_pending": True,
                "after_independent_pass": (
                    "derive source/coefficient-space pressure and mixed radius-one-ball bounds, then close "
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
            },
            "same_cycle_defect": {
                "producer": "agent3",
                "ready": True,
                "api": (
                    "evaluate_disjoint_heldin_heldout(velocity,velocity_dt,pressure,restricted_forcing,"
                    "held_in_points,held_out_points,viscosity,spatial_step)"
                ),
                "caller_supplied_residual_allowed": False,
                "caller_supplied_stress_allowed": False,
                "caller_supplied_target_allowed": False,
                "same_cycle_identity_enforced": True,
                "restricted_forcing_provider_required": True,
                "restricted_forcing_truth_must_be_independently_bound_upstream": True,
                "heldin_heldout_disjointness_enforced": True,
                "normalization_or_final_gate_assessed_by_contract": False,
            },
            "correction": {
                "producer": "agent3",
                "ready": False,
                "radial_operator_reuse_without_retuning": True,
                "next_required": (
                    "after independently admitted Agent-1 global leading/cross + matched pressure and restricted "
                    "forcing exist in one immutable cycle identity, evaluate the real full defect with the admitted "
                    "same-cycle contract; only then form theta/axial requestedStress, full finite-head debt, signed "
                    "inverse and guarded correction"
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
            "st006_same_protocol_baseline": {
                "momentum_sampled_max": ST006_MOMENTUM_MAX,
                "momentum_volume_l2": ST006_MOMENTUM_L2,
                "pde_validated": False,
            },
            "kokuno_full_same_protocol_residual_available": False,
            "current_agent3_contract_regression_is_comparable_to_st006": False,
            "reason": (
                "the Agent-3 sqrt(13) receipt is a manufactured analytic interface regression, not a Kokuno "
                "candidate residual and not the repository held-out PDE protocol"
            ),
        },
        "formal_gates_unchanged": dict(FORMAL_GATES),
        "truth_boundary": {
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "free_residual_defined_forcing_allowed": False,
            "surrogate_defect_may_promote_correction": False,
            "agent4_replay_may_substitute_for_final_independent_pde_validation": False,
            "scientific_thresholds_changed": False,
        },
    }
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    return checkpoint


def validate_checkpoint(checkpoint: dict[str, Any]) -> None:
    if checkpoint.get("task") != TASK or checkpoint.get("schema") != SCHEMA:
        raise ValueError("wrong Agent-5 checkpoint identity")
    if checkpoint.get("round") != ROUND:
        raise ValueError("wrong Agent-5 checkpoint round")
    if checkpoint.get("checkpoint_sha256") != _canonical_sha256(checkpoint):
        raise ValueError("checkpoint SHA256 mismatch")
    if checkpoint.get("formal_gates_unchanged") != FORMAL_GATES:
        raise ValueError("formal PDE gates changed")

    upstream = checkpoint.get("upstream", {})
    a3 = upstream.get("agent3_same_cycle_defect_contract", {})
    if a3.get("pr") != AGENT3_PR or a3.get("head") != AGENT3_HEAD:
        raise ValueError("wrong Agent-3 same-cycle contract identity")
    if a3.get("dedicated_workflow_run") != AGENT3_DEDICATED_RUN:
        raise ValueError("wrong Agent-3 dedicated workflow identity")
    if a3.get("dedicated_workflow_conclusion") != "success":
        raise ValueError("Agent-3 dedicated workflow is not green")
    if a3.get("artifact_id") != AGENT3_ARTIFACT_ID:
        raise ValueError("wrong Agent-3 artifact identity")
    if a3.get("artifact_zip_digest") != AGENT3_ARTIFACT_DIGEST:
        raise ValueError("wrong Agent-3 artifact digest")
    if a3.get("analytic_expected_residual") != [2.0, 0.0, 3.0]:
        raise ValueError("manufactured residual changed")
    if float(a3.get("analytic_maximum_absolute_error", math.inf)) > 2.0e-11:
        raise ValueError("manufactured FD4 regression exceeded guard")

    states = checkpoint.get("states", {})
    required_true = (
        "oscillatory_ready",
        "correction_ingest_allowed",
        "compact_radial_stress_operator_independently_admitted",
        "finite_correction_cycle_radial_operator_reuse_allowed_without_retuning",
        "same_cycle_momentum_defect_contract_ready",
        "actual_momentum_defect_formula_executable",
        "same_cycle_identity_enforced",
        "restricted_forcing_provider_required",
        "heldin_heldout_disjointness_enforced",
    )
    for key in required_true:
        if states.get(key) is not True:
            raise ValueError(f"required admitted state is not true: {key}")
    required_false = (
        "leading_ready",
        "selected_pressure_local_coordinate_independent_preflight_passed",
        "restricted_forcing_semantics_independently_validated",
        "real_full_candidate_defect_consumed",
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
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
            raise ValueError(f"premature downstream promotion: {key}")

    defect = checkpoint.get("typed_handoffs", {}).get("same_cycle_defect", {})
    if defect.get("ready") is not True:
        raise ValueError("typed same-cycle defect handoff not admitted")
    if defect.get("caller_supplied_residual_allowed") is not False:
        raise ValueError("caller residual laundering allowed")
    if defect.get("caller_supplied_stress_allowed") is not False:
        raise ValueError("caller stress laundering allowed")
    if defect.get("caller_supplied_target_allowed") is not False:
        raise ValueError("caller target laundering allowed")
    if defect.get("restricted_forcing_truth_must_be_independently_bound_upstream") is not True:
        raise ValueError("restricted forcing truth boundary weakened")
    if defect.get("normalization_or_final_gate_assessed_by_contract") is not False:
        raise ValueError("engineering contract laundered into final validation")

    baseline = checkpoint.get("baseline_vs_kokuno", {})
    if baseline.get("kokuno_full_same_protocol_residual_available") is not False:
        raise ValueError("full Kokuno residual prematurely claimed")
    if baseline.get("current_agent3_contract_regression_is_comparable_to_st006") is not False:
        raise ValueError("manufactured regression mislabelled as ST006-comparable")
    st006 = baseline.get("st006_same_protocol_baseline", {})
    if st006.get("momentum_sampled_max") != ST006_MOMENTUM_MAX:
        raise ValueError("ST006 max baseline changed")
    if st006.get("momentum_volume_l2") != ST006_MOMENTUM_L2:
        raise ValueError("ST006 L2 baseline changed")

    truth = checkpoint.get("truth_boundary", {})
    for key in (
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "free_residual_defined_forcing_allowed",
        "surrogate_defect_may_promote_correction",
        "agent4_replay_may_substitute_for_final_independent_pde_validation",
        "scientific_thresholds_changed",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"truth boundary weakened: {key}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)
    text = json.dumps(checkpoint, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
