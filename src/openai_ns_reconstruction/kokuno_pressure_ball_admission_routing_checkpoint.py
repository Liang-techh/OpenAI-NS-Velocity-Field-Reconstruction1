"""Agent-5 fail-closed admission of the independently audited PA.10 pressure-ball algebra.

This routing increment consumes Agent 4 PR #638's independent exact-arithmetic
cross-audit of Agent 1 PR #635's public conditional PA.10 pressure-ball bound
calculator.  It admits only the coefficient/operator algebra.  Diagnostic ball
inputs are not source radius-one-ball data, so this module deliberately keeps
source ball bounds, R1/R2, M/K, the global leading field, correction, export and
PDE validation closed.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_agent4_pa10_pressure_ball_audit import run_audit

TASK = "KOKUNO-A5-PA10-PRESSURE-BALL-ADMISSION-049"
SCHEMA = "kokuno-agent5-pa10-pressure-ball-admission-routing-checkpoint-v49"
ROUND = 49

AGENT1_PRESSURE_BALL_PR = 635
AGENT1_PRESSURE_BALL_HEAD = "2df7cdb5841613c0990afad188671fff252c96cd"
AGENT4_AUDIT_PR = 638
AGENT4_AUDIT_HEAD = "344fc3ac59952f82074a441dbf747876ba9d18cf"
AGENT4_DEDICATED_RUN = 35443797754
AGENT4_ARTIFACT_ID = 10585293940
AGENT4_ARTIFACT_DIGEST = "sha256:16581e2d580664d6e40dfd55918949f762f7d14110f4a4dd161824fe9acf9140"
AGENT4_REPORT_FILE_SHA256 = "e87be9e8a187c90660ce84685d3438260260ce7d30720f253fb5d5b5217303c5"
AGENT4_RECEIPT_SHA256 = "f5190441c34d12e0bca841037e121bcc6814cbaac36dd188016aef2417de1432"

PRIOR_A5_LOCAL_PRESSURE_PR = 639
PRIOR_A5_LOCAL_PRESSURE_HEAD = "8c8e140e4946b50a46e5aa68d51477e941ab0c9c"

LATEST_AGENT1_SOURCE_AXIS_PR = 644
LATEST_AGENT1_SOURCE_AXIS_HEAD = "1bd16378651f6441ec71c6cf1ae4c3abccb70cc5"
LATEST_AGENT2_BUNDLE_PR = 643
LATEST_AGENT2_BUNDLE_HEAD = "dfea7e7297d07154adbb59c30e74c81410cf0d58"
LATEST_AGENT3_LEDGER_PR = 645
LATEST_AGENT3_LEDGER_HEAD = "ce93d1681f6ad879b9ca2d33e595ca479cacbf9f"

ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622
FORMAL_GATES = {
    "held_out_normalized_momentum_max": 1.0e-3,
    "held_out_normalized_momentum_l2": 1.0e-3,
    "held_out_divergence_max": 1.0e-5,
    "held_out_divergence_l2": 1.0e-5,
}


def _canonical_sha256(payload: dict[str, Any]) -> str:
    unsigned = copy.deepcopy(payload)
    unsigned.pop("checkpoint_sha256", None)
    raw = json.dumps(
        unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _checked_a4_audit() -> dict[str, Any]:
    receipt = json.loads(json.dumps(run_audit(), sort_keys=True, allow_nan=False))
    if receipt.get("schema") != "kokuno-agent4-pa10-pressure-ball-independent-audit-v1":
        raise ValueError("unexpected Agent-4 pressure-ball audit schema")
    upstream = receipt.get("upstream", {})
    if upstream.get("agent1_head") != AGENT1_PRESSURE_BALL_HEAD:
        raise ValueError("Agent-4 audit points at the wrong Agent-1 pressure-ball head")
    if receipt.get("failed_guards") != []:
        raise ValueError(f"Agent-4 pressure-ball audit failed: {receipt.get('failed_guards')}")
    if receipt.get("pressure_ball_operator_algebra_independent_preflight_passed") is not True:
        raise ValueError("Agent-4 pressure-ball preflight did not pass")
    if receipt.get("case_count") != 37:
        raise ValueError("Agent-4 frozen case count changed")
    if receipt.get("eta_factor_underbounds") != [] or receipt.get("public_pressure_underbounds") != []:
        raise ValueError("Agent-4 found a public pressure underbound")
    if receipt.get("minimum_public_to_independent_upper_ratio", 0.0) < 1.0:
        raise ValueError("public pressure bound is not outward relative to independent oracle")
    mutation = receipt.get("downward_mutation", {})
    if mutation.get("detected") != mutation.get("total") or mutation.get("total") != 296:
        raise ValueError("one-percent downward mutation guard is not fully sensitive")
    source = receipt.get("source_weight_audit", {})
    if source.get("failures") != []:
        raise ValueError("source weight identity audit failed")
    if source.get("first_global_algebraic_guard") is not True:
        raise ValueError("source first global factor guard failed")
    if source.get("second_global_factor_guard") is not True:
        raise ValueError("source eta-I global factor guard failed")
    protocol = receipt.get("frozen_protocol", {})
    if protocol.get("final_project_gates_unchanged") != {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_L2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
    }:
        raise ValueError("Agent-4 final project gates changed")
    if receipt.get("receipt_sha256") != AGENT4_RECEIPT_SHA256:
        raise ValueError("Agent-4 deterministic receipt identity changed")

    truth = receipt.get("truth_boundary", {})
    if truth.get("pressure_ball_operator_algebra_independently_audited") is not True:
        raise ValueError("Agent-4 did not mark the algebra independently audited")
    for key in (
        "diagnostic_inputs_are_source_radius_one_ball_bounds",
        "source_rho_machine_bound",
        "source_g_coefficient_norm_machine_bound",
        "source_Phi_radius_one_ball_norm_machine_bound",
        "source_Phi_radius_one_ball_lipschitz_machine_bound",
        "source_pressure_radius_one_ball_norm_machine_bound",
        "source_pressure_radius_one_ball_lipschitz_machine_bound",
        "source_R1_R2_machine_bound",
        "source_operator_M_K_machine_bound",
        "global_pressure_matched",
        "global_leading_profile_reconstructed",
        "heldout_ns_momentum_residual_assessed",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"Agent-4 truth boundary over-promoted: {key}")
    return receipt


def build_checkpoint() -> dict[str, Any]:
    audit = _checked_a4_audit()
    checkpoint: dict[str, Any] = {
        "task": TASK,
        "schema": SCHEMA,
        "round": ROUND,
        "purpose": (
            "formally admit Agent 4 PR #638's independent PASS of Agent 1 PR #635's conditional PA.10 "
            "pressure-ball coefficient/operator algebra, without treating diagnostic inputs as source ball data"
        ),
        "upstream": {
            "agent1_pressure_ball": {
                "pr": AGENT1_PRESSURE_BALL_PR,
                "head": AGENT1_PRESSURE_BALL_HEAD,
                "role": "conditional coefficient-space PA.10 pressure-ball calculator",
            },
            "agent4_independent_pressure_ball_audit": {
                "pr": AGENT4_AUDIT_PR,
                "head": AGENT4_AUDIT_HEAD,
                "dedicated_workflow_run": AGENT4_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "artifact_id": AGENT4_ARTIFACT_ID,
                "artifact_zip_digest": AGENT4_ARTIFACT_DIGEST,
                "artifact_report_file_sha256": AGENT4_REPORT_FILE_SHA256,
                "receipt_sha256": AGENT4_RECEIPT_SHA256,
                "fresh_replay_failed_guards": audit["failed_guards"],
                "fresh_replay_passed": audit[
                    "pressure_ball_operator_algebra_independent_preflight_passed"
                ],
            },
            "prior_agent5_local_pressure_admission": {
                "pr": PRIOR_A5_LOCAL_PRESSURE_PR,
                "head": PRIOR_A5_LOCAL_PRESSURE_HEAD,
                "historical_selected_center_only": True,
                "transfers_to_agent1_644_new_sigma": False,
            },
            "newer_sibling_work": {
                "agent1_source_axis_domain": {
                    "pr": LATEST_AGENT1_SOURCE_AXIS_PR,
                    "head": LATEST_AGENT1_SOURCE_AXIS_HEAD,
                    "role": "source-compatible sigma/rho/normalized-g existence choice",
                    "independent_agent4_audit_required": True,
                    "new_sigma_selected_center_pressure_preflight_required": True,
                    "admitted_by_this_checkpoint": False,
                },
                "agent2_public_bundle": {
                    "pr": LATEST_AGENT2_BUNDLE_PR,
                    "head": LATEST_AGENT2_BUNDLE_HEAD,
                    "role": "non-retuning composability bundle for frozen oscillatory public interfaces",
                    "admitted_by_this_checkpoint": False,
                },
                "agent3_finite_cycle_ledger": {
                    "pr": LATEST_AGENT3_LEDGER_PR,
                    "head": LATEST_AGENT3_LEDGER_HEAD,
                    "role": "engineering multi-step actual-defect ledger; no real full candidate consumed",
                    "admitted_by_this_checkpoint": False,
                },
            },
        },
        "independent_pressure_ball_metrics": {
            "case_count": audit["case_count"],
            "minimum_public_to_independent_upper_ratio": audit[
                "minimum_public_to_independent_upper_ratio"
            ],
            "public_to_independent_product_constant_ratio": audit[
                "public_to_independent_product_constant_ratio"
            ],
            "source_max_first_ratio_on_lattice": audit["source_weight_audit"][
                "max_first_ratio_on_lattice"
            ],
            "source_max_rho_times_eta_I_ratio_on_lattice": audit["source_weight_audit"][
                "max_rho_times_eta_I_ratio_on_lattice"
            ],
            "downward_mutation_detected": audit["downward_mutation"]["detected"],
            "downward_mutation_total": audit["downward_mutation"]["total"],
            "public_underbounds": len(audit["public_pressure_underbounds"]),
            "eta_factor_underbounds": len(audit["eta_factor_underbounds"]),
        },
        "states": {
            "leading_ready": False,
            "selected_pressure_local_coordinate_independently_admitted": True,
            "pressure_ball_operator_algebra_independent_preflight_passed": True,
            "pressure_ball_operator_algebra_independently_admitted": True,
            "conditional_pressure_ball_operator_algebra_ready": True,
            "diagnostic_pressure_ball_inputs_are_source_radius_one_ball_bounds": False,
            "source_axis_domain_candidate_available_in_agent1_644": True,
            "source_axis_domain_independently_admitted": False,
            "source_rho_machine_bound_in_unified_pipeline": False,
            "source_g_coefficient_norm_machine_bound_in_unified_pipeline": False,
            "source_Phi_radius_one_ball_norm_ready": False,
            "source_Phi_radius_one_ball_lipschitz_ready": False,
            "source_mixed_Y_Phi_Y_Y_u_Y_bounds_ready": False,
            "source_pressure_radius_one_ball_bounds_ready": False,
            "source_R1_R2_MK_ready": False,
            "global_matched_pressure_ready": False,
            "global_leading_velocity_ready": False,
            "oscillatory_ready": True,
            "public_oscillatory_time_derivative_independently_validated": True,
            "correction_ingest_allowed": True,
            "compact_radial_stress_operator_independently_admitted": True,
            "same_cycle_momentum_defect_contract_ready": True,
            "actual_momentum_defect_formula_executable": True,
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
            "leading_pressure": {
                "producer": "agent1",
                "validator": "agent4",
                "conditional_pressure_ball_operator_algebra_ready": True,
                "source_ball_inputs_ready": False,
                "global_matched_pressure_ready": False,
                "next_required": (
                    "independently audit Agent-1 #644 source-axis-domain sigma/rho/normalized-g contract; rebuild and "
                    "preflight selected-center pressure at its new PA.8-admissible sigma; bind source Phi radius-one-ball "
                    "norm/Lipschitz and mixed Y Phi_Y / Y u_Y seams; then R1/R2 -> M,K -> B0/T_sh -> PA.16 -> global leading/pressure"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "validator": "agent4",
                "ready": True,
                "candidate_retune_allowed_by_this_checkpoint": False,
            },
            "correction": {
                "producer": "agent3",
                "ready": False,
                "radial_operator_reuse_without_retuning": True,
                "same_cycle_defect_contract_ready": True,
                "reason": "global leading/cross + matched-pressure + independently restricted-forcing inputs do not yet exist",
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
            "pressure_ball_audit_is_comparable_to_st006": False,
            "reason": "pressure-ball outward-bound audit metrics are operator-certificate metrics, not NS momentum residuals",
        },
        "formal_gates_unchanged": dict(FORMAL_GATES),
        "truth_boundary": {
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "free_residual_defined_forcing_allowed": False,
            "surrogate_defect_may_promote_correction": False,
            "conditional_operator_algebra_may_substitute_for_source_ball_inputs": False,
            "old_selected_center_pressure_pass_may_transfer_to_new_sigma": False,
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
    a4 = upstream.get("agent4_independent_pressure_ball_audit", {})
    if a4.get("pr") != AGENT4_AUDIT_PR or a4.get("head") != AGENT4_AUDIT_HEAD:
        raise ValueError("wrong Agent-4 audit identity")
    if a4.get("dedicated_workflow_run") != AGENT4_DEDICATED_RUN:
        raise ValueError("wrong Agent-4 dedicated workflow identity")
    if a4.get("dedicated_workflow_conclusion") != "success":
        raise ValueError("Agent-4 dedicated audit was not green")
    if a4.get("artifact_id") != AGENT4_ARTIFACT_ID:
        raise ValueError("wrong Agent-4 artifact identity")
    if a4.get("artifact_zip_digest") != AGENT4_ARTIFACT_DIGEST:
        raise ValueError("wrong Agent-4 artifact digest")
    if a4.get("receipt_sha256") != AGENT4_RECEIPT_SHA256:
        raise ValueError("wrong Agent-4 receipt identity")
    if a4.get("fresh_replay_failed_guards") != [] or a4.get("fresh_replay_passed") is not True:
        raise ValueError("fresh Agent-4 pressure-ball replay did not pass")

    metrics = checkpoint.get("independent_pressure_ball_metrics", {})
    if metrics.get("case_count") != 37:
        raise ValueError("wrong frozen pressure-ball case count")
    if metrics.get("minimum_public_to_independent_upper_ratio", 0.0) < 1.0:
        raise ValueError("pressure-ball public bound under independent oracle")
    if metrics.get("downward_mutation_detected") != 296 or metrics.get("downward_mutation_total") != 296:
        raise ValueError("pressure-ball mutation sensitivity changed")
    if metrics.get("public_underbounds") != 0 or metrics.get("eta_factor_underbounds") != 0:
        raise ValueError("pressure-ball underbound detected")

    states = checkpoint.get("states", {})
    required_true = (
        "selected_pressure_local_coordinate_independently_admitted",
        "pressure_ball_operator_algebra_independent_preflight_passed",
        "pressure_ball_operator_algebra_independently_admitted",
        "conditional_pressure_ball_operator_algebra_ready",
        "oscillatory_ready",
        "same_cycle_momentum_defect_contract_ready",
    )
    for key in required_true:
        if states.get(key) is not True:
            raise ValueError(f"required admitted state missing: {key}")

    required_false = (
        "leading_ready",
        "diagnostic_pressure_ball_inputs_are_source_radius_one_ball_bounds",
        "source_axis_domain_independently_admitted",
        "source_rho_machine_bound_in_unified_pipeline",
        "source_g_coefficient_norm_machine_bound_in_unified_pipeline",
        "source_Phi_radius_one_ball_norm_ready",
        "source_Phi_radius_one_ball_lipschitz_ready",
        "source_mixed_Y_Phi_Y_Y_u_Y_bounds_ready",
        "source_pressure_radius_one_ball_bounds_ready",
        "source_R1_R2_MK_ready",
        "global_matched_pressure_ready",
        "global_leading_velocity_ready",
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
            raise ValueError(f"forbidden state promotion: {key}")

    siblings = upstream.get("newer_sibling_work", {})
    if siblings.get("agent1_source_axis_domain", {}).get("admitted_by_this_checkpoint") is not False:
        raise ValueError("Agent-1 #644 must not be admitted by this checkpoint")
    if siblings.get("agent1_source_axis_domain", {}).get("new_sigma_selected_center_pressure_preflight_required") is not True:
        raise ValueError("new-sigma local-pressure preflight requirement was lost")

    truth = checkpoint.get("truth_boundary", {})
    for key in (
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "free_residual_defined_forcing_allowed",
        "surrogate_defect_may_promote_correction",
        "conditional_operator_algebra_may_substitute_for_source_ball_inputs",
        "old_selected_center_pressure_pass_may_transfer_to_new_sigma",
        "agent4_replay_may_substitute_for_final_independent_pde_validation",
        "scientific_thresholds_changed",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"truth boundary violated: {key}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "checkpoint_sha256": checkpoint["checkpoint_sha256"],
        "pressure_ball_operator_algebra_independently_admitted": checkpoint["states"]["pressure_ball_operator_algebra_independently_admitted"],
        "leading_ready": checkpoint["states"]["leading_ready"],
        "correction_ready": checkpoint["states"]["correction_ready"],
        "pde_validated": checkpoint["states"]["pde_validated"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
