"""Agent-5 fail-closed admission of the independently audited PA.10 local pressure API.

This increment consumes Agent 4 PR #628's public-value-only audit of Agent 1 PR
#624's cancellation-safe ``evaluate_local(Y, s)`` pressure interface.  It admits
only that selected-center interface.  It does not turn selected-center numerical
consistency into source radius-one-ball pressure bounds, a global matched
pressure/leading field, a correction, or a Navier--Stokes PDE pass.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_agent4_pa10_local_pressure_coordinate_audit import run_audit

TASK = "KOKUNO-A5-LOCAL-PRESSURE-ADMISSION-048"
SCHEMA = "kokuno-agent5-local-pressure-admission-routing-checkpoint-v48"
ROUND = 48

AGENT1_LOCAL_PR = 624
AGENT1_LOCAL_HEAD = "e036b057d9c9fd320e2bb87aeb0fcdb1395f2387"
AGENT1_LOCAL_CANDIDATE_SHA256 = "2c41504eeb448b840280dbb80f40f1209f0efff9d0db2ca0e8826820e73944f6"
AGENT4_AUDIT_PR = 628
AGENT4_AUDIT_HEAD = "21ec7cb03fe91dac008160419207d9cbe654eadf"
AGENT4_DEDICATED_RUN = 35441057979
AGENT4_ARTIFACT_ID = 10584116650
AGENT4_ARTIFACT_DIGEST = "sha256:48db81bf4ec3fd04cfe66079608807e9de3000cadd0ab12f565cb943361a1cbf"

PRIOR_A5_PR = 629
PRIOR_A5_HEAD = "d984cfa3f3e1c42c049e322a8bd1f0a153861a71"
PRIOR_A5_DEDICATED_RUN = 35441444326
PRIOR_A5_ARTIFACT_ID = 10583897283
PRIOR_A5_ARTIFACT_DIGEST = "sha256:e8fd6ecfdab1ecd137481cfb6bc88c4ccd56396c7a739138e85938a040e755cb"

LATEST_AGENT1_PRESSURE_BALL_PR = 635
LATEST_AGENT1_PRESSURE_BALL_HEAD = "2df7cdb5841613c0990afad188671fff252c96cd"
LATEST_AGENT3_CYCLE_PR = 637
LATEST_AGENT3_CYCLE_HEAD = "f7830a34784127a9ad7894047039b0ea1e8795eb"
LATEST_AGENT4_PRESSURE_BALL_AUDIT_PR = 638
LATEST_AGENT4_PRESSURE_BALL_AUDIT_HEAD = "344fc3ac59952f82074a441dbf747876ba9d18cf"

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
    encoded = json.dumps(
        unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _checked_a4_audit() -> dict[str, Any]:
    receipt = json.loads(json.dumps(run_audit(), sort_keys=True, allow_nan=False))
    if receipt.get("schema") != "kokuno-agent4-pa10-local-pressure-coordinate-independent-audit-v1":
        raise ValueError("unexpected Agent-4 local-pressure audit schema")
    if receipt.get("agent1_pr") != AGENT1_LOCAL_PR:
        raise ValueError("Agent-4 audit points at the wrong Agent-1 PR")
    if receipt.get("agent1_head") != AGENT1_LOCAL_HEAD:
        raise ValueError("Agent-4 audit points at the wrong Agent-1 head")
    if receipt.get("candidate_sha256") != AGENT1_LOCAL_CANDIDATE_SHA256:
        raise ValueError("Agent-1 local-pressure candidate identity changed")
    if receipt.get("failed_guards") != []:
        raise ValueError(f"Agent-4 local-pressure audit failed: {receipt.get('failed_guards')}")
    if receipt.get("local_pressure_coordinate_independent_preflight_passed") is not True:
        raise ValueError("Agent-4 local-pressure preflight did not pass")

    protocol = receipt.get("frozen_protocol", {})
    if protocol.get("final_project_gates_unchanged") != {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_L2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
    }:
        raise ValueError("Agent-4 final project gates changed")

    truth = receipt.get("truth_boundary", {})
    for key in (
        "source_pressure_radius_one_ball_bound_assessed",
        "source_pressure_radius_one_ball_lipschitz_assessed",
        "source_R1_R2_bounds_assessed",
        "source_MK_assessed",
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
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"Agent-4 truth boundary over-promoted: {key}")
    return receipt


def build_checkpoint() -> dict[str, Any]:
    audit = _checked_a4_audit()
    p_s = audit["p_s_fd8_levels"][-1]
    fto = audit["local_fundamental_theorem_levels"][-1]
    chain = audit["physical_chain_rule_from_value_fd8"]
    checkpoint: dict[str, Any] = {
        "task": TASK,
        "schema": SCHEMA,
        "round": ROUND,
        "purpose": (
            "formally admit Agent 4 PR #628's independent PASS for Agent 1 PR #624's selected-center "
            "cancellation-safe local pressure coordinate, without promoting source ball bounds, global leading, "
            "correction, artifact export, or PDE validation"
        ),
        "upstream": {
            "agent1_selected_center_local_pressure": {
                "pr": AGENT1_LOCAL_PR,
                "head": AGENT1_LOCAL_HEAD,
                "candidate_sha256": AGENT1_LOCAL_CANDIDATE_SHA256,
                "api": "KokunoPA10LocalPressureCoordinate.evaluate_local(Y,s)",
            },
            "agent4_independent_local_pressure_audit": {
                "pr": AGENT4_AUDIT_PR,
                "head": AGENT4_AUDIT_HEAD,
                "dedicated_workflow_run": AGENT4_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "artifact_id": AGENT4_ARTIFACT_ID,
                "artifact_zip_digest": AGENT4_ARTIFACT_DIGEST,
                "fresh_replay_failed_guards": audit["failed_guards"],
                "fresh_replay_passed": audit["local_pressure_coordinate_independent_preflight_passed"],
            },
            "prior_agent5_same_cycle_defect_admission": {
                "pr": PRIOR_A5_PR,
                "head": PRIOR_A5_HEAD,
                "dedicated_workflow_run": PRIOR_A5_DEDICATED_RUN,
                "artifact_id": PRIOR_A5_ARTIFACT_ID,
                "artifact_zip_digest": PRIOR_A5_ARTIFACT_DIGEST,
            },
            "newer_sibling_work": {
                "agent1_pressure_ball": {
                    "pr": LATEST_AGENT1_PRESSURE_BALL_PR,
                    "head": LATEST_AGENT1_PRESSURE_BALL_HEAD,
                    "role": "conditional coefficient-space pressure-ball algebra; source numerical inputs still missing",
                    "admitted_by_this_checkpoint": False,
                },
                "agent3_finite_cycle": {
                    "pr": LATEST_AGENT3_CYCLE_PR,
                    "head": LATEST_AGENT3_CYCLE_HEAD,
                    "role": "actual-defect decrease engineering contract; real full candidate still unavailable",
                    "admitted_by_this_checkpoint": False,
                },
                "agent4_pressure_ball_audit": {
                    "pr": LATEST_AGENT4_PRESSURE_BALL_AUDIT_PR,
                    "head": LATEST_AGENT4_PRESSURE_BALL_AUDIT_HEAD,
                    "role": "independent audit of Agent-1 conditional pressure-ball algebra",
                    "resolved_at_checkpoint_freeze": False,
                    "admitted_by_this_checkpoint": False,
                },
            },
        },
        "independent_local_pressure_metrics": {
            "fd8_finest_relative_rms": p_s["relative_rms"],
            "fd8_finest_relative_max": p_s["relative_max"],
            "fd8_refinement_ratios": audit["p_s_fd8_refinement_ratios"],
            "fundamental_theorem_finest_relative_rms": fto["relative_rms"],
            "fundamental_theorem_finest_relative_max": fto["relative_max"],
            "fundamental_theorem_refinement_ratios": audit[
                "local_fundamental_theorem_refinement_ratios"
            ],
            "physical_chain_rule_relative_rms": chain["relative_rms"],
            "physical_chain_rule_relative_max": chain["relative_max"],
            "two_term_eta_mapping_scaled_max": audit["two_term_eta_mapping"][
                "scaled_max_in_layer_widths"
            ],
            "mutation_sign_flip_relative_rms": audit["p_s_sign_flip_relative_rms"],
        },
        "states": {
            "leading_ready": False,
            "selected_pressure_local_coordinate_interface_exposed": True,
            "selected_pressure_local_coordinate_independent_preflight_passed": True,
            "selected_pressure_local_coordinate_independently_admitted": True,
            "source_pressure_radius_one_ball_bounds_ready": False,
            "source_R1_R2_MK_ready": False,
            "global_matched_pressure_ready": False,
            "global_leading_velocity_ready": False,
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
            "leading_pressure": {
                "producer": "agent1",
                "validator": "agent4",
                "selected_center_local_coordinate_ready": True,
                "source_ball_ready": False,
                "global_matched_pressure_ready": False,
                "next_required": (
                    "bind source-valid numerical rho, ||g||_rho and Phi radius-one-ball norm/Lipschitz; admit the "
                    "conditional pressure-ball algebra only after its independent audit resolves; then close mixed "
                    "Y Phi_Y / Y u_Y bounds, R1/R2 -> M,K -> B0/T_sh -> PA.16 -> global matched leading/pressure"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "validator": "agent4",
                "ready": True,
                "candidate_retune_allowed_by_this_checkpoint": False,
            },
            "same_cycle_defect": {
                "producer": "agent3",
                "ready": True,
                "real_full_candidate_inputs_available": False,
            },
            "correction": {
                "producer": "agent3",
                "ready": False,
                "radial_operator_reuse_without_retuning": True,
                "reason": "full same-cycle leading+oscillatory+pressure+restricted-forcing defect is still unavailable",
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
            "local_pressure_audit_is_comparable_to_st006": False,
            "reason": "local pressure consistency errors are interface-audit metrics, not NS momentum residuals",
        },
        "formal_gates_unchanged": dict(FORMAL_GATES),
        "truth_boundary": {
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "free_residual_defined_forcing_allowed": False,
            "surrogate_defect_may_promote_correction": False,
            "agent4_local_pressure_pass_may_promote_global_leading": False,
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
    a4 = upstream.get("agent4_independent_local_pressure_audit", {})
    if a4.get("pr") != AGENT4_AUDIT_PR or a4.get("head") != AGENT4_AUDIT_HEAD:
        raise ValueError("wrong Agent-4 audit identity")
    if a4.get("dedicated_workflow_run") != AGENT4_DEDICATED_RUN:
        raise ValueError("wrong Agent-4 workflow identity")
    if a4.get("dedicated_workflow_conclusion") != "success":
        raise ValueError("Agent-4 dedicated audit was not green")
    if a4.get("artifact_id") != AGENT4_ARTIFACT_ID:
        raise ValueError("wrong Agent-4 artifact identity")
    if a4.get("artifact_zip_digest") != AGENT4_ARTIFACT_DIGEST:
        raise ValueError("wrong Agent-4 artifact digest")
    if a4.get("fresh_replay_failed_guards") != [] or a4.get("fresh_replay_passed") is not True:
        raise ValueError("fresh Agent-4 local-pressure replay did not pass")

    states = checkpoint.get("states", {})
    for key in (
        "selected_pressure_local_coordinate_interface_exposed",
        "selected_pressure_local_coordinate_independent_preflight_passed",
        "selected_pressure_local_coordinate_independently_admitted",
        "oscillatory_ready",
        "public_oscillatory_time_derivative_independently_validated",
        "correction_ingest_allowed",
        "compact_radial_stress_operator_independently_admitted",
        "same_cycle_momentum_defect_contract_ready",
        "actual_momentum_defect_formula_executable",
        "same_cycle_identity_enforced",
        "restricted_forcing_provider_required",
        "heldin_heldout_disjointness_enforced",
    ):
        if states.get(key) is not True:
            raise ValueError(f"required admitted state is not true: {key}")
    for key in (
        "leading_ready",
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
    ):
        if states.get(key) is not False:
            raise ValueError(f"premature downstream promotion: {key}")

    newer = upstream.get("newer_sibling_work", {})
    if newer.get("agent4_pressure_ball_audit", {}).get("resolved_at_checkpoint_freeze") is not False:
        raise ValueError("unresolved pressure-ball audit was silently promoted")
    for name in ("agent1_pressure_ball", "agent3_finite_cycle", "agent4_pressure_ball_audit"):
        if newer.get(name, {}).get("admitted_by_this_checkpoint") is not False:
            raise ValueError(f"newer sibling work silently admitted: {name}")

    baseline = checkpoint.get("baseline_vs_kokuno", {})
    if baseline.get("kokuno_full_same_protocol_residual_available") is not False:
        raise ValueError("full Kokuno residual prematurely claimed")
    if baseline.get("local_pressure_audit_is_comparable_to_st006") is not False:
        raise ValueError("pressure interface metric mislabelled as NS residual")

    truth = checkpoint.get("truth_boundary", {})
    for key in (
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "free_residual_defined_forcing_allowed",
        "surrogate_defect_may_promote_correction",
        "agent4_local_pressure_pass_may_promote_global_leading",
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
