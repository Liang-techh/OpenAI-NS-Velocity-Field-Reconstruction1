"""Agent-5 routing checkpoint for the first source-scheduled multi-band seam.

This module records immutable upstream receipts and an integration failure without
promoting it to a scientific rejection.  Agent 4's first independent audit never
reached its frozen curl/divergence guards because its label/schedule mutation was
rejected by Agent 2's fail-closed band-aggregation invariant.  The numerical
mutation threshold is therefore retained unchanged for the next in-contract audit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-multiband-audit-routing-checkpoint-v28"
TASK_ID = "KOKUNO-A5-MULTIBAND-AUDIT-ROUTING-028"

SOURCE_PROVENANCE = {
    "repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "document": "navier-stokes/navier_stokes_workbench.tex",
    "corrected_reader_date": "2026-09-09",
    "zenodo_record": "22678406",
    "paper_exact": False,
}

LIVE_INTEGRATION_RECEIPT = {
    "branch": "codex/cr001-constraints",
    "head_at_routing_audit": "26b9aa874b2812837cdd598e56dfcd08b836ce2a",
    "latest_commit": "CR002: govern material-path evidence semantics (#444)",
    "active_kokuno_route_promoted_by_latest_commit": False,
}

AGENT1_RECEIPT = {
    "pr": 443,
    "head": "490c9c6d528f5bc1f6b53190eff198c958991cf6",
    "standard_run": 35374048705,
    "standard_status": "success",
    "PA10_T_sh_formula_executable": True,
    "selected_B0_C0_envelope_numerically_checked": True,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_outer_pressure_datum_bound": False,
    "actual_source_incoming_five_moment_discrepancy_bound": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "consumed_in_executable_ancestry": False,
}

AGENT2_RECEIPT = {
    "pr": 450,
    "head": "054d5b380a4ca2fcf5c2c8cf5e5a6532e82a5703",
    "dedicated_run": 35375913786,
    "standard_run": 35375913668,
    "dedicated_status": "success",
    "standard_status": "success",
    "source_scheduled_multiband_supplied_mode_family_executable": True,
    "requires_two_distinct_ell_bands": True,
    "per_band_Q_epsilon_schedule_bound": True,
    "by_beta_cartesian_columns_exposed": True,
    "by_band_cartesian_columns_exposed": True,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_source_bound_xyz_t_velocity_ready": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "consumed_in_executable_ancestry": True,
}

AGENT3_RECEIPT = {
    "pr": 449,
    "head": "383eeb622e40bc4ac46f02a369bd205ca003adc8",
    "dedicated_run": 35375639370,
    "standard_run": 35375639385,
    "agent5_reference_run": 35375639420,
    "dedicated_status": "success",
    "standard_status": "success",
    "agent5_reference_status": "success",
    "shared_selected_damping": 0.4926349922838171,
    "joint_admissible_interval": [0.02862512785452287, 0.9858531565488349],
    "joint_gain_margin_preflight_passed": True,
    "actual_multiband_family_consumed": False,
    "genuine_rank_two_passed": False,
    "finite_correction_cycle_rerun_allowed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

AGENT4_FAILED_AUDIT_RECEIPT = {
    "pr": 451,
    "head": "c314771d90a7d67c447526320ec91be2e4908bd9",
    "dedicated_run": 35376636412,
    "standard_run": 35376636366,
    "dedicated_status": "failure",
    "standard_status": "failure",
    "focused_tests_passed_before_error": 7,
    "focused_tests_errors": 3,
    "scientific_report_generated": False,
    "formal_full_domain_pde_gate_assessed": False,
    "failure_phase": "mutation_control_report_construction",
    "failure_exception": "RuntimeError: band aggregation does not reproduce the total Cartesian correction",
    "failure_interpretation": "the swapped-label mutation violated Agent-2 label/data binding and was rejected by the public fail-closed aggregation invariant before frozen curl/divergence guards ran",
    "scientific_preflight_status": "not_assessed_due_to_contract_mismatch",
    "pde_validated": False,
}

FROZEN_AGENT4_LOCAL_GUARDS = {
    "source_schedule_relative_error_max": 5.0e-14,
    "finest_cartesian_curl_relative_rms_max": 2.0e-5,
    "curl_refinement_ratio_min": 4.0,
    "finest_normalized_divergence_rms_max": 2.0e-5,
    "finest_normalized_divergence_point_max": 8.0e-5,
    "divergence_refinement_ratio_min": 3.0,
    "shared_first_band_schedule_mutation_relative_rms_min": 0.10,
    "label_schedule_swap_mutation_relative_rms_min": 0.10,
    "changed_after_failure": False,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 442,
    "head": "b67884c691647badd35f4b4b00947745a39bc332",
    "schema": "kokuno-agent5-rational-rectangle-routing-checkpoint-v27",
}

ST006_BASELINE = {
    "candidate": "ST006",
    "held_out_seed": 9172801,
    "points": 4096,
    "times": 6,
    "finest_spatial_step": 0.005,
    "momentum_sampled_max": 0.1082289305112118,
    "momentum_volume_l2": 0.10758432876230622,
    "formal_momentum_gate": 1.0e-3,
    "formal_divergence_gate": 1.0e-5,
    "pde_validated": False,
}

CANDIDATE_ARTIFACT_CONTRACT = {
    "schema_reserved": "kokuno-derived-3d-candidate-v1",
    "instantiated": False,
    "required_provenance": [
        "source repository/version/commit",
        "repository-autonomous choices",
        "Agent-1 leading receipt",
        "Agent-2 oscillatory receipt",
        "Agent-3 correction receipt",
        "Agent-4 independent validation receipt",
    ],
    "required_public_api": [
        "velocity(x,y,z,t)->[u,v,w]",
        "pressure(x,y,z,t)",
        "restricted_forcing(x,y,z,t)",
        "deterministic save/load",
    ],
    "required_exports_after_instantiation": ["Python smoke", "MATLAB smoke"],
    "forcing_policy": "fixed/restricted only; residual-defined free forcing forbidden",
}


def checkpoint_sha256(payload: dict[str, Any]) -> str:
    clean = dict(payload)
    clean.pop("checkpoint_sha256", None)
    raw = json.dumps(clean, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(raw).hexdigest()


def build_checkpoint() -> dict[str, Any]:
    states = {
        "leading_ready": False,
        "source_scheduled_multiband_supplied_mode_interface_ready": True,
        "source_multiband_independent_cartesian_audit_assessed": False,
        "source_multiband_independent_cartesian_audit_passed": False,
        "actual_positive_order_background_bound": False,
        "actual_auxiliary_torus_mode_family_bound": False,
        "public_source_bound_xyz_t_velocity_ready": False,
        "oscillatory_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "correction_ready": False,
        "finite_correction_cycle_run": False,
        "candidate_artifact_instantiated": False,
        "velocity_export_ready": False,
        "formal_full_domain_pde_gate_assessed": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    pipeline = [
        {"stage": "source_profile_ingest", "ready": False, "status": "displayed source constants/schedules and autonomous support geometry are executable; hidden positive-order/background data remain unbound"},
        {"stage": "leading_candidate", "ready": False, "status": "PA.10 T_sh formula is executable but analytic B0/source pressure datum/complete join/global pressure remain open"},
        {"stage": "oscillatory_augmentation", "ready": False, "status": "first source-scheduled >=2-band supplied-mode interface exists; independent audit #451 stopped at a label/data mutation contract mismatch before scientific guards"},
        {"stage": "mean_radial_corrections", "ready": False, "status": "joint damping margin is green on the prior one-band negative control; no actual multiband rank-two receipt exists"},
        {"stage": "finite_correction_cycle", "ready": False, "status": "blocked until an independently audited actual multiband family passes genuine rank/unit/budget/spacetime/radial/quadratic guards"},
        {"stage": "candidate_artifact", "ready": False, "status": "schema reserved, not instantiated"},
        {"stage": "independent_validation", "ready": False, "status": "Agent-4 #451 is an infrastructure/API mismatch, not a scientific pass or rejection; no global candidate exists for the formal PDE gate"},
        {"stage": "report_and_export", "ready": False, "status": "blocked until a frozen candidate artifact exists"},
    ]
    routing = {
        "new_fact": "Agent 2 #450 is the first current source-scheduled multi-band supplied-mode family with distinct dyadic ell bands and per-band Q/epsilon. Agent 4 #451 attempted an independent Cartesian audit, but its swapped-label mutation violated Agent-2 label/data binding and the public family correctly failed closed before the frozen numerical guards could be evaluated.",
        "classification": "integration_contract_mismatch_before_scientific_assessment",
        "shortest_next_closure": [
            "Agent 4: keep every #451 frozen numerical guard unchanged and replace only the invalid swapped-label mutation construction with an in-contract perturbation that returns a numeric altered field; rerun the same fresh-point FD4 audit. Do not reinterpret the present exception as satisfying the preregistered >=0.10 numeric mutation threshold.",
            "Agent 2: after that audit passes, bind provenance-labelled positive-order/background and actual auxiliary-torus mode data to the already source-scheduled multi-band interface and expose public Q-scaled by-beta/total velocity(x,y,z,t) with explicit coefficient units.",
            "Agent 3: consume the exact audited physical family and require genuine rank-two plus unchanged unit/budget/spacetime/radial/quadratic guards before materializing correction or rerunning the finite correction cycle.",
            "Agent 1: close source pressure datum, analytic/validated T_sh requirements, full inner-to-outer join and global matched pressure under fixed/restricted forcing.",
            "Agent 5: instantiate deterministic candidate artifact/save-load/Python-MATLAB smoke only when leading, oscillatory and correction lanes coexist in one executable ancestry; then hand the frozen artifact to Agent 4 for the formal held-out PDE gate and same-protocol ST006 comparison.",
        ],
        "do_not_do": [
            "do not weaken or rewrite Agent-4 frozen thresholds after the failed run",
            "do not count a fail-closed exception as the preregistered numeric mutation >=0.10 success",
            "do not call manufactured supplied-mode data recovered Kokuno source data",
            "do not infer genuine covariance rank from multiple dyadic band labels",
            "do not relax the 1e-3 momentum or 1e-5 divergence gates",
            "do not use residual-defined free forcing",
        ],
    }
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "source_provenance": SOURCE_PROVENANCE,
        "live_integration": LIVE_INTEGRATION_RECEIPT,
        "upstream": {
            "agent1": AGENT1_RECEIPT,
            "agent2": AGENT2_RECEIPT,
            "agent3": AGENT3_RECEIPT,
            "agent4_failed_audit": AGENT4_FAILED_AUDIT_RECEIPT,
            "previous_agent5": PREVIOUS_AGENT5_RECEIPT,
        },
        "states": states,
        "pipeline_frontier": pipeline,
        "agent4_frozen_local_guards": FROZEN_AGENT4_LOCAL_GUARDS,
        "candidate_artifact_contract": CANDIDATE_ARTIFACT_CONTRACT,
        "baseline_vs_kokuno": {
            "st006": ST006_BASELINE,
            "kokuno_current_comparable_full_domain_receipt": None,
            "comparison_status": "not_assessable_until_global_frozen_kokuno_candidate_exists",
        },
        "formal_gates": {
            "held_out_normalized_momentum_max": 1.0e-3,
            "held_out_normalized_momentum_l2": 1.0e-3,
            "held_out_divergence_max": 1.0e-5,
            "held_out_divergence_l2": 1.0e-5,
            "changed_this_round": False,
        },
        "routing": routing,
    }
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def validate_checkpoint(payload: dict[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or payload.get("task_id") != TASK_ID:
        raise ValueError("wrong Agent-5 checkpoint identity")
    if payload.get("checkpoint_sha256") != checkpoint_sha256(payload):
        raise ValueError("checkpoint sha256 mismatch")

    states = payload["states"]
    hard_false = (
        "leading_ready",
        "source_multiband_independent_cartesian_audit_assessed",
        "source_multiband_independent_cartesian_audit_passed",
        "actual_positive_order_background_bound",
        "actual_auxiliary_torus_mode_family_bound",
        "public_source_bound_xyz_t_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    if any(states[name] for name in hard_false):
        raise ValueError("fail-closed state was promoted without evidence")
    if not states["source_scheduled_multiband_supplied_mode_interface_ready"]:
        raise ValueError("stable Agent-2 multiband interface receipt was lost")

    a2 = payload["upstream"]["agent2"]
    if a2["dedicated_status"] != "success" or a2["standard_status"] != "success":
        raise ValueError("Agent-2 exact-head CI is not stable")
    if a2["actual_positive_order_background_bound"] or a2["actual_auxiliary_torus_mode_family_bound"]:
        raise ValueError("manufactured supplied-mode interface was promoted to source-bound data")

    a4 = payload["upstream"]["agent4_failed_audit"]
    if a4["dedicated_status"] != "failure" or a4["standard_status"] != "failure":
        raise ValueError("Agent-4 failed CI receipt changed")
    if a4["scientific_report_generated"] or a4["formal_full_domain_pde_gate_assessed"]:
        raise ValueError("failed pre-report audit was promoted to scientific assessment")
    if a4["scientific_preflight_status"] != "not_assessed_due_to_contract_mismatch":
        raise ValueError("Agent-4 failure classification changed")

    local = payload["agent4_frozen_local_guards"]
    expected_local = {
        "source_schedule_relative_error_max": 5.0e-14,
        "finest_cartesian_curl_relative_rms_max": 2.0e-5,
        "curl_refinement_ratio_min": 4.0,
        "finest_normalized_divergence_rms_max": 2.0e-5,
        "finest_normalized_divergence_point_max": 8.0e-5,
        "divergence_refinement_ratio_min": 3.0,
        "shared_first_band_schedule_mutation_relative_rms_min": 0.10,
        "label_schedule_swap_mutation_relative_rms_min": 0.10,
        "changed_after_failure": False,
    }
    if local != expected_local:
        raise ValueError("Agent-4 frozen local guards changed after failure")

    gates = payload["formal_gates"]
    if gates != {
        "held_out_normalized_momentum_max": 1.0e-3,
        "held_out_normalized_momentum_l2": 1.0e-3,
        "held_out_divergence_max": 1.0e-5,
        "held_out_divergence_l2": 1.0e-5,
        "changed_this_round": False,
    }:
        raise ValueError("formal PDE gates changed")

    baseline = payload["baseline_vs_kokuno"]["st006"]
    if baseline["momentum_sampled_max"] != 0.1082289305112118 or baseline["momentum_volume_l2"] != 0.10758432876230622:
        raise ValueError("ST006 baseline changed")
    if payload["candidate_artifact_contract"]["instantiated"]:
        raise ValueError("candidate artifact was instantiated before lane closure")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = build_checkpoint()
    validate_checkpoint(payload)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(args.output)
    print(payload["checkpoint_sha256"])


if __name__ == "__main__":
    main()
