"""Agent-5 routing checkpoint for the first source-scheduled multi-band seam.

The independent Cartesian numerical guards are now clean, but the public Agent-2
multi-band contract fails a semantically identical simultaneous beta permutation at
12/36 held-out points because its total and by-band aggregation use different
floating summation orders after large Q**(-A) scaling.  This checkpoint preserves
that structural rejection and routes the smallest upstream API fix without changing
any local or formal validation threshold.
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
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_outer_pressure_datum_bound": False,
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

AGENT4_RECEIPT = {
    "pr": 451,
    "head": "3e8170bde17cfeee5189e453f304388e65564931",
    "dedicated_run": 35377042248,
    "standard_run": 35377042241,
    "dedicated_status": "success",
    "standard_status_at_routing_audit": "in_progress",
    "artifact_id": 10560697796,
    "artifact_digest": "sha256:9accb9aa62709d6cd865574e6e01f7d199d7150c9a1ebd447a8b936c4b2b7c17",
    "scientific_report_generated": True,
    "source_schedule_max_relative_error": 1.4432899320127035e-14,
    "cartesian_curl_relative_rms": [8.601710891e-4, 5.910638257e-5, 3.789822318e-6],
    "curl_refinement_ratios": [14.55293069, 15.59608277],
    "normalized_divergence_rms": [1.174437055e-3, 8.397047230e-5, 5.450640031e-6],
    "divergence_refinement_ratios": [13.98630999, 15.40561692],
    "finest_normalized_divergence_point_max": 2.157397612e-5,
    "shared_first_band_schedule_mutation_relative_rms": 0.5546591826,
    "reversed_label_schedule_mutation_relative_rms": 1.1243576802,
    "numeric_guards_passed": True,
    "simultaneous_beta_permutation_rejections": 12,
    "simultaneous_beta_permutation_points": 36,
    "permutation_max_relative_change_where_evaluable": 1.6259252002601664e-16,
    "local_structural_preflight_passed": False,
    "failure_classification": "public_contract_summation_order_invariance_failure",
    "formal_full_domain_pde_gate_assessed": False,
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
    "changed_after_result": False,
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
        "source_multiband_independent_cartesian_numeric_guards_passed": True,
        "source_multiband_public_contract_structural_preflight_passed": False,
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
        {"stage": "source_profile_ingest", "ready": False, "status": "source constants/schedules and autonomous support geometry are executable; hidden positive-order/background data remain unbound"},
        {"stage": "leading_candidate", "ready": False, "status": "PA.10 T_sh formula executable; analytic B0/source pressure datum/complete join/global pressure remain open"},
        {"stage": "oscillatory_augmentation", "ready": False, "status": "independent curl/divergence/mutation numerics pass, but simultaneous beta permutation exposes a public aggregation-order invariant failure at 12/36 held-out points"},
        {"stage": "mean_radial_corrections", "ready": False, "status": "joint damping margin is green only on prior one-band negative control; actual audited multiband rank-two receipt absent"},
        {"stage": "finite_correction_cycle", "ready": False, "status": "blocked until actual multiband family passes structural API, genuine rank, unit/budget/spacetime/radial/quadratic guards"},
        {"stage": "candidate_artifact", "ready": False, "status": "schema reserved, not instantiated"},
        {"stage": "independent_validation", "ready": False, "status": "Agent-4 local multiband audit is a structural REJECT; no global candidate exists for formal PDE gate"},
        {"stage": "report_and_export", "ready": False, "status": "blocked until frozen candidate artifact exists"},
    ]
    routing = {
        "new_fact": "Agent 4 #451 independently clears every frozen source-schedule/curl/divergence/mutation numerical guard on Agent 2 #450, but a semantically identical simultaneous beta-data + beta-label permutation is rejected at 12/36 held-out points. Where evaluable, the field changes only 1.63e-16, so the mathematical field is permutation invariant while the public aggregation check is not.",
        "classification": "local_structural_reject_after_numeric_guards_pass",
        "shortest_next_closure": [
            "Agent 2: make cart_total and by-band aggregation use deterministic identical summation/order (or an algebraically identical ordered check) so simultaneous beta permutation is invariant; do not inflate the absolute 2e-12 tolerance after observing the failure.",
            "Agent 4: rerun the exact same #451 independent audit and frozen guards against that API-only fix; only a structural PASS clears this seam.",
            "Agent 2: after structural PASS, bind provenance-labelled positive-order/background and actual auxiliary-torus mode data and expose public Q-scaled by-beta/total velocity(x,y,z,t) with explicit coefficient units.",
            "Agent 3: consume that exact audited physical family and require genuine rank-two plus unchanged unit/budget/spacetime/radial/quadratic guards before correction materialization or finite-cycle rerun.",
            "Agent 1: close source pressure datum, validated T_sh/inner-to-outer join and global matched pressure under fixed/restricted forcing.",
            "Agent 5: instantiate deterministic candidate artifact/save-load/Python-MATLAB smoke only after leading, oscillatory and correction lanes coexist; then hand frozen artifact to Agent 4 for formal held-out PDE gate and same-protocol ST006 comparison.",
        ],
        "do_not_do": [
            "do not relax Agent-2 aggregation tolerance merely to make the permutation test green",
            "do not promote passing local curl/divergence guards to oscillatory_ready or PDE validation",
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
            "agent4": AGENT4_RECEIPT,
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
    s = payload["states"]
    hard_false = (
        "leading_ready",
        "source_multiband_public_contract_structural_preflight_passed",
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
    if any(s[name] for name in hard_false):
        raise ValueError("fail-closed state was promoted without evidence")
    if not s["source_scheduled_multiband_supplied_mode_interface_ready"] or not s["source_multiband_independent_cartesian_numeric_guards_passed"]:
        raise ValueError("stable multiband numerical receipt was lost")

    a2 = payload["upstream"]["agent2"]
    if a2["dedicated_status"] != "success" or a2["standard_status"] != "success":
        raise ValueError("Agent-2 exact-head CI is not stable")
    if a2["actual_positive_order_background_bound"] or a2["actual_auxiliary_torus_mode_family_bound"]:
        raise ValueError("manufactured supplied-mode interface was promoted to source-bound data")

    a4 = payload["upstream"]["agent4"]
    if a4["dedicated_status"] != "success" or not a4["scientific_report_generated"] or not a4["numeric_guards_passed"]:
        raise ValueError("Agent-4 numerical audit receipt is not stable")
    if a4["local_structural_preflight_passed"] or a4["simultaneous_beta_permutation_rejections"] != 12:
        raise ValueError("Agent-4 structural rejection was lost")
    if a4["formal_full_domain_pde_gate_assessed"]:
        raise ValueError("local audit was promoted to formal PDE assessment")

    expected_local = {
        "source_schedule_relative_error_max": 5.0e-14,
        "finest_cartesian_curl_relative_rms_max": 2.0e-5,
        "curl_refinement_ratio_min": 4.0,
        "finest_normalized_divergence_rms_max": 2.0e-5,
        "finest_normalized_divergence_point_max": 8.0e-5,
        "divergence_refinement_ratio_min": 3.0,
        "shared_first_band_schedule_mutation_relative_rms_min": 0.10,
        "label_schedule_swap_mutation_relative_rms_min": 0.10,
        "changed_after_result": False,
    }
    if payload["agent4_frozen_local_guards"] != expected_local:
        raise ValueError("Agent-4 frozen local guards changed after result")

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
