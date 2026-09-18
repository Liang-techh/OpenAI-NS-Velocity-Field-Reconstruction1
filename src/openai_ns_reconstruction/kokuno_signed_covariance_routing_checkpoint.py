"""Agent-5 routing checkpoint after the source signed-covariance seam.

Agent 2 #472 identifies the corrected-reader covariance axis as auxiliary-rectangle
sign sigma=+/- and makes the displayed two-sign reference inverse executable.
Agent 4 #473 independently audits that algebraic/differential seam under frozen
guards. Agent 3 #474 then makes the physical/reference epsilon unit conversion
explicit. This module records those receipts and the shortest remaining closure.

It does not infer physical complete-curl covariance rank from the reference
matrix, materialize a correction, instantiate a final candidate, or claim a
Navier--Stokes residual reduction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-signed-covariance-routing-checkpoint-v30"
TASK_ID = "KOKUNO-A5-SIGNED-COVARIANCE-ROUTING-030"

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
    "head_at_routing_audit": "60626c05b9340ba8e966702bcaefb9c52d3e5d20",
    "latest_commit": "CR002: govern cross-backbone redistribution identity (#467)",
    "active_kokuno_candidate_promoted_by_latest_commit": False,
}

AGENT1_RECEIPT = {
    "pr": 464,
    "head": "3de16f42f71994b192ba90fc4231314fb84c3f0f",
    "standard_run": 35387382491,
    "standard_status": "success",
    "source_rescaled_reference_continuation_executable": True,
    "selected_source_scale_pressure_carried_through_reference_stage": True,
    "source_fixed_point_solved": False,
    "source_contraction_threshold_verified": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "consumed_in_executable_ancestry": False,
}

AGENT2_RECEIPT = {
    "pr": 472,
    "head": "2bfb4c6f3a9314eea19b7dc0aba4fc11f2877978",
    "dedicated_run": 35386864399,
    "standard_run": 35386864444,
    "dedicated_status": "success",
    "standard_status": "success",
    "source_covariance_axis": "auxiliary_rectangle_sign_sigma_plus_minus",
    "source_covariance_axis_is_dyadic_band_contrast": False,
    "source_covariance_axis_is_fourier_harmonic_sign": False,
    "displayed_signed_reference_inverse_executable": True,
    "displayed_reference_determinant_nonzero_on_open_shell": True,
    "signed_amplitude_directional_derivative_executable": True,
    "actual_positive_order_background_bound": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_signed_auxiliary_rectangles_bound": False,
    "physical_complete_curl_signed_family_bound": False,
    "public_source_bound_xyz_t_velocity_ready": False,
    "consumed_in_executable_ancestry": True,
}

AGENT3_RECEIPT = {
    "pr": 474,
    "head": "7b189d7b59d7942aaaec5586526b2b09eabc3e00",
    "dedicated_run": 35388198317,
    "standard_run": 35388198356,
    "dedicated_status": "success",
    "standard_status": "success",
    "physical_covariance_relation": "DeltaC = epsilon * H_ref * y",
    "reference_inverse_relation": "H_ref * y = DeltaC / epsilon",
    "epsilon_division_required": True,
    "structural_calibration_only": True,
    "real_candidate_defect_consumed": False,
    "structural_squared_amplitude_examples": [[0.4, 0.6], [1.6, 2.4], [6.4, 9.6]],
    "structural_cone_margin_fraction": 0.8,
    "structural_reference_condition_number": 1.5,
    "public_velocity_correction_materialized": False,
    "finite_correction_cycle_rerun_allowed": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": True,
}

AGENT4_RECEIPT = {
    "pr": 473,
    "head": "d6e0632f499528f00b8d8a9ffed73c27fd2c1485",
    "dedicated_run": 35388072406,
    "standard_run": 35388072286,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact_id": 10564942297,
    "artifact_digest": "sha256:7de16b6a7c4cc3ceb09aee077343a9ef8a1235305571fe35f1d45db4d1746bda",
    "heldout_positive_cone_cases": 48,
    "reference_value_max_relative_error": 4.68851475280062e-16,
    "directional_derivative_relative_rms": [1.2751736273487234e-8, 7.96824446149774e-10, 4.9789788684764504e-11],
    "directional_derivative_refinement_ratios": [16.003194097649924, 16.00377240391059],
    "semantic_sigma_swap_relative_max": 0.0,
    "wrong_q_sign_mutation_relative_rms_min": 0.2821526546786201,
    "dropped_dh_sigma_mutation_relative_rms_min": 0.6187952896876538,
    "outside_cone_rejected": 48,
    "strict_direction_gap_rejected": 48,
    "local_structural_preflight_passed": True,
    "physical_complete_curl_signed_family_audited": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 465,
    "head": "a37894b67a4a7c2b2891f4d95a33b48dea3f15d1",
    "schema": "kokuno-agent5-covariance-rank-routing-checkpoint-v29",
}

ST006_BASELINE = {
    "candidate": "ST006",
    "held_out_seed": 9172801,
    "points": 4096,
    "times": 6,
    "finest_spatial_step": 0.005,
    "momentum_sampled_max": 0.1082289305112118,
    "momentum_volume_l2": 0.10758432876230622,
    "pde_validated": False,
}

FORMAL_GATES = {
    "held_out_normalized_momentum_max": 1.0e-3,
    "held_out_normalized_momentum_l2": 1.0e-3,
    "held_out_divergence_max": 1.0e-5,
    "held_out_divergence_l2": 1.0e-5,
    "changed_this_round": False,
}

CANDIDATE_ARTIFACT_CONTRACT = {
    "schema_reserved": "kokuno-derived-3d-candidate-v1",
    "instantiated": False,
    "required_provenance": [
        "source repository/version/document",
        "source-displayed versus repository-autonomous choices",
        "Agent 1-4 exact stage receipts",
        "truth-boundary flags",
    ],
    "required_public_api": [
        "velocity(x,y,z,t)->[u,v,w]",
        "pressure(x,y,z,t)",
        "restricted_forcing(x,y,z,t)",
        "deterministic save/load",
    ],
    "required_oscillatory_subcontract": [
        "velocity_osc(x,y,z,t)",
        "sign-resolved sigma=+/- contributions distinct from Fourier m=+/-1",
        "by-beta and total Q-scaled physical velocity",
        "explicit epsilon scaling and coefficient units",
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
        "displayed_source_signed_covariance_reference_executable": True,
        "signed_covariance_reference_independently_audited": True,
        "signed_covariance_physical_to_reference_unit_bridge_ready": True,
        "actual_positive_order_background_bound": False,
        "actual_source_h_sigma_pulse_integrals_bound": False,
        "actual_signed_auxiliary_rectangles_bound": False,
        "physical_complete_curl_signed_family_bound": False,
        "public_source_bound_xyz_t_velocity_ready": False,
        "oscillatory_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "correction_ready": False,
        "public_velocity_correction_materialized": False,
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
        {"stage": "source_profile_ingest", "ready": False, "status": "displayed signed covariance reference is executable and independently audited, but actual h_sigma pulse integrals, signed rectangles/modes and positive-order background remain unbound"},
        {"stage": "leading_candidate", "ready": False, "status": "source-rescaled reference continuation is executable; fixed point/contraction/global join/matched pressure remain open"},
        {"stage": "oscillatory_augmentation", "ready": False, "status": "correct source covariance axis sigma=+/- is identified; reference inverse and epsilon unit bridge are ready, but no source-bound sign-resolved complete-curl xyz,t family exists"},
        {"stage": "mean_radial_corrections", "ready": False, "status": "Agent-3 epsilon bridge is structural only; physical signed family has not yet been independently phase-mean audited or rank-screened"},
        {"stage": "finite_correction_cycle", "ready": False, "status": "blocked until physical signed complete-curl family passes genuine covariance rank plus unchanged budget/spacetime/radial/quadratic/joint-gain guards"},
        {"stage": "candidate_artifact", "ready": False, "status": "schema and signed oscillatory subcontract reserved, not instantiated"},
        {"stage": "independent_validation", "ready": False, "status": "Agent 4 cleared only the displayed signed reference inverse/derivative contract; no complete global artifact exists for the formal PDE gate"},
        {"stage": "report_and_export", "ready": False, "status": "blocked until one frozen leading+oscillatory+corrected artifact exists"},
    ]
    routing = {
        "new_fact": "The source covariance direction is auxiliary-rectangle sign sigma=+/- rather than dyadic-band contrast or Fourier m sign. Agent 4 independently clears the displayed signed reference inverse and slow-coordinate amplitude derivative, while Agent 3 freezes the required DeltaC/epsilon unit conversion.",
        "classification": "signed_reference_inverse_cleared_waiting_for_physical_complete_curl_family",
        "historical_negative_control": "Agent 3 #462 remains a valid rejection of the old manufactured dyadic-band contrast as a phase-mean covariance direction; it is not the source signed axis and is not retroactively rewritten as a pass.",
        "shortest_next_closure": [
            "Agent 2: carry sigma=+/- through provenance-labelled actual/source-motivated h_sigma pulses, signed auxiliary rectangles/modes, slow derivatives and support gradients into localized complete curls; expose public Q-scaled sign/by-beta/total velocity_osc(x,y,z,t) with explicit epsilon and coefficient units.",
            "Agent 4: independently audit the first physical signed complete-curl family, including phase-mean covariance columns and complete-curl remainder effects rather than inferring rank from det(H_ref).",
            "Agent 3: consume exactly that audited physical family and the real candidate defect through H_ref*y=DeltaC/epsilon; require genuine covariance rank two plus unchanged unit/budget/spacetime/radial +d_z sigma_1/quadratic/joint-gain guards before materializing delta_u or rerunning the finite correction cycle.",
            "Agent 1: continue the source-rescaled reference stage through the actual fixed-point/contraction/global join path and complete matched pressure under fixed/restricted forcing.",
            "Agent 5: instantiate deterministic save/load plus Python/MATLAB smoke only after leading, oscillatory and correction lanes coexist in one executable ancestry; then freeze that artifact for Agent 4's formal held-out PDE gate and same-protocol ST006 comparison.",
        ],
        "do_not_do": [
            "do not infer physical covariance rank readiness from the nonzero displayed reference determinant alone",
            "do not relabel sigma=+/- as dyadic-band contrast or Fourier m=+/-1",
            "do not omit the required division by epsilon when mapping a physical covariance target into reference units",
            "do not treat Agent-3 structural calibration as a real candidate defect or correction",
            "do not call local algebraic/differential errors Navier-Stokes residuals",
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
        "candidate_artifact_contract": CANDIDATE_ARTIFACT_CONTRACT,
        "baseline_vs_kokuno": {
            "st006": ST006_BASELINE,
            "kokuno_current_comparable_full_domain_receipt": None,
            "comparison_status": "not_assessable_until_global_frozen_kokuno_candidate_exists",
        },
        "formal_gates": FORMAL_GATES,
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
    for name in (
        "displayed_source_signed_covariance_reference_executable",
        "signed_covariance_reference_independently_audited",
        "signed_covariance_physical_to_reference_unit_bridge_ready",
    ):
        if not states[name]:
            raise ValueError(f"cleared signed-covariance seam was lost: {name}")

    for name in (
        "leading_ready",
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_signed_auxiliary_rectangles_bound",
        "physical_complete_curl_signed_family_bound",
        "public_source_bound_xyz_t_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "correction_ready",
        "public_velocity_correction_materialized",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if states[name]:
            raise ValueError(f"fail-closed state promoted without evidence: {name}")

    a2 = payload["upstream"]["agent2"]
    if a2["source_covariance_axis"] != "auxiliary_rectangle_sign_sigma_plus_minus":
        raise ValueError("source signed covariance axis changed")
    if a2["source_covariance_axis_is_dyadic_band_contrast"] or a2["source_covariance_axis_is_fourier_harmonic_sign"]:
        raise ValueError("source signed axis was conflated with another label")

    a3 = payload["upstream"]["agent3"]
    if not a3["epsilon_division_required"] or a3["real_candidate_defect_consumed"]:
        raise ValueError("Agent-3 unit/truth boundary changed")

    a4 = payload["upstream"]["agent4"]
    if not a4["local_structural_preflight_passed"]:
        raise ValueError("Agent-4 signed reference audit must remain passed")
    if a4["physical_complete_curl_signed_family_audited"]:
        raise ValueError("reference audit was promoted to physical-family audit")
    if a4["outside_cone_rejected"] != a4["heldout_positive_cone_cases"]:
        raise ValueError("outside-cone negative-control receipt changed")
    if a4["strict_direction_gap_rejected"] != a4["heldout_positive_cone_cases"]:
        raise ValueError("direction-gap negative-control receipt changed")

    if payload["formal_gates"] != FORMAL_GATES:
        raise ValueError("formal PDE gates changed")
    if payload["candidate_artifact_contract"]["instantiated"]:
        raise ValueError("candidate artifact must remain reserved only")


def write_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = write_checkpoint(args.output)
    print(payload["checkpoint_sha256"])


if __name__ == "__main__":
    main()
