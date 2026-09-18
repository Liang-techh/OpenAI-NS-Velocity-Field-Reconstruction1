"""Agent-5 routing checkpoint after the displayed source phase/frame audit.

This checkpoint integrates the independently audited source-displayed phase/frame
contract from Agents 2/4 and records the latest Agent-1 leading and Agent-3 mean
correction receipts without pretending those sibling branches are executable
ancestry.  It deliberately keeps unreleased numerical source data, correction
materialization, the global candidate artifact, and the formal Navier--Stokes
validation gate fail-closed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-source-phase-mean-routing-checkpoint-v32"
TASK_ID = "KOKUNO-A5-SOURCE-PHASE-MEAN-ROUTING-032"

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
    "head_at_routing_audit": "09682c9a4c70148113639220387344ea6c9b9538",
    "latest_commit": "CR008: replay frozen ST051-B 3-D streamline/vorticity render (#490)",
    "active_kokuno_candidate_promoted_by_latest_commit": False,
}

AGENT1_RECEIPT = {
    "pr": 491,
    "head": "f44182fe310cc87a010f8ae96a7c66d173a3229d",
    "dedicated_run": 35397679969,
    "standard_run": 35397685224,
    "dedicated_status": "success",
    "standard_status": "success",
    "selected_source_rescaled_prefix_moments_executable": True,
    "shared_C_PA15_scaling_applied": True,
    "log_space_ideal_prefix_executable": True,
    "PA16_input_tuple_formable": True,
    "selected_real_axis_C_is_autonomous_not_source_hidden_C": True,
    "source_fixed_point_solved": False,
    "source_contraction_threshold_verified": False,
    "source_complex_C_bound_verified": False,
    "source_all_PA11_C_bounds_verified": False,
    "actual_source_incoming_five_moment_discrepancy_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "consumed_in_executable_ancestry": False,
}

AGENT2_RECEIPT = {
    "pr": 492,
    "head": "b6c9a5c0c4722fb9531142c64be347f522e1e307",
    "dedicated_run": 35397926675,
    "standard_run": 35397926635,
    "dedicated_status": "success",
    "standard_status": "success",
    "source_covariance_axis": "auxiliary_rectangle_sign_sigma_plus_minus",
    "displayed_source_phase_frame_bound": True,
    "phase_and_n_Phi_derived_not_free_on_adapter_path": True,
    "deterministic_nearest_nonzero_carrier_convention_recorded": True,
    "supplied_signed_complete_curl_physical_family_executable": True,
    "actual_positive_order_background_bound": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_signed_auxiliary_rectangles_or_modes_bound": False,
    "actual_source_partition_labels_instantiated": False,
    "public_source_bound_xyz_t_oscillatory_velocity_ready": False,
    "consumed_in_executable_ancestry": True,
}

AGENT3_RECEIPT = {
    "pr": 493,
    "head": "d9c85c86f169c47f23fc9ea5df1a26e048c0059a",
    "dedicated_run": 35398408955,
    "standard_run": 35398408897,
    "dedicated_status": "success",
    "standard_status": "success",
    "requested_cross_defect_identity_admitted": True,
    "finite_head_band_admitted": True,
    "requested_cross_defect_theorem_machine_replayed": False,
    "actual_mean_cross_values_materialized": False,
    "missing_weight_values_materialized": False,
    "finite_head_mean_debt_materialized": False,
    "signed_mean_inverse_input_ready": False,
    "physical_to_reference_rule": "H_ref * y = DeltaC / epsilon",
    "real_candidate_defect_consumed": False,
    "public_velocity_correction_materialized": False,
    "finite_correction_cycle_rerun_allowed": False,
    "finite_correction_cycle_run": False,
    "consumed_in_executable_ancestry": False,
}

AGENT4_RECEIPT = {
    "pr": 494,
    "head": "31f7510545488818b3bcf979c5cac3cb5d826a49",
    "parent_agent2_head": "b6c9a5c0c4722fb9531142c64be347f522e1e307",
    "dedicated_run": 35398848549,
    "standard_run": 35398848498,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact_id": 10569677115,
    "artifact_digest": "sha256:ff66fe2d26789dfa8fd5406659ad8d2ce0f9c09cab02cda15f63d742db5583af",
    "audit_seed": 9173201,
    "case_count": 36,
    "sign_case_count": 72,
    "axis_near_sign_case_count": 18,
    "off_grid_sign_case_count": 54,
    "carrier_k_values_observed": [2, 3, 4],
    "value_relative_max": 6.661338147750939e-16,
    "frame_relative_rms": 2.7633777086542752e-17,
    "source_identity_relative_max": 6.685798127696402e-16,
    "angular_winding_relative_max": 3.3306690738754696e-16,
    "fd_steps": [0.01, 0.005, 0.0025],
    "fd_relative_rms": [1.7629841009440995e-08, 2.7665385909720154e-10, 4.3321017055747775e-12],
    "fd_refinement_ratios": [63.72526689839812, 63.86134904016421],
    "wrong_axial_normalization_relative_rms": 2.0965112461866924,
    "wrong_tilt_mutation_relative_rms": 1.5780867798741534,
    "frozen_value_relative_max_guard": 5.0e-13,
    "frozen_frame_relative_rms_guard": 5.0e-13,
    "frozen_source_identity_relative_max_guard": 5.0e-13,
    "frozen_angular_winding_relative_max_guard": 5.0e-13,
    "frozen_fd_finest_relative_rms_guard": 1.0e-9,
    "frozen_fd_min_refinement_ratio_guard": 20.0,
    "frozen_wrong_axial_normalization_relative_rms_min": 0.20,
    "frozen_wrong_tilt_mutation_relative_rms_min": 0.30,
    "oracle_public_values_only": True,
    "oracle_reuses_agent2_fd_helper": False,
    "oracle_reuses_complete_curl_helper": False,
    "pressure_or_forcing_fit_used": False,
    "training_tensor_or_loss_used": False,
    "local_structural_preflight_passed": True,
    "displayed_source_phase_frame_independently_audited": True,
    "actual_positive_order_background_bound": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_signed_auxiliary_rectangles_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "source_actual_partition_labels_instantiated": False,
    "public_source_bound_velocity_osc_ready": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": True,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 485,
    "head": "689872a76ee2dbb03164101fde9c91cf4172a152",
    "schema": "kokuno-agent5-signed-physical-covariance-routing-checkpoint-v31",
    "supplied_physical_covariance_independently_audited": True,
    "actual_source_physical_covariance_rank_two_assessed": False,
    "consumed_in_executable_ancestry": False,
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
        "rectangle sign sigma=+/- distinct from Fourier m=+/-1",
        "by-sign/by-beta/total Q-scaled physical velocity",
        "displayed source phase/frame bound to public evaluation",
        "explicit epsilon scaling and coefficient units",
    ],
    "required_exports_after_instantiation": ["Python smoke", "MATLAB smoke"],
    "forcing_policy": "fixed/restricted only; residual-defined free forcing forbidden",
}


def _copy(value: Any) -> Any:
    return json.loads(json.dumps(value))


def checkpoint_sha256(payload: dict[str, Any]) -> str:
    clean = _copy(payload)
    clean.pop("checkpoint_sha256", None)
    raw = json.dumps(clean, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(raw).hexdigest()


def build_checkpoint() -> dict[str, Any]:
    states = {
        "leading_ready": False,
        "selected_source_rescaled_prefix_moments_executable": True,
        "displayed_source_phase_frame_bound": True,
        "displayed_source_phase_frame_independently_audited": True,
        "supplied_signed_complete_curl_physical_family_executable": True,
        "supplied_signed_complete_curl_physical_covariance_independently_audited": True,
        "formal_signed_mean_defect_identity_admitted": True,
        "actual_numeric_signed_mean_defect_materialized": False,
        "actual_positive_order_background_bound": False,
        "actual_source_h_sigma_pulse_integrals_bound": False,
        "actual_signed_auxiliary_rectangles_or_modes_bound": False,
        "actual_source_partition_labels_instantiated": False,
        "actual_source_physical_covariance_rank_two_assessed": False,
        "public_source_bound_xyz_t_velocity_ready": False,
        "oscillatory_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "signed_mean_inverse_input_ready": False,
        "real_candidate_defect_consumed": False,
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
        {
            "stage": "source_profile_ingest",
            "ready": False,
            "status": "displayed signed phase/frame formulas are bound and independently audited, but numerical positive-order/background, h_sigma pulses, signed rectangles/modes and actual partition labels remain unbound",
        },
        {
            "stage": "leading_candidate",
            "ready": False,
            "status": "source-rescaled PA.15 prefix moments and PA.16 handoff tuple are executable on Agent 1 #491, but actual-source five-moment discrepancy, same-scale T_sh, join and global matched pressure remain open",
        },
        {
            "stage": "oscillatory_augmentation",
            "ready": False,
            "status": "phase/frame algebra is no longer a blocker; the remaining blocker is binding numerical source/source-motivated background, h_sigma, signed rectangles/modes and partitions into public source-bound velocity_osc(x,y,z,t)",
        },
        {
            "stage": "mean_radial_corrections",
            "ready": False,
            "status": "the formal finite-head mean-defect identity is admitted, but its actual mean-cross, missing-weight and finite-head debt values are not machine-materialized; the actual source-bound covariance family is also absent",
        },
        {
            "stage": "finite_correction_cycle",
            "ready": False,
            "status": "blocked until both the actual numeric mean debt and an independently audited source-bound signed physical covariance family exist, then the unchanged inverse/budget/spacetime/radial/quadratic/joint-gain guards must pass",
        },
        {
            "stage": "candidate_artifact",
            "ready": False,
            "status": "schema reserved only; no leading+oscillatory+corrected executable ancestry exists",
        },
        {
            "stage": "independent_validation",
            "ready": False,
            "status": "Agent 4 independently clears the displayed phase/frame seam only; no frozen global candidate exists for formal held-out NS validation",
        },
        {
            "stage": "report_and_export",
            "ready": False,
            "status": "Python/MATLAB candidate smoke remains intentionally blocked until one global deterministic candidate artifact is instantiated",
        },
    ]

    routing = {
        "new_fact": "The corrected reader's displayed carrier/phase/frame formulas are now bound into the signed complete-curl adapter and independently validated under a value-only black-box oracle; phase/frame algebra is no longer the shortest-path oscillatory blocker. Separately, Agent 3 makes explicit that the formal mean-defect identity still lacks machine-materialized numeric debt, so correction input remains fail-closed.",
        "classification": "displayed_source_phase_frame_independently_cleared_numeric_source_family_and_mean_debt_still_open",
        "shortest_next_closure": [
            "Agent 2: stop revisiting displayed phase/frame algebra; bind provenance-labelled numerical positive-order/background, h_sigma pulse integrals, signed rectangles/modes and actual/source-compatible partition labels into the audited phase-bound complete-curl path, exposing Q-scaled by-sign/by-beta/total velocity_osc(x,y,z,t) with explicit units.",
            "Agent 4: independently rerun the public-value covariance audit on that first numerically bound physical family; do not inherit the earlier supplied/manufactured rank PASS as actual-source evidence.",
            "Agent 3: in parallel machine-materialize the actual finite-head mean debt (requestedStress, missingWeight, cycle state/point) from the admitted theorem path; only after both that numeric debt and the audited physical family exist may H_ref*y=DeltaC/epsilon enter the unchanged bounded-inverse/budget/spacetime/radial +d_z sigma_1/quadratic/joint-gain guards.",
            "Agent 1: continue #491's rescaled PA.15 incoming moments through a same-scale T_sh certificate, PA.16 inner-to-outer join, I3/I4 and global matched pressure while preserving the autonomous/source-compatible truth boundary.",
            "Agent 5: instantiate deterministic save/load and Python/MATLAB smoke only after leading, oscillatory and correction lanes coexist in one executable ancestry, then freeze that artifact for Agent 4's formal held-out PDE gate and same-protocol ST006 comparison.",
        ],
        "do_not_do": [
            "do not call displayed phase/frame binding recovery of unpublished numerical background, h_sigma, rectangle/mode or partition data",
            "do not hand-fill the formal mean-defect identity and call it a real numeric candidate defect",
            "do not promote supplied/source-compatible physical covariance rank to actual-source covariance readiness",
            "do not instantiate a candidate artifact before leading, oscillatory and correction lanes coexist",
            "do not call phase/frame derivative errors or covariance diagnostics Navier-Stokes momentum residuals",
            "do not relax the 1e-3 momentum or 1e-5 divergence gates",
            "do not use residual-defined free forcing",
        ],
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "source_provenance": _copy(SOURCE_PROVENANCE),
        "live_integration": _copy(LIVE_INTEGRATION_RECEIPT),
        "upstream": {
            "agent1_sibling": _copy(AGENT1_RECEIPT),
            "agent2_ancestry": _copy(AGENT2_RECEIPT),
            "agent3_sibling": _copy(AGENT3_RECEIPT),
            "agent4_ancestry": _copy(AGENT4_RECEIPT),
            "previous_agent5_sibling": _copy(PREVIOUS_AGENT5_RECEIPT),
        },
        "states": states,
        "pipeline_frontier": pipeline,
        "candidate_artifact_contract": _copy(CANDIDATE_ARTIFACT_CONTRACT),
        "baseline_vs_kokuno": {
            "st006": _copy(ST006_BASELINE),
            "kokuno_current_comparable_full_domain_receipt": None,
            "comparison_status": "not_assessable_until_global_frozen_kokuno_candidate_exists",
        },
        "formal_gates": _copy(FORMAL_GATES),
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
        "selected_source_rescaled_prefix_moments_executable",
        "displayed_source_phase_frame_bound",
        "displayed_source_phase_frame_independently_audited",
        "supplied_signed_complete_curl_physical_family_executable",
        "supplied_signed_complete_curl_physical_covariance_independently_audited",
        "formal_signed_mean_defect_identity_admitted",
    ):
        if not states[name]:
            raise ValueError(f"cleared structural/source-display seam was lost: {name}")

    for name in (
        "leading_ready",
        "actual_numeric_signed_mean_defect_materialized",
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_signed_auxiliary_rectangles_or_modes_bound",
        "actual_source_partition_labels_instantiated",
        "actual_source_physical_covariance_rank_two_assessed",
        "public_source_bound_xyz_t_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "signed_mean_inverse_input_ready",
        "real_candidate_defect_consumed",
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

    a1 = payload["upstream"]["agent1_sibling"]
    if a1["consumed_in_executable_ancestry"]:
        raise ValueError("Agent-1 sibling receipt cannot be relabeled as executable ancestry")
    if not a1["selected_source_rescaled_prefix_moments_executable"]:
        raise ValueError("Agent-1 rescaled prefix-moment receipt was lost")
    if a1["global_leading_profile_reconstructed"] or a1["global_pressure_matched"]:
        raise ValueError("Agent-1 leading sibling was over-promoted")

    a2 = payload["upstream"]["agent2_ancestry"]
    if a2["source_covariance_axis"] != "auxiliary_rectangle_sign_sigma_plus_minus":
        raise ValueError("source signed covariance axis changed")
    if not a2["displayed_source_phase_frame_bound"]:
        raise ValueError("Agent-2 displayed source phase/frame binding was lost")
    for source_flag in (
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_signed_auxiliary_rectangles_or_modes_bound",
        "actual_source_partition_labels_instantiated",
        "public_source_bound_xyz_t_oscillatory_velocity_ready",
    ):
        if a2[source_flag]:
            raise ValueError(f"Agent-2 displayed/source-hidden truth boundary changed: {source_flag}")

    a3 = payload["upstream"]["agent3_sibling"]
    if a3["consumed_in_executable_ancestry"]:
        raise ValueError("Agent-3 sibling receipt cannot be relabeled as executable ancestry")
    if not a3["requested_cross_defect_identity_admitted"] or not a3["finite_head_band_admitted"]:
        raise ValueError("Agent-3 formal mean-defect admission was lost")
    for numeric_flag in (
        "requested_cross_defect_theorem_machine_replayed",
        "actual_mean_cross_values_materialized",
        "missing_weight_values_materialized",
        "finite_head_mean_debt_materialized",
        "signed_mean_inverse_input_ready",
        "real_candidate_defect_consumed",
        "public_velocity_correction_materialized",
        "finite_correction_cycle_rerun_allowed",
        "finite_correction_cycle_run",
    ):
        if a3[numeric_flag]:
            raise ValueError(f"Agent-3 opaque formal debt was promoted numerically: {numeric_flag}")

    a4 = payload["upstream"]["agent4_ancestry"]
    if not a4["local_structural_preflight_passed"]:
        raise ValueError("Agent-4 independent phase/frame preflight must remain passed")
    if not a4["displayed_source_phase_frame_independently_audited"]:
        raise ValueError("Agent-4 displayed source phase/frame audit was lost")
    if a4["value_relative_max"] > a4["frozen_value_relative_max_guard"]:
        raise ValueError("Agent-4 value-reference guard no longer passes")
    if a4["frame_relative_rms"] > a4["frozen_frame_relative_rms_guard"]:
        raise ValueError("Agent-4 phase-frame guard no longer passes")
    if a4["source_identity_relative_max"] > a4["frozen_source_identity_relative_max_guard"]:
        raise ValueError("Agent-4 source-identity guard no longer passes")
    if a4["angular_winding_relative_max"] > a4["frozen_angular_winding_relative_max_guard"]:
        raise ValueError("Agent-4 angular-winding guard no longer passes")
    if a4["fd_relative_rms"][-1] > a4["frozen_fd_finest_relative_rms_guard"]:
        raise ValueError("Agent-4 finest FD phase-gradient guard no longer passes")
    if min(a4["fd_refinement_ratios"]) < a4["frozen_fd_min_refinement_ratio_guard"]:
        raise ValueError("Agent-4 FD refinement guard no longer passes")
    if a4["wrong_axial_normalization_relative_rms"] < a4["frozen_wrong_axial_normalization_relative_rms_min"]:
        raise ValueError("Agent-4 wrong-axial-normalization mutation lost sensitivity")
    if a4["wrong_tilt_mutation_relative_rms"] < a4["frozen_wrong_tilt_mutation_relative_rms_min"]:
        raise ValueError("Agent-4 wrong-tilt mutation lost sensitivity")
    if a4["carrier_k_values_observed"] != [2, 3, 4]:
        raise ValueError("Agent-4 carrier-transition coverage changed")
    if a4["formal_full_domain_pde_gate_assessed"] or a4["pde_validated"]:
        raise ValueError("local Agent-4 phase audit was promoted to PDE evidence")

    if payload["formal_gates"] != FORMAL_GATES:
        raise ValueError("formal PDE gates changed")
    if payload["candidate_artifact_contract"]["instantiated"]:
        raise ValueError("candidate artifact must remain reserved only")
    if payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is not None:
        raise ValueError("no comparable full-domain Kokuno receipt exists yet")


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
