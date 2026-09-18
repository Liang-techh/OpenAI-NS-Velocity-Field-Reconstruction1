"""Agent-5 routing checkpoint after independent signed physical-covariance preflight.

This checkpoint integrates the first supplied/source-compatible sign-resolved
complete-curl physical family with Agent-3's covariance-rank screen and Agent-4's
independent black-box audit.  It deliberately does not promote manufactured
inputs to recovered Kokuno data, materialize a correction, or assess the formal
Navier--Stokes gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-signed-physical-covariance-routing-checkpoint-v31"
TASK_ID = "KOKUNO-A5-SIGNED-PHYSICAL-COVARIANCE-ROUTING-031"

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
    "head_at_routing_audit": "c0f20e0712f4347698485a2de1e3cb980b5ecfb1",
    "latest_commit": "CR002: govern delivery-readiness scope (#478)",
    "active_kokuno_candidate_promoted_by_latest_commit": False,
}

AGENT1_RECEIPT = {
    "pr": 481,
    "head": "907430357a23dd611a324cefb25cdab7054d2bf8",
    "standard_run": 35391555503,
    "standard_status": "success",
    "shared_C_to_log_X_R_binding_executable": True,
    "selected_real_axis_C_is_autonomous_not_source_hidden_C": True,
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
    "pr": 482,
    "head": "bcf983ff0b8e059b9f5d2b051b46e839a9afb2c5",
    "latest_exact_head_standard_run": 35392997968,
    "latest_exact_head_standard_status": "success",
    "source_covariance_axis": "auxiliary_rectangle_sign_sigma_plus_minus",
    "supplied_signed_complete_curl_physical_family_executable": True,
    "sign_beta_total_Q_scaled_velocity_exposed": True,
    "amplitude_gradient_terms_retained": True,
    "support_gradient_terms_retained": True,
    "actual_positive_order_background_bound": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_signed_auxiliary_rectangles_or_modes_bound": False,
    "actual_source_partition_labels_instantiated": False,
    "public_source_bound_xyz_t_oscillatory_velocity_ready": False,
    "consumed_in_executable_ancestry": True,
}

AGENT3_RECEIPT = {
    "pr": 483,
    "head": "d57cf8e94768bce1a77d51ccc8828010695f9ea8",
    "dedicated_run": 35393033617,
    "standard_run": 35393033593,
    "dedicated_status": "success",
    "standard_status": "success",
    "supplied_source_compatible_structural_screen_only": True,
    "raw_sign_column_norms": [68.48634461905719, 48.771382276448676],
    "raw_velocity_singular_values": [68.8507850799386, 48.25553359265506],
    "raw_velocity_rank_ratio": 0.7008712179044639,
    "radial_cells": 5,
    "rank_two_cells": 5,
    "minimum_covariance_s_min": 17.945011875018206,
    "minimum_covariance_rank_ratio": 0.3203750745181774,
    "minimum_signed_response_novelty": 0.9072837449778101,
    "maximum_total_reconstruction_error": 0.0,
    "maximum_sign_tangent_sum_error": 8.881784197001252e-16,
    "maximum_quadratic_homogeneity_error": 1.6765755792695767e-16,
    "physical_to_reference_rule": "H_ref * y = DeltaC / epsilon",
    "real_candidate_defect_consumed": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "public_velocity_correction_materialized": False,
    "finite_correction_cycle_rerun_allowed": False,
    "consumed_in_executable_ancestry": False,
}

AGENT4_RECEIPT = {
    "pr": 484,
    "head": "44793b7f632bb0344bbd24d2a6074828b263b04b",
    "parent_agent2_head": "bcf983ff0b8e059b9f5d2b051b46e839a9afb2c5",
    "dedicated_run": 35393207116,
    "standard_run": 35393180096,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact_id": 10567070589,
    "artifact_digest": "sha256:7a196556ee07baa315f7a02f5943a3932b7ee4ddde173263fc6c64c90e35a7dd",
    "audit_seed": 9173191,
    "phase_resolutions": [12, 24, 48],
    "radii": [0.08, 0.35, 0.90],
    "target_step": 2.0e-4,
    "oracle_consumes_public_cartesian_total_only": True,
    "oracle_reuses_complete_curl_helper": False,
    "oracle_reads_returned_tangents_or_amplitudes": False,
    "minimum_physical_rank_ratio": 0.35573170032291024,
    "maximum_resolution_relative_drift": 4.5397547533405476e-13,
    "minimum_velocity_rms": 21.242258025207185,
    "maximum_duplicated_sign_mutation_rank_ratio": 4.565665180524497e-13,
    "frozen_minimum_rank_ratio_guard": 2.0e-2,
    "frozen_maximum_resolution_drift_guard": 1.0e-8,
    "frozen_minimum_velocity_rms_guard": 1.0e-4,
    "frozen_maximum_duplicated_sign_rank_ratio_guard": 1.0e-8,
    "local_physical_covariance_preflight_passed": True,
    "source_compatible_manufactured_input": True,
    "actual_source_physical_covariance_rank_two_assessed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": True,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 477,
    "head": "3950b5cc445cfe3ffe34e8284d6984268b231ff5",
    "schema": "kokuno-agent5-signed-covariance-routing-checkpoint-v30",
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
        "displayed_source_signed_covariance_reference_executable": True,
        "signed_covariance_reference_independently_audited": True,
        "signed_covariance_physical_to_reference_unit_bridge_ready": True,
        "supplied_signed_complete_curl_physical_family_executable": True,
        "supplied_signed_complete_curl_physical_covariance_rank_two": True,
        "supplied_signed_complete_curl_physical_covariance_independently_audited": True,
        "actual_positive_order_background_bound": False,
        "actual_source_h_sigma_pulse_integrals_bound": False,
        "actual_signed_auxiliary_rectangles_or_modes_bound": False,
        "actual_source_partition_labels_instantiated": False,
        "actual_source_physical_covariance_rank_two_assessed": False,
        "public_source_bound_xyz_t_velocity_ready": False,
        "oscillatory_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
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
            "status": "signed formulas are executable, but actual positive-order/background, h_sigma integrals, signed rectangles/modes and partition labels remain unbound",
        },
        {
            "stage": "leading_candidate",
            "ready": False,
            "status": "rescaled C/log-X_R bridge is executable; actual incoming moments, T_sh, PA.16 join and global matched pressure remain open",
        },
        {
            "stage": "oscillatory_augmentation",
            "ready": False,
            "status": "a supplied/source-compatible signed complete-curl physical family now passes independent covariance-rank preflight, but it is not the source-bound xyz,t family",
        },
        {
            "stage": "mean_radial_corrections",
            "ready": False,
            "status": "signed physical covariance direction survives supplied preflight, but no real candidate defect has been consumed and source binding is missing",
        },
        {
            "stage": "finite_correction_cycle",
            "ready": False,
            "status": "blocked until an independently audited source-bound signed family passes real-defect inverse plus unchanged unit/budget/spacetime/radial/quadratic/joint-gain guards",
        },
        {
            "stage": "candidate_artifact",
            "ready": False,
            "status": "schema reserved; no leading+oscillatory+corrected executable ancestry exists",
        },
        {
            "stage": "independent_validation",
            "ready": False,
            "status": "Agent 4 independently clears only the supplied physical covariance seam; formal full-domain PDE validation remains unassessed",
        },
        {
            "stage": "report_and_export",
            "ready": False,
            "status": "blocked until one frozen global candidate artifact exists",
        },
    ]
    routing = {
        "new_fact": "The source-correct sigma=+/- axis survives Q scaling and complete-curl materialization as a rank-two phase-mean covariance response on the supplied source-compatible physical family, and Agent 4 independently confirms that fact through a value-only Cartesian black-box oracle.",
        "classification": "supplied_signed_physical_covariance_seam_independently_cleared_actual_source_binding_still_open",
        "shortest_next_closure": [
            "Agent 2: stop repeating manufactured covariance-rank tests; bind provenance-labelled actual/source-motivated positive-order/background, h_sigma pulse integrals, signed rectangles/modes and source partition labels into the same complete-curl path, exposing a source-bound Q-scaled by-sign/by-beta/total velocity_osc(x,y,z,t) with explicit units.",
            "Agent 4: independently rerun the black-box covariance audit on that first source-bound/source-motivated physical family; the supplied manufactured PASS cannot be inherited as actual-source evidence.",
            "Agent 3: only after that source-bound audit, consume the real candidate defect through H_ref*y=DeltaC/epsilon and require the existing bounded inverse, coefficient budget, spacetime, radial +d_z sigma_1, quadratic and joint-gain guards before materializing delta_u and rerunning the finite correction cycle.",
            "Agent 1: continue the shared-C/log-X_R path into rescaled PA.15 incoming moments, a same-scale numerical T_sh certificate, PA.16 inner-to-outer join, I3/I4 and global matched pressure under fixed/restricted forcing.",
            "Agent 5: instantiate deterministic candidate save/load plus Python/MATLAB smoke only after leading, oscillatory and correction lanes coexist in one executable ancestry; then freeze it for Agent 4's formal held-out PDE gate and same-protocol ST006 comparison.",
        ],
        "do_not_do": [
            "do not promote supplied/source-compatible manufactured prototypes to actual Kokuno h_sigma, rectangles, modes, partitions or background",
            "do not promote the supplied physical covariance PASS to genuinely_independent_second_covariance_column_ready for the real candidate",
            "do not consume a structural calibration as the real candidate defect",
            "do not call covariance-rank, derivative or divergence implementation errors Navier-Stokes momentum residuals",
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
            "agent1": _copy(AGENT1_RECEIPT),
            "agent2": _copy(AGENT2_RECEIPT),
            "agent3": _copy(AGENT3_RECEIPT),
            "agent4": _copy(AGENT4_RECEIPT),
            "previous_agent5": _copy(PREVIOUS_AGENT5_RECEIPT),
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
        "displayed_source_signed_covariance_reference_executable",
        "signed_covariance_reference_independently_audited",
        "signed_covariance_physical_to_reference_unit_bridge_ready",
        "supplied_signed_complete_curl_physical_family_executable",
        "supplied_signed_complete_curl_physical_covariance_rank_two",
        "supplied_signed_complete_curl_physical_covariance_independently_audited",
    ):
        if not states[name]:
            raise ValueError(f"cleared supplied signed-physical seam was lost: {name}")

    for name in (
        "leading_ready",
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_signed_auxiliary_rectangles_or_modes_bound",
        "actual_source_partition_labels_instantiated",
        "actual_source_physical_covariance_rank_two_assessed",
        "public_source_bound_xyz_t_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
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

    a2 = payload["upstream"]["agent2"]
    if a2["source_covariance_axis"] != "auxiliary_rectangle_sign_sigma_plus_minus":
        raise ValueError("source signed covariance axis changed")
    if not a2["supplied_signed_complete_curl_physical_family_executable"]:
        raise ValueError("Agent-2 physical handoff was lost")
    for source_flag in (
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_signed_auxiliary_rectangles_or_modes_bound",
        "actual_source_partition_labels_instantiated",
        "public_source_bound_xyz_t_oscillatory_velocity_ready",
    ):
        if a2[source_flag]:
            raise ValueError(f"Agent-2 manufactured/source truth boundary changed: {source_flag}")

    a3 = payload["upstream"]["agent3"]
    if a3["rank_two_cells"] != a3["radial_cells"]:
        raise ValueError("Agent-3 supplied physical covariance rank receipt changed")
    if a3["minimum_covariance_rank_ratio"] <= 0.0:
        raise ValueError("Agent-3 supplied covariance rank receipt lost")
    if a3["real_candidate_defect_consumed"] or a3["genuinely_independent_second_covariance_column_ready"]:
        raise ValueError("Agent-3 supplied screen was promoted to real-candidate readiness")

    a4 = payload["upstream"]["agent4"]
    if not a4["local_physical_covariance_preflight_passed"]:
        raise ValueError("Agent-4 independent physical covariance preflight must remain passed")
    if a4["minimum_physical_rank_ratio"] < a4["frozen_minimum_rank_ratio_guard"]:
        raise ValueError("Agent-4 physical rank guard no longer passes")
    if a4["maximum_resolution_relative_drift"] > a4["frozen_maximum_resolution_drift_guard"]:
        raise ValueError("Agent-4 resolution guard no longer passes")
    if a4["minimum_velocity_rms"] < a4["frozen_minimum_velocity_rms_guard"]:
        raise ValueError("Agent-4 nontriviality guard no longer passes")
    if a4["maximum_duplicated_sign_mutation_rank_ratio"] > a4["frozen_maximum_duplicated_sign_rank_ratio_guard"]:
        raise ValueError("Agent-4 duplicated-sign negative control no longer passes")
    if a4["actual_source_physical_covariance_rank_two_assessed"]:
        raise ValueError("supplied Agent-4 audit was promoted to actual-source evidence")

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
