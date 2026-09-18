"""Agent-5 integration checkpoint for the shortest Kokuno closure.

This file binds immutable lane receipts and validation semantics only. It does
not reconstruct missing mathematics, instantiate a candidate, or promote local
diagnostics to the formal Navier--Stokes gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-disjoint-validation-routing-checkpoint-v26"
TASK_ID = "KOKUNO-A5-DISJOINT-VALIDATION-ROUTING-026"
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
    "head_at_routing_audit": "938145416760477308103f3f251edaa201471d7c",
    "latest_commit": "CR002: govern ST048-S temporal Piola semantics (#426)",
    "active_route_promoted_by_latest_commit": False,
}
AGENT1_RECEIPT = {
    "pr": 429,
    "head": "6d65d563db03b5672eebfb0bb5476e3cb52799ca",
    "standard_run": 35363221632,
    "standard_status": "success",
    "appendix_B_Xi_profile_executable": True,
    "selected_realization_boundary_values_executable": True,
    "source_hidden_numeric_choices_recovered": False,
    "incoming_five_moment_discrepancy_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "inner_to_outer_join_completed": False,
    "i3_ready": False,
    "i4_ready": False,
    "global_matched_pressure_ready": False,
    "global_leading_ready": False,
    "consumed_in_executable_ancestry": False,
}
AGENT2_RECEIPT = {
    "pr": 430,
    "head": "167d8e750d7ba5e4abaeaec7fc93ccf70145b7fd",
    "dedicated_run": 35363035071,
    "standard_run": 35363035013,
    "dedicated_status": "success",
    "standard_status": "success",
    "source_band_covering_schedule_executable": True,
    "displayed_source_torus_constants_bound": True,
    "source_rectangle_centers_recovered": False,
    "source_rectangle_radius_r0_recovered": False,
    "source_partition_grid_origin_recovered": False,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_q_scaled_by_beta_total_xyz_t_velocity_ready": False,
    "consumed_in_executable_ancestry": False,
}
AGENT3_RECEIPT = {
    "pr": 431,
    "head": "9d2c234234e46b53bf141d08d458ae206409ca4d",
    "dedicated_run": 35364119783,
    "standard_run": 35364119580,
    "agent5_reference_run": 35364119732,
    "dedicated_status": "success",
    "standard_status": "success",
    "agent5_reference_status": "success",
    "artifact_name": "kokuno-agent3-spacetime-damping-holdout-transfer",
    "artifact_id": 10555722075,
    "artifact_digest": "sha256:cc069df53d630fbb66aae8af3a374fc9d61fe148071673d45fb66e763d612b00",
    "calibration_times": [0.375, 0.5, 0.625],
    "calibration_z": [0.06, 0.08, 0.10],
    "validation_times": [0.3125, 0.6875],
    "validation_z": [0.07, 0.09],
    "damping_calibration_validation_disjoint": True,
    "heldout_damping_reoptimized": False,
    "validation_objective_used_for_selection": False,
    "shared_selected_damping": 0.4926349922838171,
    "calibration_relative_stress_residual_rms": 0.6223507487797298,
    "validation_fixed_damping_relative_stress_residual_rms": 0.6223757026222282,
    "validation_full_step_relative_stress_residual_rms": 0.9861371010802802,
    "validation_transfer_preflight_passed": True,
    "required_rank_two_nodes": 25,
    "rank_two_nodes": 0,
    "one_band_algebraic_relative_stress_residual": 0.2101945021227428,
    "upstream_bounded_inverse_passed": False,
    "actual_source_multiband_family_consumed": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "finite_correction_cycle_rerun_allowed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": True,
}
PREVIOUS_AGENT4_RECEIPT = {
    "pr": 423,
    "head": "a0258945aef7e24690b03baab996edead17b6716",
    "displayed_source_torus_binding_independently_audited": True,
    "worst_finest_primary_relative_rms": 6.131909380022552e-10,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}
AGENT4_RECEIPT = {
    "pr": 432,
    "head": "ada7e63fe8cd95acbb5c1b8bb2d646a1eefc644d",
    "base_agent2_head": AGENT2_RECEIPT["head"],
    "dedicated_run": 35364815755,
    "standard_run": 35364815707,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact_id": 10556286420,
    "artifact_digest": "sha256:eba4bb13efc2a9c19fe406621baaa681e6d921d83ee8d27592b1f404153dc5a3",
    "source_band_covering_independently_audited": True,
    "local_structural_preflight_passed": True,
    "broad_cases": 54,
    "worst_value_relative_error": 2.2898921537313096e-14,
    "covering_level_mismatches": 0,
    "transition_near_cases": 16,
    "transition_near_mismatches": 0,
    "interacting_band_pairs": 1245,
    "interacting_band_violations": 0,
    "maximum_observed_level_difference": 2,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}
PREVIOUS_AGENT5_RECEIPT = {
    "pr": 424,
    "head": "dc795ac0729e8d0e1e75ccfb5d629a5f461d2a4e",
    "dedicated_run": 35359194541,
    "standard_run": 35359194831,
    "artifact_id": 10553412016,
    "artifact_digest": "sha256:09968d17a6928ab1cf7a3812f151c510c8743d0a607cbc653024994d11d0432e",
    "checkpoint_payload_sha256": "bb3717915db2ec07d86ccda1b5f4a82b56864c7dac5059ab61db1cfc6e0eb840",
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
        "Agent-3 correction receipt with calibration/validation split",
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


def _pairs(times: list[float], zs: list[float]) -> set[tuple[float, float]]:
    return {(float(t), float(z)) for z in zs for t in times}


def build_checkpoint() -> dict[str, Any]:
    states = {
        "leading_ready": False,
        "appendix_B_Xi_profile_executable": True,
        "source_band_covering_schedule_executable": True,
        "source_band_covering_independently_audited": True,
        "displayed_source_torus_binding_independently_audited": True,
        "correction_calibration_validation_disjoint": True,
        "local_damping_transfer_preflight_passed": True,
        "actual_positive_order_background_bound": False,
        "actual_auxiliary_torus_mode_family_bound": False,
        "public_q_scaled_by_beta_total_velocity_ready": False,
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
        {"stage": "source_profile_ingest", "ready": False, "status": "Xi_boundary_and_band_schedule_audited_missing_incoming_five_moment_Tsh_and_actual_positive_order_multiband_modes"},
        {"stage": "leading_candidate", "ready": False, "status": "selected_Xi_boundary_ready_missing_incoming_discrepancy_Tsh_inner_join_I3_I4_global_pressure"},
        {"stage": "oscillatory_augmentation", "ready": False, "status": "source_torus_and_band_covering_audited_missing_actual_positive_order_background_multiband_modes_and_public_by_beta_total_velocity"},
        {"stage": "mean_radial_corrections", "ready": False, "status": "disjoint_damping_transfer_passes_but_one_band_rank_two_is_0_of_25"},
        {"stage": "finite_correction_cycle", "ready": False, "status": "blocked_until_actual_multiband_family_rank_budget_spacetime_radial_quadratic_guards_pass"},
        {"stage": "candidate_artifact", "ready": False, "status": "schema_reserved_not_instantiated"},
        {"stage": "independent_validation", "ready": False, "status": "A4_source_torus_and_band_local_audits_exist_A3_holdout_is_not_formal_A4_full_domain_PDE_gate"},
        {"stage": "report_and_export", "ready": False, "status": "blocked_until_frozen_candidate_artifact_exists"},
    ]
    routing = {
        "new_fact": "A3 #431 now uses disjoint calibration/validation for frozen damping, while latest A4 #432 independently clears A2 #430's deterministic source band-covering schedule with high-precision arithmetic. Neither result supplies the missing actual multi-band physical velocity family.",
        "interpretation": "Validation reuse and source-band arithmetic are no longer shortest-path blockers. The main cross-lane blocker is still A2's actual multi-band physical family and genuinely independent covariance direction; A1 global leading closure proceeds in parallel.",
        "shortest_next_closure": [
            "Agent 2: bind actual/source-compatible positive-order background plus a real multi-band auxiliary-torus mode family to the cleared partition/curl/physical-evaluation/band-covering chain; expose Q-scaled by-beta and total velocity(x,y,z,t), explicit units, and at least two active slow bands.",
            "Agent 3: consume exactly that family, require genuine rank-two coverage plus unchanged unit/budget/spacetime/radial/quadratic guards, then repeat the disjoint damping-transfer protocol before materialization and finite-cycle rerun.",
            "Agent 1: propagate the executable Xi boundary through incoming five-moment discrepancy and numerical T_sh certificate, complete inner-to-outer join, I3/I4, and global matched pressure under fixed/restricted forcing.",
            "Agent 4: stop re-auditing source-band arithmetic; independently audit the first actual multi-band physical family and later the frozen full-domain composite with same-protocol ST006 comparison.",
            "Agent 5: instantiate/save/load/export only after leading, oscillatory and correction stages coexist in one executable ancestry.",
        ],
        "do_not_do": [
            "do not call A3 local stress-transfer residual an NS residual",
            "do not treat disjoint A3 validation as Agent-4 independent PDE validation",
            "do not fabricate rank two with duplicate columns or per-label free amplitudes",
            "do not call autonomous/caller-supplied hidden source choices recovered data",
            "do not promote the finite cycle before all frozen guards pass",
            "do not weaken the 1e-3 momentum or 1e-5 divergence gates",
            "do not use residual-defined free forcing",
        ],
    }
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "source_provenance": SOURCE_PROVENANCE,
        "live_integration": LIVE_INTEGRATION_RECEIPT,
        "upstream": {"agent1": AGENT1_RECEIPT, "agent2": AGENT2_RECEIPT, "agent3": AGENT3_RECEIPT, "agent4": AGENT4_RECEIPT, "previous_agent4": PREVIOUS_AGENT4_RECEIPT, "previous_agent5": PREVIOUS_AGENT5_RECEIPT},
        "states": states,
        "pipeline_frontier": pipeline,
        "candidate_artifact_contract": CANDIDATE_ARTIFACT_CONTRACT,
        "baseline_vs_kokuno": {"st006": ST006_BASELINE, "kokuno_current_comparable_full_domain_receipt": None, "comparison_status": "not_assessable_until_global_frozen_kokuno_candidate_exists"},
        "formal_gates": {"held_out_normalized_momentum_max": 1.0e-3, "held_out_normalized_momentum_l2": 1.0e-3, "held_out_divergence_max": 1.0e-5, "held_out_divergence_l2": 1.0e-5, "changed_this_round": False},
        "routing": routing,
    }
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def validate_checkpoint(payload: dict[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or payload.get("task_id") != TASK_ID:
        raise ValueError("wrong Agent-5 checkpoint identity")
    if payload.get("checkpoint_sha256") != checkpoint_sha256(payload):
        raise ValueError("checkpoint sha256 mismatch")
    a3 = payload["upstream"]["agent3"]
    if _pairs(a3["calibration_times"], a3["calibration_z"]) & _pairs(a3["validation_times"], a3["validation_z"]):
        raise ValueError("damping calibration and validation grids overlap")
    if not a3["damping_calibration_validation_disjoint"]:
        raise ValueError("damping calibration/validation split must remain disjoint")
    if a3["heldout_damping_reoptimized"] or a3["validation_objective_used_for_selection"]:
        raise ValueError("held-out damping validation must not select or retune damping")
    if not a3["validation_transfer_preflight_passed"]:
        raise ValueError("recorded disjoint damping transfer must remain passed")
    if a3["upstream_bounded_inverse_passed"] or a3["finite_correction_cycle_rerun_allowed"]:
        raise ValueError("local damping transfer cannot bypass the rank/bounded-inverse blocker")
    a4 = payload["upstream"]["agent4"]
    if a4["base_agent2_head"] != payload["upstream"]["agent2"]["head"] or not a4["source_band_covering_independently_audited"]:
        raise ValueError("latest Agent-4 source-band audit is not bound to Agent-2 exact head")
    if a4["covering_level_mismatches"] or a4["transition_near_mismatches"] or a4["interacting_band_violations"]:
        raise ValueError("recorded Agent-4 source-band audit is not clean")
    states = payload["states"]
    for key in ("leading_ready", "actual_positive_order_background_bound", "actual_auxiliary_torus_mode_family_bound", "public_q_scaled_by_beta_total_velocity_ready", "oscillatory_ready", "genuinely_independent_second_covariance_column_ready", "correction_ready", "finite_correction_cycle_run", "candidate_artifact_instantiated", "velocity_export_ready", "formal_full_domain_pde_gate_assessed", "pde_validated", "paper_exact"):
        if states[key]:
            raise ValueError(f"truth boundary violated: {key}=true")
    if payload["formal_gates"]["held_out_normalized_momentum_max"] != 1.0e-3:
        raise ValueError("formal momentum gate changed")
    if payload["formal_gates"]["held_out_divergence_max"] != 1.0e-5:
        raise ValueError("formal divergence gate changed")
    if payload["candidate_artifact_contract"]["instantiated"]:
        raise ValueError("candidate artifact must remain reserved, not instantiated")
    if a4["formal_full_domain_pde_gate_assessed"] or a4["pde_validated"]:
        raise ValueError("local Agent-4 audit cannot be promoted to formal PDE validation")


def write_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/kokuno_agent5/disjoint_validation_routing_checkpoint_v26/checkpoint.json")
    args = parser.parse_args(argv)
    payload = write_checkpoint(args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
