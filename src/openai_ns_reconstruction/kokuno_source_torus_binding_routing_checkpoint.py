"""Agent-5 routing checkpoint after the displayed source-torus binding audit.

This module is integration-only. It binds the exact Agent-2 source-constant
handoff and Agent-4 independent audit, while recording the latest Agent-1 and
Agent-3 sibling receipts without pretending they are executable ancestry.
It does not invent the still-missing positive-order/background or auxiliary-
torus mode family and it does not instantiate a final candidate.

Formal project gates remain unchanged: held-out normalized momentum max/L2
<= 1e-3 and divergence max/L2 <= 1e-5. The Agent-4 values recorded here are
local source-binding/derivative diagnostics, not Navier--Stokes residuals.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-source-torus-binding-routing-checkpoint-v25"
TASK_ID = "KOKUNO-A5-SOURCE-TORUS-BINDING-ROUTING-025"

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
    "head_at_routing_audit": "7749a3e4d8d42b3cfc74670a713f279152ec5f39",
    "note": "live integration provenance only; this stacked PR executes on Agent-4 #423 ancestry",
}

AGENT1_RECEIPT = {
    "pr": 421,
    "head": "cac0fe6d6172b544f800d7330944c4ad074d2ea1",
    "standard_run": 35358052461,
    "standard_status": "success",
    "source_later_inner_join_formula_executable": True,
    "actual_upstream_appendix_B_boundary_data_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "inner_to_outer_join_completed": False,
    "i3_ready": False,
    "i4_ready": False,
    "global_matched_pressure_ready": False,
    "global_leading_ready": False,
    "consumed_in_executable_ancestry": False,
}

AGENT2_RECEIPT = {
    "pr": 420,
    "head": "1c77b07a0e8286b8ffa487917f48d54f6017c7e0",
    "dedicated_run": 35356699998,
    "standard_run": 35356699947,
    "dedicated_status": "success",
    "standard_status": "success",
    "displayed_source_torus_constants_bound": True,
    "source_torus_vectors_displayed_and_bound": True,
    "source_group_constants_displayed_and_bound": True,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_q_scaled_by_beta_total_xyz_t_velocity_ready": False,
    "consumed_in_executable_ancestry": True,
}

AGENT3_RECEIPT = {
    "pr": 422,
    "head": "96c250baedc6ad7af12acf954ed9166b045cb916",
    "dedicated_run": 35357590985,
    "standard_run": 35357590966,
    "agent5_reference_run": 35357591217,
    "artifact_id": 10552751682,
    "artifact_digest": "sha256:d87288a44f34c6a32b9f279181c2cfa12cfcb75836a0fba56cbe6484ff81e668",
    "shared_spacetime_quadratic_damping_ready": True,
    "negative_shared_lambda": 0.4926349922838171,
    "negative_shared_damped_l1_update": 0.049314953534270234,
    "negative_aggregate_exact_relative_stress_residual": 0.6223507487797298,
    "required_rank_two_nodes": 25,
    "rank_two_nodes": 0,
    "one_band_algebraic_relative_stress_residual": 0.2101945021227428,
    "actual_source_multiband_family_consumed": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "finite_correction_cycle_rerun_allowed": False,
    "consumed_in_executable_ancestry": False,
}

AGENT4_RECEIPT = {
    "pr": 423,
    "head": "a0258945aef7e24690b03baab996edead17b6716",
    "base_agent2_head": AGENT2_RECEIPT["head"],
    "dedicated_run": 35358659588,
    "standard_run": 35358659698,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact_name": "kokuno-agent4-source-torus-binding-independent-audit",
    "artifact_id": 10554315174,
    "artifact_digest": "sha256:d8e6b3873461eb708e350cab280b1c0383d7023c5a65bf00c0753e4008f4b5ba",
    "seed": 9173151,
    "fd6_steps": [0.006, 0.003, 0.0015],
    "worst_finest_primary_relative_rms": 6.131909380022552e-10,
    "worst_phase_embedding_max_abs_error": 1.1379786002407855e-15,
    "worst_refinement_ratio": 68.37941008216457,
    "worst_source_binding_relative_error": 0.0,
    "weakest_wrong_sign_mutation_relative_rms": 1.2578253545072033,
    "weakest_exponent_mutation_relative_rms": 0.016006724435756885,
    "local_guards": {
        "phase_embedding_max_abs": 2.0e-13,
        "source_binding_relative": 2.0e-14,
        "finest_relative_rms": 1.0e-7,
        "minimum_refinement_ratio": 20.0,
        "wrong_sign_mutation_relative_rms_floor": 0.5,
        "exponent_mutation_relative_rms_floor": 1.0e-2,
    },
    "displayed_source_torus_binding_independently_audited": True,
    "local_structural_preflight_passed": True,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": True,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 413,
    "head": "9bfaba92a9b70d733b3fc95f011ef9dfa1236430",
    "dedicated_run": 35353272545,
    "standard_run": 35353272583,
    "reference_run": 35353272556,
    "artifact_id": 10550797884,
    "artifact_digest": "sha256:460355d02c112c725e915a5754222a023910e4ea0cf6c305522f11698ff1c621",
    "checkpoint_payload_sha256": "c371f6417a91ebcfd0544714d6b8e5c745bfa12c78502ccc9098563513670e31",
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
    "formal_momentum_gate": 1.0e-3,
    "formal_divergence_gate": 1.0e-5,
    "pde_validated": False,
    "comparison_to_agent4_local_source_binding_audit": "not directly comparable",
}

CANDIDATE_ARTIFACT_CONTRACT = {
    "schema_reserved": "kokuno-derived-3d-candidate-v1",
    "instantiated": False,
    "required_provenance": [
        "source repository/version/commit",
        "repository-autonomous parameter choices",
        "Agent-1 leading receipt",
        "Agent-2 oscillatory receipt",
        "Agent-3 correction receipt",
        "Agent-4 independent-validation receipt",
    ],
    "required_public_api": [
        "velocity(x,y,z,t)->[u,v,w]",
        "pressure(x,y,z,t)",
        "restricted_forcing(x,y,z,t)",
        "save/load deterministic parameter artifact",
    ],
    "required_exports_after_instantiation": ["Python smoke", "MATLAB smoke"],
    "forcing_policy": "fixed/restricted only; residual-defined free forcing forbidden",
}


def checkpoint_sha256(payload: dict[str, Any]) -> str:
    clean = dict(payload)
    clean.pop("checkpoint_sha256", None)
    raw = json.dumps(clean, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def build_checkpoint() -> dict[str, Any]:
    states = {
        "leading_ready": False,
        "oscillatory_machinery_ready": True,
        "auxiliary_torus_physical_evaluation_independently_audited": True,
        "displayed_source_torus_constants_bound": True,
        "displayed_source_torus_binding_independently_audited": True,
        "actual_positive_order_background_bound": False,
        "actual_auxiliary_torus_mode_family_bound": False,
        "public_q_scaled_by_beta_total_velocity_ready": False,
        "oscillatory_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "shared_spacetime_quadratic_damping_ready": True,
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
        {
            "stage": "source_profile_ingest",
            "status": "displayed_auxiliary_torus_constants_bound_and_independently_audited_missing_positive_order_background_and_mode_family",
            "ready": False,
        },
        {
            "stage": "leading_candidate",
            "status": "later_inner_join_formula_executable_actual_upstream_appendix_B_data_i3_i4_and_global_pressure_missing",
            "ready": False,
        },
        {
            "stage": "oscillatory_augmentation",
            "status": "torus_constants_and_pullback_cleared_waiting_for_positive_order_background_multiband_modes_and_public_xyz_t_velocity",
            "ready": False,
        },
        {
            "stage": "mean_radial_corrections",
            "status": "shared_damping_ready_waiting_for_actual_multiband_rank_and_bounded_inverse",
            "ready": False,
        },
        {
            "stage": "finite_correction_cycle",
            "status": "blocked_until_real_second_direction_and_all_frozen_guards_pass",
            "ready": False,
        },
        {
            "stage": "candidate_artifact",
            "status": "schema_reserved_not_instantiated",
            "ready": False,
        },
        {
            "stage": "independent_validation",
            "status": "displayed_source_torus_binding_preflight_passed_formal_full_domain_gate_not_assessed",
            "ready": False,
        },
        {
            "stage": "report_and_export",
            "status": "blocked_until_candidate_artifact_exists",
            "ready": False,
        },
    ]
    routing = {
        "new_fact": (
            "Agent-2 #420 binds the corrected reader's displayed auxiliary-torus constants and Agent-4 #423 independently verifies that binding and derivative pullback under frozen FD6 guards."
        ),
        "interpretation": (
            "Numerical torus directions/group constants are no longer an oscillatory blocker. The shortest missing Agent-2 object is now the positive-order/background plus actual multi-band auxiliary-torus mode family needed to expose public Q-scaled by-beta/total xyz,t velocity and a genuinely independent covariance direction."
        ),
        "shortest_next_closure": [
            "Agent 2: bind the positive-order/background data and actual multi-band auxiliary-torus mode family into the already-cleared partition/localized-curl/physical-evaluation chain; expose Q-scaled by-beta plus total velocity(x,y,z,t), explicit coefficient units, and at least two active slow bands.",
            "Agent 3: consume that exact physical family through cross-label covariance, unit conversion, bounded inverse, frozen spacetime/radial guards, nonlinear self-covariance and the already-audited shared damping; require genuine rank-two coverage before materializing a correction or running the finite cycle.",
            "Agent 1: bind the still-missing upstream Appendix-B boundary data/T_sh contract, complete the inner-to-outer join, I3/I4 and global matched pressure under fixed/restricted forcing; only then expose leading_ready=true.",
            "Agent 4: stop re-auditing the displayed torus constants unless that interface changes; independently audit the first actual multi-band physical family and later own the formal full-domain normalized NS gate with same-protocol ST006 comparison.",
            "Agent 5: instantiate the deterministic candidate artifact and Python/MATLAB smoke only after leading, oscillatory and correction stages coexist in one executable ancestry.",
        ],
        "do_not_do": [
            "do not call unpublished pulse/background/mode choices recovered Kokuno data",
            "do not call Agent-4 local FD6 source-binding errors Navier-Stokes residuals",
            "do not fabricate a second covariance direction with duplicate or per-label free columns",
            "do not rerun the finite correction cycle before the frozen rank/budget/spacetime/radial/quadratic guards pass",
            "do not instantiate or export a final Kokuno candidate before all three mathematical lanes coexist",
            "do not weaken the fixed 1e-3 momentum or 1e-5 divergence gates",
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
    if not states["displayed_source_torus_constants_bound"]:
        raise ValueError("displayed source torus constants must remain bound")
    if not states["displayed_source_torus_binding_independently_audited"]:
        raise ValueError("displayed source torus binding must remain independently audited")
    if payload["upstream"]["agent4"]["base_agent2_head"] != payload["upstream"]["agent2"]["head"]:
        raise ValueError("Agent-4 source-binding audit is not tied to the recorded Agent-2 head")
    for key in (
        "leading_ready",
        "actual_positive_order_background_bound",
        "actual_auxiliary_torus_mode_family_bound",
        "public_q_scaled_by_beta_total_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
        "paper_exact",
    ):
        if states[key]:
            raise ValueError(f"truth boundary violated: {key}=true")
    if payload["formal_gates"]["held_out_normalized_momentum_max"] != 1.0e-3:
        raise ValueError("formal momentum gate changed")
    if payload["formal_gates"]["held_out_divergence_max"] != 1.0e-5:
        raise ValueError("formal divergence gate changed")
    if payload["candidate_artifact_contract"]["instantiated"]:
        raise ValueError("candidate artifact must remain reserved, not instantiated")
    a4 = payload["upstream"]["agent4"]
    if not a4["local_structural_preflight_passed"]:
        raise ValueError("recorded Agent-4 source-binding preflight is not passed")
    if a4["formal_full_domain_pde_gate_assessed"] or a4["pde_validated"]:
        raise ValueError("local Agent-4 audit cannot be promoted to the formal PDE gate")


def write_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent5/source_torus_binding_routing_checkpoint_v25/checkpoint.json",
    )
    args = parser.parse_args(argv)
    payload = write_checkpoint(args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
