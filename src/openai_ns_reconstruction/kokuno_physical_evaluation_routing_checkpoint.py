"""Agent-5 routing checkpoint after the auxiliary-torus physical-evaluation audit.

This module is intentionally integration-only.  It binds the exact Agent-2/4
physical-evaluation handoff plus the latest Agent-1/3 sibling receipts into one
fail-closed checkpoint.  It does not synthesize the still-missing numerical
torus/background/mode data and it does not instantiate a final candidate.

Formal project gates remain unchanged: held-out normalized momentum max/L2
<= 1e-3 and divergence max/L2 <= 1e-5.  The Agent-4 numbers recorded here are
local derivative-contract diagnostics, not Navier--Stokes residuals.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-physical-evaluation-routing-checkpoint-v24"
TASK_ID = "KOKUNO-A5-PHYSICAL-EVALUATION-ROUTING-024"

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
    "head_at_routing_audit": "a2b215790fc483766be712233a134cd48ad7c149",
    "note": "live integration provenance only; this stacked PR executes on Agent-4 #412 ancestry",
}

AGENT1_RECEIPT = {
    "pr": 409,
    "head": "d4804cba6dfa962fee4f7a8a637e5b498f9b8086",
    "standard_run": 35350985320,
    "standard_status": "success",
    "selected_outer_velocity_assembly_available": True,
    "continuous_outer_route_logX": [322.48732140326854, 520.3167947396956],
    "inner_reference_to_outer_join_ready": False,
    "i3_ready": False,
    "i4_ready": False,
    "global_matched_pressure_ready": False,
    "global_leading_ready": False,
    "consumed_in_executable_ancestry": False,
}

AGENT2_RECEIPT = {
    "pr": 410,
    "head": "482c86c9e8cf10dc87e2c857c97112e74c28809f",
    "dedicated_run": 35351242764,
    "standard_run": 35351242746,
    "dedicated_status": "success",
    "standard_status": "success",
    "physical_evaluation_map_executable": True,
    "source_torus_vectors_or_group_constants_recovered": False,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_q_scaled_by_beta_total_xyz_t_velocity_ready": False,
    "consumed_in_executable_ancestry": True,
}

AGENT3_RECEIPT = {
    "pr": 411,
    "head": "196030f75952befff199635d5225b6c4121e27df",
    "dedicated_run": 35351943779,
    "standard_run": 35351943712,
    "agent5_reference_run": 35351943874,
    "artifact_id": 10550075968,
    "artifact_digest": "sha256:cdb0a167b2192287305907530ec9f7002a7db6c028569ea64ee68b44b913133c",
    "quadratic_covariance_gain_guard_ready": True,
    "fractional_common_budget_at_reference_amplitude": 0.10010444711945854,
    "nonlinear_to_linear_covariance_ratio_at_budget_boundary": 0.05005222355972927,
    "actual_source_multiband_family_consumed": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "finite_correction_cycle_rerun_allowed": False,
    "consumed_in_executable_ancestry": False,
}

AGENT4_RECEIPT = {
    "pr": 412,
    "head": "654fec264c87ad5e6d59f1841206e7335be1331a",
    "base_agent2_head": AGENT2_RECEIPT["head"],
    "dedicated_run": 35352515086,
    "standard_run": 35352515085,
    "artifact_name": "kokuno-agent4-physical-evaluation-independent-audit",
    "artifact_id": 10550406550,
    "artifact_digest": "sha256:6078379d32070284bffa6b2f7e75f794937751d2b9919bcba186b05c892fbebc",
    "seed": 9173141,
    "fd6_steps": [0.008, 0.004, 0.002],
    "case_primary_relative_rms_ladders": [
        [5.2625412326405895e-06, 5.090287477700525e-08, 7.22842318951028e-10],
        [3.6390506234329835e-05, 3.837121458902741e-07, 5.526301830533695e-09],
        [3.911401581764177e-06, 6.205579135196832e-08, 9.733354410114198e-10],
    ],
    "worst_finest_primary_relative_rms": 5.526301830533695e-09,
    "worst_refinement_ratio": 63.03040371493245,
    "weakest_omit_torus_mutation_relative_rms": 0.9872565807423554,
    "local_guards": {
        "finest_primary_relative_rms_max": 2.0e-7,
        "minimum_refinement_ratio": 6.0,
        "omit_torus_mutation_relative_rms_min": 5.0e-3,
    },
    "designed_mod1_seam_points_exercised_per_case": 4,
    "auxiliary_torus_pullback_independently_preflighted": True,
    "local_structural_preflight_passed": True,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": True,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 403,
    "head": "d5f6e0aa1d19fbd971bc6c9c9b5e204ed79c3b34",
    "dedicated_run": 35347637558,
    "standard_run": 35347637825,
    "reference_run": 35347637556,
    "artifact_id": 10547428458,
    "artifact_digest": "sha256:da37599b95a78af7c0450fbcd246726949d05960ba5b1ed14c6855fe923fd98a",
    "checkpoint_payload_sha256": "7dd1d629022c14b3a56e659fb6d4676777cb7bc2e9dfaa4706ff2040f3d54572",
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
    "comparison_to_agent4_local_fd6": "not directly comparable",
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
        "actual_source_or_autonomous_torus_data_bound": False,
        "actual_positive_order_background_bound": False,
        "actual_auxiliary_torus_mode_family_bound": False,
        "public_q_scaled_by_beta_total_velocity_ready": False,
        "oscillatory_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "quadratic_covariance_gain_guard_ready": True,
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
            "status": "source_formulas_and_versions_bound_autonomous_numeric_choices_still_required",
            "ready": False,
        },
        {
            "stage": "leading_candidate",
            "status": "selected_outer_assembly_available_inner_join_i3_i4_and_global_pressure_missing",
            "ready": False,
        },
        {
            "stage": "oscillatory_augmentation",
            "status": "physical_evaluation_contract_independently_cleared_waiting_for_numeric_torus_background_and_modes",
            "ready": False,
        },
        {
            "stage": "mean_radial_corrections",
            "status": "quadratic_guard_ready_waiting_for_actual_multiband_rank_and_bounded_inverse",
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
            "status": "local_physical_evaluation_preflight_passed_formal_full_domain_gate_not_assessed",
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
            "Agent-4 #412 independently clears Agent-2 #410's auxiliary-torus physical-evaluation differential pullback under frozen FD6 guards."
        ),
        "interpretation": (
            "The shortest oscillatory blocker is no longer chain-rule algebra. It is the missing numerical source-compatible/autonomous torus directions and group constants together with the positive-order/background and actual multi-band mode family needed to expose a public Q-scaled by-beta/total xyz,t velocity."
        ),
        "shortest_next_closure": [
            "Agent 2: bind provenance-labelled numerical torus directions/group constants, positive-order/background data and actual mode family into the already-cleared physical-evaluation/localized-curl chain; expose Q-scaled by-beta plus total velocity(x,y,z,t) with at least two active slow bands and explicit coefficient units.",
            "Agent 3: consume that exact physical family through cross-label covariance, unit conversion, bounded inverse and quadratic self-covariance; require rank, aggregate budget, spacetime and retained radial d_z sigma_1 guards before materializing a correction or running the finite cycle.",
            "Agent 1: close the inner-reference to RF40 join, I3/I4 and global matched pressure with only fixed/restricted forcing; then expose the first global leading velocity/pressure pair.",
            "Agent 4: independently audit the first actual-family derivative/tangent/rank receipt and, only after a complete frozen composite exists, run the formal held-out normalized NS gate with same-protocol ST006 comparison.",
            "Agent 5: instantiate the deterministic candidate artifact and Python/MATLAB smoke only after leading, oscillatory and correction stages coexist in one executable ancestry.",
        ],
        "do_not_do": [
            "do not repeat the cleared auxiliary-torus chain-rule audit unless that interface changes",
            "do not call repository-autonomous torus values recovered Kokuno hidden data",
            "do not call Agent-4 local FD6 derivative errors Navier-Stokes residuals",
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
    if not states["auxiliary_torus_physical_evaluation_independently_audited"]:
        raise ValueError("independent physical-evaluation audit must remain bound")
    for key in (
        "leading_ready",
        "actual_source_or_autonomous_torus_data_bound",
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
    if payload["upstream"]["agent4"]["formal_full_domain_pde_gate_assessed"]:
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
        default="artifacts/kokuno_agent5/physical_evaluation_routing_checkpoint_v24/checkpoint.json",
    )
    args = parser.parse_args(argv)
    payload = write_checkpoint(args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
