"""Agent-5 routing checkpoint for Kokuno source-realization closure.

This checkpoint records immutable lane receipts and fail-closed integration state.
It does not reconstruct missing mathematics, instantiate a candidate, or promote
local source-geometry checks to the formal Navier--Stokes validation gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-rational-rectangle-routing-checkpoint-v27"
TASK_ID = "KOKUNO-A5-RATIONAL-RECTANGLE-ROUTING-027"

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
    "head_at_routing_audit": "9bca40a2ce79c7f9438023a6c63a98b19bfc0fb1",
    "latest_commit": "CR002: govern post-Piola swirl-gain semantics (#434)",
    "active_kokuno_route_promoted_by_latest_commit": False,
}

AGENT1_RECEIPT = {
    "pr": 438,
    "head": "417a43e9259a8597549ed3345a763b0b9021691e",
    "standard_run": 35369433940,
    "standard_status": "success",
    "selected_appendix_B_prefix_moments_executable": True,
    "selected_PA15_scaled_incoming_discrepancy_executable": True,
    "source_outer_pressure_datum_bound": False,
    "actual_source_incoming_five_moment_discrepancy_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "consumed_in_executable_ancestry": False,
}

AGENT2_RECEIPT = {
    "pr": 439,
    "head": "e5bc36691d16be09cab602885fc49cb99f3d2f71",
    "dedicated_run": 35369821905,
    "standard_run": 35369821897,
    "dedicated_status": "success",
    "standard_status": "success",
    "source_rational_rectangle_separation_contract_executable": True,
    "autonomous_rational_witness_executable": True,
    "strict_r0_guards_executable": True,
    "band_pulse_length_compatibility_executable": True,
    "source_rectangle_centers_recovered": False,
    "source_rectangle_radius_r0_recovered": False,
    "source_color_count_recovered": False,
    "source_center_denominator_recovered": False,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_q_scaled_by_beta_total_xyz_t_velocity_ready": False,
    "consumed_in_executable_ancestry": True,
}

AGENT3_STABLE_RECEIPT = {
    "pr": 431,
    "head": "9d2c234234e46b53bf141d08d458ae206409ca4d",
    "dedicated_run": 35364119783,
    "standard_run": 35364119580,
    "dedicated_status": "success",
    "standard_status": "success",
    "damping_calibration_validation_disjoint": True,
    "shared_selected_damping": 0.4926349922838171,
    "validation_fixed_damping_relative_stress_residual_rms": 0.6223757026222282,
    "required_rank_two_nodes": 25,
    "rank_two_nodes": 0,
    "one_band_algebraic_relative_stress_residual": 0.2101945021227428,
    "upstream_bounded_inverse_passed": False,
    "finite_correction_cycle_rerun_allowed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

AGENT3_LATEST_PENDING = {
    "pr": 440,
    "head_at_routing_audit": "ea593e2f5ed128a4f31a8da3419755e49fa8ad3f",
    "purpose": "retain source-required radial +d_z sigma_1 in held-out damping guard",
    "dedicated_run": 35370309668,
    "standard_run": 35370309720,
    "dedicated_status_at_routing_audit": "in_progress",
    "standard_status_at_routing_audit": "in_progress",
    "promoted_to_stable_receipt": False,
}

AGENT4_RECEIPT = {
    "pr": 432,
    "head": "ada7e63fe8cd95acbb5c1b8bb2d646a1eefc644d",
    "dedicated_run": 35364815755,
    "standard_run": 35364815707,
    "dedicated_status": "success",
    "standard_status": "success",
    "source_band_covering_independently_audited": True,
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
    "pr": 433,
    "head": "4eb000cdc891ba03e00f1dc8cc34bc5b6c9251ca",
    "artifact_id": 10556815188,
    "artifact_digest": "sha256:494eb7861804331b83baa98bbe5bf01e9575830902c24accc4da5911a3245ec2",
    "checkpoint_payload_sha256": "5fdf1fe21ec8c13165cfb80575e132276d2df5b843bef47fd936553e3e3820f2",
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
        "repository-autonomous choices including rational rectangle witness",
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
        "source_rational_rectangle_separation_contract_executable": True,
        "autonomous_rational_rectangle_witness_ready": True,
        "source_rectangle_parameters_recovered": False,
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
        {"stage": "source_profile_ingest", "ready": False, "status": "source constants/band schedule/rational separation executable; hidden rectangle tuple remains autonomous and leading pressure datum remains unbound"},
        {"stage": "leading_candidate", "ready": False, "status": "selected incoming moments executable; source pressure datum, T_sh, complete join, I3/I4, global pressure remain open"},
        {"stage": "oscillatory_augmentation", "ready": False, "status": "rational rectangle existence seam executable; actual positive-order/background multi-band physical mode family and public by-beta/total velocity still absent"},
        {"stage": "mean_radial_corrections", "ready": False, "status": "stable disjoint damping exists but current one-band rank-two coverage remains 0_of_25; latest radial-aware extension still pending CI at audit"},
        {"stage": "finite_correction_cycle", "ready": False, "status": "blocked until actual multi-band family passes rank/unit/budget/spacetime/radial/quadratic guards"},
        {"stage": "candidate_artifact", "ready": False, "status": "schema reserved not instantiated"},
        {"stage": "independent_validation", "ready": False, "status": "local A4 audits exist; no frozen global Kokuno candidate for full-domain Agent-4 PDE gate"},
        {"stage": "report_and_export", "ready": False, "status": "blocked until frozen candidate artifact exists"},
    ]
    routing = {
        "new_fact": "Agent 2 #439 makes the corrected-reader finite rational auxiliary-torus center separation and strict r0 existence guards executable with exact Fraction arithmetic, while preserving all hidden rectangle/grid choices as repository-autonomous rather than recovered source data.",
        "interpretation": "Rational rectangle separation is no longer a software/geometry blocker. It does not create the missing physical oscillatory field; the shortest cross-lane blocker remains Agent 2's actual positive-order/background multi-band family with a genuine second covariance direction.",
        "shortest_next_closure": [
            "Agent 2: use the audited band schedule and autonomous rational separation witness to materialize at least two active slow bands of a provenance-labelled positive-order/background auxiliary-torus physical family; expose Q-scaled by-beta and total velocity(x,y,z,t) plus explicit coefficient units.",
            "Agent 3: consume exactly that family through genuine rank-two, unchanged unit/budget/spacetime/radial/quadratic guards; after latest radial-aware holdout CI is stable, keep +d_z sigma_1 fail-closed and only then materialize a correction and rerun the finite cycle.",
            "Agent 1: resolve the newly exposed source outer-pressure-datum binding, verify T_sh, complete inner-to-outer join, I3/I4 and global matched pressure under fixed/restricted forcing.",
            "Agent 4: independently audit the first actual multi-band physical family rather than repeating already-cleared source arithmetic; own the final full-domain normalized PDE gate and same-protocol ST006 comparison.",
            "Agent 5: instantiate deterministic save/load and Python/MATLAB export only after leading, oscillatory and correction stages coexist in one executable ancestry.",
        ],
        "do_not_do": [
            "do not call autonomous rational centers/r0 recovered source parameters",
            "do not treat exact center-separation arithmetic as a Navier-Stokes residual",
            "do not fabricate rank two with duplicate columns or unconstrained per-label amplitudes",
            "do not promote pending sibling CI to stable evidence",
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
            "agent3_stable": AGENT3_STABLE_RECEIPT,
            "agent3_latest_pending": AGENT3_LATEST_PENDING,
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
    hard_false = (
        "leading_ready",
        "source_rectangle_parameters_recovered",
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
        "openai_field_identified",
        "blowup_proved",
    )
    if any(states[name] for name in hard_false):
        raise ValueError("fail-closed state was promoted without evidence")
    if not states["source_rational_rectangle_separation_contract_executable"]:
        raise ValueError("Agent-2 rational separation receipt was lost")
    if not states["autonomous_rational_rectangle_witness_ready"]:
        raise ValueError("autonomous witness receipt was lost")

    a2 = payload["upstream"]["agent2"]
    if a2["standard_status"] != "success" or a2["dedicated_status"] != "success":
        raise ValueError("Agent-2 exact-head CI is not stable")
    if any(a2[key] for key in ("source_rectangle_centers_recovered", "source_rectangle_radius_r0_recovered", "source_color_count_recovered", "source_center_denominator_recovered")):
        raise ValueError("hidden source rectangle choices were mislabelled as recovered")

    a3_pending = payload["upstream"]["agent3_latest_pending"]
    if a3_pending["promoted_to_stable_receipt"]:
        raise ValueError("pending Agent-3 sibling was promoted")

    gates = payload["formal_gates"]
    if gates["held_out_normalized_momentum_max"] != 1.0e-3 or gates["held_out_normalized_momentum_l2"] != 1.0e-3:
        raise ValueError("momentum gate changed")
    if gates["held_out_divergence_max"] != 1.0e-5 or gates["held_out_divergence_l2"] != 1.0e-5:
        raise ValueError("divergence gate changed")
    if gates["changed_this_round"]:
        raise ValueError("formal gate mutation forbidden")

    contract = payload["candidate_artifact_contract"]
    if contract["instantiated"]:
        raise ValueError("candidate artifact cannot be instantiated at this frontier")
    if "residual-defined free forcing forbidden" not in contract["forcing_policy"]:
        raise ValueError("forcing policy weakened")


def write_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="artifacts/constrained/kokuno_agent5_rational_rectangle_routing_checkpoint_v27.json")
    args = parser.parse_args(argv)
    payload = write_checkpoint(args.output)
    print(payload["checkpoint_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
