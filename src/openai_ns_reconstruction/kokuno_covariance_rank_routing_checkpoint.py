"""Agent-5 routing checkpoint after the repaired multi-band structural audit.

Agent 4 #463 independently re-ran the frozen #451 protocol on Agent 2 #461 and
cleared the deterministic aggregation/permutation seam without changing any local
or formal threshold.  Agent 3 #462 then sharpened the remaining blocker: its
predeclared manufactured two-band family has genuinely rank-two *raw velocity*
columns, but its declared phase/angle-mean covariance response vanishes, so the
mean-correction inverse still has no second direction.

This module records those receipts and the shortest closure.  It does not promote
manufactured data to recovered Kokuno source truth, materialize a correction, or
claim a Navier--Stokes residual reduction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-covariance-rank-routing-checkpoint-v29"
TASK_ID = "KOKUNO-A5-COVARIANCE-RANK-ROUTING-029"

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
    "active_kokuno_candidate_promoted_by_latest_commit": False,
}

AGENT1_RECEIPT = {
    "pr": 459,
    "head": "95fa1b815b2f8b04bf8ae64fee6e1166f76a1497",
    "standard_run": 35380452919,
    "standard_status": "success",
    "source_rescaled_core_seed_executable": True,
    "selected_real_axis_center_velocity_executable": True,
    "source_fixed_point_solved": False,
    "source_contraction_threshold_verified": False,
    "source_pressure_datum_applied_to_selected_global_path": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "consumed_in_executable_ancestry": False,
}

AGENT2_RECEIPT = {
    "pr": 461,
    "head": "aafb67af44c40794701ed0d66dbf0a08c1e45fb3",
    "dedicated_run": 35381196493,
    "standard_run": 35381196511,
    "dedicated_status": "success",
    "standard_status": "success",
    "deterministic_multiband_aggregation_fixed": True,
    "source_scheduled_supplied_mode_family_executable": True,
    "by_beta_and_by_band_columns_exposed": True,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_source_bound_xyz_t_velocity_ready": False,
    "consumed_in_executable_ancestry": True,
}

AGENT3_RECEIPT = {
    "pr": 462,
    "head": "fe8378b997d7c1375df9be4d44e489d6f2c1d2c5",
    "dedicated_run": 35381889846,
    "standard_run": 35381889973,
    "dedicated_status": "success",
    "standard_status": "success",
    "supplied_family_provenance": "repository_manufactured_source_compatible",
    "band_column_norms": [19.382103877741642, 10.48764389964904],
    "raw_velocity_column_singular_values": [19.56145103115169, 10.14919991177965],
    "raw_velocity_rank_ratio": 0.5188367619363721,
    "radial_cells": 5,
    "phase_mean_covariance_rank_two_cells": 0,
    "phase_mean_covariance_rank_two_fraction": 0.0,
    "local_supplied_family_covariance_rank_two": False,
    "actual_source_mode_family_bound": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "finite_correction_cycle_rerun_allowed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

AGENT4_RECEIPT = {
    "pr": 463,
    "head": "89849ade6928c27c955794e55c61c397770021d0",
    "dedicated_run": 35382296059,
    "standard_run": 35382296045,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact_id": 10562362645,
    "artifact_digest": "sha256:71c91a85891fc52e464ad606665958eaa121f6c60f8c5aad538cbd426dcdc639",
    "simultaneous_beta_permutation_points": 36,
    "simultaneous_beta_permutation_failures": 0,
    "simultaneous_beta_permutation_relative_max": 0.0,
    "bitwise_identical_under_semantic_permutation": True,
    "finest_cartesian_curl_relative_rms": 3.7898223181764806e-6,
    "finest_normalized_divergence_rms": 5.450640030559649e-6,
    "finest_normalized_divergence_point_max": 2.1573976121470144e-5,
    "shared_first_band_schedule_mutation_relative_rms": 0.554659182591723,
    "reversed_label_schedule_mutation_relative_rms": 1.1243576802006303,
    "local_structural_preflight_passed": True,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": True,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 452,
    "head": "f59164ae39ef16abef78acd392b72de501a5e3b1",
    "schema": "kokuno-agent5-multiband-audit-routing-checkpoint-v28",
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
        "source_multiband_public_contract_structural_preflight_passed": True,
        "source_multiband_independent_cartesian_audit_passed": True,
        "supplied_multiband_raw_velocity_rank_two": True,
        "supplied_multiband_phase_mean_covariance_rank_two": False,
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
        {"stage": "source_profile_ingest", "ready": False, "status": "displayed constants/schedules plus autonomous support geometry are executable; actual positive-order/background and mode data remain unbound"},
        {"stage": "leading_candidate", "ready": False, "status": "source-rescaled Appendix-B core seed is executable, but fixed point/contraction threshold/global pressure/global leading path remain open"},
        {"stage": "oscillatory_augmentation", "ready": False, "status": "deterministic multi-band public aggregation and independent local curl/divergence audit now pass; supplied modes are still manufactured rather than source-bound"},
        {"stage": "mean_radial_corrections", "ready": False, "status": "manufactured two-band raw velocity columns are rank two but phase/angle-mean covariance response is zero on 5/5 cells, so no genuine second correction direction exists"},
        {"stage": "finite_correction_cycle", "ready": False, "status": "blocked until a source-bound/source-motivated family survives phase-mean covariance rank plus unchanged unit/budget/spacetime/radial/quadratic guards"},
        {"stage": "candidate_artifact", "ready": False, "status": "schema reserved, not instantiated"},
        {"stage": "independent_validation", "ready": False, "status": "Agent 4 cleared the local multiband structural seam; no complete global candidate exists for the formal held-out PDE gate"},
        {"stage": "report_and_export", "ready": False, "status": "blocked until one frozen leading+oscillatory+corrected artifact exists"},
    ]
    routing = {
        "new_fact": "Agent 4 #463 closes the deterministic aggregation/permutation blocker under the unchanged #451 audit, but Agent 3 #462 shows that robust raw velocity rank two is insufficient: the declared phase/angle mean annihilates the covariance response on all 5 radial cells.",
        "classification": "structural_multiband_pass_then_phase_mean_covariance_rank_reject",
        "shortest_next_closure": [
            "Agent 2: stop revisiting deterministic aggregation and generic local curl numerics; bind provenance-labelled positive-order/background plus a source-bound or explicitly source-motivated multi-band auxiliary-torus family with public Q-scaled by-beta/total velocity(x,y,z,t), explicit coefficient units, and nonzero phase-mean covariance response.",
            "Agent 4: independently audit that first source-bound/source-motivated physical family, including its phase-mean covariance response rather than only raw velocity-column rank.",
            "Agent 3: consume the exact audited family; require genuine covariance rank two plus unchanged unit/budget/spacetime/radial d_z sigma_1/quadratic/joint-gain guards before materializing delta_u or rerunning the finite correction cycle.",
            "Agent 1: propagate the source-rescaled core seed through the actual fixed-point/contraction/global join path and complete matched pressure under fixed/restricted forcing.",
            "Agent 5: instantiate deterministic save/load plus Python/MATLAB smoke only after leading, oscillatory and correction lanes coexist in one executable ancestry; then freeze the artifact for Agent 4's formal held-out PDE gate and same-protocol ST006 comparison.",
        ],
        "do_not_do": [
            "do not infer covariance-rank readiness from multiple bands or raw velocity rank alone",
            "do not tune manufactured phases/polarizations after seeing the zero covariance response merely to force a pass",
            "do not promote local divergence 5.45e-6 to a Navier-Stokes momentum result",
            "do not call autonomous/source-compatible hidden choices recovered Kokuno data",
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
    if not states["source_multiband_public_contract_structural_preflight_passed"]:
        raise ValueError("Agent-4 structural pass receipt was lost")
    if not states["source_multiband_independent_cartesian_audit_passed"]:
        raise ValueError("Agent-4 independent local audit pass receipt was lost")
    if not states["supplied_multiband_raw_velocity_rank_two"]:
        raise ValueError("Agent-3 raw velocity rank-two receipt was lost")
    if states["supplied_multiband_phase_mean_covariance_rank_two"]:
        raise ValueError("phase-mean covariance rank was promoted despite the Agent-3 rejection")

    for name in (
        "leading_ready",
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
    ):
        if states[name]:
            raise ValueError(f"fail-closed state promoted without evidence: {name}")

    a4 = payload["upstream"]["agent4"]
    if a4["simultaneous_beta_permutation_failures"] != 0:
        raise ValueError("deterministic aggregation seam is not closed")
    if not a4["local_structural_preflight_passed"]:
        raise ValueError("Agent-4 local structural preflight must remain passed")

    a3 = payload["upstream"]["agent3"]
    if a3["phase_mean_covariance_rank_two_cells"] != 0 or a3["radial_cells"] != 5:
        raise ValueError("Agent-3 covariance-rank rejection receipt changed")
    if a3["genuinely_independent_second_covariance_column_ready"]:
        raise ValueError("second covariance direction cannot be promoted")

    gates = payload["formal_gates"]
    if gates != FORMAL_GATES:
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
