"""Kokuno Agent-5 v38 routing after the independently audited axial-support repair.

Agent 2 now supplies a public provenance-labelled ``velocity_osc(x,y,z,t)`` whose
registered axial support is independently confirmed by Agent 4.  The same frozen
black-box protocol still rejects the realization on divergence/resolution and
phase-mean covariance rank.  This checkpoint records that narrower frontier and
keeps correction/composite PDE validation fail-closed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-post-axial-oscillatory-routing-checkpoint-v38"
TASK_ID = "KOKUNO-A5-POST-AXIAL-OSCILLATORY-ROUTING-038"

FORMAL_GATES = {
    "held_out_normalized_momentum_max": 1.0e-3,
    "held_out_normalized_momentum_l2": 1.0e-3,
    "held_out_divergence_max": 1.0e-5,
    "held_out_divergence_l2": 1.0e-5,
}

ST006_BASELINE = {
    "candidate": "ST006",
    "momentum_sampled_max": 0.1082289305112118,
    "momentum_volume_l2": 0.10758432876230622,
    "same_protocol_comparison_available_for_kokuno": False,
    "pde_validated": False,
}

PREVIOUS_AGENT5 = {
    "pr": 532,
    "head": "c65e7265a14e9d7e6e394f8b3cacfef5b005d9d0",
    "schema": "kokuno-agent5-public-oscillatory-rejection-routing-checkpoint-v37",
    "dedicated_run": 35414204283,
    "standard_run": 35414204396,
    "dedicated_status": "success",
    "standard_status": "success",
}

AGENT1 = {
    "pr": 537,
    "head": "7375fc961bbbc3ee469936104fb1eea3e742a979",
    "dedicated_run": 35415781315,
    "standard_run": 35415790337,
    "dedicated_status": "success",
    "standard_status": "success",
    "selected_lambda": 2.503192875997364e27,
    "x0": 1.5979591658138821e-27,
    "l0": 66.40155058292201,
    "l0_machine_bound": True,
    "remaining_source_b0_dependencies": [
        "||log phi_* + log Phi(4,.)||_{C^0_eta}",
        "M_0",
    ],
    "source_B0_dependencies_machine_bound": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
}

AGENT2 = {
    "pr": 539,
    "head": "985d1c3fc43081d2d87feb1f0c27643cc4ebe75b",
    "dedicated_run": 35416211440,
    "standard_run": 35416211479,
    "dedicated_status": "success",
    "standard_status": "success",
    "provider": "openai_ns_reconstruction.kokuno_public_candidate_velocity:velocity_osc",
    "serialization_schema": "v2",
    "public_xyz_t_provider_executable": True,
    "vector_potential_axial_support_localization": True,
    "analytic_D_z_coefficient_derivative_supplied": True,
    "registered_axial_support": [-2.0, 2.0],
    "support_repair_claims_divergence_fixed": False,
    "support_repair_claims_covariance_rank_fixed": False,
    "actual_positive_order_background_bound": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}

AGENT3 = {
    "pr": 541,
    "head": "e59cdd518460983bea2073f59524d4082769d94b",
    "dedicated_run": 35416312767,
    "standard_run": 35416312769,
    "dedicated_status": "success",
    "standard_status": "success",
    "correction_admission_fail_closed": True,
    "embedded_protocol_origin_parent_agent2_head": "46c1699c7fc7ab2dc103d5cf430ead3a39bbbc1e",
    "actual_audited_candidate_identity_bound_in_admission": False,
    "actual_agent4_artifact_identity_bound_in_admission": False,
    "correction_ingest_allowed": False,
    "same_cycle_requested_stress_materialized": False,
    "candidate_numeric_finite_head_mean_debt_materialized": False,
    "finite_correction_cycle_rerun_allowed": False,
}

A4_GUARDS = {
    "finest_relative_divergence_rms": 2.0e-5,
    "finest_relative_divergence_max": 1.0e-4,
    "minimum_divergence_refinement_ratio": 3.0,
    "minimum_covariance_rank_ratio": 2.0e-2,
    "maximum_covariance_resolution_drift": 1.0e-6,
    "axis_near_absolute_max": 1.0e-14,
    "radial_exterior_absolute_max": 1.0e-14,
    "project_axial_exterior_absolute_max": 1.0e-12,
    "minimum_nontrivial_velocity_rms": 1.0e-8,
    "minimum_parameter_perturbation_relative_change": 1.0e-4,
    "minimum_divergence_mutation_relative_rms": 5.0e-2,
    "maximum_duplicated_covariance_rank_ratio": 1.0e-8,
}

AGENT4 = {
    "pr": 543,
    "head": "9498209376fd9dc5abf228248967958759e9f64a",
    "actual_audited_agent2_head": "985d1c3fc43081d2d87feb1f0c27643cc4ebe75b",
    "protocol_origin_parent_agent2_head": "46c1699c7fc7ab2dc103d5cf430ead3a39bbbc1e",
    "dedicated_run": 35416733758,
    "standard_run": 35416733748,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact_id": 10576031714,
    "artifact_digest": "sha256:8d7dbdc7654fc74820fb6916c5726984b42116055df146936b76607e6d1a0981",
    "seed": 9173241,
    "fd4_steps": [0.02, 0.01, 0.005],
    "phase_resolutions": [32, 64, 128],
    "guards": A4_GUARDS,
    "metrics": {
        "relative_divergence_rms_by_step": [
            0.7274175799848378,
            0.5291197810511111,
            0.15167269397404395,
        ],
        "finest_relative_divergence_max": 0.19755153174056717,
        "divergence_refinement_ratios": [1.3747692035625707, 3.488563215878927],
        "minimum_covariance_rank_ratio": 0.01873858406014804,
        "maximum_covariance_resolution_drift": 3.215333113729077e-11,
        "axis_near_absolute_max": 0.0,
        "radial_exterior_absolute_max": 0.0,
        "project_axial_exterior_absolute_max": 0.0,
        "heldout_velocity_rms": 204022.57816331778,
        "h_minus_10pct_relative_change": 1.3356399641375716,
        "h_plus_10pct_relative_change": 0.00815993371683721,
        "divergence_mutation_relative_rms": 0.2445506270671321,
        "duplicated_covariance_rank_ratio": 7.923184306254537e-18,
    },
    "local_divergence_passed": False,
    "covariance_rank_preflight_passed": False,
    "radial_axis_support_passed": True,
    "project_axial_support_passed": True,
    "project_support_preflight_passed": True,
    "public_oscillatory_preflight_passed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
}

STATES = {
    "leading_ready": False,
    "public_provenance_labelled_xyz_t_oscillatory_velocity_ready": True,
    "public_oscillatory_independent_audit_assessed": True,
    "public_oscillatory_project_support_passed": True,
    "public_oscillatory_local_divergence_passed": False,
    "public_oscillatory_covariance_rank_passed": False,
    "oscillatory_ready": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "same_cycle_requested_stress_materialized": False,
    "candidate_numeric_finite_head_mean_debt_materialized": False,
    "correction_ingest_allowed": False,
    "correction_ready": False,
    "finite_correction_cycle_run": False,
    "candidate_artifact_instantiated": False,
    "velocity_export_ready": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}

TYPED_HANDOFFS = {
    "leading": {
        "status": "blocked_before_PA16_after_L0_binding",
        "provider": "Agent-1 source-ordered PA.10 dependency path",
        "next_required": [
            "machine-bind or independently upper-bound ||log phi_* + log Phi(4,.)||_{C^0_eta}",
            "machine-bind or independently upper-bound M_0",
            "evaluate source-ordered B_0/T_sh; only then reconsider PA.16 and global pressure matching",
        ],
        "ready": False,
    },
    "oscillatory": {
        "status": "public_provider_support_passed_but_divergence_and_covariance_rejected",
        "provider": AGENT2["provider"],
        "public_output": "velocity_osc(x,y,z,t)->[...,3]",
        "serialization": "SHA-bound JSON save/load, schema v2",
        "closed_seam": "registered axial support is independently confirmed: exterior max 0 under the unchanged Agent-4 protocol",
        "blocking_failures": [
            "fixed-ladder FD4 divergence: finest relative RMS 1.5167269397404395e-1 versus 2e-5 and max 1.9755153174056717e-1 versus 1e-4",
            "fixed-ladder FD4 refinement: first ratio 1.3747692035625707 versus >=3",
            "phase-mean covariance rank: minimum s_min/s_max 1.873858406014804e-2 versus >=2e-2",
        ],
        "next_required": [
            "preserve vector-potential support localization and public API while exposing/quantifying the effective carrier wavelength and derivative scales seen by the frozen FD4 ladder",
            "repair or reparameterize the public complete-curl realization so the unchanged .02/.01/.005 black-box divergence protocol resolves it and passes the frozen guards",
            "strengthen the signed physical covariance direction above the unchanged 0.02 rank gate without provenance laundering",
            "rerun the byte-identical Agent-4 protocol; only its PASS may set oscillatory_ready",
        ],
        "ready": False,
    },
    "correction": {
        "status": "blocked_on_rejected_oscillatory_handoff_and_receipt_identity_upgrade",
        "provenance_fix_required": "bind the actual audited Agent-2 head plus Agent-4 workflow/artifact identity; the embedded legacy parent_agent2_head is protocol-origin metadata only",
        "required_actual_identity": {
            "agent2_head": AGENT4["actual_audited_agent2_head"],
            "agent4_head": AGENT4["head"],
            "agent4_artifact_id": AGENT4["artifact_id"],
            "agent4_artifact_digest": AGENT4["artifact_digest"],
        },
        "required_scientific_inputs_after_A4_pass": [
            "same-cycle theta/axial requestedStress",
            "candidate-specific finite-head mean debt",
            "H_ref y = DeltaC/epsilon with existing bounded inverse/budget/spacetime/radial/quadratic/joint-gain guards",
        ],
        "ready": False,
    },
    "independent_validator": {
        "status": "support_passed_component_rejected_before_composite_PDE_gate",
        "current_protocol": "public velocity_osc only; held-out Cartesian FD4 + phase-mean covariance + support checks",
        "must_keep_frozen": A4_GUARDS,
        "final_gate": FORMAL_GATES,
        "ready_for_final_gate": False,
    },
}


def _canonical(payload: dict[str, Any]) -> bytes:
    body = {key: value for key, value in payload.items() if key != "checkpoint_sha256"}
    return json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def checkpoint_sha256(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical(payload)).hexdigest()


def build_checkpoint() -> dict[str, Any]:
    payload = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "formal_gates": FORMAL_GATES,
        "baseline_vs_kokuno": {
            "st006": ST006_BASELINE,
            "kokuno_current_comparable_full_domain_receipt": None,
        },
        "upstream": {
            "previous_agent5": PREVIOUS_AGENT5,
            "agent1_latest": AGENT1,
            "agent2_support_repair": AGENT2,
            "agent3_admission": AGENT3,
            "agent4_frozen_rerun": AGENT4,
        },
        "states": STATES,
        "typed_handoffs": TYPED_HANDOFFS,
        "routing": {
            "new_fact": "Agent 4 independently confirms Agent 2's vector-potential axial-support repair under the byte-identical frozen protocol; support is no longer the oscillatory blocker, while divergence/resolution and covariance rank remain rejected.",
            "shortest_next_closure": [
                "Agent 2: diagnose and repair fixed-ladder resolvability/divergence first, preserving support localization and the public provider; do not alter Agent-4 guards.",
                "Agent 2: then strengthen the signed physical covariance direction above the unchanged 0.02 rank gate.",
                "Agent 4: rerun the byte-identical public black-box protocol after each real candidate change; only an independent PASS opens oscillatory_ready.",
                "Agent 3: keep correction closed and upgrade admission provenance to bind the actual Agent-2 candidate head and Agent-4 artifact/workflow identity before any future PASS can unlock same-cycle defect materialization.",
                "Agent 1: bind the two remaining source-ordered B0 dependencies after L0; PA.16/global pressure remain closed.",
                "Agent 5: instantiate the composite artifact/export only after leading + independently accepted oscillatory + guarded correction coexist in one executable ancestry, then hand the frozen composite to Agent 4 for the unchanged full NS gate.",
            ],
        },
    }
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def validate_checkpoint(payload: dict[str, Any]) -> None:
    expected = build_checkpoint()
    if payload.get("checkpoint_sha256") != checkpoint_sha256(payload):
        raise ValueError("checkpoint sha256 mismatch")
    if payload.get("schema") != SCHEMA or payload.get("task_id") != TASK_ID:
        raise ValueError("checkpoint schema/task changed")
    if payload.get("formal_gates") != FORMAL_GATES:
        raise ValueError("formal gates changed")
    if payload.get("upstream") != expected["upstream"]:
        raise ValueError("upstream receipt changed")
    if payload.get("states") != STATES:
        raise ValueError("routing states changed")
    if payload.get("typed_handoffs") != TYPED_HANDOFFS:
        raise ValueError("typed handoffs changed")
    if payload.get("baseline_vs_kokuno") != expected["baseline_vs_kokuno"]:
        raise ValueError("baseline/comparison semantics changed")

    a4 = payload["upstream"]["agent4_frozen_rerun"]
    metrics = a4["metrics"]
    guards = a4["guards"]
    if a4["actual_audited_agent2_head"] == a4["protocol_origin_parent_agent2_head"]:
        raise ValueError("actual candidate identity collapsed into protocol-origin metadata")
    if not a4["project_support_preflight_passed"] or not a4["project_axial_support_passed"]:
        raise ValueError("independently confirmed support pass lost")
    if metrics["project_axial_exterior_absolute_max"] > guards["project_axial_exterior_absolute_max"]:
        raise ValueError("axial support no longer passes frozen guard")
    if metrics["relative_divergence_rms_by_step"][-1] <= guards["finest_relative_divergence_rms"]:
        raise ValueError("frozen divergence rejection unexpectedly changed")
    if metrics["finest_relative_divergence_max"] <= guards["finest_relative_divergence_max"]:
        raise ValueError("frozen divergence max rejection unexpectedly changed")
    if min(metrics["divergence_refinement_ratios"]) >= guards["minimum_divergence_refinement_ratio"]:
        raise ValueError("frozen refinement rejection unexpectedly changed")
    if metrics["minimum_covariance_rank_ratio"] >= guards["minimum_covariance_rank_ratio"]:
        raise ValueError("frozen covariance rejection unexpectedly changed")
    if a4["public_oscillatory_preflight_passed"] or a4["pde_validated"]:
        raise ValueError("rejected component promoted")
    if payload["states"]["correction_ingest_allowed"] or payload["states"]["pde_validated"]:
        raise ValueError("downstream gate promoted")


def write_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = write_checkpoint(args.output)
    print(json.dumps({"output": str(args.output), "checkpoint_sha256": payload["checkpoint_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
