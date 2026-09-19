"""Kokuno Agent-5 v37 routing checkpoint after public oscillatory black-box rejection.

The first deterministic public ``velocity_osc(x,y,z,t)`` now exists, but the
latest exact-head Agent-4 black-box audit rejects it under guards frozen before
Actions.  This module records that scientific REJECT as an integration receipt;
it never weakens the validator or promotes a candidate-only realization into
recovered Kokuno/OpenAI source data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-public-oscillatory-rejection-routing-checkpoint-v37"
TASK_ID = "KOKUNO-A5-PUBLIC-OSCILLATORY-REJECTION-ROUTING-037"

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
    "pr": 524,
    "head": "8cae8de187938b2d2af7f77d65c050e4d6fbc36e",
    "schema": "kokuno-agent5-coupled-c-phase-curl-routing-checkpoint-v36",
    "dedicated_run": 35411006970,
    "standard_run": 35411006890,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact": "kokuno-agent5-coupled-c-phase-curl-routing-checkpoint-v36",
    "artifact_id": 10574672159,
    "artifact_digest": "sha256:c7a75037bc812d8f25c67186750657623583f0dba953becf5c7af51a0bcfed7a",
    "checkpoint_sha256": "2a2a5b2aae0b5b489c978378748b13af6b7387d49fc5a7ee814e00f4034710e2",
}

AGENT1 = {
    "pr": 529,
    "head": "61b24c1ced2b9246df26e1ee4148392d2ff521d9",
    "dedicated_run": 35412761193,
    "standard_run": 35412762675,
    "dedicated_status": "success",
    "standard_status": "success",
    "source_choice_order_guard_executable": True,
    "fresh_coupled_C_pointwise_screen_passed": True,
    "fresh_coupled_C_pass_is_source_ordered_certificate": False,
    "source_B0_dependencies_machine_bound": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
}

AGENT2 = {
    "pr": 530,
    "head": "46c1699c7fc7ab2dc103d5cf430ead3a39bbbc1e",
    "dedicated_run": 35413533008,
    "standard_run": 35413532982,
    "dedicated_status": "success",
    "standard_status": "success",
    "provider": "openai_ns_reconstruction.kokuno_public_candidate_velocity:velocity_osc",
    "public_xyz_t_provider_executable": True,
    "save_load_executable": True,
    "by_sign_by_beta_total_access": True,
    "coordinate_contract": "dimensionless repository wave chart: R=hypot(x,y), theta=atan2(y,x), Z=z",
    "radial_support": [0.15, 1.35],
    "registered_time_interval": [0.25, 0.75],
    "actual_positive_order_background_bound": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_signed_auxiliary_rectangles_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "source_actual_partition_labels_instantiated": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}

AGENT3 = {
    "pr": 514,
    "head": "50d1b8ca965b35d5e0a7210c22e2302cd30c20e8",
    "dedicated_run": 35406792841,
    "standard_run": 35406792929,
    "dedicated_status": "success",
    "standard_status": "success",
    "autonomous_finite_head_factor_executable": True,
    "formal_theorem_missingWeight_replaced": False,
    "same_cycle_requested_stress_materialized": False,
    "candidate_numeric_finite_head_mean_debt_materialized": False,
    "real_candidate_defect_consumed": False,
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

# Latest fully-green Agent-4 exact-head receipt.  Small last-bit variation in
# independent FD/quadrature runs is expected across hosted runners, so the v37
# CI re-run checks the same frozen inequalities rather than bitwise float values.
AGENT4 = {
    "pr": 531,
    "head": "9edab1f1091d8d73c6ce1c5a8a262a08b33725ad",
    "parent_agent2_head": "46c1699c7fc7ab2dc103d5cf430ead3a39bbbc1e",
    "dedicated_run": 35413900071,
    "dedicated_status": "success",
    "standard_run": 35413900074,
    "standard_status": "success",
    "artifact": "k4osc-r0-6.9164e-01-r1-4.3300e-01-r2-1.0461e-01-rr-1.597e+00-rank-1.8748e-02-drift-2.061e-11-axz-7.6926e+05-vrms-3.1766e+05-hchg-9.102e-03-mut-2.376e-01-div-0-cov-0-rad-1-axial-0-overall-0",
    "artifact_id": 10575321691,
    "artifact_digest": "sha256:cc7032987c32ce38a32e9f7206f1e58052c15620ad9fd4ba0856b4ef93058883",
    "seed": 9173241,
    "fd4_steps": [0.02, 0.01, 0.005],
    "phase_resolutions": [32, 64, 128],
    "guards": A4_GUARDS,
    "metrics": {
        "finest_relative_divergence_rms": 0.10460795155383176,
        "finest_relative_divergence_max": 0.09811443995058577,
        "divergence_refinement_ratios": [1.5973332059304777, 4.139246627791994],
        "minimum_covariance_rank_ratio": 0.01874819346346177,
        "maximum_covariance_resolution_drift": 2.060561561211433e-11,
        "axis_near_absolute_max": 0.0,
        "radial_exterior_absolute_max": 0.0,
        "project_axial_exterior_absolute_max": 769255.3191046526,
        "heldout_velocity_rms": 317659.12117004656,
        "h_minus_10pct_relative_change": 1.2614067667595157,
        "h_plus_10pct_relative_change": 0.009101848910306524,
        "divergence_mutation_relative_rms": 0.23755538605232182,
        "duplicated_covariance_rank_ratio": 9.669379904398983e-19,
    },
    "local_divergence_passed": False,
    "covariance_rank_preflight_passed": False,
    "radial_axis_support_passed": True,
    "project_axial_support_passed": False,
    "local_curl_covariance_preflight_passed": False,
    "project_support_preflight_passed": False,
    "public_oscillatory_preflight_passed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
}

STATES = {
    "leading_ready": False,
    "public_provenance_labelled_xyz_t_oscillatory_velocity_ready": True,
    "public_oscillatory_independent_audit_assessed": True,
    "public_oscillatory_local_divergence_passed": False,
    "public_oscillatory_covariance_rank_passed": False,
    "public_oscillatory_project_support_passed": False,
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
        "status": "blocked_before_PA16",
        "provider": "Agent-1 coupled-C Appendix-B candidate",
        "next_required": [
            "machine-bind or independently upper-bound the three pre-C B0 dependencies",
            "evaluate source-ordered B0/T_sh without later-C laundering",
            "only then reconsider PA.16, I3/I4 and global matched pressure",
        ],
        "ready": False,
    },
    "oscillatory": {
        "status": "public_provider_exists_but_independent_preflight_rejected",
        "provider": AGENT2["provider"],
        "public_output": "velocity_osc(x,y,z,t)->[...,3]",
        "serialization": "SHA-bound JSON save/load",
        "blocking_failures": [
            "project axial support leak: |z|>2 absolute max 7.692553191046526e5 versus 1e-12 guard",
            "FD4 divergence: finest relative RMS 1.0460795155383176e-1 versus 2e-5 guard and max 9.811443995058577e-2 versus 1e-4 guard",
            "FD4 refinement: first ratio 1.5973332059304777 versus >=3 guard",
            "phase-mean covariance rank: minimum s_min/s_max 1.874819346346177e-2 versus >=2e-2 guard",
        ],
        "preserved_passes": [
            "axis-near exact zero",
            "radial exterior exact zero",
            "nontriviality",
            "h-parameter perturbation sensitivity",
            "divergence-detector mutation sensitivity",
            "covariance-resolution stability",
            "duplicated-column rank-collapse negative control",
        ],
        "ready": False,
    },
    "correction": {
        "status": "blocked_on_rejected_oscillatory_handoff_and_same_cycle_defect",
        "required_inputs": [
            "independently accepted public oscillatory covariance",
            "same-cycle theta/axial requestedStress",
            "candidate-specific finite-head mean debt",
        ],
        "unchanged_guards": [
            "H_ref y = DeltaC/epsilon",
            "bounded inverse",
            "unit/budget",
            "spacetime",
            "radial +d_z sigma_1",
            "quadratic damping",
            "joint gain",
        ],
        "ready": False,
    },
    "independent_validator": {
        "status": "oscillatory_component_rejected_before_composite_PDE_gate",
        "current_protocol": "public velocity_osc only; held-out FD4 + phase-mean covariance + support checks",
        "must_keep_frozen": A4_GUARDS,
        "final_gate": FORMAL_GATES,
        "ready_for_final_gate": False,
    },
}


def _canonical(payload: dict[str, Any]) -> bytes:
    body = {k: v for k, v in payload.items() if k != "checkpoint_sha256"}
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
            "agent2_public_provider": AGENT2,
            "agent3_latest": AGENT3,
            "agent4_public_provider_audit": AGENT4,
        },
        "states": STATES,
        "typed_handoffs": TYPED_HANDOFFS,
        "routing": {
            "new_fact": "The first deterministic public Kokuno-structured velocity_osc provider now exists and is independently assessed, but Agent 4 rejects it under the frozen local divergence/covariance/support preflight. It must not enter correction or the composite artifact yet.",
            "shortest_next_closure": [
                "Agent 2: add project-compatible axial compact support at vector-potential level with the required analytic D_z coefficient derivative; do not post-multiply velocity by a cutoff.",
                "Agent 2: make the public complete-curl realization resolve under unchanged Agent-4 FD4 steps .02/.01/.005 so divergence RMS/max meet 2e-5/1e-4 and both refinement ratios meet >=3; do not relax the validator.",
                "Agent 2: after support/divergence repair, improve the autonomous signed covariance realization enough to clear the unchanged 0.02 rank-ratio guard while preserving provenance and the same public API.",
                "Agent 4: rerun exactly the same black-box protocol and frozen guards on the repaired public provider; only an independent PASS can open oscillatory_ready.",
                "Agent 3: remain fail-closed until Agent 4 accepts the public oscillatory handoff; then materialize same-cycle requestedStress and candidate-specific finite-head mean debt before any DeltaC/epsilon inverse or finite cycle.",
                "Agent 1: continue the source-ordered pre-C B0/T_sh route; the coupled-C pointwise pass remains candidate-only and does not authorize PA.16.",
                "Agent 5: only after non-obstructed leading + independently accepted oscillatory + guarded correction coexist in one executable ancestry, instantiate the composite velocity/pressure/restricted-forcing artifact, save/load and Python/MATLAB smoke, then freeze for held-out NS validation.",
            ],
        },
    }
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def validate_checkpoint(payload: dict[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or payload.get("task_id") != TASK_ID:
        raise ValueError("checkpoint schema/task mismatch")
    if payload.get("checkpoint_sha256") != checkpoint_sha256(payload):
        raise ValueError("checkpoint sha256 mismatch")
    if payload.get("formal_gates") != FORMAL_GATES:
        raise ValueError("formal PDE gates changed")
    if payload.get("baseline_vs_kokuno", {}).get("st006") != ST006_BASELINE:
        raise ValueError("ST006 baseline changed")
    if payload["baseline_vs_kokuno"].get("kokuno_current_comparable_full_domain_receipt") is not None:
        raise ValueError("no comparable full-domain Kokuno receipt exists yet")
    if payload.get("states") != STATES:
        raise ValueError("routing states changed")
    if payload.get("typed_handoffs") != TYPED_HANDOFFS:
        raise ValueError("typed handoffs changed")
    expected_upstream = {
        "previous_agent5": PREVIOUS_AGENT5,
        "agent1_latest": AGENT1,
        "agent2_public_provider": AGENT2,
        "agent3_latest": AGENT3,
        "agent4_public_provider_audit": AGENT4,
    }
    if payload.get("upstream") != expected_upstream:
        raise ValueError("upstream receipt changed")

    a4 = payload["upstream"]["agent4_public_provider_audit"]
    m, g = a4["metrics"], a4["guards"]
    if a4["dedicated_status"] != "success" or a4["standard_status"] != "success":
        raise ValueError("Agent-4 exact-head software CI is not green")
    if a4["public_oscillatory_preflight_passed"]:
        raise ValueError("rejected oscillatory provider was promoted")
    if a4["local_divergence_passed"] or a4["covariance_rank_preflight_passed"] or a4["project_axial_support_passed"]:
        raise ValueError("failed Agent-4 scientific guard was rewritten")
    if not m["finest_relative_divergence_rms"] > g["finest_relative_divergence_rms"]:
        raise ValueError("divergence RMS rejection no longer holds")
    if not m["finest_relative_divergence_max"] > g["finest_relative_divergence_max"]:
        raise ValueError("divergence max rejection no longer holds")
    if not min(m["divergence_refinement_ratios"]) < g["minimum_divergence_refinement_ratio"]:
        raise ValueError("divergence refinement rejection no longer holds")
    if not m["minimum_covariance_rank_ratio"] < g["minimum_covariance_rank_ratio"]:
        raise ValueError("covariance rank rejection no longer holds")
    if not m["project_axial_exterior_absolute_max"] > g["project_axial_exterior_absolute_max"]:
        raise ValueError("axial support rejection no longer holds")
    if a4["formal_full_domain_pde_gate_assessed"] or a4["heldout_ns_residual_assessed"] or a4["pde_validated"]:
        raise ValueError("oscillatory preflight was laundered into full PDE validation")


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
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
