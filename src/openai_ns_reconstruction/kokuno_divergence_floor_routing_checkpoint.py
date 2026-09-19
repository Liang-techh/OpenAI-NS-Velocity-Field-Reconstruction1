"""Agent-5 deterministic routing checkpoint after the carrier-resolved A4 audit.

This module records one integration fact only: Agent 4's unchanged public-field
black-box audit closes support and physical covariance for Agent 2's
carrier-resolved oscillatory realization, while independent FD4 divergence still
fails the frozen RMS/refinement guards.  It deliberately does not instantiate a
composite candidate, run a correction cycle, or assess the full Navier--Stokes
PDE gate.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-divergence-floor-routing-checkpoint-v39"

AGENT1_HEAD = "3e19706475b62467984eeb5bc34425f69a39309e"
AGENT2_HEAD = "92852046e6e1ec0a5889b53f989df3ea5d79cb3b"
AGENT3_HEAD = "010286320cc407c7973df5dc28726a8b3089e1f0"
AGENT4_HEAD = "2a18a87db0ef8c6682f10fc5a6355c90b4168c7c"
AGENT4_WORKFLOW_RUN = 35419518424
AGENT4_STANDARD_RUN = 35419518366
AGENT4_ARTIFACT_ID = 10576976589
AGENT4_ARTIFACT_DIGEST = (
    "sha256:3ee35376dd1fa8cb1417639ebe40513b74cc6d8fe686d3fe380fa669ac6239cf"
)

ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622

FROZEN_GUARDS = {
    "finest_relative_divergence_rms": 2.0e-5,
    "finest_relative_divergence_max": 1.0e-4,
    "minimum_divergence_refinement_ratio": 3.0,
    "minimum_covariance_rank_ratio": 2.0e-2,
    "maximum_covariance_resolution_drift": 1.0e-6,
    "axis_near_absolute_max": 1.0e-14,
    "radial_exterior_absolute_max": 1.0e-14,
    "project_axial_exterior_absolute_max": 1.0e-12,
}

CURRENT_A4_METRICS = {
    "relative_divergence_rms_by_step": {
        "0.02": 8.479626402157787e-05,
        "0.01": 5.578704316484997e-05,
        "0.005": 5.451453500506105e-05,
    },
    "finest_relative_divergence_max": 6.930909589100082e-05,
    "relative_rms_refinement_ratios": [1.5199992545044203, 1.0233425481785872],
    "covariance_rank_ratio_by_phase_resolution": {
        "32": 0.09026569603652136,
        "64": 0.09026569603658607,
        "128": 0.09026569603651097,
    },
    "covariance_relative_drift_to_finest": [
        3.288617254673393e-13,
        1.0962057515577977e-13,
    ],
    "duplicated_first_column_rank_ratio": 7.289794049607572e-17,
    "axis_near_absolute_max": 0.0,
    "radial_exterior_absolute_max": 0.0,
    "project_axial_exterior_absolute_max": 0.0,
    "heldout_velocity_rms": 5126.150075408773,
    "h_minus_10pct_relative_change": 0.0036840930390566007,
    "h_plus_10pct_relative_change": 0.0036976476083083083,
    "external_divergence_mutation_relative_rms": 0.2500090453773952,
}

PRIOR_A4_543_METRICS = {
    "finest_relative_divergence_rms": 0.15167269397404395,
    "finest_relative_divergence_max": 0.19755153174056717,
    "minimum_covariance_rank_ratio": 0.01873858406014804,
}

FORMAL_GATES = {
    "held_out_normalized_momentum_max": 1.0e-3,
    "held_out_normalized_momentum_l2": 1.0e-3,
    "held_out_divergence_max": 1.0e-5,
    "held_out_divergence_l2": 1.0e-5,
}


def _canonical_sha256(payload: dict[str, Any]) -> str:
    unsigned = copy.deepcopy(payload)
    unsigned.pop("checkpoint_sha256", None)
    encoded = json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def build_checkpoint() -> dict[str, Any]:
    """Build the deterministic v39 routing checkpoint."""
    finest_rms = CURRENT_A4_METRICS["relative_divergence_rms_by_step"]["0.005"]
    finest_max = CURRENT_A4_METRICS["finest_relative_divergence_max"]
    min_refinement = min(CURRENT_A4_METRICS["relative_rms_refinement_ratios"])
    min_cov_rank = min(CURRENT_A4_METRICS["covariance_rank_ratio_by_phase_resolution"].values())
    max_cov_drift = max(CURRENT_A4_METRICS["covariance_relative_drift_to_finest"])

    checkpoint: dict[str, Any] = {
        "schema": SCHEMA,
        "round": 39,
        "purpose": (
            "freeze the independently audited carrier-resolved handoff, close support/covariance, "
            "and route the remaining oscillatory blocker to the public complete-curl divergence floor"
        ),
        "upstream": {
            "agent1_latest_sibling": {
                "pr": 550,
                "head": AGENT1_HEAD,
                "claim": "candidate_side_selected_PA10_M0_engineering_envelope_only",
                "source_M0_machine_bound": False,
                "source_B0_dependencies_machine_bound": False,
                "selected_pa16_handoff_allowed": False,
                "ci_used_as_dependency": False,
            },
            "agent2_carrier_resolved_candidate": {
                "pr": 551,
                "head": AGENT2_HEAD,
                "public_provider": (
                    "openai_ns_reconstruction.kokuno_public_carrier_resolved_velocity."
                    "KokunoCarrierResolvedCandidateOscillatoryVelocity.velocity_osc"
                ),
                "actual_positive_order_background_bound": False,
                "paper_exact": False,
            },
            "agent3_identity_wrapper": {
                "pr": 552,
                "head": AGENT3_HEAD,
                "identity_binding_ready_for_new_receipt": True,
                "correction_ingest_allowed": False,
            },
            "agent4_carrier_resolved_independent_audit": {
                "pr": 553,
                "head": AGENT4_HEAD,
                "audited_agent2_head": AGENT2_HEAD,
                "dedicated_workflow_run": AGENT4_WORKFLOW_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT4_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "artifact_id": AGENT4_ARTIFACT_ID,
                "artifact_zip_digest": AGENT4_ARTIFACT_DIGEST,
                "scientific_guards_changed_from_agent4_543": False,
                "metrics": copy.deepcopy(CURRENT_A4_METRICS),
                "scientific_verdict": "REJECT_divergence_only",
            },
        },
        "states": {
            "leading_ready": False,
            "public_provenance_labelled_xyz_t_oscillatory_velocity_ready": True,
            "public_oscillatory_project_support_passed": True,
            "public_oscillatory_covariance_rank_passed": True,
            "candidate_independent_second_covariance_column_ready": True,
            "public_oscillatory_local_divergence_max_passed": True,
            "public_oscillatory_local_divergence_rms_passed": False,
            "public_oscillatory_divergence_refinement_passed": False,
            "public_oscillatory_local_divergence_passed": False,
            "public_oscillatory_preflight_passed": False,
            "oscillatory_ready": False,
            "correction_ingest_allowed": False,
            "same_cycle_requested_stress_materialized": False,
            "candidate_numeric_finite_head_mean_debt_materialized": False,
            "correction_ready": False,
            "finite_correction_cycle_run": False,
            "candidate_artifact_instantiated": False,
            "velocity_export_ready": False,
            "formal_full_domain_pde_gate_assessed": False,
            "heldout_ns_residual_assessed": False,
            "pde_validated": False,
        },
        "closed_seams": {
            "registered_project_support": True,
            "physical_phase_mean_covariance_rank": True,
            "finest_relative_divergence_max_guard": True,
        },
        "remaining_blocker": {
            "name": "independent_public_complete_curl_divergence_floor",
            "finest_relative_rms": finest_rms,
            "rms_guard": FROZEN_GUARDS["finest_relative_divergence_rms"],
            "finest_relative_max": finest_max,
            "max_guard": FROZEN_GUARDS["finest_relative_divergence_max"],
            "minimum_refinement_ratio": min_refinement,
            "refinement_guard": FROZEN_GUARDS["minimum_divergence_refinement_ratio"],
            "plateau_observed": True,
            "plateau_steps": [0.01, 0.005],
            "routing": (
                "isolate coefficient/remainder/public-coordinate/cancellation seam in the existing "
                "complete-curl provider; preserve support, covariance response, public API, and every frozen A4 guard"
            ),
        },
        "typed_handoffs": {
            "leading": {
                "producer": "agent1",
                "ready": False,
                "next_required": (
                    "machine-bind or rigorously upper-bound the remaining source-ordered B0 dependencies, "
                    "then verify B0/T_sh before PA.16/global matched pressure"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "validator": "agent4",
                "public_api": "velocity_osc(x,y,z,t)->[...,3]",
                "candidate_head": AGENT2_HEAD,
                "validator_head": AGENT4_HEAD,
                "validator_artifact_id": AGENT4_ARTIFACT_ID,
                "validator_artifact_zip_digest": AGENT4_ARTIFACT_DIGEST,
                "support_passed": True,
                "covariance_rank_passed": True,
                "local_divergence_passed": False,
                "admitted": False,
            },
            "correction": {
                "producer": "agent3",
                "ingest_allowed": False,
                "covariance_input_structurally_available": True,
                "blocked_reason": "overall_independent_oscillatory_preflight_is_false",
                "required_actual_identity": {
                    "agent2_head": AGENT2_HEAD,
                    "agent4_head": AGENT4_HEAD,
                    "agent4_workflow_run": AGENT4_WORKFLOW_RUN,
                    "agent4_artifact_id": AGENT4_ARTIFACT_ID,
                    "agent4_artifact_zip_digest": AGENT4_ARTIFACT_DIGEST,
                },
                "after_admission": (
                    "materialize same-cycle theta/axial requestedStress and candidate-specific finite-head mean debt, "
                    "then run DeltaC/epsilon inverse plus existing budget/spacetime/radial/quadratic/joint-gain guards"
                ),
            },
            "final_candidate": {
                "producer": "agent5",
                "ready": False,
                "required_inputs": [
                    "non_obstructed_global_leading_velocity_and_pressure",
                    "agent4_admitted_public_oscillatory_velocity",
                    "guarded_materialized_finite_correction_cycle",
                ],
                "planned_api": ["velocity(x,y,z,t)", "pressure(x,y,z,t)", "restricted_forcing(x,y,z,t)"],
                "save_load_required": True,
                "python_matlab_smoke_required": True,
            },
        },
        "baseline_vs_kokuno": {
            "st006": {
                "normalized_momentum_sampled_max": ST006_MOMENTUM_MAX,
                "normalized_momentum_volume_l2": ST006_MOMENTUM_L2,
                "pde_validated": False,
            },
            "kokuno_current_comparable_full_domain_receipt": None,
            "component_preflight_only": {
                "prior_agent4_543": copy.deepcopy(PRIOR_A4_543_METRICS),
                "current_agent4_553": {
                    "finest_relative_divergence_rms": finest_rms,
                    "finest_relative_divergence_max": finest_max,
                    "minimum_covariance_rank_ratio": min_cov_rank,
                    "maximum_covariance_resolution_drift": max_cov_drift,
                },
                "finest_rms_improvement_factor": (
                    PRIOR_A4_543_METRICS["finest_relative_divergence_rms"] / finest_rms
                ),
                "finest_max_improvement_factor": (
                    PRIOR_A4_543_METRICS["finest_relative_divergence_max"] / finest_max
                ),
                "is_full_ns_comparison": False,
            },
        },
        "formal_gates": copy.deepcopy(FORMAL_GATES),
        "truth_boundary": {
            "kokuno_hidden_background_recovered": False,
            "paper_exact": False,
            "free_forcing_used_to_cancel_residual": False,
            "thresholds_relaxed_after_result": False,
            "agent4_replay_counted_as_full_pde_validation": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    validate_checkpoint(checkpoint)
    return checkpoint


def validate_checkpoint(checkpoint: dict[str, Any]) -> None:
    """Fail closed if a routing receipt overstates the audited scientific state."""
    if checkpoint.get("schema") != SCHEMA:
        raise ValueError("unexpected checkpoint schema")
    if checkpoint.get("checkpoint_sha256") != _canonical_sha256(checkpoint):
        raise ValueError("checkpoint SHA256 mismatch")

    states = checkpoint["states"]
    if states["leading_ready"]:
        raise ValueError("leading must remain closed")
    if not states["public_provenance_labelled_xyz_t_oscillatory_velocity_ready"]:
        raise ValueError("public oscillatory provider must remain recorded")
    if not states["public_oscillatory_project_support_passed"]:
        raise ValueError("independently accepted support seam was lost")
    if not states["public_oscillatory_covariance_rank_passed"]:
        raise ValueError("independently accepted covariance seam was lost")
    if not states["candidate_independent_second_covariance_column_ready"]:
        raise ValueError("candidate covariance second direction must remain recorded")
    if not states["public_oscillatory_local_divergence_max_passed"]:
        raise ValueError("finest divergence max now passes and must remain recorded")
    if states["public_oscillatory_local_divergence_rms_passed"]:
        raise ValueError("divergence RMS guard is still failed")
    if states["public_oscillatory_divergence_refinement_passed"]:
        raise ValueError("divergence refinement guard is still failed")
    if states["public_oscillatory_local_divergence_passed"]:
        raise ValueError("overall local divergence must remain failed")
    if states["public_oscillatory_preflight_passed"] or states["oscillatory_ready"]:
        raise ValueError("oscillatory admission cannot open before divergence passes")
    if states["correction_ingest_allowed"] or states["correction_ready"]:
        raise ValueError("correction must remain fail-closed")
    if states["candidate_artifact_instantiated"] or states["velocity_export_ready"]:
        raise ValueError("composite artifact/export cannot exist yet")
    if states["formal_full_domain_pde_gate_assessed"] or states["heldout_ns_residual_assessed"]:
        raise ValueError("full-domain PDE gate has not been assessed")
    if states["pde_validated"]:
        raise ValueError("PDE validation must remain false")

    a4 = checkpoint["upstream"]["agent4_carrier_resolved_independent_audit"]
    if a4["head"] != AGENT4_HEAD or a4["audited_agent2_head"] != AGENT2_HEAD:
        raise ValueError("actual A2/A4 audit identity mismatch")
    if a4["artifact_id"] != AGENT4_ARTIFACT_ID or a4["artifact_zip_digest"] != AGENT4_ARTIFACT_DIGEST:
        raise ValueError("Agent-4 artifact identity mismatch")
    if a4["scientific_guards_changed_from_agent4_543"]:
        raise ValueError("scientific guards must remain frozen")

    metrics = a4["metrics"]
    finest_rms = metrics["relative_divergence_rms_by_step"]["0.005"]
    finest_max = metrics["finest_relative_divergence_max"]
    if not finest_rms > FROZEN_GUARDS["finest_relative_divergence_rms"]:
        raise ValueError("divergence RMS blocker unexpectedly cleared")
    if not finest_max <= FROZEN_GUARDS["finest_relative_divergence_max"]:
        raise ValueError("finest divergence max should be a closed seam")
    if not min(metrics["relative_rms_refinement_ratios"]) < FROZEN_GUARDS["minimum_divergence_refinement_ratio"]:
        raise ValueError("divergence refinement blocker unexpectedly cleared")
    if not min(metrics["covariance_rank_ratio_by_phase_resolution"].values()) >= FROZEN_GUARDS["minimum_covariance_rank_ratio"]:
        raise ValueError("covariance rank must remain independently passed")
    if not max(metrics["covariance_relative_drift_to_finest"]) <= FROZEN_GUARDS["maximum_covariance_resolution_drift"]:
        raise ValueError("covariance resolution stability must remain passed")
    if metrics["axis_near_absolute_max"] > FROZEN_GUARDS["axis_near_absolute_max"]:
        raise ValueError("axis support guard failed")
    if metrics["radial_exterior_absolute_max"] > FROZEN_GUARDS["radial_exterior_absolute_max"]:
        raise ValueError("radial support guard failed")
    if metrics["project_axial_exterior_absolute_max"] > FROZEN_GUARDS["project_axial_exterior_absolute_max"]:
        raise ValueError("axial support guard failed")

    formal = checkpoint["formal_gates"]
    if formal != FORMAL_GATES:
        raise ValueError("formal PDE gates changed")
    truth = checkpoint["truth_boundary"]
    if truth["free_forcing_used_to_cancel_residual"]:
        raise ValueError("free forcing is forbidden")
    if truth["thresholds_relaxed_after_result"]:
        raise ValueError("post-result threshold relaxation is forbidden")
    if truth["paper_exact"] or truth["pde_validated"]:
        raise ValueError("truth boundary overstated")
    if checkpoint["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is not None:
        raise ValueError("no comparable Kokuno full-domain receipt exists yet")


def write_checkpoint(path: Path) -> dict[str, Any]:
    checkpoint = build_checkpoint()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n")
    return checkpoint


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    checkpoint = write_checkpoint(args.output)
    print("checkpoint_sha256=", checkpoint["checkpoint_sha256"])
    print("finest_relative_divergence_rms=", checkpoint["remaining_blocker"]["finest_relative_rms"])
    print("minimum_covariance_rank_ratio=", min(CURRENT_A4_METRICS["covariance_rank_ratio_by_phase_resolution"].values()))
    print("oscillatory_ready=", checkpoint["states"]["oscillatory_ready"])
    print("pde_validated=", checkpoint["states"]["pde_validated"])


if __name__ == "__main__":
    main()
