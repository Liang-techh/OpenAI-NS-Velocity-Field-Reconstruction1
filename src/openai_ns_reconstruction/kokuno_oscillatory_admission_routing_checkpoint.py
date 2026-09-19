"""Agent-5 deterministic routing checkpoint after independent public-z oscillatory PASS.

This module records a single integration increment.  Agent 4 independently
re-ran the frozen public black-box protocol on Agent 2's public-z pullback
``velocity_osc(x,y,z,t)`` and every registered local oscillatory preflight guard
passed.  The oscillatory component may therefore be admitted as an independently
validated *component* handoff.  This does not create a global leading profile,
materialize the mean/radial correction, assess the full Navier--Stokes residual,
or promote the candidate to paper-exact/source-exact status.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-oscillatory-admission-routing-checkpoint-v40"

AGENT1_PR = 560
AGENT1_HEAD = "7a549c9b063a1ca36145811da56070faf0e9f10f"
AGENT2_PR = 561
AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
AGENT3_PR = 562
AGENT3_HEAD = "81986135c2b54866def945f2ff6343438407a95c"
AGENT4_PR = 563
AGENT4_HEAD = "9b0f86012c53fa8e32a19f766dbc150931870425"
AGENT4_WORKFLOW_RUN = 35422203621
AGENT4_STANDARD_RUN = 35422203626
AGENT4_ARTIFACT_ID = 10577644750
AGENT4_ARTIFACT_DIGEST = (
    "sha256:7205e041c31d2a68b345235bbd7b0f3316a50ec6c41fa3012c107789c6a07029"
)

ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622

FROZEN_LOCAL_GUARDS = {
    "finest_relative_divergence_rms": 2.0e-5,
    "finest_relative_divergence_max": 1.0e-4,
    "minimum_divergence_refinement_ratio": 3.0,
    "minimum_covariance_rank_ratio": 2.0e-2,
    "maximum_covariance_resolution_drift": 1.0e-6,
    "axis_near_absolute_max": 1.0e-14,
    "radial_exterior_absolute_max": 1.0e-14,
    "project_axial_exterior_absolute_max": 1.0e-12,
}

AGENT4_PASS_METRICS = {
    "relative_divergence_rms_by_step": {
        "0.02": 4.462053597423048e-05,
        "0.01": 2.8197350552150228e-06,
        "0.005": 1.767246633658725e-07,
    },
    "finest_relative_divergence_max": 2.0538436904463344e-07,
    "relative_rms_refinement_ratios": [15.824371829440508, 15.955526532124916],
    "covariance_rank_ratio_by_phase_resolution": {
        "32": 0.09276278453864097,
        "64": 0.09276278453878781,
        "128": 0.09276278453868625,
    },
    "covariance_relative_drift_to_finest": [
        1.1554379812996678e-13,
        2.0012772886058432e-13,
    ],
    "duplicated_first_column_rank_ratio": 7.034620062075487e-17,
    "axis_near_absolute_max": 0.0,
    "radial_exterior_absolute_max": 0.0,
    "project_axial_exterior_absolute_max": 0.0,
    "heldout_velocity_rms": 5087.870320629883,
    "h_minus_10pct_relative_change": 0.0027257496677792563,
    "h_plus_10pct_relative_change": 0.002731779028717074,
    "external_divergence_mutation_relative_rms": 0.25000004186441394,
}

PRIOR_AGENT4_553_METRICS = {
    "finest_relative_divergence_rms": 5.451453500506105e-05,
    "finest_relative_divergence_max": 6.930909589100082e-05,
    "minimum_covariance_rank_ratio": 0.09026569603651097,
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
    encoded = json.dumps(
        unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _derived_passes() -> dict[str, bool]:
    metrics = AGENT4_PASS_METRICS
    finest_rms = metrics["relative_divergence_rms_by_step"]["0.005"]
    finest_max = metrics["finest_relative_divergence_max"]
    min_refinement = min(metrics["relative_rms_refinement_ratios"])
    min_rank = min(metrics["covariance_rank_ratio_by_phase_resolution"].values())
    max_drift = max(metrics["covariance_relative_drift_to_finest"])
    return {
        "divergence_rms": finest_rms <= FROZEN_LOCAL_GUARDS["finest_relative_divergence_rms"],
        "divergence_max": finest_max <= FROZEN_LOCAL_GUARDS["finest_relative_divergence_max"],
        "divergence_refinement": min_refinement
        >= FROZEN_LOCAL_GUARDS["minimum_divergence_refinement_ratio"],
        "covariance_rank": min_rank >= FROZEN_LOCAL_GUARDS["minimum_covariance_rank_ratio"],
        "covariance_resolution": max_drift
        <= FROZEN_LOCAL_GUARDS["maximum_covariance_resolution_drift"],
        "support": (
            metrics["axis_near_absolute_max"] <= FROZEN_LOCAL_GUARDS["axis_near_absolute_max"]
            and metrics["radial_exterior_absolute_max"]
            <= FROZEN_LOCAL_GUARDS["radial_exterior_absolute_max"]
            and metrics["project_axial_exterior_absolute_max"]
            <= FROZEN_LOCAL_GUARDS["project_axial_exterior_absolute_max"]
        ),
    }


def build_checkpoint() -> dict[str, Any]:
    """Build the deterministic v40 integration checkpoint."""
    passes = _derived_passes()
    if not all(passes.values()):
        raise ValueError("frozen Agent-4 public-z receipt no longer satisfies its frozen guards")

    checkpoint: dict[str, Any] = {
        "schema": SCHEMA,
        "round": 40,
        "purpose": (
            "promote only the independently validated public oscillatory component to an admitted typed handoff, "
            "then route correction identity binding and leading closure without fabricating a composite PDE result"
        ),
        "upstream": {
            "agent1_latest_sibling": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "claim": "conditional_rescaled_core_fixed_point_radius_formula_bridge_only",
                "source_operator_constant_M_machine_bound": False,
                "source_operator_constant_K_machine_bound": False,
                "selected_pa16_handoff_allowed": False,
                "global_leading_profile_reconstructed": False,
                "ci_used_as_dependency": False,
            },
            "agent2_public_z_pullback_candidate": {
                "pr": AGENT2_PR,
                "head": AGENT2_HEAD,
                "public_provider": (
                    "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc"
                ),
                "coordinate_contract": "bandwise Z_beta=epsilon_beta*z so public partial_z equals source D_z",
                "repository_autonomous_source_compatible": True,
                "actual_positive_order_background_bound": False,
                "paper_exact": False,
            },
            "agent3_latest_identity_gate": {
                "pr": AGENT3_PR,
                "head": AGENT3_HEAD,
                "pins_agent4_553_reject_receipt": True,
                "pins_agent4_563_pass_receipt": False,
                "correction_ingest_allowed": False,
                "next_required": "pin exact Agent-4 #563 PASS execution identity before correction ingest",
            },
            "agent4_public_z_independent_audit": {
                "pr": AGENT4_PR,
                "head": AGENT4_HEAD,
                "audited_agent2_head": AGENT2_HEAD,
                "dedicated_workflow_run": AGENT4_WORKFLOW_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT4_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "artifact_id": AGENT4_ARTIFACT_ID,
                "artifact_zip_digest": AGENT4_ARTIFACT_DIGEST,
                "scientific_guards_changed_from_agent4_543": False,
                "metrics": copy.deepcopy(AGENT4_PASS_METRICS),
                "derived_guard_passes": passes,
                "scientific_verdict": "PASS_local_oscillatory_preflight_only",
            },
        },
        "states": {
            "leading_ready": False,
            "public_provenance_labelled_xyz_t_oscillatory_velocity_ready": True,
            "public_oscillatory_project_support_passed": True,
            "public_oscillatory_covariance_rank_passed": True,
            "candidate_independent_second_covariance_column_ready": True,
            "public_oscillatory_local_divergence_passed": True,
            "public_oscillatory_preflight_passed": True,
            "oscillatory_ready": True,
            "correction_receipt_identity_pinned_to_current_pass": False,
            "correction_ingest_allowed": False,
            "same_cycle_requested_stress_materialized": False,
            "candidate_numeric_finite_head_mean_debt_materialized": False,
            "real_candidate_defect_consumed": False,
            "correction_ready": False,
            "finite_correction_cycle_run": False,
            "candidate_artifact_instantiated": False,
            "velocity_export_ready": False,
            "formal_full_domain_pde_gate_assessed": False,
            "heldout_ns_residual_assessed": False,
            "pde_validated": False,
        },
        "typed_handoffs": {
            "leading": {
                "producer": "agent1",
                "ready": False,
                "next_required": (
                    "machine-bind source operator constants M,K and remaining source-ordered B0/T_sh prerequisites; "
                    "then complete PA.16/global matched pressure"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "validator": "agent4",
                "public_api": "velocity_osc(x,y,z,t)->[...,3]",
                "candidate_head": AGENT2_HEAD,
                "validator_head": AGENT4_HEAD,
                "validator_workflow_run": AGENT4_WORKFLOW_RUN,
                "validator_artifact_id": AGENT4_ARTIFACT_ID,
                "validator_artifact_zip_digest": AGENT4_ARTIFACT_DIGEST,
                "support_passed": True,
                "covariance_rank_passed": True,
                "local_divergence_passed": True,
                "admitted": True,
                "scope": "oscillatory_component_preflight_only",
            },
            "correction": {
                "producer": "agent3",
                "ingest_allowed": False,
                "blocked_reason": "agent3_has_not_yet_pinned_exact_agent4_563_PASS_execution_identity",
                "required_actual_identity": {
                    "agent2_head": AGENT2_HEAD,
                    "agent4_head": AGENT4_HEAD,
                    "agent4_workflow_run": AGENT4_WORKFLOW_RUN,
                    "agent4_artifact_id": AGENT4_ARTIFACT_ID,
                    "agent4_artifact_zip_digest": AGENT4_ARTIFACT_DIGEST,
                },
                "after_identity_pin": (
                    "materialize same-cycle theta/axial requestedStress and candidate-specific finite-head mean debt; "
                    "then run DeltaC/epsilon inverse plus existing bounded-inverse, budget, spacetime, radial, quadratic, and joint-gain guards"
                ),
            },
            "final_candidate": {
                "producer": "agent5",
                "ready": False,
                "required_inputs": [
                    "non_obstructed_global_leading_velocity_and_matched_pressure",
                    "current_agent4_admitted_public_oscillatory_velocity",
                    "agent3_guarded_materialized_finite_correction_cycle",
                ],
                "planned_api": [
                    "velocity(x,y,z,t)",
                    "pressure(x,y,z,t)",
                    "restricted_forcing(x,y,z,t)",
                ],
                "artifact_schema_required": True,
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
            "oscillatory_component_preflight": {
                "prior_agent4_553_finest_relative_divergence_rms": PRIOR_AGENT4_553_METRICS[
                    "finest_relative_divergence_rms"
                ],
                "current_agent4_563_finest_relative_divergence_rms": AGENT4_PASS_METRICS[
                    "relative_divergence_rms_by_step"
                ]["0.005"],
                "rms_improvement_factor": PRIOR_AGENT4_553_METRICS[
                    "finest_relative_divergence_rms"
                ]
                / AGENT4_PASS_METRICS["relative_divergence_rms_by_step"]["0.005"],
                "is_full_ns_comparison": False,
            },
        },
        "formal_gates": copy.deepcopy(FORMAL_GATES),
        "truth_boundary": {
            "oscillatory_component_preflight_only": True,
            "kokuno_hidden_background_recovered": False,
            "paper_exact": False,
            "free_forcing_used_to_cancel_residual": False,
            "thresholds_relaxed_after_result": False,
            "agent4_component_preflight_counted_as_full_pde_validation": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    validate_checkpoint(checkpoint)
    return checkpoint


def validate_checkpoint(checkpoint: dict[str, Any]) -> None:
    """Fail closed on any over-promotion beyond the independently audited component PASS."""
    if checkpoint.get("schema") != SCHEMA:
        raise ValueError("unexpected checkpoint schema")
    if checkpoint.get("checkpoint_sha256") != _canonical_sha256(checkpoint):
        raise ValueError("checkpoint SHA256 mismatch")

    audit = checkpoint["upstream"]["agent4_public_z_independent_audit"]
    if audit["audited_agent2_head"] != AGENT2_HEAD or audit["head"] != AGENT4_HEAD:
        raise ValueError("oscillatory audit identity mismatch")
    if audit["artifact_id"] != AGENT4_ARTIFACT_ID:
        raise ValueError("oscillatory audit artifact mismatch")
    if audit["artifact_zip_digest"] != AGENT4_ARTIFACT_DIGEST:
        raise ValueError("oscillatory audit digest mismatch")
    if audit["scientific_guards_changed_from_agent4_543"]:
        raise ValueError("scientific guards must remain frozen")
    if audit["derived_guard_passes"] != _derived_passes() or not all(
        audit["derived_guard_passes"].values()
    ):
        raise ValueError("Agent-4 PASS metrics do not satisfy frozen guards")

    states = checkpoint["states"]
    if states["leading_ready"]:
        raise ValueError("leading must remain closed")
    for key in (
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready",
        "public_oscillatory_project_support_passed",
        "public_oscillatory_covariance_rank_passed",
        "candidate_independent_second_covariance_column_ready",
        "public_oscillatory_local_divergence_passed",
        "public_oscillatory_preflight_passed",
        "oscillatory_ready",
    ):
        if not states[key]:
            raise ValueError(f"independently admitted oscillatory state lost: {key}")

    if states["correction_receipt_identity_pinned_to_current_pass"]:
        raise ValueError("Agent 3 has not yet pinned the current PASS receipt")
    if states["correction_ingest_allowed"]:
        raise ValueError("correction ingest cannot open before Agent-3 identity pin")
    for key in (
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "real_candidate_defect_consumed",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ):
        if states[key]:
            raise ValueError(f"downstream state must remain closed: {key}")

    truth = checkpoint["truth_boundary"]
    if truth["paper_exact"] or truth["kokuno_hidden_background_recovered"]:
        raise ValueError("source truth boundary violated")
    if truth["free_forcing_used_to_cancel_residual"]:
        raise ValueError("free forcing is forbidden")
    if truth["thresholds_relaxed_after_result"]:
        raise ValueError("threshold relaxation is forbidden")
    if truth["agent4_component_preflight_counted_as_full_pde_validation"]:
        raise ValueError("component preflight is not full PDE validation")
    if truth["formal_full_domain_pde_gate_assessed"] or truth["pde_validated"]:
        raise ValueError("full-domain PDE gate remains unassessed")
    if checkpoint["formal_gates"] != FORMAL_GATES:
        raise ValueError("formal project gates changed")


def save_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = build_checkpoint()
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_checkpoint(payload)
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = save_checkpoint(args.output)
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
