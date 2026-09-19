"""Agent-5 routing checkpoint after Agent-3 pins the accepted oscillatory receipt.

This is one narrow integration increment.  Agent 3 PR #571 binds the exact
Agent-2 #561 / Agent-4 #563 independent oscillatory PASS and opens only the
same-cycle requested-stress / candidate finite-head mean-debt materialization
stage.  Agent 4 PR #572 additionally supplies a disjoint-seed FD6 robustness
PASS for the *same unchanged* public oscillatory field.

Neither receipt materializes a correction, produces a global leading field, or
assesses the full Navier--Stokes momentum residual.  This module therefore
promotes correction *ingest permission* only and keeps every downstream PDE
state fail closed.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_public_z_pass_admission import evaluate_public_z_pass_admission

TASK = "KOKUNO-A5-CORRECTION-INGEST-ROUTING-041"
SCHEMA = "kokuno-agent5-correction-ingest-routing-checkpoint-v41"

# Retained repository PDE baseline.
ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622

# Agent 1 latest leading-side sibling.  The universal arithmetic is executable,
# but the actual source-dependent inputs needed to certify M,K remain absent.
AGENT1_PR = 569
AGENT1_HEAD = "8836af19f77518f4e160fd1e8d29b2eb57dc61e5"
AGENT1_DEDICATED_RUN = 35424195154
AGENT1_ARTIFACT_ID = 10579052096
AGENT1_ARTIFACT_DIGEST = (
    "sha256:3f8a12512b49121a7ed788bf9b0e0bc165baec54d57c88f0543e789101f93e45"
)

# Admitted public oscillatory candidate identity.  Agent 2 #570 is diagnostics
# only and explicitly does not change these candidate bytes.
AGENT2_ADMITTED_PR = 561
AGENT2_ADMITTED_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
AGENT2_MORPHOLOGY_PR = 570
AGENT2_MORPHOLOGY_HEAD = "134095115fe88e4e7c285a02d8592fd81451b44c"
AGENT2_MORPHOLOGY_DEDICATED_RUN = 35424316352
AGENT2_MORPHOLOGY_STANDARD_RUN = 35424316335
AGENT2_MORPHOLOGY_ARTIFACT_ID = 10578161907
AGENT2_MORPHOLOGY_ARTIFACT_DIGEST = (
    "sha256:6f04cf0a229b33e7e9772c902f7a6bc7ae4808120e0f4f17856cf6e1f4cd9c6e"
)

# Agent 3 exact admission identity.  Its checked-in audit receipt is re-evaluated
# rather than trusting a reported PASS bit.
AGENT3_PR = 571
AGENT3_HEAD = "51e50f20712b51c14e96a0f5276db9f9a3a3da9f"
AGENT3_DEDICATED_RUN = 35424680539
AGENT3_STANDARD_RUN = 35424680559
AGENT3_ARTIFACT_ID = 10578072258
AGENT3_ARTIFACT_DIGEST = (
    "sha256:757c9d09e6ed6e3b30cfe7a5c53497ce2c19e442503c7abdf9a5b37d34e50eb2"
)

# Agent 4 #572 is supplemental generalization evidence for the exact same A2
# field.  It does NOT replace the #561/#563 identity admitted by Agent 3.
AGENT4_GENERALIZATION_PR = 572
AGENT4_GENERALIZATION_HEAD = "9385af81ea4f4e3b2b2dfdf3862534a15a8b72e8"
AGENT4_GENERALIZATION_DEDICATED_RUN = 35424884191
AGENT4_GENERALIZATION_STANDARD_RUN = 35424884153
AGENT4_GENERALIZATION_ARTIFACT_ID = 10578353800
AGENT4_GENERALIZATION_ARTIFACT_DIGEST = (
    "sha256:6f9e03d8359f22ef4794d32d12991c7e0258b5627d2ba1256e4ec6551a2dccf0"
)

A4_GENERALIZATION_GUARDS = {
    "finest_relative_divergence_rms": 2.0e-5,
    "finest_relative_divergence_max": 1.0e-4,
    "minimum_divergence_refinement_ratio": 4.0,
    "maximum_stratum_finest_relative_divergence_rms": 5.0e-5,
    "support_exterior_absolute_max": 1.0e-12,
    "minimum_nontrivial_velocity_rms": 1.0e-8,
    "minimum_parameter_perturbation_relative_change": 1.0e-4,
    "minimum_divergence_mutation_relative_rms": 5.0e-2,
}

A4_GENERALIZATION_METRICS = {
    "seed": 9173251,
    "point_count": 48,
    "points_per_stratum": 16,
    "fd6_steps": [0.016, 0.008, 0.004],
    "relative_divergence_rms_by_step": [
        4.6147875751247766e-7,
        7.317464739749743e-9,
        1.1476081461351159e-10,
    ],
    "finest_relative_divergence_max": 1.4435022950460222e-10,
    "refinement_ratios": [63.0653886182252, 63.76274658203939],
    "finest_stratum_relative_rms": {
        "axial_collar": 4.392573150126145e-11,
        "core": 1.1373750688305845e-10,
        "radial_collar": 1.9080356468531223e-9,
    },
    "support_exterior_absolute_max": 0.0,
    "heldout_velocity_rms": 3803.527071088177,
    "h_minus_10pct_relative_change": 0.0015337575050953206,
    "h_plus_10pct_relative_change": 0.00153553766859264,
    "divergence_mutation_relative_rms": 0.2500000000200818,
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


def _supplemental_generalization_passes(
    metrics: dict[str, Any] | None = None,
) -> dict[str, bool]:
    m = A4_GENERALIZATION_METRICS if metrics is None else metrics
    return {
        "divergence_rms": (
            m["relative_divergence_rms_by_step"][-1]
            <= A4_GENERALIZATION_GUARDS["finest_relative_divergence_rms"]
        ),
        "divergence_max": (
            m["finest_relative_divergence_max"]
            <= A4_GENERALIZATION_GUARDS["finest_relative_divergence_max"]
        ),
        "refinement": (
            min(m["refinement_ratios"])
            >= A4_GENERALIZATION_GUARDS["minimum_divergence_refinement_ratio"]
        ),
        "strata": (
            max(m["finest_stratum_relative_rms"].values())
            <= A4_GENERALIZATION_GUARDS[
                "maximum_stratum_finest_relative_divergence_rms"
            ]
        ),
        "support": (
            m["support_exterior_absolute_max"]
            <= A4_GENERALIZATION_GUARDS["support_exterior_absolute_max"]
        ),
        "nontriviality": (
            m["heldout_velocity_rms"]
            >= A4_GENERALIZATION_GUARDS["minimum_nontrivial_velocity_rms"]
        ),
        "parameter_response": (
            min(
                m["h_minus_10pct_relative_change"],
                m["h_plus_10pct_relative_change"],
            )
            >= A4_GENERALIZATION_GUARDS[
                "minimum_parameter_perturbation_relative_change"
            ]
        ),
        "mutation_detection": (
            m["divergence_mutation_relative_rms"]
            >= A4_GENERALIZATION_GUARDS[
                "minimum_divergence_mutation_relative_rms"
            ]
        ),
    }


def build_checkpoint(raw_agent4_563_receipt: bytes) -> dict[str, Any]:
    """Build v41 by independently replaying Agent 3's exact admission gate."""
    admission = evaluate_public_z_pass_admission(raw_agent4_563_receipt)
    if admission["correction_receipt_identity_pinned_to_current_pass"] is not True:
        raise ValueError("Agent-3 exact PASS identity is not pinned")
    if admission["correction_ingest_allowed"] is not True:
        raise ValueError("Agent-3 exact PASS does not authorize correction ingest")
    for key in (
        "requested_stress_actual_state_values_materialized",
        "finite_head_mean_debt_materialized",
        "real_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "public_velocity_correction_materialized",
        "finite_correction_cycle_run",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ):
        if admission[key] is not False:
            raise ValueError(f"Agent-3 admission unexpectedly promoted {key}")

    generalization_passes = _supplemental_generalization_passes()
    if not all(generalization_passes.values()):
        raise ValueError("frozen Agent-4 #572 supplemental receipt no longer passes")

    checkpoint: dict[str, Any] = {
        "task": TASK,
        "schema": SCHEMA,
        "round": 41,
        "purpose": (
            "admit Agent-3's exact oscillatory PASS identity into the integration route, "
            "open only same-cycle correction-input materialization, and record Agent-4 "
            "cross-seed FD6 generalization without fabricating a correction or PDE result"
        ),
        "upstream": {
            "agent1_leading_sibling": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "dedicated_workflow_run": AGENT1_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "artifact_id": AGENT1_ARTIFACT_ID,
                "artifact_zip_digest": AGENT1_ARTIFACT_DIGEST,
                "universal_operator_primitives_executable": True,
                "source_rho_machine_bound": False,
                "source_M_chi_machine_bound": False,
                "source_R1_R2_ball_and_lipschitz_bounds_machine_bound": False,
                "source_operator_M_K_machine_bound": False,
                "source_B0_T_sh_verified": False,
                "selected_pa16_handoff_allowed": False,
                "global_leading_velocity_pressure_ready": False,
            },
            "agent2_admitted_public_oscillatory": {
                "pr": AGENT2_ADMITTED_PR,
                "head": AGENT2_ADMITTED_HEAD,
                "public_api": (
                    "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc"
                ),
                "coordinate_contract": "bandwise Z_beta=epsilon_beta*z",
                "velocity_candidate_frozen": True,
                "repository_autonomous_source_compatible": True,
                "paper_exact": False,
            },
            "agent2_morphology_sibling": {
                "pr": AGENT2_MORPHOLOGY_PR,
                "head": AGENT2_MORPHOLOGY_HEAD,
                "dedicated_workflow_run": AGENT2_MORPHOLOGY_DEDICATED_RUN,
                "standard_workflow_run": AGENT2_MORPHOLOGY_STANDARD_RUN,
                "artifact_id": AGENT2_MORPHOLOGY_ARTIFACT_ID,
                "artifact_zip_digest": AGENT2_MORPHOLOGY_ARTIFACT_DIGEST,
                "velocity_candidate_changed": False,
                "descriptive_morphology_only": True,
            },
            "agent3_current_pass_admission": {
                "pr": AGENT3_PR,
                "head": AGENT3_HEAD,
                "dedicated_workflow_run": AGENT3_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT3_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "artifact_id": AGENT3_ARTIFACT_ID,
                "artifact_zip_digest": AGENT3_ARTIFACT_DIGEST,
                "recomputed_admission": copy.deepcopy(admission),
            },
            "agent4_supplemental_generalization": {
                "pr": AGENT4_GENERALIZATION_PR,
                "head": AGENT4_GENERALIZATION_HEAD,
                "audited_agent2_head": AGENT2_ADMITTED_HEAD,
                "dedicated_workflow_run": AGENT4_GENERALIZATION_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT4_GENERALIZATION_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "artifact_id": AGENT4_GENERALIZATION_ARTIFACT_ID,
                "artifact_zip_digest": AGENT4_GENERALIZATION_ARTIFACT_DIGEST,
                "guards": copy.deepcopy(A4_GENERALIZATION_GUARDS),
                "metrics": copy.deepcopy(A4_GENERALIZATION_METRICS),
                "derived_guard_passes": generalization_passes,
                "scientific_verdict": "PASS_supplemental_oscillatory_generalization_only",
                "replaces_agent3_admission_identity": False,
            },
        },
        "states": {
            "leading_ready": False,
            "public_provenance_labelled_xyz_t_oscillatory_velocity_ready": True,
            "public_oscillatory_preflight_passed": True,
            "oscillatory_generalization_independently_passed": True,
            "oscillatory_ready": True,
            "correction_receipt_identity_pinned_to_current_pass": True,
            "correction_ingest_allowed": True,
            "same_cycle_requested_stress_materialization_allowed": True,
            "candidate_finite_head_mean_debt_materialization_allowed": True,
            "same_cycle_requested_stress_materialized": False,
            "candidate_numeric_finite_head_mean_debt_materialized": False,
            "real_candidate_defect_consumed": False,
            "signed_mean_inverse_input_ready": False,
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
                    "bind source-dependent rho/M_chi/R1/R2 ball and Lipschitz inputs, "
                    "then certify source M,K -> B0/T_sh -> PA.16 -> global matched pressure"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "validator": "agent4",
                "candidate_head": AGENT2_ADMITTED_HEAD,
                "public_api": "velocity_osc(x,y,z,t)->[...,3]",
                "admitted": True,
                "supplemental_cross_seed_fd6_generalization_passed": True,
                "frozen_do_not_retune": True,
            },
            "correction": {
                "producer": "agent3",
                "ingest_allowed": True,
                "stage": "materialize_same_cycle_physical_defect_inputs_only",
                "required_next_outputs": [
                    "theta_requestedStress_from_same_admitted_candidate_cycle_state",
                    "axial_requestedStress_from_same_admitted_candidate_cycle_state",
                    "candidate_specific_finite_head_mean_debt",
                ],
                "required_after_materialization": [
                    "DeltaC_over_epsilon_signed_inverse",
                    "bounded_inverse_guard",
                    "coefficient_budget_guard",
                    "spacetime_guard",
                    "radial_plus_dz_sigma1_guard",
                    "quadratic_damping_guard",
                    "joint_gain_guard",
                ],
                "correction_velocity_materialized": False,
                "finite_cycle_run": False,
            },
            "final_candidate": {
                "producer": "agent5",
                "ready": False,
                "required_inputs": [
                    "global_leading_velocity_and_matched_pressure",
                    "admitted_public_oscillatory_velocity",
                    "guarded_materialized_mean_radial_correction",
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
            "component_only_note": (
                "Agent-4 #563/#572 numbers are oscillatory divergence/component diagnostics, "
                "not full Navier-Stokes momentum residuals and not comparable to ST006 momentum"
            ),
        },
        "formal_gates": copy.deepcopy(FORMAL_GATES),
        "truth_boundary": {
            "agent3_pass_identity_is_correction_input_permission_only": True,
            "agent4_572_is_supplemental_not_admission_replacement": True,
            "same_cycle_requested_stress_values_exist": False,
            "candidate_mean_debt_exists": False,
            "correction_velocity_exists": False,
            "global_leading_velocity_pressure_exists": False,
            "composite_candidate_exists": False,
            "free_forcing_used_to_cancel_residual": False,
            "thresholds_relaxed_after_result": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
        },
    }
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    validate_checkpoint(checkpoint)
    return checkpoint


def validate_checkpoint(checkpoint: dict[str, Any]) -> None:
    """Fail closed on identity drift or any promotion beyond correction ingest."""
    if checkpoint.get("schema") != SCHEMA:
        raise ValueError("unexpected checkpoint schema")
    if checkpoint.get("checkpoint_sha256") != _canonical_sha256(checkpoint):
        raise ValueError("checkpoint SHA256 mismatch")

    upstream = checkpoint["upstream"]
    a3 = upstream["agent3_current_pass_admission"]
    if a3["head"] != AGENT3_HEAD or a3["artifact_id"] != AGENT3_ARTIFACT_ID:
        raise ValueError("Agent-3 execution identity mismatch")
    if a3["artifact_zip_digest"] != AGENT3_ARTIFACT_DIGEST:
        raise ValueError("Agent-3 artifact digest mismatch")
    admission = a3["recomputed_admission"]
    if admission["actual_audited_agent2_head"] != AGENT2_ADMITTED_HEAD:
        raise ValueError("Agent-3 admitted candidate identity mismatch")
    if admission["correction_receipt_identity_pinned_to_current_pass"] is not True:
        raise ValueError("Agent-3 PASS identity must remain pinned")
    if admission["correction_ingest_allowed"] is not True:
        raise ValueError("correction ingest must follow exact Agent-3 admission")

    a4 = upstream["agent4_supplemental_generalization"]
    if a4["head"] != AGENT4_GENERALIZATION_HEAD:
        raise ValueError("Agent-4 supplemental execution identity mismatch")
    if a4["artifact_id"] != AGENT4_GENERALIZATION_ARTIFACT_ID:
        raise ValueError("Agent-4 supplemental artifact mismatch")
    if a4["artifact_zip_digest"] != AGENT4_GENERALIZATION_ARTIFACT_DIGEST:
        raise ValueError("Agent-4 supplemental digest mismatch")
    if a4["metrics"] != A4_GENERALIZATION_METRICS:
        raise ValueError("Agent-4 supplemental frozen metrics changed")
    if a4["guards"] != A4_GENERALIZATION_GUARDS:
        raise ValueError("Agent-4 supplemental guards changed")
    expected_generalization = _supplemental_generalization_passes(a4["metrics"])
    if a4["derived_guard_passes"] != expected_generalization or not all(
        expected_generalization.values()
    ):
        raise ValueError("Agent-4 supplemental generalization lost a frozen guard")
    if a4["replaces_agent3_admission_identity"] is not False:
        raise ValueError("supplemental Agent-4 audit cannot replace Agent-3 admission identity")

    states = checkpoint["states"]
    for key in (
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready",
        "public_oscillatory_preflight_passed",
        "oscillatory_generalization_independently_passed",
        "oscillatory_ready",
        "correction_receipt_identity_pinned_to_current_pass",
        "correction_ingest_allowed",
        "same_cycle_requested_stress_materialization_allowed",
        "candidate_finite_head_mean_debt_materialization_allowed",
    ):
        if states.get(key) is not True:
            raise ValueError(f"required admitted state lost: {key}")

    for key in (
        "leading_ready",
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "real_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ):
        if states.get(key) is not False:
            raise ValueError(f"forbidden downstream promotion: {key}")

    if checkpoint["formal_gates"] != FORMAL_GATES:
        raise ValueError("formal PDE gates changed")
    if checkpoint["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is not None:
        raise ValueError("no comparable Kokuno full-domain receipt exists yet")
    truth = checkpoint["truth_boundary"]
    for key in (
        "same_cycle_requested_stress_values_exist",
        "candidate_mean_debt_exists",
        "correction_velocity_exists",
        "global_leading_velocity_pressure_exists",
        "composite_candidate_exists",
        "free_forcing_used_to_cancel_residual",
        "thresholds_relaxed_after_result",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
        "paper_exact",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"truth boundary promotion is forbidden: {key}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    checkpoint = build_checkpoint(args.audit.read_bytes())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(checkpoint, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
