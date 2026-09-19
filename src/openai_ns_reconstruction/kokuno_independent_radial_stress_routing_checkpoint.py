"""Kokuno Agent-5 routing after independent compact radial-stress admission.

This checkpoint performs one narrow integration step: it consumes Agent 3 PR #598's
exact admission of Agent 4 PR #595's independent compact/radial-stress audit and
marks that operator seam reusable *when* a real full same-cycle composite defect
exists.  It does not create that missing defect, materialize a signed inverse or
correction velocity, instantiate a composite candidate, or assess the full-domain
Navier--Stokes gate.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_independent_radial_stress_admission import (
    evaluate_independent_radial_stress_admission,
)

TASK = "KOKUNO-A5-INDEPENDENT-RADIAL-STRESS-ROUTING-044"
SCHEMA = "kokuno-agent5-independent-radial-stress-routing-checkpoint-v44"

ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_MOMENTUM_L2 = 0.10758432876230622

# Already admitted public oscillatory interfaces.  These exact identities are
# unchanged in this round.
AGENT2_VELOCITY_PR = 561
AGENT2_VELOCITY_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
AGENT2_TIME_DERIVATIVE_PR = 579
AGENT2_TIME_DERIVATIVE_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"

# Agent 3's independently audited radial-stress admission being integrated here.
AGENT3_RADIAL_ADMISSION_PR = 598
AGENT3_RADIAL_ADMISSION_HEAD = "a97f0882ab0b1fcc0d4b9575b5366fb44b35c233"
AGENT3_RADIAL_ADMISSION_DEDICATED_RUN = 35432615191
AGENT3_RADIAL_ADMISSION_STANDARD_RUN = 35432615180
AGENT3_RADIAL_ADMISSION_ARTIFACT_ID = 10581435351
AGENT3_RADIAL_ADMISSION_ARTIFACT_DIGEST = (
    "sha256:be1c19d3257086622f90ef5a8b7be9181dcbe9628d568eb8a4f9bfeac213ec6c"
)

# Exact Agent-4 audit transitively pinned by Agent 3 #598.
AGENT4_RADIAL_AUDIT_PR = 595
AGENT4_RADIAL_AUDIT_HEAD = "19296acad4f84ba05d3b095bb8c01bf5a2c95891"
AGENT4_RADIAL_AUDIT_DEDICATED_RUN = 35431599990
AGENT4_RADIAL_AUDIT_STANDARD_RUN = 35431599977
AGENT4_RADIAL_AUDIT_ARTIFACT_ID = 10579939940
AGENT4_RADIAL_AUDIT_ARTIFACT_DIGEST = (
    "sha256:295980b23ff31877f341bbf3d7866718678dcd716fc4bd5d4efb3d0ff8c3845a"
)
A4_AUDIT_RECEIPT_PATH = Path(
    "artifacts/constrained/kokuno_agent3/agent4_595_independent_radial_stress_audit.json"
)

# Latest Agent-1 sibling at the v44 freeze.  Its exact R1/R2 algebra is now
# executable, but the primitive radius-one-ball bounds remain diagnostic rather
# than source-certified, so no leading-state promotion is allowed here.
AGENT1_PR = 597
AGENT1_HEAD = "d652364ff02ffa90b25c6cc19d89d08d88b5dbf4"
AGENT1_DEDICATED_RUN = 35432188984
AGENT1_STANDARD_RUN = 35432207095
AGENT1_STANDARD_STATUS_AT_FREEZE = "in_progress"

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


def _json_normalized(value: Any) -> Any:
    """Normalize tuple/list distinctions before signing a saved checkpoint."""
    return json.loads(json.dumps(value, sort_keys=True, allow_nan=False))


def _admit_radial_stress(receipt_path: Path) -> dict[str, Any]:
    result = evaluate_independent_radial_stress_admission(receipt_path.read_bytes())
    normalized = _json_normalized(result)
    if not isinstance(normalized, dict):
        raise ValueError("Agent-3 radial-stress admission did not return an object")
    required_true = (
        "execution_identity_bound",
        "independent_compact_radial_stress_cross_audit_passed",
        "compact_radial_stress_operator_independently_admitted",
        "full_composite_radial_stress_execution_allowed_when_actual_defect_available",
        "real_oscillatory_component_defect_independently_consumed",
    )
    for key in required_true:
        if normalized.get(key) is not True:
            raise ValueError(f"required radial-stress admission state is false: {key}")
    required_false = (
        "full_same_cycle_composite_requested_stress_materialized",
        "agent1_leading_cross_terms_included",
        "matched_pressure_included",
        "restricted_forcing_included",
        "candidate_finite_head_mean_debt_materialized",
        "real_full_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "public_velocity_correction_materialized",
        "finite_correction_cycle_rerun_allowed",
        "finite_correction_cycle_run",
        "heldout_ns_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
        "surrogate_defect_used",
    )
    for key in required_false:
        if normalized.get(key) is not False:
            raise ValueError(f"radial-stress admission over-promoted downstream state: {key}")
    return normalized


def build_checkpoint(receipt_path: Path = A4_AUDIT_RECEIPT_PATH) -> dict[str, Any]:
    radial = _admit_radial_stress(receipt_path)
    theta = radial["theta_e2"]
    axial = radial["axial_e1"]

    checkpoint: dict[str, Any] = {
        "task": TASK,
        "schema": SCHEMA,
        "round": 44,
        "purpose": (
            "admit the independently cross-audited compact radial-stress operator as a reusable "
            "correction-stage handoff while keeping the missing full same-cycle defect and every "
            "downstream correction/composite/PDE gate fail closed"
        ),
        "upstream": {
            "agent1_latest_leading": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "dedicated_workflow_run": AGENT1_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT1_STANDARD_RUN,
                "standard_workflow_status_at_checkpoint_freeze": AGENT1_STANDARD_STATUS_AT_FREEZE,
                "exact_displayed_R1_R2_monomial_bookkeeping_executable": True,
                "one_factor_at_a_time_lipschitz_bookkeeping_executable": True,
                "source_R1_radius_one_ball_norm_machine_bound": False,
                "source_R2_radius_one_ball_norm_machine_bound": False,
                "source_R1_radius_one_ball_lipschitz_machine_bound": False,
                "source_R2_radius_one_ball_lipschitz_machine_bound": False,
                "source_operator_M_K_machine_bound": False,
                "source_B0_T_sh_machine_bound": False,
                "global_leading_velocity_pressure_ready": False,
            },
            "agent2_admitted_oscillatory": {
                "velocity_pr": AGENT2_VELOCITY_PR,
                "velocity_head": AGENT2_VELOCITY_HEAD,
                "time_derivative_pr": AGENT2_TIME_DERIVATIVE_PR,
                "time_derivative_head": AGENT2_TIME_DERIVATIVE_HEAD,
                "public_velocity_api": "velocity_osc(x,y,z,t)->[...,3]",
                "public_time_derivative_api": "velocity_osc_dt(x,y,z,t)->[...,3]",
                "velocity_candidate_changed_in_this_round": False,
            },
            "agent3_independent_radial_stress_admission": {
                "pr": AGENT3_RADIAL_ADMISSION_PR,
                "head": AGENT3_RADIAL_ADMISSION_HEAD,
                "dedicated_workflow_run": AGENT3_RADIAL_ADMISSION_DEDICATED_RUN,
                "dedicated_workflow_conclusion": "success",
                "standard_workflow_run": AGENT3_RADIAL_ADMISSION_STANDARD_RUN,
                "standard_workflow_conclusion": "success",
                "artifact_id": AGENT3_RADIAL_ADMISSION_ARTIFACT_ID,
                "artifact_zip_digest": AGENT3_RADIAL_ADMISSION_ARTIFACT_DIGEST,
                "agent4_audit_pr": AGENT4_RADIAL_AUDIT_PR,
                "agent4_audit_head": AGENT4_RADIAL_AUDIT_HEAD,
                "agent4_dedicated_workflow_run": AGENT4_RADIAL_AUDIT_DEDICATED_RUN,
                "agent4_standard_workflow_run": AGENT4_RADIAL_AUDIT_STANDARD_RUN,
                "agent4_artifact_id": AGENT4_RADIAL_AUDIT_ARTIFACT_ID,
                "agent4_artifact_zip_digest": AGENT4_RADIAL_AUDIT_ARTIFACT_DIGEST,
                "failed_guards": radial["failed_guards"],
                "theta_e2": theta,
                "axial_e1": axial,
                "independent_compact_radial_stress_cross_audit_passed": True,
                "compact_radial_stress_operator_independently_admitted": True,
                "full_composite_radial_stress_execution_allowed_when_actual_defect_available": True,
            },
        },
        "states": {
            "leading_ready": False,
            "public_provenance_labelled_xyz_t_oscillatory_velocity_ready": True,
            "public_oscillatory_preflight_passed": True,
            "oscillatory_ready": True,
            "public_oscillatory_time_derivative_handoff_ready": True,
            "public_oscillatory_time_derivative_independently_validated": True,
            "correction_receipt_identity_pinned_to_current_pass": True,
            "correction_ingest_allowed": True,
            "oscillatory_requested_stress_component_materialized": True,
            "oscillatory_component_finite_head_mean_debt_materialized": True,
            "independent_compact_radial_stress_cross_audit_passed": True,
            "compact_radial_stress_operator_independently_admitted": True,
            "full_composite_radial_stress_execution_allowed_when_actual_defect_available": True,
            "same_cycle_requested_stress_materialized": False,
            "candidate_numeric_finite_head_mean_debt_materialized": False,
            "real_full_candidate_defect_consumed": False,
            "signed_mean_inverse_input_ready": False,
            "public_velocity_correction_materialized": False,
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
                    "bind the explicit PA.10 primitive radius-one-ball inputs, especially the pressure "
                    "primitive and mixed radial/eta derivative seams, then promote genuine source R1/R2 "
                    "norm/Lipschitz -> M,K -> B0/T_sh -> PA.16 -> global matched pressure"
                ),
            },
            "oscillatory": {
                "producer": "agent2",
                "validator": "agent4",
                "ready": True,
                "velocity_api": "velocity_osc(x,y,z,t)->[...,3]",
                "time_derivative_api": "velocity_osc_dt(x,y,z,t)->[...,3]",
                "retune_allowed_by_this_checkpoint": False,
            },
            "correction": {
                "producer": "agent3",
                "ingest_allowed": True,
                "oscillatory_component_requested_stress_materialized": True,
                "oscillatory_component_mean_debt_materialized": True,
                "compact_radial_stress_operator_independently_admitted": True,
                "full_composite_radial_stress_execution_allowed_when_actual_defect_available": True,
                "full_same_cycle_requested_stress_materialized": False,
                "full_candidate_mean_debt_materialized": False,
                "signed_mean_inverse_input_ready": False,
                "ready": False,
                "next_required": (
                    "after Agent-1 supplies the actual same-cycle leading/cross + matched-pressure + "
                    "restricted-forcing handoff, form the full theta/axial defect, run the now-independently-"
                    "admitted compact radial-stress reconstruction on that full defect, materialize the full "
                    "finite-head target, and only then enter DeltaC/epsilon -> signed inverse -> bounded "
                    "inverse -> budget -> spacetime -> radial correction -> damping/joint-gain guards"
                ),
            },
            "final_candidate": {
                "ready": False,
                "artifact_schema_ready": False,
                "save_load_ready": False,
                "python_visualization_smoke_ready": False,
                "matlab_visualization_smoke_ready": False,
                "independent_pde_validation_ready": False,
            },
        },
        "baseline_vs_kokuno": {
            "st006": {
                "same_protocol_full_domain_baseline": True,
                "momentum_sampled_max": ST006_MOMENTUM_MAX,
                "momentum_volume_l2": ST006_MOMENTUM_L2,
                "pde_validated": False,
            },
            "kokuno": {
                "comparable_full_domain_momentum_receipt_exists": False,
                "radial_stress_metrics_are_component_operator_diagnostics_not_ns_residuals": True,
                "theta_finest_operator_relative_rms": theta["finest_relative_rms"],
                "axial_finest_operator_relative_rms": axial["finest_relative_rms"],
            },
        },
        "formal_gates_unchanged": copy.deepcopy(FORMAL_GATES),
        "truth_boundary": {
            "kokuno_reconstruction_is_paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "free_residual_defined_forcing_used": False,
            "surrogate_defect_used_for_promotion": False,
            "validation_thresholds_relaxed": False,
            "component_radial_stress_laundered_as_full_same_cycle_correction": False,
            "component_operator_metric_laundered_as_pde_residual": False,
            "agent4_replay_laundered_as_final_independent_pde_validation": False,
        },
        "blockers": [
            "Agent-1 global leading/cross terms, matched pressure and restricted-forcing handoff are not ready",
            "full same-cycle composite theta/axial defect has not been materialized",
            "candidate full finite-head mean debt has not been materialized",
            "signed mean inverse and public correction velocity remain closed",
            "no composite velocity/pressure/restricted-forcing artifact exists",
            "no formal held-out full-domain NS residual has been assessed",
        ],
    }
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    return checkpoint


def validate_checkpoint(checkpoint: dict[str, Any]) -> None:
    if checkpoint.get("schema") != SCHEMA or checkpoint.get("task") != TASK:
        raise ValueError("wrong Agent-5 v44 schema/task")
    if checkpoint.get("checkpoint_sha256") != _canonical_sha256(checkpoint):
        raise ValueError("checkpoint SHA256 mismatch")

    states = checkpoint["states"]
    required_true = (
        "oscillatory_ready",
        "public_oscillatory_time_derivative_independently_validated",
        "correction_ingest_allowed",
        "oscillatory_requested_stress_component_materialized",
        "oscillatory_component_finite_head_mean_debt_materialized",
        "independent_compact_radial_stress_cross_audit_passed",
        "compact_radial_stress_operator_independently_admitted",
        "full_composite_radial_stress_execution_allowed_when_actual_defect_available",
    )
    for key in required_true:
        if states.get(key) is not True:
            raise ValueError(f"required admitted state is false: {key}")

    required_false = (
        "leading_ready",
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "real_full_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "public_velocity_correction_materialized",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    )
    for key in required_false:
        if states.get(key) is not False:
            raise ValueError(f"downstream state promoted too early: {key}")

    if checkpoint["formal_gates_unchanged"] != FORMAL_GATES:
        raise ValueError("formal final gates changed")

    leading = checkpoint["upstream"]["agent1_latest_leading"]
    for key in (
        "source_R1_radius_one_ball_norm_machine_bound",
        "source_R2_radius_one_ball_norm_machine_bound",
        "source_R1_radius_one_ball_lipschitz_machine_bound",
        "source_R2_radius_one_ball_lipschitz_machine_bound",
        "source_operator_M_K_machine_bound",
        "source_B0_T_sh_machine_bound",
        "global_leading_velocity_pressure_ready",
    ):
        if leading.get(key) is not False:
            raise ValueError(f"Agent-1 diagnostic algebra laundered into source certificate: {key}")

    radial = checkpoint["upstream"]["agent3_independent_radial_stress_admission"]
    if radial.get("head") != AGENT3_RADIAL_ADMISSION_HEAD:
        raise ValueError("wrong Agent-3 radial-stress admission head")
    if radial.get("artifact_id") != AGENT3_RADIAL_ADMISSION_ARTIFACT_ID:
        raise ValueError("wrong Agent-3 radial-stress admission artifact")
    if radial.get("artifact_zip_digest") != AGENT3_RADIAL_ADMISSION_ARTIFACT_DIGEST:
        raise ValueError("wrong Agent-3 radial-stress admission artifact digest")
    if radial.get("agent4_audit_head") != AGENT4_RADIAL_AUDIT_HEAD:
        raise ValueError("wrong Agent-4 radial-stress audit head")
    if radial.get("agent4_artifact_id") != AGENT4_RADIAL_AUDIT_ARTIFACT_ID:
        raise ValueError("wrong Agent-4 radial-stress audit artifact")
    if radial.get("failed_guards") != []:
        raise ValueError("independent radial-stress admission no longer passes all frozen guards")
    for channel in ("theta_e2", "axial_e1"):
        row = radial[channel]
        if row.get("passed") is not True or row.get("failed_guards") != []:
            raise ValueError(f"independent radial-stress channel no longer passes: {channel}")

    truth = checkpoint["truth_boundary"]
    if truth.get("free_residual_defined_forcing_used") is not False:
        raise ValueError("free residual-defined forcing is prohibited")
    if truth.get("surrogate_defect_used_for_promotion") is not False:
        raise ValueError("surrogate defect cannot be used for promotion")
    if truth.get("validation_thresholds_relaxed") is not False:
        raise ValueError("validation thresholds cannot be relaxed")
    if truth.get("component_radial_stress_laundered_as_full_same_cycle_correction") is not False:
        raise ValueError("component radial stress cannot be laundered as full correction")
    if truth.get("component_operator_metric_laundered_as_pde_residual") is not False:
        raise ValueError("component operator metric cannot be laundered as PDE residual")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, default=A4_AUDIT_RECEIPT_PATH)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    checkpoint = build_checkpoint(args.audit)
    validate_checkpoint(checkpoint)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(checkpoint, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
