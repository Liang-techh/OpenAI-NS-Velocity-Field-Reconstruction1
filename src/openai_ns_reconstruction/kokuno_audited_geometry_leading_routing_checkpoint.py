"""Agent-5 checkpoint for audited autonomous geometry and the PA.10 leading obstruction.

This is one minimal integration increment on top of the v33 autonomous-provenance
checkpoint.  It records two independently useful routing facts without changing
candidate/PDE readiness:

* Agent 4 independently audited Agent 2's repository-autonomous/source-compatible
  signed-rectangle geometry handoff.  The geometry bookkeeping/separation/schedule
  seam is therefore closed for the labelled candidate route, but remains non-source-exact.
* Agent 1's currently selected autonomous shared-C Appendix-B normalization is
  PA.10-obstructed by a necessary pointwise T_sh lower bound exceeding the geometric
  ceiling.  That selected leading path may not proceed to PA.16; a new upstream
  normalization/scale choice must be propagated and screened instead.

Neither fact instantiates a full 3D Kokuno candidate or assesses the formal NS gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .kokuno_autonomous_provenance_routing_checkpoint import (
    AUTONOMOUS_REALIZATION_POLICY,
    FORMAL_GATES,
    SCHEMA as PARENT_SCHEMA,
    TASK_ID as PARENT_TASK_ID,
    build_checkpoint as build_parent_checkpoint,
    checkpoint_sha256,
    validate_checkpoint as validate_parent_checkpoint,
)

SCHEMA = "kokuno-agent5-audited-geometry-leading-routing-checkpoint-v34"
TASK_ID = "KOKUNO-A5-AUDITED-GEOMETRY-LEADING-ROUTING-034"

AGENT1_RECEIPT = {
    "pr": 500,
    "head": "a0037e5cff60be97aa01cc7273eb7af18389ed7a",
    "evidence_class": "candidate_derived",
    "dedicated_run": 35404218836,
    "dedicated_status": "success",
    "standard_run_at_freeze": 35404221549,
    "standard_status_at_freeze": "in_progress",
    "artifact": "kokuno-agent1-rescaled-tsh-certificate-v2",
    "artifact_id": 10572046239,
    "artifact_digest": "sha256:981a858c766780efcce2d261a96bc506bf35da6f942b1ee474ca8570364792f2",
    "selected_shared_C_normalization_is_repository_autonomous": True,
    "selected_shared_C_pa10_obstructed": True,
    "selected_pa16_handoff_allowed": False,
    "eta_w": -0.004449378508888404,
    "ell_i_at_eta_w": 4.451490800486385e23,
    "pointwise_B0_lower_bound": 4.451490800486385e23,
    "necessary_T_sh_lower_bound": 7.122385280778216e25,
    "geometric_T_sh_ceiling": 4.451490800486385e24,
    "separation_margin_upper_bound": -6.677236200729578e25,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

AGENT4_RECEIPT = {
    "pr": 504,
    "head": "aa9d2bb0a6c196ea86daffd5a6f539785a9b0363",
    "evidence_class": "independent_validation",
    "dedicated_run": 35403144495,
    "standard_run": 35403144460,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact": "kokuno-agent4-autonomous-signed-rectangle-independent-audit-v1",
    "artifact_id": 10570869244,
    "artifact_digest": "sha256:0cefa22292ca082d4bff50a9cbae3c9f8eca6f84160aa15701ea13033092226c",
    "autonomous_signed_rectangle_geometry_independently_audited": True,
    "local_structural_preflight_passed": True,
    "schedule_relative_max": 0.0,
    "center_float_relative_max": 0.0,
    "permutation_relative_max": 0.0,
    "strict_geometry_guard_violations": 0,
    "forged_squared_partition_rejections": 3,
    "caller_order_assignment_mutation_mismatch_fraction": 0.8333333333333334,
    "wrong_Ls_formula_mutation_relative_rms": 0.999998125541105,
    "fixed_label_sigma_swap_mismatch_fraction": 1.0,
    "source_rectangle_centers_recovered": False,
    "source_rectangle_radius_r0_recovered": False,
    "source_actual_partition_labels_instantiated": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 506,
    "head": "96e9514766ca5c512dc30409fbfa70b732fb56e7",
    "schema": "kokuno-agent5-autonomous-provenance-routing-checkpoint-v33",
    "dedicated_run": 35404043854,
    "standard_run": 35404043821,
    "dedicated_status": "success",
    "standard_status": "success",
    "consumed_in_executable_ancestry": True,
}


def _copy(value: Any) -> Any:
    return json.loads(json.dumps(value))


def _set_pipeline_status(payload: dict[str, Any], stage: str, status: str) -> None:
    for item in payload["pipeline"]:
        if item["stage"] == stage:
            item["ready"] = False
            item["status"] = status
            return
    raise ValueError(f"missing pipeline stage: {stage}")


def build_checkpoint() -> dict[str, Any]:
    payload = _copy(build_parent_checkpoint())
    payload.pop("checkpoint_sha256", None)
    payload["schema"] = SCHEMA
    payload["task_id"] = TASK_ID
    payload["upstream"]["agent1_latest_sibling"] = _copy(AGENT1_RECEIPT)
    payload["upstream"]["agent4_sibling"] = _copy(AGENT4_RECEIPT)
    payload["upstream"]["previous_agent5_ancestry"] = _copy(PREVIOUS_AGENT5_RECEIPT)

    states = payload["states"]
    states["autonomous_signed_rectangle_geometry_independently_audited"] = True
    states["selected_shared_C_pa10_path_obstructed"] = True
    states["selected_shared_C_pa16_handoff_allowed"] = False
    states["leading_reparameterization_required"] = True

    _set_pipeline_status(
        payload,
        "leading_candidate",
        "the currently selected autonomous shared-C normalization is PA.10-obstructed; PA.16 handoff is forbidden and a new upstream normalization/scale choice must be propagated and screened before join/I3/I4/matched-pressure work resumes",
    )
    _set_pipeline_status(
        payload,
        "oscillatory_augmentation",
        "the repository-autonomous signed-rectangle geometry handoff is now independently audited; do not repeat rectangle arithmetic, and next bind provenance-labelled numerical h_sigma/background/modes into the audited phase/complete-curl path to emit public Q-scaled by-sign/by-beta/total velocity_osc(x,y,z,t)",
    )
    _set_pipeline_status(
        payload,
        "independent_validation",
        "Agent 4 has independently closed the autonomous signed-rectangle geometry seam; the next valuable audit is the first fully numerical provenance-labelled physical oscillatory velocity/covariance interface, while full-domain NS validation still waits for a frozen global candidate",
    )

    payload["routing"] = {
        "new_fact": "The autonomous signed-rectangle geometry seam is independently closed, while the currently selected autonomous shared-C leading normalization is conclusively PA.10-obstructed. Geometry arithmetic should not be repeated, and the rejected leading normalization must not be pushed through PA.16.",
        "shortest_next_closure": [
            "Agent 1: retire only the obstructed selected shared-C normalization, choose a new explicitly provenance-labelled source-compatible/autonomous normalization or scale upstream, propagate it through the Appendix-B path, and rerun the same PA.10/T_sh necessary screen before any PA.16 join; do not bypass the obstruction post hoc.",
            "Agent 2: consume #501 geometry now independently audited by #504 and bind the first provenance-labelled numerical positive-order/background, h_sigma pulses and auxiliary modes into the audited source-phase/complete-curl path; emit Q-scaled by-sign/by-beta/total velocity_osc(x,y,z,t) with source-recovery flags remaining false.",
            "Agent 3: on the candidate route, serialize one explicitly autonomous numeric bump/profile and compute requestedStress from the same candidate cycle state so a candidate-specific finite-head mean debt can be materialized without claiming the formal source missingWeight was numerically recovered.",
            "Agent 4: black-box audit the first fully numerical provenance-labelled physical oscillatory velocity, complete-curl/divergence behavior, and phase-mean covariance rank; do not spend another round on already-audited autonomous rectangle bookkeeping.",
            "Agent 5: once a non-obstructed leading path plus independently audited oscillatory family plus guarded correction coexist, instantiate the deterministic candidate artifact and Python/MATLAB smoke, then freeze it for Agent 4 held-out normalized NS validation and same-protocol ST006 comparison.",
        ],
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
    if payload.get("autonomous_realization_policy") != AUTONOMOUS_REALIZATION_POLICY:
        raise ValueError("autonomous realization policy changed")

    # Reuse all v33 fail-closed semantics after translating only schema/task/hash.
    parent_view = _copy(payload)
    parent_view["schema"] = PARENT_SCHEMA
    parent_view["task_id"] = PARENT_TASK_ID
    parent_view["checkpoint_sha256"] = checkpoint_sha256(parent_view)
    validate_parent_checkpoint(parent_view)

    states = payload["states"]
    if states.get("autonomous_signed_rectangle_geometry_independently_audited") is not True:
        raise ValueError("independent autonomous geometry audit was removed")
    if states.get("selected_shared_C_pa10_path_obstructed") is not True:
        raise ValueError("selected shared-C PA.10 obstruction was removed")
    if states.get("selected_shared_C_pa16_handoff_allowed") is not False:
        raise ValueError("obstructed shared-C path was incorrectly allowed into PA.16")
    if states.get("leading_reparameterization_required") is not True:
        raise ValueError("leading reparameterization requirement was removed")
    if states.get("leading_ready") is not False:
        raise ValueError("leading readiness promoted despite PA.10 obstruction")

    a1 = payload["upstream"]["agent1_latest_sibling"]
    if a1 != AGENT1_RECEIPT:
        raise ValueError("Agent 1 obstruction receipt changed")
    if not (a1["necessary_T_sh_lower_bound"] > a1["geometric_T_sh_ceiling"]):
        raise ValueError("Agent 1 obstruction inequality no longer holds")
    if a1["selected_pa16_handoff_allowed"]:
        raise ValueError("Agent 1 obstructed handoff was promoted")
    if a1["source_B0_analytic_bound_proved"] or a1["source_T_sh_lower_bound_verified"]:
        raise ValueError("selected candidate obstruction was laundered into source theorem truth")

    a4 = payload["upstream"]["agent4_sibling"]
    if a4 != AGENT4_RECEIPT:
        raise ValueError("Agent 4 independent geometry receipt changed")
    if not a4["local_structural_preflight_passed"]:
        raise ValueError("independent autonomous geometry preflight lost")
    if a4["source_rectangle_centers_recovered"] or a4["source_rectangle_radius_r0_recovered"]:
        raise ValueError("autonomous geometry was promoted to recovered source data")
    if a4["formal_full_domain_pde_gate_assessed"] or a4["pde_validated"]:
        raise ValueError("local geometry audit was promoted to PDE validation")

    previous = payload["upstream"].get("previous_agent5_ancestry")
    if previous != PREVIOUS_AGENT5_RECEIPT:
        raise ValueError("previous Agent 5 receipt changed")
    if payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is not None:
        raise ValueError("no comparable full-domain Kokuno receipt exists yet")


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
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
