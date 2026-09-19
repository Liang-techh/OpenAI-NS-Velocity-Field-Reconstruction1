"""Fail-closed Agent-5 admission firewall for the Kokuno source-leading lane.

This integration module freezes only evidence that already has an independent
Agent-4 admission receipt and keeps unresolved Agent-1 / Agent-4 seams false.
It exists to prevent an integration branch from promoting ``leading_ready`` or
full R2 merely because an Agent-1 self-certificate or a queued audit exists.

It does not evaluate a Navier--Stokes residual, does not alter any candidate,
and does not expose a path for residual-defined forcing.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA = "kokuno-a5-source-leading-admission-firewall-v1"
TASK = "KOKUNO-A5-SOURCE-LEADING-ADMISSION-FIREWALL-055"
PARENT_A5_PR = 688
PARENT_A5_HEAD = "0233d530e5e24de02195c5b4bd3edf06c2759c19"

FINAL_PROJECT_GATES = {
    "normalized_momentum_max": 1.0e-3,
    "normalized_momentum_L2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_L2": 1.0e-5,
}


@dataclass(frozen=True)
class AdmissionEvidence:
    """One frozen scientific-admission datum."""

    key: str
    producer: str
    pr: int
    head: str
    workflow_id: int | None
    workflow_conclusion_at_freeze: str
    independently_admitted: bool
    scope: str
    note: str


FROZEN_EVIDENCE: dict[str, AdmissionEvidence] = {
    "source_axis": AdmissionEvidence(
        key="source_axis",
        producer="Agent 4",
        pr=671,
        head="6d469766f7cf54730fcd0af3ccc597a2a300d4d4",
        workflow_id=35456226819,
        workflow_conclusion_at_freeze="success",
        independently_admitted=True,
        scope="PA.8 / complex tube / rho / normalized-g source-compatible prerequisite",
        note="Independent PASS only for the source-axis prerequisite.",
    ),
    "pressure_r2_ordinary": AdmissionEvidence(
        key="pressure_r2_ordinary",
        producer="Agent 4",
        pr=678,
        head="a15697b7043a67d9addd90aa5a767c0ff8b5c967",
        workflow_id=35459240658,
        workflow_conclusion_at_freeze="success",
        independently_admitted=True,
        scope="three pressure ordinary R2 slots after J1",
        note="Does not admit the other ordinary R2 slots or full R2.",
    ),
    "source_u_ball": AdmissionEvidence(
        key="source_u_ball",
        producer="Agent 4",
        pr=687,
        head="dcb7450623bedd849efa642b77df610693a1b9d7",
        workflow_id=35462565728,
        workflow_conclusion_at_freeze="success",
        independently_admitted=True,
        scope="source-compatible u0 formula and radius-one u ball",
        note="Independent PASS conditional on the already-admitted source-axis prerequisite.",
    ),
    "algebraic_r2_ordinary": AdmissionEvidence(
        key="algebraic_r2_ordinary",
        producer="Agent 4",
        pr=695,
        head="c38112a8d0172825119ef4f174b3b57af3a369ea",
        workflow_id=35465298813,
        workflow_conclusion_at_freeze="queued",
        independently_admitted=False,
        scope="four algebraic ordinary R2 slots from Agent-1 PR #684",
        note="Queued independent audit is not an admission receipt.",
    ),
    "derivative_r2_ordinary": AdmissionEvidence(
        key="derivative_r2_ordinary",
        producer="Agent 1 / pending Agent 4",
        pr=693,
        head="f2f79dc27259d35f0ca341d0a5fa03fe97d6c51e",
        workflow_id=35464666115,
        workflow_conclusion_at_freeze="queued",
        independently_admitted=False,
        scope="four derivative-bearing ordinary R2 slots after J1",
        note=(
            "Agent-1 code claims an 11/11 ordinary R2 self-certificate, but its exact-head "
            "dedicated workflow is queued and no independent Agent-4 audit of this subset exists."
        ),
    ),
}


REQUIRED_R2_ORDINARY_EVIDENCE = (
    "source_axis",
    "source_u_ball",
    "pressure_r2_ordinary",
    "algebraic_r2_ordinary",
    "derivative_r2_ordinary",
)


def derive_state(
    evidence: Mapping[str, AdmissionEvidence] = FROZEN_EVIDENCE,
) -> dict[str, bool]:
    """Derive only states justified by the supplied independent evidence.

    This intentionally stops before full-R2 promotion: the ordinary ledger is
    only one prerequisite and the source mixed seam / global leading handoff are
    separate scientific objects.
    """

    missing = [key for key in REQUIRED_R2_ORDINARY_EVIDENCE if key not in evidence]
    if missing:
        raise ValueError(f"missing required admission evidence: {missing}")

    ordinary_r2_ready = all(
        bool(evidence[key].independently_admitted)
        for key in REQUIRED_R2_ORDINARY_EVIDENCE
    )

    return {
        "source_axis_independently_admitted": evidence["source_axis"].independently_admitted,
        "source_u_ball_independently_admitted": evidence["source_u_ball"].independently_admitted,
        "source_r2_pressure_ordinary_slots_independently_admitted": evidence[
            "pressure_r2_ordinary"
        ].independently_admitted,
        "source_r2_algebraic_ordinary_slots_independently_admitted": evidence[
            "algebraic_r2_ordinary"
        ].independently_admitted,
        "source_r2_derivative_ordinary_slots_independently_admitted": evidence[
            "derivative_r2_ordinary"
        ].independently_admitted,
        "source_r2_all_ordinary_slots_independently_admitted": ordinary_r2_ready,
        "source_full_post_J1_R2_independently_admitted": False,
        "all_R1_ordinary_slots_source_bound": False,
        "source_full_post_J2_R1_machine_bound": False,
        "source_operator_MK_machine_bound": False,
        "source_B0_Tsh_machine_bound": False,
        "selected_PA16_handoff_allowed": False,
        "global_matched_pressure_materialized": False,
        "global_leading_velocity_materialized": False,
        "restricted_forcing_semantics_independently_validated": False,
        "leading_ready": False,
        "oscillatory_ready": True,
        "same_cycle_actual_defect_chain_present_on_a5_stack": True,
        "actual_mean_radial_finite_head_chain_present_on_a5_stack": True,
        "real_full_candidate_defect_consumed": False,
        "correction_ready": False,
        "kokuno_candidate_artifact_materialized": False,
        "velocity_export_ready": False,
        "formal_heldout_pde_gate_assessed": False,
        "pde_validated": False,
    }


def _canonical(payload: dict[str, Any]) -> bytes:
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    ).encode("utf-8")


def deterministic_receipt(*, exact_head: str | None = None) -> dict[str, Any]:
    state = derive_state()
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "evidence": {key: asdict(value) for key, value in sorted(FROZEN_EVIDENCE.items())},
        "derived_state": state,
        "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        "truth_boundary": {
            "integration_state_registry_only": True,
            "queued_workflow_treated_as_success": False,
            "agent1_self_certificate_treated_as_independent_admission": False,
            "surrogate_data_promoted_to_real_candidate": False,
            "free_residual_defined_forcing_allowed": False,
            "threshold_relaxed": False,
            "kokuno_replay_used_as_final_independent_pde_validation": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "leading_ready": state["leading_ready"],
            "correction_ready": state["correction_ready"],
            "velocity_export_ready": state["velocity_export_ready"],
            "pde_validated": state["pde_validated"],
        },
        "shortest_legal_next_chain": [
            "resolve Agent-4 #695 exact-head algebraic-R2 audit",
            "obtain independent Agent-4 audit for Agent-1 #693 derivative-R2 subset",
            "bind and independently admit all 11 R1 ordinary slots",
            "materialize source full post-J R1/R2 then M/K, B0/T_sh and PA.16",
            "materialize global leading velocity + matched pressure + independently restricted forcing",
            "consume real composite through existing A5 actual-defect/mean/radial/finite-head chain",
            "materialize signed correction velocity and guarded finite correction cycle",
            "serialize candidate and run independent Agent-4 held-out normalized NS validation",
        ],
    }
    receipt["receipt_sha256"] = hashlib.sha256(_canonical(receipt)).hexdigest()
    return receipt


def validate_receipt(receipt: Mapping[str, Any]) -> None:
    """Fail closed on accidental state/gate promotion in the frozen receipt."""

    if receipt.get("schema") != SCHEMA:
        raise ValueError("unexpected admission-firewall schema")
    evidence = receipt["evidence"]
    for key in ("source_axis", "pressure_r2_ordinary", "source_u_ball"):
        if evidence[key]["independently_admitted"] is not True:
            raise RuntimeError(f"green prerequisite lost admission: {key}")
    for key in ("algebraic_r2_ordinary", "derivative_r2_ordinary"):
        if evidence[key]["independently_admitted"] is not False:
            raise RuntimeError(f"pending prerequisite was promoted without receipt: {key}")

    state = receipt["derived_state"]
    forbidden_true = (
        "source_r2_all_ordinary_slots_independently_admitted",
        "source_full_post_J1_R2_independently_admitted",
        "all_R1_ordinary_slots_source_bound",
        "source_operator_MK_machine_bound",
        "global_matched_pressure_materialized",
        "global_leading_velocity_materialized",
        "restricted_forcing_semantics_independently_validated",
        "leading_ready",
        "real_full_candidate_defect_consumed",
        "correction_ready",
        "kokuno_candidate_artifact_materialized",
        "velocity_export_ready",
        "formal_heldout_pde_gate_assessed",
        "pde_validated",
    )
    promoted = [name for name in forbidden_true if state[name] is True]
    if promoted:
        raise RuntimeError(f"fail-closed states were promoted: {promoted}")

    gates = receipt["final_project_gates_unchanged"]
    if gates != FINAL_PROJECT_GATES:
        raise RuntimeError("final project gates changed")
    truth = receipt["truth_boundary"]
    if truth["free_residual_defined_forcing_allowed"] is not False:
        raise RuntimeError("free residual-defined forcing must remain forbidden")
    if truth["threshold_relaxed"] is not False:
        raise RuntimeError("project thresholds must not be relaxed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--exact-head")
    args = parser.parse_args()
    receipt = deterministic_receipt(exact_head=args.exact_head)
    validate_receipt(receipt)
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
