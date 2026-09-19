"""Agent-5 source-leading admission firewall refresh after Agent-4 PR #704.

This module performs one narrow integration-state promotion: the four
source-compatible derivative-bearing ordinary R2 slots are admitted because
Agent 4's exact-head independent audit completed successfully.  Every broader
state remains fail-closed.  In particular, full ordinary R2 still depends on
Agent-4 PR #695, and full post-J1 R2, R1, global leading fields, correction,
export and PDE validation remain false.

This is an integration-state receipt only.  It does not evaluate a
Navier--Stokes residual, alter a candidate, relax a threshold, or permit a
residual-defined forcing.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_source_leading_admission_firewall import (
    AdmissionEvidence,
    FINAL_PROJECT_GATES,
    FROZEN_EVIDENCE as V1_FROZEN_EVIDENCE,
    derive_state as _derive_state_v1,
)

SCHEMA = "kokuno-a5-source-leading-admission-firewall-v2"
TASK = "KOKUNO-A5-ADMIT-INDEPENDENT-DERIVATIVE-R2-ORDINARY-060"
PARENT_A5_PR = 705
PARENT_A5_HEAD = "f151bbea2937871d4bba30e2b58b5d460ac90f8c"

A4_DERIVATIVE_R2_RECEIPT = {
    "pr": 704,
    "head": "a96bef02c0cb5c2fdc292d4dc340e5eb464f8d42",
    "workflow_id": 35468474402,
    "workflow_name": "kokuno-agent4-pa10-source-derivative-slots-independent-audit",
    "workflow_conclusion": "success",
    "artifact_name": "kokuno-agent4-pa10-source-derivative-slots-independent-audit-v1",
    "artifact_id": 10592174398,
    "artifact_digest": "sha256:deec2ddd2b04f7a6d8dcfd4eba0d743482c90055d5062149891f43d64eed3dfb",
}

FROZEN_EVIDENCE: dict[str, AdmissionEvidence] = dict(V1_FROZEN_EVIDENCE)
FROZEN_EVIDENCE["derivative_r2_ordinary"] = AdmissionEvidence(
    key="derivative_r2_ordinary",
    producer="Agent 4",
    pr=704,
    head=A4_DERIVATIVE_R2_RECEIPT["head"],
    workflow_id=A4_DERIVATIVE_R2_RECEIPT["workflow_id"],
    workflow_conclusion_at_freeze="success",
    independently_admitted=True,
    scope="four derivative-bearing ordinary R2 slots from Agent-1 PR #693",
    note=(
        "Independent Agent-4 exact-head PASS only for the derivative-bearing "
        "ordinary R2 subset; full ordinary R2 remains blocked by PR #695."
    ),
)


def derive_state(
    evidence: Mapping[str, AdmissionEvidence] = FROZEN_EVIDENCE,
) -> dict[str, bool]:
    """Reuse the v1 fail-closed state derivation with refreshed evidence."""

    return _derive_state_v1(evidence)


def _canonical(payload: dict[str, Any]) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


def deterministic_receipt(*, exact_head: str | None = None) -> dict[str, Any]:
    state = derive_state()
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "new_independent_receipt": dict(A4_DERIVATIVE_R2_RECEIPT),
        "evidence": {
            key: asdict(value) for key, value in sorted(FROZEN_EVIDENCE.items())
        },
        "derived_state": state,
        "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        "truth_boundary": {
            "single_promotion": "source_r2_derivative_ordinary_slots_independently_admitted",
            "agent4_independent_receipt_required": True,
            "queued_workflow_treated_as_success": False,
            "agent1_self_certificate_treated_as_independent_admission": False,
            "surrogate_data_promoted_to_real_candidate": False,
            "free_residual_defined_forcing_allowed": False,
            "threshold_relaxed": False,
            "kokuno_replay_used_as_final_independent_pde_validation": False,
            "paper_exact": False,
            "leading_ready": state["leading_ready"],
            "correction_ready": state["correction_ready"],
            "velocity_export_ready": state["velocity_export_ready"],
            "pde_validated": state["pde_validated"],
        },
        "shortest_legal_next_chain": [
            "resolve Agent-4 PR #695 exact-head algebraic-R2 audit",
            "only then admit all ordinary R2 slots",
            "bind and independently admit all 11 R1 ordinary slots",
            "materialize source full post-J R1/R2 then M/K, B0/T_sh and PA.16",
            "materialize global leading velocity + matched pressure + independently restricted forcing",
            "consume the real composite through the existing A5 defect/mean/radial/debt/inverse spine",
            "materialize signed correction velocity and guarded finite correction cycle",
            "serialize the candidate and run Agent-4 held-out normalized NS validation",
        ],
    }
    receipt["receipt_sha256"] = hashlib.sha256(_canonical(receipt)).hexdigest()
    return receipt


def validate_receipt(receipt: Mapping[str, Any]) -> None:
    """Fail closed unless the one allowed promotion is exactly evidenced."""

    if receipt.get("schema") != SCHEMA:
        raise ValueError("unexpected admission-firewall schema")

    independent = receipt["new_independent_receipt"]
    if independent != A4_DERIVATIVE_R2_RECEIPT:
        raise RuntimeError("Agent-4 derivative-R2 receipt changed")
    if independent["workflow_conclusion"] != "success":
        raise RuntimeError("Agent-4 derivative-R2 workflow is not successful")

    evidence = receipt["evidence"]
    derivative = evidence["derivative_r2_ordinary"]
    if derivative["producer"] != "Agent 4" or derivative["pr"] != 704:
        raise RuntimeError("derivative-R2 evidence is not the frozen Agent-4 audit")
    if derivative["workflow_id"] != 35468474402:
        raise RuntimeError("unexpected derivative-R2 workflow id")
    if derivative["workflow_conclusion_at_freeze"] != "success":
        raise RuntimeError("derivative-R2 audit was not frozen as success")
    if derivative["independently_admitted"] is not True:
        raise RuntimeError("green derivative-R2 audit lost admission")

    algebraic = evidence["algebraic_r2_ordinary"]
    if algebraic["independently_admitted"] is not False:
        raise RuntimeError("pending algebraic-R2 audit was promoted without a receipt")

    state = receipt["derived_state"]
    if state["source_r2_derivative_ordinary_slots_independently_admitted"] is not True:
        raise RuntimeError("derivative-R2 state was not promoted")
    forbidden_true = (
        "source_r2_all_ordinary_slots_independently_admitted",
        "source_full_post_J1_R2_independently_admitted",
        "all_R1_ordinary_slots_source_bound",
        "source_full_post_J2_R1_machine_bound",
        "source_operator_MK_machine_bound",
        "source_B0_Tsh_machine_bound",
        "selected_PA16_handoff_allowed",
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

    if receipt["final_project_gates_unchanged"] != FINAL_PROJECT_GATES:
        raise RuntimeError("final project gates changed")
    truth = receipt["truth_boundary"]
    if truth["free_residual_defined_forcing_allowed"] is not False:
        raise RuntimeError("free residual-defined forcing must remain forbidden")
    if truth["threshold_relaxed"] is not False:
        raise RuntimeError("project thresholds must not be relaxed")
    if truth["pde_validated"] is not False:
        raise RuntimeError("this admission refresh cannot validate the PDE")


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
