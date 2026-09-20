"""Agent-5 admission refresh for the independently audited averaged-eta R1 slot.

This module performs one narrow integration-state promotion: Agent 4 PR #718
completed an exact-head push audit of Agent 1 PR #715's
``Lambda^-1 d (partial_eta A(u)) Phi`` ordinary R1 post-J2 slot.

Only that slot is newly admitted here.  The previously admitted zeta slot stays
admitted.  Wstar-logradial, the remaining/prior ordinary R1 slots, full post-J2
R1, full post-J1 R2, M/K, global leading fields, correction readiness, export,
and PDE validation remain fail-closed.

No Navier--Stokes residual is evaluated here.  This receipt does not alter a
candidate, pressure, forcing, or any scientific threshold.
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
)
from .kokuno_a5_source_leading_admission_firewall_v3 import (
    FROZEN_EVIDENCE as V3_FROZEN_EVIDENCE,
    ST006_BASELINE,
    derive_state as _derive_state_v3,
)

SCHEMA = "kokuno-a5-source-leading-admission-firewall-v4"
TASK = "KOKUNO-A5-ADMIT-INDEPENDENT-R1-AVERAGED-ETA-SLOT-063"
PARENT_A5_PR = 728
PARENT_A5_HEAD = "ce4ebb28089546819074cf6e32431f550cebce27"

A4_AVERAGED_ETA_R1_RECEIPT = {
    "pr": 718,
    "head": "c4515a3c7803b7fa28d06f5b546cb943804dcd35",
    "workflow_id": 35474229697,
    "workflow_name": "kokuno-agent4-pa10-source-r1-averaged-eta-independent-audit",
    "workflow_event": "push",
    "workflow_conclusion": "success",
    "artifact_name": "kokuno-agent4-pa10-source-r1-averaged-eta-independent-audit-v1",
    "artifact_id": 10595470641,
    "artifact_digest": "sha256:65907488f0cab2fe431e6e0786152ee13e5656a227009e7b1d05df2bee6df158",
    "receipt_sha256": "1cd9906fb18cee52f95d658a690ff1a1fc51a25343662ca03b3091b0694046e1",
}

AGENT1_R1_PROGRESS = {
    "latest_pr": 733,
    "latest_head": "e2333474234e0f509d6cfab04b94b01f5a952b04",
    "ordinary_r1_source_bound_slots": 9,
    "ordinary_r1_total_slots": 11,
    "ordinary_r1_remaining_slots": 2,
    "remaining_slot_names": [
        "H_* Phi_eta",
        "Lambda^-1 d u Phi_eta",
    ],
    "agent1_self_certificate_is_not_independent_admission": True,
}

PENDING_A4_R1_AUDITS = {
    "wstar_logradial": {
        "pr": 727,
        "head": "279de4c794935bcb3331627046be3c96b865bac7",
        "workflow_id": 35477137066,
        "workflow_conclusion_at_freeze": "queued",
    },
    "averaged_logradial": {
        "agent1_pr": 733,
        "agent1_head": "e2333474234e0f509d6cfab04b94b01f5a952b04",
        "independent_agent4_audit_available_at_freeze": False,
    },
}

FROZEN_EVIDENCE: dict[str, AdmissionEvidence] = dict(V3_FROZEN_EVIDENCE)
FROZEN_EVIDENCE["r1_averaged_eta_ordinary"] = AdmissionEvidence(
    key="r1_averaged_eta_ordinary",
    producer="Agent 4",
    pr=718,
    head=A4_AVERAGED_ETA_R1_RECEIPT["head"],
    workflow_id=A4_AVERAGED_ETA_R1_RECEIPT["workflow_id"],
    workflow_conclusion_at_freeze="success",
    independently_admitted=True,
    scope=(
        "Agent-1 PR #715 Lambda^-1 d (partial_eta A(u)) Phi "
        "ordinary R1 post-J2 slot"
    ),
    note=(
        "Exact-head Agent-4 push workflow PASS only for this averaged-eta "
        "ordinary-R1 slot; it does not retroactively admit other R1 slots "
        "or full R1."
    ),
)


def derive_state(
    evidence: Mapping[str, AdmissionEvidence] = FROZEN_EVIDENCE,
) -> dict[str, bool]:
    """Extend the v3 fail-closed state with one narrow R1 admission."""

    state = dict(_derive_state_v3(evidence))
    averaged_eta = evidence.get("r1_averaged_eta_ordinary")
    averaged_eta_admitted = bool(
        averaged_eta is not None and averaged_eta.independently_admitted
    )
    state.update(
        {
            "source_r1_averaged_eta_ordinary_slot_independently_admitted":
                averaged_eta_admitted,
            "source_r1_wstar_logradial_ordinary_slot_independently_admitted":
                False,
            "source_r1_averaged_logradial_ordinary_slot_independently_admitted":
                False,
            "source_r1_all_ordinary_slots_independently_admitted": False,
        }
    )
    return state


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
        "new_independent_receipt": dict(A4_AVERAGED_ETA_R1_RECEIPT),
        "evidence": {
            key: asdict(value) for key, value in sorted(FROZEN_EVIDENCE.items())
        },
        "upstream_agent1_r1_progress": dict(AGENT1_R1_PROGRESS),
        "pending_agent4_r1_audits": {
            key: dict(value) for key, value in PENDING_A4_R1_AUDITS.items()
        },
        "derived_state": state,
        "baseline_vs_kokuno": {
            "retained_repository_baseline": dict(ST006_BASELINE),
            "kokuno_same_protocol_full_candidate_residual_available": False,
            "kokuno_same_protocol_comparison_performed": False,
            "reason": (
                "This increment admits one source/operator slot only; no complete "
                "global Kokuno velocity/pressure/forcing candidate exists yet."
            ),
        },
        "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        "truth_boundary": {
            "single_integration_increment":
                "admit independently audited R1 averaged-eta ordinary slot",
            "agent4_exact_head_push_receipt_required": True,
            "queued_workflow_treated_as_success": False,
            "agent1_self_certificate_treated_as_independent_admission": False,
            "prior_r1_slots_retroactively_admitted": False,
            "full_r1_promoted": False,
            "full_r2_promoted": False,
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
            "resolve Agent-4 #727 Wstar-logradial exact-head audit",
            "independently audit Agent-1 #733 averaged-logradial R1 slot",
            "bind Agent-1's remaining 2/11 ordinary R1 slots and independently audit them",
            "obtain independent admission for earlier ordinary-R1 slots not yet covered",
            "materialize source full post-J R1/R2 then M/K, B0/T_sh and PA.16",
            "materialize global leading velocity + matched pressure + independently restricted forcing",
            "consume the real composite through the existing A5 defect/mean/radial/debt/inverse spine",
            "materialize signed correction velocity and guarded finite correction cycle",
            "serialize candidate and run Agent-4 held-out normalized NS validation",
        ],
    }
    receipt["receipt_sha256"] = hashlib.sha256(_canonical(receipt)).hexdigest()
    return receipt


def validate_receipt(receipt: Mapping[str, Any]) -> None:
    """Fail closed unless the exact narrow Agent-4 admission is preserved."""

    if receipt.get("schema") != SCHEMA:
        raise ValueError("unexpected admission-firewall schema")

    independent = receipt["new_independent_receipt"]
    if independent != A4_AVERAGED_ETA_R1_RECEIPT:
        raise RuntimeError("Agent-4 R1 averaged-eta receipt changed")
    if independent["workflow_event"] != "push":
        raise RuntimeError("R1 averaged-eta evidence is not an exact-head push workflow")
    if independent["head"] != "c4515a3c7803b7fa28d06f5b546cb943804dcd35":
        raise RuntimeError("unexpected Agent-4 R1 averaged-eta head")
    if independent["workflow_conclusion"] != "success":
        raise RuntimeError("Agent-4 R1 averaged-eta workflow is not successful")

    evidence = receipt["evidence"]["r1_averaged_eta_ordinary"]
    if evidence["producer"] != "Agent 4" or evidence["pr"] != 718:
        raise RuntimeError("R1 averaged-eta evidence is not the frozen Agent-4 audit")
    if evidence["workflow_id"] != 35474229697:
        raise RuntimeError("unexpected R1 averaged-eta workflow id")
    if evidence["workflow_conclusion_at_freeze"] != "success":
        raise RuntimeError("R1 averaged-eta audit was not frozen as success")
    if evidence["independently_admitted"] is not True:
        raise RuntimeError("green R1 averaged-eta audit lost admission")

    state = receipt["derived_state"]
    if state["source_r2_all_ordinary_slots_independently_admitted"] is not True:
        raise RuntimeError("previously admitted ordinary-R2 state regressed")
    if state["source_r1_zeta_ordinary_slot_independently_admitted"] is not True:
        raise RuntimeError("previously admitted R1-zeta state regressed")
    if state[
        "source_r1_averaged_eta_ordinary_slot_independently_admitted"
    ] is not True:
        raise RuntimeError("R1 averaged-eta admission missing")

    forbidden_true = (
        "source_r1_wstar_logradial_ordinary_slot_independently_admitted",
        "source_r1_averaged_logradial_ordinary_slot_independently_admitted",
        "source_r1_all_ordinary_slots_independently_admitted",
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
    baseline = receipt["baseline_vs_kokuno"]
    if baseline["retained_repository_baseline"] != ST006_BASELINE:
        raise RuntimeError("ST006 baseline changed")
    if baseline["kokuno_same_protocol_full_candidate_residual_available"] is not False:
        raise RuntimeError("a nonexistent full Kokuno residual was promoted")

    progress = receipt["upstream_agent1_r1_progress"]
    if progress["ordinary_r1_source_bound_slots"] != 9:
        raise RuntimeError("Agent-1 R1 progress snapshot changed")
    if progress["ordinary_r1_remaining_slots"] != 2:
        raise RuntimeError("Agent-1 R1 remaining-slot count changed")
    if progress["agent1_self_certificate_is_not_independent_admission"] is not True:
        raise RuntimeError("Agent-1 self-certificate was promoted to admission")

    truth = receipt["truth_boundary"]
    if truth["prior_r1_slots_retroactively_admitted"] is not False:
        raise RuntimeError("prior R1 slots cannot be retroactively admitted")
    if truth["free_residual_defined_forcing_allowed"] is not False:
        raise RuntimeError("free residual-defined forcing must remain forbidden")
    if truth["threshold_relaxed"] is not False:
        raise RuntimeError("project thresholds must remain frozen")
    if truth["pde_validated"] is not False:
        raise RuntimeError("this source admission cannot validate the PDE")


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
