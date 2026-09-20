"""Agent-5 admission refresh for the independently audited Wstar-logradial R1 slot.

This module performs one narrow integration-state promotion: Agent 4 PR #727
completed an exact-PR-head independent audit of Agent 1 PR #724's
``W_* Y Phi_Y`` ordinary R1 post-J2 slot.

Only that slot is newly admitted here. Previously admitted ordinary-R2, zeta-R1
and averaged-eta-R1 state stays admitted. Averaged-logradial, Hstar-Phi_eta,
d-u-Phi_eta, earlier Agent-1-only R1 slots, full post-J2 R1, full post-J1 R2,
M/K, global leading fields, correction readiness, export, and PDE validation
remain fail-closed.

No Navier--Stokes residual is evaluated here. This receipt changes no candidate,
pressure, forcing, or scientific threshold.
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
from .kokuno_a5_source_leading_admission_firewall_v4 import (
    FROZEN_EVIDENCE as V4_FROZEN_EVIDENCE,
    ST006_BASELINE,
    derive_state as _derive_state_v4,
)

SCHEMA = "kokuno-a5-source-leading-admission-firewall-v5"
TASK = "KOKUNO-A5-ADMIT-INDEPENDENT-R1-WSTAR-LOGRADIAL-SLOT-065"
PARENT_A5_PR = 745
PARENT_A5_HEAD = "069fdd19c7c6a7d63117c84f8db9f486a74c213d"

A4_WSTAR_R1_RECEIPT = {
    "pr": 727,
    "head": "279de4c794935bcb3331627046be3c96b865bac7",
    "workflow_id": 35477137066,
    "workflow_name": "kokuno-agent4-pa10-source-r1-wstar-logradial-independent-audit",
    "workflow_event": "pull_request",
    "workflow_conclusion": "success",
    "artifact_name": "kokuno-agent4-pa10-source-r1-wstar-logradial-independent-audit-v1",
    "artifact_id": 10596193830,
    "artifact_digest": "sha256:08b514794efa085297ce9ab14d6e8eebc0e9d737b1a753e3ffd2dfd03baaa2b7",
    "report_sha256": "76488aa1440e9651ded65480a529eac7b1fed86ea299e698eb5e6d8c4a315c6e",
    "receipt_sha256": "f0816da8a9c1264ef6822ed99034bfed9a15e6a5d1200fae1d9b378426da60ec",
    "failed_guards": [],
    "independent_preflight_passed": True,
}

AGENT1_R1_PROGRESS = {
    "latest_pr": 751,
    "latest_head": "080090f32214a9da102542ce4e927ec62d028baf",
    "ordinary_r1_source_bound_slots": 11,
    "ordinary_r1_total_slots": 11,
    "ordinary_r1_remaining_slots": 0,
    "remaining_slot_names": [],
    "agent1_full_post_J2_R1_self_envelope_materialized": True,
    "agent1_self_certificate_is_not_independent_admission": True,
    "agent1_dedicated_workflow_id": 35484022718,
    "agent1_dedicated_workflow_conclusion_at_freeze": "queued",
}

PENDING_A4_R1_AUDITS = {
    "averaged_logradial": {
        "pr": 736,
        "head": "1e5b66dba633a868fc0b898178d712f3d640f022",
        "workflow_id": 35479780213,
        "workflow_conclusion_at_freeze": "queued",
    },
    "hstar_phi_eta": {
        "pr": 744,
        "head": "a45b9964f8bb36355c0a3bca999aeb2bee163c40",
        "workflow_id": 35482164500,
        "workflow_conclusion_at_freeze": "queued",
    },
    "du_phi_eta": {
        "pr": 755,
        "head": "fe8e12f194770f0751ad10ca9f9edf8f3588cb52",
        "workflow_id": 35484899251,
        "workflow_conclusion_at_freeze": "queued",
    },
    "earlier_algebraic_r1_slots": {
        "independent_agent4_coverage_complete_at_freeze": False,
        "note": (
            "Existing slot-scoped Agent-4 receipts do not independently admit every "
            "earlier Agent-1-only ordinary R1 slot."
        ),
    },
}

PARENT_A5_CI_AT_FREEZE = {
    "pr": 745,
    "head": PARENT_A5_HEAD,
    "dedicated_workflow_id": 35482609891,
    "dedicated_workflow_conclusion_at_freeze": "queued",
    "repository_tests_workflow_id": 35482609955,
    "repository_tests_conclusion_at_freeze": "queued",
    "parent_queued_ci_treated_as_success": False,
}

FROZEN_EVIDENCE: dict[str, AdmissionEvidence] = dict(V4_FROZEN_EVIDENCE)
FROZEN_EVIDENCE["r1_wstar_logradial_ordinary"] = AdmissionEvidence(
    key="r1_wstar_logradial_ordinary",
    producer="Agent 4",
    pr=727,
    head=A4_WSTAR_R1_RECEIPT["head"],
    workflow_id=A4_WSTAR_R1_RECEIPT["workflow_id"],
    workflow_conclusion_at_freeze="success",
    independently_admitted=True,
    scope="Agent-1 PR #724 W_* Y Phi_Y ordinary R1 post-J2 slot",
    note=(
        "Exact PR-head Agent-4 dedicated PASS for this Wstar-logradial "
        "ordinary-R1 slot only; it does not admit other R1 slots or full R1."
    ),
)


def derive_state(
    evidence: Mapping[str, AdmissionEvidence] = FROZEN_EVIDENCE,
) -> dict[str, bool]:
    """Extend v4 with exactly one independently admitted R1 slot."""

    state = dict(_derive_state_v4(evidence))
    wstar = evidence.get("r1_wstar_logradial_ordinary")
    wstar_admitted = bool(wstar is not None and wstar.independently_admitted)
    state.update(
        {
            "source_r1_wstar_logradial_ordinary_slot_independently_admitted":
                wstar_admitted,
            "source_r1_averaged_logradial_ordinary_slot_independently_admitted":
                False,
            "source_r1_hstar_phi_eta_ordinary_slot_independently_admitted":
                False,
            "source_r1_du_phi_eta_ordinary_slot_independently_admitted":
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
        "parent_a5_ci_at_freeze": dict(PARENT_A5_CI_AT_FREEZE),
        "new_independent_receipt": dict(A4_WSTAR_R1_RECEIPT),
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
                "global Kokuno velocity/pressure/restricted-forcing candidate exists."
            ),
        },
        "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        "truth_boundary": {
            "single_integration_increment":
                "admit independently audited R1 Wstar-logradial ordinary slot",
            "agent4_exact_pr_head_receipt_required": True,
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
            "resolve Agent-4 #736 averaged-logradial R1 exact-head audit",
            "resolve Agent-4 #744 Hstar-Phi_eta R1 exact-head audit",
            "resolve Agent-4 #755 d-u-Phi_eta R1 exact-head audit",
            "independently audit/admit earlier Agent-1-only ordinary R1 slots still uncovered",
            "independently admit full post-J2 R1 and full post-J1 R2 without laundering slot receipts",
            "materialize M/K, B0/T_sh and PA.16",
            "materialize global leading velocity + matched pressure + independently restricted forcing",
            "consume the real composite through the existing A5 defect/mean/radial/debt/inverse spine",
            "materialize signed correction velocity and run guarded finite correction cycles",
            "serialize candidate, smoke-test Python/MATLAB load/use, then run Agent-4 held-out normalized NS validation",
        ],
    }
    receipt["receipt_sha256"] = hashlib.sha256(_canonical(receipt)).hexdigest()
    return receipt


def validate_receipt(receipt: Mapping[str, Any]) -> None:
    """Fail closed unless only the Wstar-logradial slot was newly admitted."""

    if receipt.get("schema") != SCHEMA:
        raise ValueError("unexpected admission-firewall schema")

    independent = receipt["new_independent_receipt"]
    if independent != A4_WSTAR_R1_RECEIPT:
        raise RuntimeError("Agent-4 Wstar-logradial receipt changed")
    if independent["head"] != "279de4c794935bcb3331627046be3c96b865bac7":
        raise RuntimeError("unexpected Agent-4 Wstar-logradial head")
    if independent["workflow_id"] != 35477137066:
        raise RuntimeError("unexpected Agent-4 Wstar-logradial workflow id")
    if independent["workflow_conclusion"] != "success":
        raise RuntimeError("Agent-4 Wstar-logradial workflow is not successful")
    if independent["artifact_id"] != 10596193830:
        raise RuntimeError("unexpected Agent-4 Wstar-logradial artifact id")
    if independent["failed_guards"] != []:
        raise RuntimeError("Agent-4 Wstar-logradial receipt has failed guards")
    if independent["independent_preflight_passed"] is not True:
        raise RuntimeError("Agent-4 Wstar-logradial preflight did not pass")

    evidence = receipt["evidence"]["r1_wstar_logradial_ordinary"]
    if evidence["producer"] != "Agent 4" or evidence["pr"] != 727:
        raise RuntimeError("Wstar-logradial evidence is not the frozen Agent-4 audit")
    if evidence["workflow_id"] != 35477137066:
        raise RuntimeError("unexpected Wstar-logradial workflow id")
    if evidence["workflow_conclusion_at_freeze"] != "success":
        raise RuntimeError("Wstar-logradial audit was not frozen as success")
    if evidence["independently_admitted"] is not True:
        raise RuntimeError("green Wstar-logradial audit lost admission")

    state = receipt["derived_state"]
    required_true = (
        "source_r2_all_ordinary_slots_independently_admitted",
        "source_r1_zeta_ordinary_slot_independently_admitted",
        "source_r1_averaged_eta_ordinary_slot_independently_admitted",
        "source_r1_wstar_logradial_ordinary_slot_independently_admitted",
    )
    missing = [name for name in required_true if state[name] is not True]
    if missing:
        raise RuntimeError(f"previously/newly admitted state missing: {missing}")

    forbidden_true = (
        "source_r1_averaged_logradial_ordinary_slot_independently_admitted",
        "source_r1_hstar_phi_eta_ordinary_slot_independently_admitted",
        "source_r1_du_phi_eta_ordinary_slot_independently_admitted",
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
    if progress["ordinary_r1_source_bound_slots"] != 11:
        raise RuntimeError("Agent-1 R1 self-ledger snapshot is not 11/11")
    if progress["ordinary_r1_remaining_slots"] != 0:
        raise RuntimeError("Agent-1 R1 self-ledger remaining count changed")
    if progress["agent1_full_post_J2_R1_self_envelope_materialized"] is not True:
        raise RuntimeError("Agent-1 full post-J2 self-envelope snapshot missing")
    if progress["agent1_self_certificate_is_not_independent_admission"] is not True:
        raise RuntimeError("Agent-1 self-certificate was promoted to admission")
    if progress["agent1_dedicated_workflow_conclusion_at_freeze"] != "queued":
        raise RuntimeError("Agent-1 queued CI freeze fact changed")

    parent_ci = receipt["parent_a5_ci_at_freeze"]
    if parent_ci["parent_queued_ci_treated_as_success"] is not False:
        raise RuntimeError("queued parent A5 CI was promoted")

    pending = receipt["pending_agent4_r1_audits"]
    for key in ("averaged_logradial", "hstar_phi_eta", "du_phi_eta"):
        if pending[key]["workflow_conclusion_at_freeze"] != "queued":
            raise RuntimeError(f"pending Agent-4 audit {key} changed freeze state")
    if pending["earlier_algebraic_r1_slots"][
        "independent_agent4_coverage_complete_at_freeze"
    ] is not False:
        raise RuntimeError("earlier R1 independent coverage was fabricated")

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
