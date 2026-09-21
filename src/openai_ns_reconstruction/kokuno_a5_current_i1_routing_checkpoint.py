"""Agent-5 routing/artifact checkpoint for the current modulated Kokuno I1 lineage.

This module is integration glue only.  It consumes the implementation-distinct
Agent-4 #1053 receipt for the exact Agent-1 #1051 current I1 candidate, records
the latest Agent-2 #1052 and Agent-3 #1054 deliveries as *non-consumed sibling
state*, and fails closed against cross-lineage evidence transfer from the older
unmodulated X4 correction route (#1045/#1046/#1047).

It does not create pressure, forcing, an NS defect, a correction velocity, or a
finite correction cycle.  A scoped divergence PASS is not a complete momentum
residual and cannot promote PDE validity.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a4_current_i1_leading_divergence_independent_audit as a4

SCHEMA = "kokuno-a5-current-i1-routing-checkpoint-v1"
TASK = "KOKUNO-A5-CURRENT-I1-ROUTING-CHECKPOINT-106"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")

AGENT1_CURRENT_I1 = {
    "pr": 1051,
    "head": "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5",
    "branch": "kokuno-agent1/freeze-i1-closure-gate-112",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_i1_frozen_gate.py",
    "source_blob": "fe14dab8dccdbca828d722bdfe2f75f287aa9ea8",
    "frozen_i1_closure_gate": 5.0e-7,
    "public_velocity_api": "velocity(x,y,z,t)",
    "through_i1": True,
    "through_i2": False,
    "global_leading": False,
}

AGENT4_CURRENT_I1_AUDIT = {
    "pr": 1053,
    "head": "ad9f8722b6ec89df975f9b2402e63066dc14422d",
    "branch": "agent-kokuno-4/current-i1-leading-divergence-108",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_i1_leading_divergence_independent_audit.py",
    "source_blob": "cc8c8f9ed0e8f727035c14232bb7f80d2703ac76",
    "audited_agent1_pr": 1051,
    "audited_agent1_head": "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5",
    "operator": "centered Cartesian FD2 from public velocity only",
    "spatial_steps": [0.02, 0.01, 0.005],
    "seed": 9173801,
    "implementation_distinct": True,
    "scope": "leading-only divergence through current I1",
    "complete_ns_residual_evidence": False,
}

AGENT2_LATEST_SIBLING = {
    "pr": 1052,
    "head": "d58cee2bf38bdca8da13a796215a079a569ad7d7",
    "branch": "codex/kokuno-a2-source-amplitude-sensitivity-094",
    "source_path": "src/openai_ns_reconstruction/kokuno_source_zero_data_amplitude_sensitivity.py",
    "source_blob": "eba3c00703bd5763b12e6f48a8c34eb918dec82e",
    "delivers": "zero-data t_m, D_r t_m, D_z t_m sensitivity route into existing complete-curl path",
    "caller_supplied_background_and_forcing_directional_jets_remain": True,
    "self_contained_cartesian_composite_on_agent1_1051": False,
    "consumed_by_this_checkpoint": False,
}

AGENT3_LATEST_SIBLING = {
    "pr": 1054,
    "head": "8f2f485e2000955bab47fc08f478b512d3b49e5d",
    "branch": "agent3/complete-ns-defect-cycle-admission-20260922",
    "source_path": "src/openai_ns_reconstruction/kokuno_complete_ns_defect_cycle_admission.py",
    "source_blob": "d62ffa13f910e508e7b387b2e553fdc6b40c6ed8",
    "parent_unmodulated_x4_pr": 1045,
    "delivers": "fail-closed typed prerequisites for complete NS defect and real finite-cycle admission",
    "finite_cycle_admitted": False,
    "consumed_by_current_i1_identity": False,
}

PRIOR_A5_UNMODULATED_X4 = {
    "pr": 1047,
    "head": "1f6cf97c386080ea3547568a2743e708562a4b01",
    "branch": "codex/kokuno-a5-rf40-power-law-radial-force-ingest",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_rf40_power_law_radial_force_ingest_contract.py",
    "source_blob": "b8991b6c381620123e3f82c8eb78c16a3e8ced13",
    "lineage": "older unmodulated X4 correction-side mean->stress->radial-force route",
    "evidence_transfer_to_current_i1_allowed": False,
}

FINAL_GATES = {
    "normalized_momentum_sampled_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "normalized_divergence_sampled_max": 1.0e-5,
    "normalized_divergence_volume_l2": 1.0e-5,
    "canonical_quadrature": [24, 48, 96],
    "residual_defined_free_forcing_allowed": False,
}

ST006_BASELINE = {
    "role": "same-protocol complete-PDE baseline only",
    "momentum_sampled_max": 0.1082289305112118,
    "momentum_volume_l2": 0.10758432876230622,
    "pde_validated": False,
}

READINESS = {
    "leading_ready": False,
    "oscillatory_ready": True,
    "correction_ready": False,
    "velocity_export_ready": False,
    "pde_validated": False,
}

TRUTH_BOUNDARY = {
    "current_modulated_leading_through_i1_materialized": True,
    "current_i1_configuration_save_load_available": True,
    "current_i1_a4_independent_divergence_audit_registered": True,
    "current_i1_a4_audit_is_scoped_divergence_only": True,
    "current_i1_matching_leading_plus_oscillatory_composite_materialized": False,
    "agent2_1052_consumed_into_current_i1_composite": False,
    "agent3_1054_consumed_into_current_i1_correction": False,
    "older_unmodulated_x4_correction_evidence_transferred": False,
    "current_i2_applied": False,
    "current_i3_applied": False,
    "current_i4_applied": False,
    "outer_global_leading_velocity_materialized": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "preregistered_restricted_forcing_materialized": False,
    "forcing_proven_not_residual_defined": False,
    "complete_identity_bound_ns_defect_materialized": False,
    "current_cartesian_correction_velocity_materialized": False,
    "finite_correction_cycle_admitted": False,
    "real_candidate_finite_correction_cycle_run": False,
    "heldout_normalized_full_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "scientific_admission": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

PIPELINE_POSITION = {
    "current_stage": "exact current autonomous/modulated leading through I1 plus implementation-distinct scoped divergence audit",
    "current_identity": "A1 #1051 -> A4 #1053",
    "non_consumed_siblings": "A2 #1052 source-amplitude sensitivities; A3 #1054 complete-defect/cycle admission contract; A5 #1047 older unmodulated X4 correction route",
    "next_shortest_blocker": "advance exact #1051 leading through I2/I3/I4/global completion and recompose a matching A2 oscillatory Cartesian candidate on that same identity; only then bind matched pressure/restricted forcing, complete defect, and a current A3 correction cycle",
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _require_hex40(value: Any, label: str) -> None:
    if not isinstance(value, str) or _HEX40.fullmatch(value) is None:
        raise ValueError(f"{label} is not an exact 40-hex git identity")


def _require_hex64(value: Any, label: str) -> None:
    if not isinstance(value, str) or _HEX64.fullmatch(value) is None:
        raise ValueError(f"{label} is not a 64-hex semantic digest")


def validate_upstream_a4_receipt(receipt: Mapping[str, Any]) -> None:
    """Fail closed unless this is the exact preregistered #1053 scientific receipt."""
    if receipt.get("schema") != a4.SCHEMA:
        raise ValueError("unexpected A4 current-I1 audit schema")
    if receipt.get("upstream_pr") != AGENT1_CURRENT_I1["pr"]:
        raise ValueError("A4 receipt is detached from exact A1 #1051")
    if receipt.get("upstream_head") != AGENT1_CURRENT_I1["head"]:
        raise ValueError("A4 receipt upstream head drifted")
    _require_hex64(receipt.get("candidate_semantic_sha256"), "candidate_semantic_sha256")

    protocol = receipt.get("protocol")
    gates = receipt.get("gates")
    truth = receipt.get("truth_boundary")
    resolutions = receipt.get("resolutions")
    if not isinstance(protocol, Mapping) or not isinstance(gates, Mapping) or not isinstance(truth, Mapping):
        raise ValueError("A4 receipt is missing protocol/gates/truth boundary")
    if int(protocol.get("seed", -1)) != AGENT4_CURRENT_I1_AUDIT["seed"]:
        raise ValueError("A4 held-out seed drifted")
    if list(protocol.get("spatial_steps", [])) != AGENT4_CURRENT_I1_AUDIT["spatial_steps"]:
        raise ValueError("A4 spatial step ladder drifted")
    if protocol.get("derivative_operator") != AGENT4_CURRENT_I1_AUDIT["operator"]:
        raise ValueError("A4 independent derivative operator drifted")
    if float(protocol.get("frozen_i1_closure_tolerance", float("nan"))) != AGENT1_CURRENT_I1["frozen_i1_closure_gate"]:
        raise ValueError("A1 frozen I1 closure gate drifted in A4 receipt")
    if float(gates.get("final_project_momentum_gate_unchanged", float("nan"))) != 1.0e-3:
        raise ValueError("final momentum gate drifted")
    if float(gates.get("final_project_divergence_gate_unchanged", float("nan"))) != 1.0e-5:
        raise ValueError("final divergence gate drifted")
    if truth.get("current_i1_leading_velocity_consumed") is not True:
        raise ValueError("A4 did not consume current I1 leading velocity")
    for key in (
        "leading_only_ns_residual_assessed",
        "leading_plus_oscillatory_ns_residual_assessed",
        "after_correction_ns_residual_assessed",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "canonical_24_48_96_volume_admission_assessed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"scoped A4 receipt illegally promoted {key}")
    if not isinstance(resolutions, list) or len(resolutions) != 3:
        raise ValueError("A4 receipt must preserve the three-resolution ladder")
    if [float(row.get("step")) for row in resolutions] != AGENT4_CURRENT_I1_AUDIT["spatial_steps"]:
        raise ValueError("A4 receipt resolution steps drifted")
    if not isinstance(receipt.get("audit_pass"), bool):
        raise ValueError("A4 audit_pass must be a boolean")


def _selected_a4_evidence(receipt: Mapping[str, Any]) -> dict[str, Any]:
    fine = receipt["resolutions"][-1]
    return {
        "audit_schema": receipt["schema"],
        "candidate_semantic_sha256": receipt["candidate_semantic_sha256"],
        "audit_pass": bool(receipt["audit_pass"]),
        "receipt_sha256": _sha256(receipt),
        "finest_step": float(fine["step"]),
        "finest_sampled_max": float(fine["sampled_max"]),
        "finest_volume_l2_estimate": float(fine["pooled_volume_l2_estimate"]),
        "finest_weighted_rms": float(fine["pooled_weighted_rms"]),
        "finest_i1_entry_seam_sampled_max": float(fine["i1_entry_seam_sampled_max"]),
        "finest_axis_axis_near_sampled_max": float(fine["axis_axis_near_sampled_max"]),
        "finest_speed_rms": float(fine["speed_rms"]),
        "scientific_scope": "scoped leading-only divergence through current I1; not a complete NS momentum residual",
    }


def build_checkpoint(a4_receipt: Mapping[str, Any]) -> dict[str, Any]:
    validate_upstream_a4_receipt(a4_receipt)
    payload = {
        "agent1_current_i1": copy.deepcopy(AGENT1_CURRENT_I1),
        "agent4_current_i1_audit": copy.deepcopy(AGENT4_CURRENT_I1_AUDIT),
        "agent2_latest_sibling": copy.deepcopy(AGENT2_LATEST_SIBLING),
        "agent3_latest_sibling": copy.deepcopy(AGENT3_LATEST_SIBLING),
        "prior_a5_unmodulated_x4": copy.deepcopy(PRIOR_A5_UNMODULATED_X4),
        "a4_scoped_evidence": _selected_a4_evidence(a4_receipt),
        "final_gates": copy.deepcopy(FINAL_GATES),
        "st006_baseline": copy.deepcopy(ST006_BASELINE),
        "readiness": copy.deepcopy(READINESS),
        "truth_boundary": copy.deepcopy(TRUTH_BOUNDARY),
        "pipeline_position": copy.deepcopy(PIPELINE_POSITION),
    }
    checkpoint = {"schema": SCHEMA, "task": TASK, **payload}
    checkpoint["digest"] = _sha256(payload)
    validate_checkpoint(checkpoint)
    return checkpoint


def validate_checkpoint(checkpoint: Mapping[str, Any]) -> None:
    if checkpoint.get("schema") != SCHEMA or checkpoint.get("task") != TASK:
        raise ValueError("unexpected A5 current-I1 checkpoint identity")
    for name in (
        "agent1_current_i1",
        "agent4_current_i1_audit",
        "agent2_latest_sibling",
        "agent3_latest_sibling",
        "prior_a5_unmodulated_x4",
    ):
        obj = checkpoint.get(name)
        if not isinstance(obj, Mapping):
            raise ValueError(f"missing {name}")
        _require_hex40(obj.get("head"), f"{name}.head")
        _require_hex40(obj.get("source_blob"), f"{name}.source_blob")

    if checkpoint["agent1_current_i1"] != AGENT1_CURRENT_I1:
        raise ValueError("A1 current-I1 identity drifted")
    if checkpoint["agent4_current_i1_audit"] != AGENT4_CURRENT_I1_AUDIT:
        raise ValueError("A4 current-I1 audit identity drifted")
    if checkpoint["agent2_latest_sibling"] != AGENT2_LATEST_SIBLING:
        raise ValueError("A2 sibling identity/truth drifted")
    if checkpoint["agent3_latest_sibling"] != AGENT3_LATEST_SIBLING:
        raise ValueError("A3 sibling identity/truth drifted")
    if checkpoint["prior_a5_unmodulated_x4"] != PRIOR_A5_UNMODULATED_X4:
        raise ValueError("prior A5 X4 lineage identity/truth drifted")
    if checkpoint["agent4_current_i1_audit"]["audited_agent1_head"] != checkpoint["agent1_current_i1"]["head"]:
        raise ValueError("A4 is not bound to the exact current A1 identity")
    if checkpoint["agent2_latest_sibling"]["self_contained_cartesian_composite_on_agent1_1051"] is not False:
        raise ValueError("A2 #1052 cannot be laundered into a matching current-I1 composite")
    if checkpoint["agent3_latest_sibling"]["consumed_by_current_i1_identity"] is not False:
        raise ValueError("older unmodulated A3 correction evidence cannot transfer to current I1")
    if checkpoint["prior_a5_unmodulated_x4"]["evidence_transfer_to_current_i1_allowed"] is not False:
        raise ValueError("prior unmodulated X4 A5 evidence transfer must remain forbidden")

    if checkpoint.get("final_gates") != FINAL_GATES:
        raise ValueError("fixed project gates drifted")
    if checkpoint.get("st006_baseline") != ST006_BASELINE:
        raise ValueError("ST006 baseline drifted")
    if checkpoint.get("readiness") != READINESS:
        raise ValueError("core readiness state drifted")
    if checkpoint.get("truth_boundary") != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drifted")
    if checkpoint.get("pipeline_position") != PIPELINE_POSITION:
        raise ValueError("pipeline routing drifted")

    evidence = checkpoint.get("a4_scoped_evidence")
    if not isinstance(evidence, Mapping):
        raise ValueError("missing A4 scoped evidence")
    _require_hex64(evidence.get("candidate_semantic_sha256"), "a4_scoped_evidence.candidate_semantic_sha256")
    _require_hex64(evidence.get("receipt_sha256"), "a4_scoped_evidence.receipt_sha256")
    if not isinstance(evidence.get("audit_pass"), bool):
        raise ValueError("A4 scoped audit status must be boolean")
    if float(evidence.get("finest_step", float("nan"))) != 0.005:
        raise ValueError("A4 finest derivative step drifted")

    payload = {k: copy.deepcopy(v) for k, v in checkpoint.items() if k not in {"schema", "task", "digest"}}
    if checkpoint.get("digest") != _sha256(payload):
        raise ValueError("A5 routing checkpoint digest mismatch")


def enforce_scoped_a4_gate(checkpoint: Mapping[str, Any]) -> None:
    """Require the scoped A4 divergence gate, without promoting PDE validity."""
    validate_checkpoint(checkpoint)
    if checkpoint["a4_scoped_evidence"]["audit_pass"] is not True:
        raise AssertionError("exact current-I1 A4 scoped divergence audit did not pass")
    if checkpoint["truth_boundary"]["pde_validated"] is not False:
        raise AssertionError("scoped A4 evidence must not promote pde_validated")


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-input")
    parser.add_argument("--output", required=True)
    parser.add_argument("--enforce-existing", action="store_true")
    args = parser.parse_args()
    target = Path(args.output)
    if args.enforce_existing:
        checkpoint = json.loads(target.read_text(encoding="utf-8"))
        enforce_scoped_a4_gate(checkpoint)
        return
    if not args.audit_input:
        raise SystemExit("--audit-input is required when materializing the checkpoint")
    audit = json.loads(Path(args.audit_input).read_text(encoding="utf-8"))
    checkpoint = build_checkpoint(audit)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    _main()
