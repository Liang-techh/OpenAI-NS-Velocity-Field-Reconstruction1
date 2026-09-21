"""Agent-5 routing checkpoint for the current I4 velocity / I2 correction split.

This is integration glue only.  The current callable velocity frontier is exact
A1 #1079 -> A2 #1080 through reserved I4, with A4 #1082 providing an
implementation-distinct scoped divergence audit of that save/reloaded public
Cartesian velocity.  The correction-side frontier is exact A3 #1081, which
routes the older current-I2 nonlinear mean into compact theta/axial stress.

Those frontiers are deliberately *not* treated as one candidate identity: the
I2 stress may not be transplanted onto the I4 candidate, and the I4 divergence
audit does not validate A3 stress, momentum, a complete NS defect, or a
correction cycle.  Source positive-order I3 correction, source I4 mean
correction, terminal/global completion, matched pressure, restricted forcing,
radial force, Cartesian correction velocity and PDE admission remain absent.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_current_i2_composite_mean_ingest as parent_a5

SCHEMA = "kokuno-a5-current-i4-routing-i2-stress-ingest-v1"
TASK = "KOKUNO-A5-CURRENT-I4-ROUTING-I2-STRESS-INGEST-109"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")

PARENT_A5 = {
    "pr": 1075,
    "head": "bb9beb985a44c6868fe36daa025704742508e9d4",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_i2_composite_mean_ingest.py",
    "source_blob": "bd34059f6bb2d1c4996f226b023588b4b6806680",
    "task": "KOKUNO-A5-CURRENT-I2-COMPOSITE-MEAN-INGEST-108",
    "schema": "kokuno-a5-current-i2-composite-mean-ingest-v1",
}

AGENT3_CURRENT_I2_STRESS = {
    "pr": 1081,
    "head": "cb9ee25d0c466d42ec7ca30ee89cc10517d4b5c4",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i2_nonlinear_radial_stress.py",
    "source_blob": "876cdd2cbe4060b9d249cfc4960c26da0d577ca9",
    "schema": "kokuno-a3-current-i2-nonlinear-radial-stress-v1",
    "task": "KOKUNO-A3-CURRENT-I2-NONLINEAR-RADIAL-STRESS-117",
    "parent_agent3_pr": 1073,
    "parent_agent3_head": "a798056d7ddabe0bf4020082861807bc5ac1efbd",
    "candidate_stage": "current-I2",
    "theta_e2_stress_materialized": True,
    "axial_e1_stress_materialized": True,
    "radial_force_partial_z_sigma1_materialized": False,
    "authorized_as_complete_ns_correction_target": False,
    "consumed_by_current_i4_candidate": False,
}

AGENT1_CURRENT_I4_LEADING = {
    "pr": 1079,
    "head": "b06742ca6e189499192ede3cce40f62cdc1e35ca",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_i4_leading_preservation.py",
    "source_blob": "6f04ce0a856b44430402576dad88438da90d1ebb",
    "schema": "kokuno-pa16-current-cartesian-i4-leading-preservation-v1",
    "parent_agent1_pr": 1072,
    "parent_agent1_head": "5fb7b583062b4a991db86e39fdcdb231022b70c9",
    "leading_velocity_through_i4_materialized": True,
    "source_positive_order_i3_correction_materialized": False,
    "source_i4_mean_correction_materialized": False,
    "terminal_global_leading_completion_materialized": False,
}

AGENT2_CURRENT_I4_COMPOSITE = {
    "pr": 1080,
    "head": "c40d8ddecd2971544a6e07dab093436b423cf326",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i4_leading_oscillatory_identity.py",
    "source_blob": "2a0a5aa5966b02da856bdcf51940f3c186042802",
    "schema": "kokuno-a2-current-i4-leading-oscillatory-identity-v1",
    "task": "K2-OSC-097",
    "consumes_agent1_pr": 1079,
    "consumes_agent1_head": "b06742ca6e189499192ede3cce40f62cdc1e35ca",
    "candidate_stage": "current-I4",
    "composition": "u_current_I4 = u_lead_A1_1079 + u_osc_frozen_complete_curl",
    "identity_preserving_save_load": True,
    "public_velocity_api": "velocity(x,y,z,t)",
    "source_positive_order_i3_correction_materialized": False,
    "source_i4_mean_correction_materialized": False,
    "terminal_global_completion_materialized": False,
}

AGENT4_CURRENT_I4_COMPOSITE_AUDIT = {
    "pr": 1082,
    "head": "0ff663fcabaa52af984d273ba87356a72481c85c",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_i4_composite_divergence_independent_audit.py",
    "source_blob": "518208d1feac5b0f09d08e9989b2eba9275fbd82",
    "schema": "kokuno-a4-current-i4-composite-divergence-audit-v1",
    "audited_agent2_pr": 1080,
    "audited_agent2_head": "c40d8ddecd2971544a6e07dab093436b423cf326",
    "audited_agent1_pr": 1079,
    "audited_agent1_head": "b06742ca6e189499192ede3cce40f62cdc1e35ca",
    "operator": "centered Cartesian FD2 from public save/reloaded velocity only",
    "spatial_steps": [0.02, 0.01, 0.005],
    "seed": 9173831,
    "scope": "leading, leading+oscillatory, and oscillatory-increment divergence through current I4",
    "implementation_distinct": True,
    "audits_agent3_i2_stress": False,
    "scientific_admission": False,
    "complete_ns_residual_evidence": False,
}

FINAL_GATES = copy.deepcopy(parent_a5.FINAL_GATES)
ST006_BASELINE = copy.deepcopy(parent_a5.ST006_BASELINE)
READINESS = copy.deepcopy(parent_a5.READINESS)

TRUTH_BOUNDARY = {
    "parent_current_i2_checkpoint_preserved": True,
    "velocity_frontier_stage": "current-I4",
    "correction_frontier_stage": "current-I2",
    "current_i4_leading_materialized": True,
    "current_i4_leading_plus_oscillatory_composite_materialized": True,
    "current_i4_identity_preserving_save_load_available": True,
    "agent4_matching_current_i4_composite_audit_present": True,
    "agent4_matching_current_i4_composite_audit_registered": True,
    "agent4_matching_current_i4_composite_audit_admitted": False,
    "current_i2_compact_radial_stress_materialized": True,
    "current_i2_radial_force_materialized": False,
    "current_i2_stress_consumed_by_current_i4_candidate": False,
    "current_i2_stress_authorized_as_complete_ns_correction_target": False,
    "agent4_current_i4_audit_validates_current_i2_stress": False,
    "source_positive_order_i3_correction_materialized": False,
    "source_i4_mean_correction_materialized": False,
    "terminal_global_leading_completion_materialized": False,
    "global_leading_plus_oscillatory_velocity_materialized": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "preregistered_restricted_forcing_materialized": False,
    "forcing_proven_not_residual_defined": False,
    "complete_identity_bound_ns_defect_materialized": False,
    "current_cartesian_correction_velocity_materialized": False,
    "finite_correction_cycle_admitted": False,
    "real_candidate_finite_correction_cycle_run": False,
    "heldout_normalized_full_ns_residual_assessed": False,
    "canonical_24_48_96_whole_domain_admission_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "scientific_admission": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

PIPELINE_POSITION = {
    "velocity_frontier": "A1 #1079 -> A2 #1080 current-I4 Cartesian composite; A4 #1082 independently audits exact #1080 scoped divergence",
    "correction_frontier": "A1 #1061 -> A2 #1071 -> A3 #1073 -> A3 #1081 current-I2 nonlinear mean to compact theta/axial stress",
    "stage_firewall": "current-I2 stress is not consumed by current-I4 candidate; current-I4 A4 divergence evidence does not validate current-I2 stress or complete NS momentum",
    "next_shortest_blocker": "close the stage mismatch before correction authorization: carry correction-side mean/stress/radial-force machinery onto the exact current velocity identity only after required source I3/I4 roles and terminal/global velocity are materialized; matched pressure and preregistered restricted forcing are still required before a complete NS defect exists",
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
        raise ValueError(f"{label} is not a 64-hex digest")


def _validate_parent_contract() -> None:
    if parent_a5.SCHEMA != PARENT_A5["schema"] or parent_a5.TASK != PARENT_A5["task"]:
        raise ValueError("parent A5 schema/task drifted")
    if parent_a5.AGENT3_CURRENT_I2_MEAN["head"] != AGENT3_CURRENT_I2_STRESS["parent_agent3_head"]:
        raise ValueError("A3 #1081 is detached from parent A5 current-I2 mean")
    if parent_a5.FINAL_GATES != FINAL_GATES or parent_a5.ST006_BASELINE != ST006_BASELINE:
        raise ValueError("fixed project gates or ST006 baseline drifted")
    if parent_a5.READINESS != READINESS:
        raise ValueError("core readiness drifted")


def build_registration() -> dict[str, Any]:
    _validate_parent_contract()
    payload = {
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent3_current_i2_stress": copy.deepcopy(AGENT3_CURRENT_I2_STRESS),
        "agent1_current_i4_leading": copy.deepcopy(AGENT1_CURRENT_I4_LEADING),
        "agent2_current_i4_composite": copy.deepcopy(AGENT2_CURRENT_I4_COMPOSITE),
        "agent4_current_i4_composite_audit": copy.deepcopy(AGENT4_CURRENT_I4_COMPOSITE_AUDIT),
        "final_gates": copy.deepcopy(FINAL_GATES),
        "st006_baseline": copy.deepcopy(ST006_BASELINE),
        "readiness": copy.deepcopy(READINESS),
        "truth_boundary": copy.deepcopy(TRUTH_BOUNDARY),
        "pipeline_position": copy.deepcopy(PIPELINE_POSITION),
    }
    registration = {"schema": SCHEMA, "task": TASK, **payload}
    registration["digest"] = _sha256(payload)
    validate_registration(registration)
    return registration


def validate_registration(registration: Mapping[str, Any]) -> None:
    _validate_parent_contract()
    if registration.get("schema") != SCHEMA or registration.get("task") != TASK:
        raise ValueError("unexpected A5 I4/I2 routing registration identity")
    expected = {
        "parent_a5": PARENT_A5,
        "agent3_current_i2_stress": AGENT3_CURRENT_I2_STRESS,
        "agent1_current_i4_leading": AGENT1_CURRENT_I4_LEADING,
        "agent2_current_i4_composite": AGENT2_CURRENT_I4_COMPOSITE,
        "agent4_current_i4_composite_audit": AGENT4_CURRENT_I4_COMPOSITE_AUDIT,
    }
    for name, block_expected in expected.items():
        block = registration.get(name)
        if not isinstance(block, Mapping):
            raise ValueError(f"missing {name}")
        _require_hex40(block.get("head"), f"{name}.head")
        _require_hex40(block.get("source_blob"), f"{name}.source_blob")
        if dict(block) != block_expected:
            raise ValueError(f"{name} identity/truth drifted")

    if registration.get("final_gates") != FINAL_GATES:
        raise ValueError("fixed project gates drifted")
    if registration.get("st006_baseline") != ST006_BASELINE:
        raise ValueError("ST006 baseline drifted")
    if registration.get("readiness") != READINESS:
        raise ValueError("core readiness drifted")
    if registration.get("truth_boundary") != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drifted")
    if registration.get("pipeline_position") != PIPELINE_POSITION:
        raise ValueError("pipeline position drifted")

    a1 = registration["agent1_current_i4_leading"]
    a2 = registration["agent2_current_i4_composite"]
    a3 = registration["agent3_current_i2_stress"]
    a4 = registration["agent4_current_i4_composite_audit"]
    if a2["consumes_agent1_head"] != a1["head"]:
        raise ValueError("A2 current-I4 composite detached from exact A1 #1079")
    if a4["audited_agent2_head"] != a2["head"] or a4["audited_agent1_head"] != a1["head"]:
        raise ValueError("A4 current-I4 audit detached from exact A1/A2 identity")
    if a3["candidate_stage"] != "current-I2" or a2["candidate_stage"] != "current-I4":
        raise ValueError("velocity/correction stage split was erased")
    if a3["consumed_by_current_i4_candidate"] is not False:
        raise ValueError("current-I2 correction evidence was transplanted to current I4")
    if a3["radial_force_partial_z_sigma1_materialized"] is not False:
        raise ValueError("current-I2 radial force was invented")
    if a3["authorized_as_complete_ns_correction_target"] is not False:
        raise ValueError("scoped stress was promoted to complete NS correction target")
    if a4["audits_agent3_i2_stress"] is not False:
        raise ValueError("I4 divergence audit was laundered into an A3 stress audit")
    if a4["scientific_admission"] is not False or a4["complete_ns_residual_evidence"] is not False:
        raise ValueError("scoped I4 A4 audit was promoted beyond its evidence")

    if registration["readiness"] != {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }:
        raise ValueError("core readiness was promoted")
    if FINAL_GATES["normalized_momentum_sampled_max"] != 1.0e-3 or FINAL_GATES["normalized_momentum_volume_l2"] != 1.0e-3:
        raise ValueError("fixed momentum gate drifted")
    if FINAL_GATES["normalized_divergence_sampled_max"] != 1.0e-5 or FINAL_GATES["normalized_divergence_volume_l2"] != 1.0e-5:
        raise ValueError("fixed divergence gate drifted")
    if FINAL_GATES["canonical_quadrature"] != [24, 48, 96]:
        raise ValueError("canonical quadrature drifted")
    if FINAL_GATES["residual_defined_free_forcing_allowed"] is not False:
        raise ValueError("residual-defined free forcing was enabled")

    digest = registration.get("digest")
    _require_hex64(digest, "digest")
    payload = {k: copy.deepcopy(v) for k, v in registration.items() if k not in {"schema", "task", "digest"}}
    if digest != _sha256(payload):
        raise ValueError("registration digest mismatch")


def write_registration(path: str | Path) -> dict[str, Any]:
    registration = build_registration()
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(registration, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return registration


def load_registration(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("registration JSON must contain one object")
    validate_registration(payload)
    return payload


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--enforce-existing", action="store_true")
    args = parser.parse_args()
    if args.enforce_existing:
        load_registration(args.output)
    else:
        write_registration(args.output)


if __name__ == "__main__":
    _main()
