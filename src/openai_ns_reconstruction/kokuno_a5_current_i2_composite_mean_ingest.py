"""Agent-5 integration checkpoint for the exact current-I2 composite/mean seam.

This module is glue only. It advances A5 routing after Agent-2 #1071
materialized the exact Agent-1 #1061 current-I2 leading plus frozen complete-
curl oscillation and Agent-3 #1073 consumed that same semantic identity to
materialize the current-I2 nonlinear m=0 mean attribution.

It deliberately does not copy Agent-1/2/3/4 mathematics. Agent-4 #1064 remains
implementation-distinct scoped divergence evidence for the older current-I1
composite only; it is not transplanted to current I2. At this checkpoint no
dedicated current-I2 Agent-4 audit exists. Agent-1 #1072 is recorded as a newer
leading-only sibling through reserved I3; it is not consumed here and the
source positive-order I3 correction remains unmaterialized.

No pressure, forcing, complete Navier--Stokes defect, correction velocity,
finite correction cycle, held-out full-NS residual, or PDE admission is created
here. The high-precision I2 repair may remain sub-ulp at the public float64
velocity boundary, so public_float64_I2_delta_verified stays false. Kokuno
source provenance is structural/reconstruction provenance, not paper-exact or
OpenAI-field identity.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_current_i1_composite_mean_ingest as parent_a5

SCHEMA = "kokuno-a5-current-i2-composite-mean-ingest-v1"
TASK = "KOKUNO-A5-CURRENT-I2-COMPOSITE-MEAN-INGEST-108"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")

PARENT_A5 = {
    "pr": 1065,
    "head": "2d2cd726107a221b88b3bdfeaa0adcce83c3d57f",
    "branch": "kokuno-agent5/current-i1-composite-mean-ingest-107",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_i1_composite_mean_ingest.py",
    "source_blob": "c4287eea0c49510e5f1be0b2f001a6c3cce9f108",
    "task": "KOKUNO-A5-CURRENT-I1-COMPOSITE-MEAN-INGEST-107",
    "schema": "kokuno-a5-current-i1-composite-mean-ingest-v1",
}

AGENT1_CURRENT_I2 = {
    "pr": 1061,
    "head": "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3",
    "branch": "kokuno-agent1/current-i2-heat-repair-113",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_i2_heat_repair.py",
    "source_blob": "bfeff163d304c4cf1fb84eb8ff1cdfcf430870b4",
    "leading_velocity_through_i2_materialized": True,
    "public_velocity_api": "velocity(x,y,z,t)",
    "public_float64_i2_delta_verified": False,
    "through_i3": False,
    "global_velocity": False,
    "paper_exact": False,
}

AGENT2_CURRENT_I2_COMPOSITE = {
    "pr": 1071,
    "head": "48d69e37e78e7f7f0e4e9936f28ff7974719288d",
    "branch": "kokuno-agent2/current-i2-leading-oscillatory-096",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i2_leading_oscillatory_identity.py",
    "source_blob": "be68248aca97173a15474324f63efff2c9ffce56",
    "schema": "kokuno-a2-current-i2-leading-oscillatory-identity-v1",
    "task": "K2-OSC-096",
    "consumes_agent1_pr": 1061,
    "consumes_agent1_head": "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3",
    "composition": "u_current_I2 = u_lead_A1_1061 + u_osc_frozen_complete_curl",
    "full_concrete_oscillatory_runtime_digest_bound": True,
    "identity_preserving_save_load": True,
    "public_velocity_api": "velocity(x,y,z,t)",
    "through_i2": True,
    "through_i3": False,
    "public_float64_i2_delta_verified": False,
    "global_velocity": False,
    "paper_exact": False,
}

AGENT3_CURRENT_I2_MEAN = {
    "pr": 1073,
    "head": "a798056d7ddabe0bf4020082861807bc5ac1efbd",
    "branch": "agent3/current-i2-nonlinear-mean-20260921",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i2_nonlinear_mean_attribution.py",
    "source_blob": "a9288dc3fe5e3f5e37b33c10762d42b93dfca529",
    "schema": "kokuno-a3-current-i2-nonlinear-mean-v1",
    "task": "KOKUNO-A3-CURRENT-I2-NONLINEAR-MEAN-116",
    "consumes_agent2_pr": 1071,
    "consumes_agent2_head": "48d69e37e78e7f7f0e4e9936f28ff7974719288d",
    "consumes_agent1_pr": 1061,
    "consumes_agent1_head": "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3",
    "decomposition": "A=(u_lead.grad)u_osc; B=(u_osc.grad)u_lead; Q=(u_osc.grad)u_osc; N=A+B+Q; cylindrical m=0 projection",
    "public_float64_i2_delta_verified": False,
    "radial_stress_materialized": False,
    "radial_force_materialized": False,
    "complete_ns_defect": False,
    "correction_velocity": False,
}

AGENT1_CURRENT_I3_SIBLING = {
    "pr": 1072,
    "head": "5fb7b583062b4a991db86e39fdcdb231022b70c9",
    "branch": "kokuno-agent1/current-leading-preserve-i3-115",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_i3_leading_preservation.py",
    "source_blob": "fd1d7dcbd68ce42a715fdf18d7c4d05daaef0042",
    "consumes_agent1_i2_head": "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3",
    "leading_velocity_through_i3_materialized": True,
    "source_positive_order_i3_profiles_materialized": False,
    "i3_positive_order_correction_materialized": False,
    "current_i4_mean_correction_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "matching_agent2_i3_composite_materialized": False,
    "consumed_by_this_i2_checkpoint": False,
}

AGENT4_CURRENT_I1_ONLY_EVIDENCE = {
    "pr": 1064,
    "head": "33626e55f301bbb5a28f75684de49d020bc6039d",
    "branch": "kokuno-agent4/current-i1-leading-oscillatory-divergence-109",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_i1_composite_divergence_independent_audit.py",
    "source_blob": "22b9947684f318e892e17b0a455efac8da63fdb8",
    "scope": "leading, leading+oscillatory, and oscillatory-increment divergence through current I1 only",
    "implementation_distinct": True,
    "audits_agent2_1062_current_i1_composite": True,
    "audits_agent2_1071_current_i2_composite": False,
    "audits_agent3_1073_current_i2_mean": False,
    "scientific_admission": False,
    "complete_ns_residual_evidence": False,
}

AGENT4_CURRENT_I2_AUDIT_STATUS = {
    "dedicated_composite_audit_present": False,
    "independent_composite_audit_registered": False,
    "independent_composite_audit_admitted": False,
    "independent_nonlinear_mean_audit_present": False,
    "complete_ns_residual_evidence": False,
}

FINAL_GATES = copy.deepcopy(parent_a5.FINAL_GATES)
ST006_BASELINE = copy.deepcopy(parent_a5.ST006_BASELINE)
READINESS = copy.deepcopy(parent_a5.READINESS)

TRUTH_BOUNDARY = {
    "current_i1_composite_and_mean_parent_checkpoint_preserved": True,
    "current_i2_leading_materialized": True,
    "current_i2_matching_leading_plus_oscillatory_composite_materialized": True,
    "current_i2_composite_identity_preserving_save_load_available": True,
    "current_i2_full_concrete_oscillatory_runtime_digest_bound": True,
    "current_i2_nonlinear_m0_mean_attribution_materialized": True,
    "current_i2_public_float64_delta_verified": False,
    "agent4_current_i1_only_evidence_retained": True,
    "agent4_matching_current_i2_composite_audit_present": False,
    "agent4_matching_current_i2_composite_audit_registered": False,
    "agent4_matching_current_i2_composite_audit_admitted": False,
    "agent4_matching_current_i2_nonlinear_mean_audit_present": False,
    "current_i2_radial_stress_materialized": False,
    "current_i2_radial_force_materialized": False,
    "current_i2_mean_authorized_as_complete_ns_correction_target": False,
    "agent1_current_i3_leading_sibling_materialized": True,
    "agent1_current_i3_sibling_consumed_into_this_i2_checkpoint": False,
    "source_positive_order_i3_profiles_materialized": False,
    "current_i3_matching_leading_plus_oscillatory_composite_materialized": False,
    "current_i4_applied": False,
    "outer_global_leading_velocity_materialized": False,
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
    "current_stage": "exact current-I2 leading+oscillatory composite plus matching nonlinear m=0 mean attribution; no current-I2 A4 audit yet",
    "current_identity": "A1 #1061 -> A2 #1071 -> A3 #1073",
    "independent_evidence": "A4 #1064 remains implementation-distinct scoped divergence evidence for exact current-I1 A2 #1062 only and cannot be transplanted to current I2",
    "newer_non_consumed_sibling": "A1 #1072 preserves leading velocity through I3; source positive-order I3 correction and matching A2 I3 composite remain absent",
    "next_shortest_blocker": "obtain a matching implementation-distinct A4 audit of exact A2 #1071 current-I2 composite, then carry current-I2 mean to stress/force on the same identity; in parallel recompose the oscillation only when the exact I3 identity is ready without laundering I2 evidence",
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
    if parent_a5.AGENT1_CURRENT_I2_SIBLING["head"] != AGENT1_CURRENT_I2["head"]:
        raise ValueError("parent A5 current-I2 Agent-1 identity drifted")
    if parent_a5.AGENT4_CURRENT_I1_COMPOSITE_AUDIT["head"] != AGENT4_CURRENT_I1_ONLY_EVIDENCE["head"]:
        raise ValueError("parent A5 current-I1 Agent-4 evidence drifted")
    if parent_a5.FINAL_GATES != FINAL_GATES:
        raise ValueError("fixed project gates drifted from parent A5")
    if parent_a5.ST006_BASELINE != ST006_BASELINE:
        raise ValueError("ST006 baseline drifted from parent A5")
    if parent_a5.READINESS != READINESS:
        raise ValueError("core readiness drifted from parent A5")


def build_registration() -> dict[str, Any]:
    """Build one deterministic provenance/truth checkpoint; accepts no tuning inputs."""
    _validate_parent_contract()
    payload = {
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent1_current_i2": copy.deepcopy(AGENT1_CURRENT_I2),
        "agent2_current_i2_composite": copy.deepcopy(AGENT2_CURRENT_I2_COMPOSITE),
        "agent3_current_i2_mean": copy.deepcopy(AGENT3_CURRENT_I2_MEAN),
        "agent1_current_i3_sibling": copy.deepcopy(AGENT1_CURRENT_I3_SIBLING),
        "agent4_current_i1_only_evidence": copy.deepcopy(AGENT4_CURRENT_I1_ONLY_EVIDENCE),
        "agent4_current_i2_audit_status": copy.deepcopy(AGENT4_CURRENT_I2_AUDIT_STATUS),
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
        raise ValueError("unexpected A5 current-I2 composite/mean registration identity")

    expected_blocks = {
        "parent_a5": PARENT_A5,
        "agent1_current_i2": AGENT1_CURRENT_I2,
        "agent2_current_i2_composite": AGENT2_CURRENT_I2_COMPOSITE,
        "agent3_current_i2_mean": AGENT3_CURRENT_I2_MEAN,
        "agent1_current_i3_sibling": AGENT1_CURRENT_I3_SIBLING,
        "agent4_current_i1_only_evidence": AGENT4_CURRENT_I1_ONLY_EVIDENCE,
    }
    for name, expected in expected_blocks.items():
        block = registration.get(name)
        if not isinstance(block, Mapping):
            raise ValueError(f"missing {name}")
        _require_hex40(block.get("head"), f"{name}.head")
        _require_hex40(block.get("source_blob"), f"{name}.source_blob")
        if dict(block) != expected:
            raise ValueError(f"{name} identity/truth drifted")

    if registration.get("agent4_current_i2_audit_status") != AGENT4_CURRENT_I2_AUDIT_STATUS:
        raise ValueError("current-I2 A4 audit status drifted")
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

    a1 = registration["agent1_current_i2"]
    a2 = registration["agent2_current_i2_composite"]
    a3 = registration["agent3_current_i2_mean"]
    i3 = registration["agent1_current_i3_sibling"]
    old_a4 = registration["agent4_current_i1_only_evidence"]
    a4_i2 = registration["agent4_current_i2_audit_status"]

    if a2["consumes_agent1_head"] != a1["head"]:
        raise ValueError("A2 current-I2 composite detached from exact A1 #1061 identity")
    if a3["consumes_agent2_head"] != a2["head"] or a3["consumes_agent1_head"] != a1["head"]:
        raise ValueError("A3 current-I2 mean detached from exact A1/A2 identity")
    if a1["public_float64_i2_delta_verified"] is not False or a2["public_float64_i2_delta_verified"] is not False or a3["public_float64_i2_delta_verified"] is not False:
        raise ValueError("sub-ulp I2 evidence was promoted to a verified public float64 delta")
    if i3["consumes_agent1_i2_head"] != a1["head"]:
        raise ValueError("newer I3 leading sibling detached from exact I2 parent")
    if i3["consumed_by_this_i2_checkpoint"] is not False or i3["matching_agent2_i3_composite_materialized"] is not False:
        raise ValueError("newer I3 leading sibling was illegally consumed")
    if i3["source_positive_order_i3_profiles_materialized"] is not False or i3["i3_positive_order_correction_materialized"] is not False:
        raise ValueError("source positive-order I3 correction was invented")
    if old_a4["audits_agent2_1071_current_i2_composite"] is not False or old_a4["audits_agent3_1073_current_i2_mean"] is not False:
        raise ValueError("current-I1 A4 evidence was laundered into current-I2 validation")
    if old_a4["scientific_admission"] is not False or old_a4["complete_ns_residual_evidence"] is not False:
        raise ValueError("scoped current-I1 A4 evidence was promoted beyond its scope")
    if any(bool(a4_i2[key]) for key in a4_i2):
        raise ValueError("missing current-I2 A4 evidence was promoted")

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
    payload = {
        key: copy.deepcopy(value)
        for key, value in registration.items()
        if key not in {"schema", "task", "digest"}
    }
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
