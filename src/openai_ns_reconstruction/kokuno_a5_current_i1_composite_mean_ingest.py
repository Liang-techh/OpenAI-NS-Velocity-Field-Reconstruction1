"""Agent-5 integration checkpoint for the exact current-I1 composite/mean seam.

This module is glue only. It advances the A5 routing truth after Agent-2 #1062
materialized the exact Agent-1 #1051 current-I1 leading plus the frozen
complete-curl oscillation, and Agent-3 #1063 consumed that exact composite to
materialize the current-I1 nonlinear m=0 mean attribution.

It deliberately does not copy Agent-1/2/3 mathematics. Agent-4 #1053 remains a
leading-only divergence audit and is not laundered into an audit of the #1062
composite or #1063 nonlinear mean. Agent-1 #1061 is recorded as a newer I2
leading-only sibling; no matching I2 oscillatory composite exists at this
checkpoint.

No pressure, forcing, complete Navier--Stokes defect, correction velocity,
finite correction cycle, held-out full-NS residual, or PDE admission is created
here. Kokuno source provenance remains structural/reconstruction provenance,
not a paper-exact field identity.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_current_i1_routing_checkpoint as parent_a5

SCHEMA = "kokuno-a5-current-i1-composite-mean-ingest-v1"
TASK = "KOKUNO-A5-CURRENT-I1-COMPOSITE-MEAN-INGEST-107"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")

PARENT_A5 = {
    "pr": 1055,
    "head": "ee6e1909e06edabbdf44e204fab1bc072c37d4b6",
    "branch": "kokuno-agent5/current-i1-routing-checkpoint-106",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_i1_routing_checkpoint.py",
    "source_blob": "56b51c7ac8accb5c0f846bc89af8ad72b4eb6cc2",
    "task": "KOKUNO-A5-CURRENT-I1-ROUTING-CHECKPOINT-106",
    "schema": "kokuno-a5-current-i1-routing-checkpoint-v1",
}

AGENT1_CURRENT_I1 = {
    "pr": 1051,
    "head": "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_i1_frozen_gate.py",
    "source_blob": "fe14dab8dccdbca828d722bdfe2f75f287aa9ea8",
    "frozen_i1_closure_gate": 5.0e-7,
    "public_velocity_api": "velocity(x,y,z,t)",
}

AGENT2_CURRENT_I1_COMPOSITE = {
    "pr": 1062,
    "head": "109527f520abb29bbe10372b0517eda44bcad0b6",
    "branch": "kokuno-agent2/current-i1-leading-oscillatory-095",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i1_leading_oscillatory_identity.py",
    "source_blob": "3a8db3e12f57050fd6579b36b9e20ce0e6b95af7",
    "schema": "kokuno-a2-current-i1-leading-oscillatory-identity-v1",
    "task": "K2-OSC-095",
    "consumes_agent1_pr": 1051,
    "consumes_agent1_head": "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5",
    "composition": "u_current_I1 = u_lead_A1_1051 + u_osc_frozen_complete_curl",
    "full_concrete_oscillatory_runtime_digest_bound": True,
    "identity_preserving_save_load": True,
    "public_velocity_api": "velocity(x,y,z,t)",
    "through_i1": True,
    "through_i2": False,
    "global_velocity": False,
    "paper_exact": False,
}

AGENT3_CURRENT_I1_MEAN = {
    "pr": 1063,
    "head": "e6ac8bed0ef775232e18e81e87f6a7572dcba8fb",
    "branch": "agent3/current-i1-nonlinear-mean-20260921",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i1_nonlinear_mean_attribution.py",
    "source_blob": "33f69761d5ba014bf0dab2c74c08e4e72166538b",
    "schema": "kokuno-a3-current-i1-nonlinear-mean-v1",
    "task": "KOKUNO-A3-CURRENT-I1-NONLINEAR-MEAN-115",
    "consumes_agent2_pr": 1062,
    "consumes_agent2_head": "109527f520abb29bbe10372b0517eda44bcad0b6",
    "consumes_agent1_pr": 1051,
    "consumes_agent1_head": "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5",
    "decomposition": "A=(u_lead.grad)u_osc; B=(u_osc.grad)u_lead; Q=(u_osc.grad)u_osc; N=A+B+Q; cylindrical m=0 projection",
    "radial_stress_materialized": False,
    "radial_force_materialized": False,
    "complete_ns_defect": False,
    "correction_velocity": False,
}

AGENT1_CURRENT_I2_SIBLING = {
    "pr": 1061,
    "head": "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3",
    "branch": "kokuno-agent1/current-i2-heat-repair-113",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_i2_heat_repair.py",
    "source_blob": "bfeff163d304c4cf1fb84eb8ff1cdfcf430870b4",
    "consumes_agent1_i1_head": "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5",
    "leading_velocity_through_i2_materialized": True,
    "matching_agent2_i2_composite_materialized": False,
    "consumed_by_this_i1_checkpoint": False,
}

AGENT4_CURRENT_I1_LEADING_AUDIT = {
    "pr": 1053,
    "head": "ad9f8722b6ec89df975f9b2402e63066dc14422d",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_i1_leading_divergence_independent_audit.py",
    "source_blob": "cc8c8f9ed0e8f727035c14232bb7f80d2703ac76",
    "scope": "leading-only divergence through current I1",
    "implementation_distinct": True,
    "audits_agent1_1051_leading": True,
    "audits_agent2_1062_composite": False,
    "audits_agent3_1063_nonlinear_mean": False,
    "complete_ns_residual_evidence": False,
}

FINAL_GATES = copy.deepcopy(parent_a5.FINAL_GATES)
ST006_BASELINE = copy.deepcopy(parent_a5.ST006_BASELINE)
READINESS = copy.deepcopy(parent_a5.READINESS)

TRUTH_BOUNDARY = {
    "current_modulated_leading_through_i1_materialized": True,
    "current_i1_leading_plus_frozen_complete_curl_composite_materialized": True,
    "current_i1_composite_identity_preserving_save_load_available": True,
    "current_i1_full_concrete_oscillatory_runtime_digest_bound": True,
    "current_i1_nonlinear_m0_mean_attribution_materialized": True,
    "agent4_leading_only_audit_registered": True,
    "agent4_matching_current_i1_composite_audit_present": False,
    "agent4_matching_current_i1_composite_audit_registered": False,
    "agent4_matching_current_i1_composite_audit_admitted": False,
    "agent4_matching_current_i1_nonlinear_mean_audit_present": False,
    "agent4_matching_current_i1_nonlinear_mean_audit_registered": False,
    "agent1_current_i2_leading_sibling_materialized": True,
    "agent1_current_i2_sibling_consumed_into_this_i1_checkpoint": False,
    "current_i2_matching_leading_plus_oscillatory_composite_materialized": False,
    "current_i2_nonlinear_mean_materialized": False,
    "current_i1_radial_stress_materialized": False,
    "current_i1_radial_force_materialized": False,
    "current_i1_mean_authorized_as_complete_ns_correction_target": False,
    "current_i3_applied": False,
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
    "current_stage": "exact current-I1 leading+oscillatory composite plus nonlinear m=0 mean attribution",
    "current_identity": "A1 #1051 -> A2 #1062 -> A3 #1063",
    "independent_evidence": "A4 #1053 audits only A1 #1051 leading divergence; it does not audit A2 #1062 composite or A3 #1063 mean",
    "newer_non_consumed_sibling": "A1 #1061 leading-only through I2; no matching A2 I2 composite yet",
    "next_shortest_blocker": "recompose frozen oscillation onto exact A1 #1061 I2 identity and obtain implementation-distinct A4 audit of the matching current composite; do not transfer the #1053 leading-only audit",
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
    if parent_a5.AGENT1_CURRENT_I1["head"] != AGENT1_CURRENT_I1["head"]:
        raise ValueError("parent A5 current-I1 Agent-1 identity drifted")
    if parent_a5.AGENT4_CURRENT_I1_AUDIT["head"] != AGENT4_CURRENT_I1_LEADING_AUDIT["head"]:
        raise ValueError("parent A5 current-I1 Agent-4 identity drifted")
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
        "agent1_current_i1": copy.deepcopy(AGENT1_CURRENT_I1),
        "agent2_current_i1_composite": copy.deepcopy(AGENT2_CURRENT_I1_COMPOSITE),
        "agent3_current_i1_mean": copy.deepcopy(AGENT3_CURRENT_I1_MEAN),
        "agent1_current_i2_sibling": copy.deepcopy(AGENT1_CURRENT_I2_SIBLING),
        "agent4_current_i1_leading_audit": copy.deepcopy(AGENT4_CURRENT_I1_LEADING_AUDIT),
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
        raise ValueError("unexpected A5 composite/mean registration identity")

    expected_blocks = {
        "parent_a5": PARENT_A5,
        "agent1_current_i1": AGENT1_CURRENT_I1,
        "agent2_current_i1_composite": AGENT2_CURRENT_I1_COMPOSITE,
        "agent3_current_i1_mean": AGENT3_CURRENT_I1_MEAN,
        "agent1_current_i2_sibling": AGENT1_CURRENT_I2_SIBLING,
        "agent4_current_i1_leading_audit": AGENT4_CURRENT_I1_LEADING_AUDIT,
    }
    for name, expected in expected_blocks.items():
        block = registration.get(name)
        if not isinstance(block, Mapping):
            raise ValueError(f"missing {name}")
        _require_hex40(block.get("head"), f"{name}.head")
        _require_hex40(block.get("source_blob"), f"{name}.source_blob")
        if dict(block) != expected:
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

    a2 = registration["agent2_current_i1_composite"]
    a3 = registration["agent3_current_i1_mean"]
    i2 = registration["agent1_current_i2_sibling"]
    a4 = registration["agent4_current_i1_leading_audit"]
    if a2["consumes_agent1_head"] != AGENT1_CURRENT_I1["head"]:
        raise ValueError("A2 composite detached from exact current-I1 leading")
    if a3["consumes_agent2_head"] != a2["head"] or a3["consumes_agent1_head"] != AGENT1_CURRENT_I1["head"]:
        raise ValueError("A3 nonlinear mean detached from exact current-I1 composite")
    if i2["consumed_by_this_i1_checkpoint"] is not False or i2["matching_agent2_i2_composite_materialized"] is not False:
        raise ValueError("newer I2 leading sibling was illegally consumed")
    if a4["audits_agent2_1062_composite"] is not False or a4["audits_agent3_1063_nonlinear_mean"] is not False:
        raise ValueError("leading-only A4 evidence was laundered into composite/mean validation")

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
