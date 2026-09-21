"""Typed/provenance ingest seam for the current RF40 X4 radial-force stage.

Agent-5 integration glue only.  This binds the existing unmodulated production
lineage A1 #1005 -> A2 #1010 -> A3 #1037 -> A3 #1045 and the matching
implementation-distinct A4 #1046 FD4 audit.  It does not create new mathematics,
pressure, forcing, an NS defect, a correction velocity, or PDE admission.
"""
from __future__ import annotations

import copy
import hashlib
import json
import re
from typing import Any, Mapping

from openai_ns_reconstruction import kokuno_a5_rf40_power_law_radial_stress_ingest_contract as parent

SCHEMA_NAME = "kokuno-agent5-rf40-power-law-radial-force-ingest-v1"
TASK_ID = "KOKUNO-A5-RF40-POWER-LAW-RADIAL-FORCE-INGEST-105"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 1039,
    "head": "6d88e8a83b51704f65aa125006b7671546c77c43",
    "branch": "codex/kokuno-a5-rf40-power-law-radial-stress-ingest-104",
}

AGENT3_POWER_LAW_FORCE = {
    "pr": 1045,
    "head": "3ec20a77b551be819a71e3308e9ecab2a2226406",
    "branch": "codex/kokuno-a3-rf40-power-law-radial-force-114",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_power_law_nonlinear_radial_force.py",
    "source_blob": "39c98d3f4df05846848b897c240b3b56fb6fbe2c",
    "test_path": "tests/test_constrained_kokuno_current_rf40_power_law_nonlinear_radial_force.py",
    "workflow_path": ".github/workflows/kokuno-agent3-current-rf40-power-law-nonlinear-radial-force.yml",
    "parent_stress_pr": 1037,
    "parent_stress_head": "9a5cdcdf78b7862d5ebe171efb9ca56d53664612",
    "agent2_composite_head": "e36d9da4b7f037e998f5b1658f8c0ea291a76b80",
    "agent1_leading_head": "2c76ebdc41d6c566f43a1305034ba2ff9dce410b",
    "formula": "(div T)_r = partial_z sigma_1",
    "production_derivative": "centered_fd2",
    "z_step_ladder": [0.02, 0.01, 0.005],
    "piece_closure_gate": 5.0e-10,
    "fine_pair_stability_gate": 5.0e-2,
    "complete_ns_defect_evidence": False,
    "authorized_as_ns_correction_target": False,
}

AGENT4_POWER_LAW_FORCE_AUDIT = {
    "pr": 1046,
    "head": "df5af9c503d69d665cea93d66fa30bf6f00863a8",
    "branch": "codex/kokuno-a4-rf40-power-law-radial-force-independent-audit-107",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_rf40_power_law_radial_force_independent_audit.py",
    "test_path": "tests/test_constrained_kokuno_a4_current_rf40_power_law_radial_force_independent_audit.py",
    "workflow_path": ".github/workflows/kokuno-agent4-current-rf40-power-law-radial-force-independent-audit.yml",
    "audited_agent3_pr": 1045,
    "audited_agent3_head": "3ec20a77b551be819a71e3308e9ecab2a2226406",
    "audited_agent3_source_blob": "39c98d3f4df05846848b897c240b3b56fb6fbe2c",
    "production_operator": "centered_fd2",
    "independent_operator": "five_point_centered_fd4",
    "frozen_seed": 9173791,
    "frozen_time": 0.39,
    "heldout_z_centers": [-0.143, -0.067, 0.049],
    "z_step_ladder": [0.02, 0.01, 0.005],
    "radial_count": 169,
    "radial_interval": [0.005, 0.44],
    "fine_relative_rms_gate": 5.0e-2,
    "fine_relative_max_gate": 1.5e-1,
    "fine_relative_integral_l2_gate": 5.0e-2,
    "axis_near_normalized_error_gate": 1.5e-1,
    "independent_fine_pair_relative_rms_gate": 5.0e-2,
    "piece_closure_gate": 5.0e-10,
    "nontrivial_force_rms_floor": 1.0e-12,
    "implementation_distinct": True,
    "complete_ns_residual_evidence": False,
    "authorized_as_ns_correction_target": False,
}

FROZEN_SCIENCE = copy.deepcopy(parent.FROZEN_SCIENCE)
FINAL_GATE = copy.deepcopy(parent.FINAL_GATE)
ST006_BASELINE = copy.deepcopy(parent.ST006_BASELINE)
READINESS = copy.deepcopy(parent.READINESS)

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5
CANONICAL_QUADRATURE = (24, 48, 96)

TRUTH_BOUNDARY = copy.deepcopy(parent.TRUTH_BOUNDARY)
TRUTH_BOUNDARY.update({
    "current_nonlinear_radial_stress_through_power_law_materialized": True,
    "current_nonlinear_radial_force_through_power_law_materialized": True,
    "agent4_dedicated_power_law_radial_force_audit_present": True,
    "agent4_independent_power_law_radial_force_audit_registered": True,
    "agent4_independent_power_law_radial_force_audit_admitted": False,
    "scoped_power_law_radial_force_authorized_as_ns_correction_target": False,
    "newer_modulated_a1_lineage_consumed": False,
    "evidence_transferred_from_newer_modulated_lineage": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_defect_materialized": False,
    "real_agent3_ns_correction_velocity_materialized": False,
    "real_candidate_finite_correction_cycle_run": False,
    "complete_candidate_api_ready": False,
    "heldout_normalized_ns_residual_assessed": False,
    "scientific_admission": False,
    "pde_validated": False,
})

PIPELINE_POSITION = {
    "stage": "unmodulated correction-side X4 radial force plus matching independent A4 FD4 audit registered",
    "production_velocity": "A1 #1005 -> A2 #1010 through RF40 power-law X4",
    "correction_side": "A3 #1028 mean -> A3 #1037 compact stress -> A3 #1045 radial force partial_z sigma_1",
    "independent_evidence": "A4 #1046 independently reconstructs the exact #1045 force from public stress receipts with centered five-point FD4",
    "not_output": "pressure, restricted forcing, complete NS defect, correction authorization/velocity, finite correction cycle, export-ready candidate, held-out full-NS assessment, or PDE validation",
    "next_shortest_blocker": "post-X4/global velocity plus matched pressure and preregistered restricted forcing; only then form the complete NS defect before correction authorization",
}

_EXPECTED = {
    "parent_a5": PARENT_A5,
    "agent3_power_law_force": AGENT3_POWER_LAW_FORCE,
    "agent4_power_law_force_audit": AGENT4_POWER_LAW_FORCE_AUDIT,
    "frozen_science": FROZEN_SCIENCE,
    "final_gate": FINAL_GATE,
    "st006_baseline": ST006_BASELINE,
    "readiness": READINESS,
    "truth_boundary": TRUTH_BOUNDARY,
    "pipeline_position": PIPELINE_POSITION,
    "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
    "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    "canonical_quadrature": list(CANONICAL_QUADRATURE),
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _require_hex40(value: str, label: str) -> None:
    if not _HEX40.fullmatch(value):
        raise ValueError(f"{label} is not an exact 40-hex git identity")


def _enforce_internal_relations(registration: Mapping[str, Any]) -> None:
    p = registration["parent_a5"]
    a3 = registration["agent3_power_law_force"]
    a4 = registration["agent4_power_law_force_audit"]
    for label, obj in (("parent_a5", p), ("agent3_power_law_force", a3), ("agent4_power_law_force_audit", a4)):
        _require_hex40(obj["head"], f"{label}.head")
    _require_hex40(a3["source_blob"], "agent3_power_law_force.source_blob")
    _require_hex40(a4["audited_agent3_source_blob"], "agent4_power_law_force_audit.audited_agent3_source_blob")

    if a4["audited_agent3_pr"] != a3["pr"] or a4["audited_agent3_head"] != a3["head"]:
        raise ValueError("A4 #1046 must audit the exact registered A3 #1045 force artifact")
    if a4["audited_agent3_source_blob"] != a3["source_blob"]:
        raise ValueError("A4/A3 force source identity drifted")
    if not a4["implementation_distinct"] or a4["independent_operator"] == a4["production_operator"]:
        raise ValueError("A4 #1046 must remain implementation-distinct")
    if a3["formula"] != "(div T)_r = partial_z sigma_1":
        raise ValueError("source-required radial-force formula drifted")
    if a3["complete_ns_defect_evidence"] or a3["authorized_as_ns_correction_target"]:
        raise ValueError("scoped A3 radial force cannot become complete-NS evidence")
    if a4["complete_ns_residual_evidence"] or a4["authorized_as_ns_correction_target"]:
        raise ValueError("scoped A4 force audit cannot authorize an NS correction")

    truth = registration["truth_boundary"]
    required_true = (
        "current_nonlinear_radial_stress_through_power_law_materialized",
        "current_nonlinear_radial_force_through_power_law_materialized",
        "agent4_dedicated_power_law_radial_force_audit_present",
        "agent4_independent_power_law_radial_force_audit_registered",
    )
    for key in required_true:
        if truth.get(key) is not True:
            raise ValueError(f"required X4 radial-force truth regressed: {key}")
    required_false = (
        "agent4_independent_power_law_radial_force_audit_admitted",
        "scoped_power_law_radial_force_authorized_as_ns_correction_target",
        "newer_modulated_a1_lineage_consumed",
        "evidence_transferred_from_newer_modulated_lineage",
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "real_agent3_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "complete_candidate_api_ready",
        "heldout_normalized_ns_residual_assessed",
        "scientific_admission",
        "pde_validated",
    )
    for key in required_false:
        if truth.get(key) is not False:
            raise ValueError(f"truth boundary illegally promoted: {key}")

    if registration["readiness"] != parent.READINESS:
        raise ValueError("readiness state drifted")
    if registration["final_gate"] != parent.FINAL_GATE:
        raise ValueError("final scientific gate drifted")
    if registration["st006_baseline"] != parent.ST006_BASELINE:
        raise ValueError("ST006 baseline drifted")
    if registration["frozen_science"] != parent.FROZEN_SCIENCE:
        raise ValueError("frozen science drifted")
    if float(registration["final_normalized_momentum_gate"]) != 1.0e-3:
        raise ValueError("final normalized momentum gate drifted")
    if float(registration["final_normalized_divergence_gate"]) != 1.0e-5:
        raise ValueError("final normalized divergence gate drifted")
    if tuple(registration["canonical_quadrature"]) != CANONICAL_QUADRATURE:
        raise ValueError("canonical quadrature drifted")


def validate_registration(registration: Mapping[str, Any]) -> None:
    if registration.get("schema") != SCHEMA_NAME:
        raise ValueError("registration schema drifted")
    if registration.get("task_id") != TASK_ID:
        raise ValueError("registration task id drifted")
    payload = {k: copy.deepcopy(v) for k, v in registration.items() if k not in {"schema", "task_id", "digest"}}
    if payload != _EXPECTED:
        raise ValueError("registration payload drifted from exact frozen identities/truth boundary")
    if registration.get("digest") != _sha256(payload):
        raise ValueError("registration digest mismatch")
    _enforce_internal_relations(registration)


def build_registration() -> dict[str, Any]:
    payload = copy.deepcopy(_EXPECTED)
    out: dict[str, Any] = {"schema": SCHEMA_NAME, "task_id": TASK_ID, **payload}
    out["digest"] = _sha256(payload)
    validate_registration(out)
    return out


def registration_json() -> str:
    return json.dumps(build_registration(), sort_keys=True, indent=2, allow_nan=False) + "\n"
