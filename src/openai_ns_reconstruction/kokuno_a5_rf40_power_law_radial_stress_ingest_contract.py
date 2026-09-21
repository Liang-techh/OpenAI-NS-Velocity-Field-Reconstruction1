"""Kokuno Agent-5 ingest contract for the current RF40 power-law radial-stress seam.

Integration/provenance glue only. Exact current chain:
A1 #1005 -> A2 #1010 -> A3 #1028 X4 nonlinear m=0 mean -> A3 #1037
compact theta(e=2)/axial(e=1) stress -> A4 #1038 implementation-distinct
stress audit. A4 #1029 remains the independent audit of the parent X4 mean.
No radial force, complete NS defect, correction authorization, or PDE promotion.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from openai_ns_reconstruction import kokuno_a5_rf40_axial_shutdown_radial_force_audit_registration as parent

SCHEMA_NAME = "kokuno-agent5-rf40-power-law-radial-stress-ingest-v1"
TASK_ID = "KOKUNO-A5-RF40-POWER-LAW-RADIAL-STRESS-INGEST-104"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 1030,
    "head": "2dab469f3d48268bfb7edac18ae454b8265a99c0",
    "branch": "codex/kokuno-a5-rf40-axial-shutdown-radial-force-audit-reg-103",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_rf40_axial_shutdown_radial_force_audit_registration.py",
    "source_blob": "4cab62d7cab7b88cbfb87b474321e3e1c04b1673",
    "test_path": "tests/test_kokuno_a5_rf40_axial_shutdown_radial_force_audit_registration.py",
    "test_blob": "4c2a1ea17d740710572aa1582a12f3963e912af7",
    "workflow_path": ".github/workflows/kokuno-agent5-rf40-axial-shutdown-radial-force-audit-registration.yml",
    "workflow_blob": "25b8b7d037dd4afdf3364d5eb7c0a55c02d047a5",
}

AGENT3_POWER_LAW_STRESS = {
    "pr": 1037,
    "head": "9a5cdcdf78b7862d5ebe171efb9ca56d53664612",
    "branch": "codex/kokuno-a3-rf40-power-law-radial-stress-113",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_power_law_nonlinear_radial_stress.py",
    "source_blob": "44dc6d61bc965cef08b119aa5df3026175abd5de",
    "test_path": "tests/test_constrained_kokuno_current_rf40_power_law_nonlinear_radial_stress.py",
    "test_blob": "40ee7a854e8f7649b4e5dba6fba191d066a67805",
    "workflow_path": ".github/workflows/kokuno-agent3-current-rf40-power-law-nonlinear-radial-stress.yml",
    "workflow_blob": "0c604cc3e7ae84f58ad712a25145ebbd1dc17296",
    "parent_agent3_pr": 1028,
    "parent_agent3_head": "4b7a4400dcf1ee0be6bb07f23c2f8af6f893de88",
    "parent_agent3_source_blob": "4ec98a6adbb7afd16d6595cc24072fe18e62fa87",
    "agent2_composite_pr": 1010,
    "agent2_composite_head": "e36d9da4b7f037e998f5b1658f8c0ea291a76b80",
    "agent1_leading_pr": 1005,
    "agent1_leading_head": "2c76ebdc41d6c566f43a1305034ba2ff9dce410b",
    "theta_exponent": 2,
    "axial_exponent": 1,
    "radial_force_requires_later_partial_z_sigma_1": True,
    "complete_ns_defect_evidence": False,
    "authorized_as_ns_correction_target": False,
}

AGENT4_POWER_LAW_MEAN_AUDIT = {
    "pr": 1029,
    "head": "8a34d11605d095334175df9618f8f5b18f5c210e",
    "branch": "codex/kokuno-a4-rf40-power-law-nonlinear-mean-independent-audit-105",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_rf40_power_law_nonlinear_mean_independent_audit.py",
    "source_blob": "87ddbe38d7b442847fe06423b3dd38fcc139253e",
    "test_path": "tests/test_constrained_kokuno_a4_current_rf40_power_law_nonlinear_mean_independent_audit.py",
    "test_blob": "ee8f0362ccd3f50dd9c91f03665f6c4122877cb7",
    "workflow_path": ".github/workflows/kokuno-agent4-current-rf40-power-law-nonlinear-mean-independent-audit.yml",
    "workflow_blob": "fb868e992488704b15f3d9c936f4eb89715ead76",
    "audited_agent3_pr": 1028,
    "audited_agent3_head": "4b7a4400dcf1ee0be6bb07f23c2f8af6f893de88",
    "audit_scope": "parent_x4_nonlinear_mean_only",
    "audits_agent3_1037_radial_stress": False,
    "implementation_distinct": True,
}

AGENT4_POWER_LAW_STRESS_AUDIT = {
    "pr": 1038,
    "head": "a4546cacfa224c82398ef7245ba41cca41159800",
    "branch": "codex/kokuno-a4-rf40-power-law-radial-stress-audit-106",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_rf40_power_law_radial_stress_independent_audit.py",
    "source_blob": "17f98a0c1054132d251dfe7310120e8d41f60244",
    "test_path": "tests/test_constrained_kokuno_a4_current_rf40_power_law_radial_stress_independent_audit.py",
    "test_blob": "2d476a7dc71dfa90351b1c3de2be1c2bc50cb8d7",
    "workflow_path": ".github/workflows/kokuno-agent4-current-rf40-power-law-radial-stress-independent-audit.yml",
    "workflow_blob": "016a80bdb311e547db4fb4cde85e2a9accb2b241",
    "audited_agent3_pr": 1037,
    "audited_agent3_head": "9a5cdcdf78b7862d5ebe171efb9ca56d53664612",
    "audited_agent3_source_blob": "44dc6d61bc965cef08b119aa5df3026175abd5de",
    "production_operator": "a3_inherited_compact_moment_complement_first_cell_trapezoidal",
    "independent_operator": "local_piecewise_cubic_plus_order8_gauss_legendre_on_json_public_receipts",
    "frozen_seed": 9173781,
    "frozen_time": 0.39,
    "frozen_z": -0.08,
    "radial_interval": [0.005, 0.44],
    "radial_counts": [43, 85, 169],
    "offgrid_count": 169,
    "fine_relative_rms_gate": 5.0e-2,
    "fine_relative_max_gate": 1.5e-1,
    "weighted_moment_relative_error_gate": 3.0e-2,
    "axis_near_normalized_error_gate": 1.5e-1,
    "outer_edge_normalized_stress_gate": 1.0e-8,
    "nontrivial_stress_rms_floor": 1.0e-12,
    "implementation_distinct": True,
    "complete_ns_residual_evidence": False,
    "authorized_as_ns_correction_target": False,
}

OBSERVED_CI_AT_FREEZE = {
    "parent_a5_1030": {"repository_tests_run": 35635405082, "dedicated_run": 35635405080, "status": "queued", "conclusion": None},
    "agent3_1037_power_law_stress": {"repository_tests_run": 35640260781, "dedicated_run": 35640260778, "status": "queued", "conclusion": None},
    "agent4_1029_power_law_mean_audit": {"repository_tests_run": 35635278480, "dedicated_run": 35635278504, "status": "queued", "conclusion": None},
    "agent4_1038_power_law_stress_audit": {"repository_tests_run": 35641711974, "dedicated_run": 35641712034, "status": "queued", "conclusion": None},
}

FROZEN_SCIENCE = copy.deepcopy(parent.FROZEN_SCIENCE)
FINAL_GATE = copy.deepcopy(parent.FINAL_GATE)
ST006_BASELINE = copy.deepcopy(parent.ST006_BASELINE)
READINESS = copy.deepcopy(parent.READINESS)

TRUTH_BOUNDARY = copy.deepcopy(parent.TRUTH_BOUNDARY)
TRUTH_BOUNDARY.update({
    "current_nonlinear_mean_through_power_law_materialized": True,
    "agent4_1029_power_law_mean_audit_registered": True,
    "agent3_1037_power_law_radial_stress_registered": True,
    "current_nonlinear_radial_stress_through_power_law_materialized": True,
    "agent4_dedicated_power_law_radial_stress_audit_present": True,
    "agent4_independent_power_law_radial_stress_audit_registered": True,
    "agent4_independent_power_law_radial_stress_audit_admitted": False,
    "current_nonlinear_radial_force_through_power_law_materialized": False,
    "scoped_power_law_radial_stress_authorized_as_ns_correction_target": False,
})

PIPELINE_POSITION = {
    "stage": "production velocity is identity-bound through X_4; correction-side X_4 mean reaches compact stress with a matching independent A4 stress audit registered",
    "production_velocity": "A1 #1005 -> A2 #1010 through RF40 power-law X_4",
    "correction_side": "A3 #1028 X_4 nonlinear m=0 mean -> A3 #1037 compact theta(e=2)/axial(e=1) stress",
    "independent_evidence": "A4 #1029 audits the parent X_4 nonlinear mean; A4 #1038 independently audits exact A3 #1037 stress via JSON-round-tripped public receipts and a distinct cubic+GL8 operator",
    "new_output": "checksum-bound A5 ingest/registration receipt for exact A3 #1037 plus matching A4 #1038 stress audit",
    "not_output": "scientific admission while CI is unresolved, X_4 partial_z sigma_1 force, complete NS defect, correction authorization/velocity, finite cycle, export-ready candidate, or PDE validation",
    "next_shortest_blocker": "materialize X_4 partial_z sigma_1 radial force on exact #1037 and obtain matching implementation-distinct A4 force audit while A1/A2 continue post-X4/global completion",
}

_EXPECTED = {
    "parent_a5": PARENT_A5,
    "agent3_power_law_stress": AGENT3_POWER_LAW_STRESS,
    "agent4_power_law_mean_audit": AGENT4_POWER_LAW_MEAN_AUDIT,
    "agent4_power_law_stress_audit": AGENT4_POWER_LAW_STRESS_AUDIT,
    "observed_ci_at_freeze": OBSERVED_CI_AT_FREEZE,
    "frozen_science": FROZEN_SCIENCE,
    "final_gate": FINAL_GATE,
    "st006_baseline": ST006_BASELINE,
    "readiness": READINESS,
    "truth_boundary": TRUTH_BOUNDARY,
    "pipeline_position": PIPELINE_POSITION,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _require_hex40(value: str, label: str) -> None:
    if not _HEX40.fullmatch(value):
        raise ValueError(f"{label} is not an exact 40-hex git identity")


def _enforce_internal_relations(registration: Mapping[str, Any]) -> None:
    for label in ("parent_a5", "agent3_power_law_stress", "agent4_power_law_mean_audit", "agent4_power_law_stress_audit"):
        obj = registration[label]
        _require_hex40(obj["head"], f"{label}.head")
        for key in ("source_blob", "test_blob", "workflow_blob"):
            _require_hex40(obj[key], f"{label}.{key}")

    stress = registration["agent3_power_law_stress"]
    mean_audit = registration["agent4_power_law_mean_audit"]
    stress_audit = registration["agent4_power_law_stress_audit"]
    if stress["parent_agent3_pr"] != mean_audit["audited_agent3_pr"] or stress["parent_agent3_head"] != mean_audit["audited_agent3_head"]:
        raise ValueError("A3 stress must consume the exact X4 mean artifact audited by A4 #1029")
    if mean_audit["audits_agent3_1037_radial_stress"]:
        raise ValueError("A4 #1029 must not be relabeled as a stress audit")
    if stress_audit["audited_agent3_pr"] != stress["pr"] or stress_audit["audited_agent3_head"] != stress["head"] or stress_audit["audited_agent3_source_blob"] != stress["source_blob"]:
        raise ValueError("A4 #1038 must target the exact registered A3 #1037 stress artifact")
    if not stress_audit["implementation_distinct"] or stress_audit["independent_operator"] == stress_audit["production_operator"]:
        raise ValueError("A4 #1038 stress audit must remain implementation-distinct")
    if stress["complete_ns_defect_evidence"] or stress["authorized_as_ns_correction_target"] or stress_audit["complete_ns_residual_evidence"] or stress_audit["authorized_as_ns_correction_target"]:
        raise ValueError("scoped X4 radial stress/audit cannot be promoted to complete-NS evidence")

    truth = registration["truth_boundary"]
    for key in (
        "current_nonlinear_mean_through_power_law_materialized",
        "agent4_1029_power_law_mean_audit_registered",
        "agent3_1037_power_law_radial_stress_registered",
        "current_nonlinear_radial_stress_through_power_law_materialized",
        "agent4_dedicated_power_law_radial_stress_audit_present",
        "agent4_independent_power_law_radial_stress_audit_registered",
    ):
        if truth.get(key) is not True:
            raise ValueError(f"required X4 stress truth regressed: {key}")
    for key in (
        "agent4_independent_power_law_radial_stress_audit_admitted",
        "current_nonlinear_radial_force_through_power_law_materialized",
        "scoped_power_law_radial_stress_authorized_as_ns_correction_target",
        "velocity_after_X_4_materialized",
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "real_agent3_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "complete_candidate_api_ready",
        "heldout_normalized_ns_residual_assessed",
        "same_protocol_st006_comparison_available_now",
        "residual_reduction_claimed",
        "scientific_admission",
        "pde_validated",
    ):
        if truth.get(key) is True:
            raise ValueError(f"truth boundary illegally promoted: {key}")

    if registration["readiness"] != parent.READINESS:
        raise ValueError("readiness state drifted")
    if registration["final_gate"] != parent.FINAL_GATE:
        raise ValueError("final scientific gate drifted")
    if registration["st006_baseline"] != parent.ST006_BASELINE:
        raise ValueError("ST006 baseline drifted")
    if registration["frozen_science"] != parent.FROZEN_SCIENCE:
        raise ValueError("frozen science contract drifted")
    if not registration["frozen_science"]["residual_defined_free_forcing_forbidden"]:
        raise ValueError("residual-defined free forcing must remain forbidden")


def materialize_registration_receipt() -> dict[str, Any]:
    registration = copy.deepcopy(_EXPECTED)
    _enforce_internal_relations(registration)
    return {"schema": SCHEMA_NAME, "task_id": TASK_ID, "registration": registration, "registration_sha256": _sha256(registration)}


def enforce_registration_receipt(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA_NAME or receipt.get("task_id") != TASK_ID:
        raise ValueError("registration schema/task drifted")
    registration = receipt.get("registration")
    if not isinstance(registration, Mapping):
        raise ValueError("registration payload missing")
    if receipt.get("registration_sha256") != _sha256(registration):
        raise ValueError("registration SHA-256 mismatch")
    if registration != _EXPECTED:
        raise ValueError("registration payload drifted from frozen contract")
    _enforce_internal_relations(registration)


def write_registration_receipt(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    receipt = materialize_registration_receipt()
    enforce_registration_receipt(receipt)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return output


def _main() -> int:
    parser = argparse.ArgumentParser(description="Materialize frozen Kokuno A5 X4 radial-stress ingest receipt.")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    print(write_registration_receipt(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
