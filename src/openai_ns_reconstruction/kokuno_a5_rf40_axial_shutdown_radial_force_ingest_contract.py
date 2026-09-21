"""Kokuno Agent-5 registration for the RF40 axial-shutdown radial-force seam.

Integration/provenance glue only.  The production velocity lineage remains the
A5 #1013 X_4 leading+oscillatory composite.  This increment advances only the
lagging correction-side lane from the A5 #1019 registered X_2 compact stress to
A3 #1022's source-required radial force

    (div T)_r = partial_z sigma_1.

No matching Agent-4 independent radial-force audit exists at freeze time.  The
A4 #1017 stress audit remains valid only for its stress object and must not be
laundered into force evidence.  The radial force is not an authorized NS
correction target.  No pressure, forcing, complete NS defect, Cartesian
correction velocity, full candidate export, or PDE admission is created here.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from openai_ns_reconstruction import kokuno_a5_rf40_axial_shutdown_radial_stress_ingest_contract as parent

SCHEMA_NAME = "kokuno-agent5-rf40-axial-shutdown-radial-force-ingest-v1"
TASK_ID = "KOKUNO-A5-RF40-AXIAL-SHUTDOWN-RADIAL-FORCE-INGEST-102"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 1019,
    "head": "b5b9b30f6a209e13f22eb8aa3b7f0a251637cc42",
    "branch": "codex/kokuno-a5-rf40-axial-shutdown-radial-stress-ingest-101",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_rf40_axial_shutdown_radial_stress_ingest_contract.py",
    "source_blob": "b5d3c845dd8daffc7cd3dc2634d7532acd1435b6",
    "test_path": "tests/test_kokuno_a5_rf40_axial_shutdown_radial_stress_ingest_contract.py",
    "test_blob": "8d9e8f581efa79f7a7c93b4f586155ade960f2f9",
    "workflow_path": ".github/workflows/kokuno-agent5-rf40-axial-shutdown-radial-stress-ingest.yml",
    "workflow_blob": "57ed5d6e5169a332041c66136879f2ee08948777",
}

AGENT3_AXIAL_SHUTDOWN_RADIAL_FORCE = {
    "pr": 1022,
    "head": "93200f28deae4372c2338e0620b785543af892d4",
    "branch": "codex/kokuno-a3-rf40-axial-shutdown-radial-force-111",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_axial_shutdown_nonlinear_radial_force.py",
    "source_blob": "222e75ab5ecb1f0625ed3b25d4e2defe6f780b8f",
    "test_path": "tests/test_constrained_kokuno_current_rf40_axial_shutdown_nonlinear_radial_force.py",
    "test_blob": "5c4aeb46ad229e94c43a3a5ebc94746363eea01a",
    "workflow_path": ".github/workflows/kokuno-agent3-current-rf40-axial-shutdown-nonlinear-radial-force.yml",
    "workflow_blob": "a70cf9b5d04d368ad886cce3fe4a35f0521a4c30",
    "parent_stress_pr": 1015,
    "parent_stress_head": "2e7683f06bea810881afbe39d8f0cf64815ddd56",
    "parent_stress_source_blob": "9c6082e6dc5edacf7b6a64e87549cbec142b3d4a",
    "agent2_composite_pr": 999,
    "agent2_composite_head": "93b99292fcf141f24b7c6d7e4fbf95775e6a07e5",
    "agent1_leading_pr": 993,
    "agent1_leading_head": "2ac6460b483efb1c07f2fa65e7fed781a32f2718",
    "source_formula": "(div T)_r = partial_z sigma_1",
    "z_derivative_step_ladder": [0.02, 0.01, 0.005],
    "piece_closure_relative_gate": 5.0e-10,
    "fine_pair_relative_stability_gate": 5.0e-2,
    "radial_force_through_X_2_materialized": True,
    "independent_third_radial_moment_inverse_introduced": False,
    "recorded_radial_mean_equated_to_radial_force": False,
    "authorized_as_ns_correction_target": False,
    "consumes_parent_x4_composite": False,
    "candidate_residual_evidence": False,
}

AGENT4_STRESS_AUDIT_ONLY = copy.deepcopy(parent.AGENT4_AXIAL_SHUTDOWN_RADIAL_STRESS_AUDIT)
AGENT4_STRESS_AUDIT_ONLY.update({
    "scope": "independent audit of A3 #1015 compact radial stress only",
    "audits_agent3_1022_radial_force": False,
    "radial_force_independent_audit_authority": False,
})

OBSERVED_CI_AT_FREEZE = {
    "parent_a5_1019": {
        "repository_tests_run": 35623422702,
        "dedicated_run": 35623422484,
        "status": "queued",
        "conclusion": None,
    },
    "agent3_1022_axial_shutdown_radial_force": {
        "repository_tests_run": 35627421903,
        "dedicated_run": 35627422067,
        "status": "queued",
        "conclusion": None,
    },
    "agent4_1017_stress_audit_only": {
        "repository_tests_run": 35622184081,
        "dedicated_run": 35622184156,
        "status": "queued",
        "conclusion": None,
    },
}

FROZEN_SCIENCE = copy.deepcopy(parent.FROZEN_SCIENCE)
FINAL_GATE = copy.deepcopy(parent.FINAL_GATE)
ST006_BASELINE = copy.deepcopy(parent.ST006_BASELINE)
READINESS = copy.deepcopy(parent.READINESS)

TRUTH_BOUNDARY = copy.deepcopy(parent.TRUTH_BOUNDARY)
TRUTH_BOUNDARY.update({
    "current_nonlinear_radial_force_through_axial_shutdown_materialized": True,
    "agent3_1022_axial_shutdown_radial_force_registered": True,
    "agent4_1017_stress_audit_does_not_audit_radial_force": True,
    "agent4_dedicated_axial_shutdown_radial_force_audit_present": False,
    "agent4_independent_axial_shutdown_radial_force_audit_registered": False,
    "agent4_independent_axial_shutdown_radial_force_audit_admitted": False,
    "scoped_axial_shutdown_radial_force_authorized_as_ns_correction_target": False,
})

PIPELINE_POSITION = {
    "stage": "X_4 identity-bound production velocity preserved; correction-side X_2 mean->stress->radial-force registered",
    "production_velocity": "A1 #1005 -> A2 #1010 through RF40 power-law X_4; A4 #1012 remains the scoped velocity-divergence validator",
    "correction_side": "A3 #1022 consumes the exact A3 #1015 X_2 compact axial stress and materializes partial_z sigma_1; A4 #1017 audits only the parent stress",
    "new_output": "checksum-bound X_2 radial-force registration with an explicit no-independent-A4-force-audit firewall",
    "not_output": "X_4 A3 correction-side coverage, independent force audit, pressure, forcing, complete NS defect, correction authorization, Cartesian correction velocity, finite cycle, full candidate, or PDE validation",
    "next_shortest_blocker": "obtain an implementation-distinct Agent-4 audit of the X_2 radial force and carry A3 mean/stress/force plumbing toward the exact X_4 production identity while A1/A2 continue outer/global completion",
}

_EXPECTED = {
    "parent_a5": PARENT_A5,
    "agent3_axial_shutdown_radial_force": AGENT3_AXIAL_SHUTDOWN_RADIAL_FORCE,
    "agent4_stress_audit_only": AGENT4_STRESS_AUDIT_ONLY,
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
    parent_a5 = registration["parent_a5"]
    a3 = registration["agent3_axial_shutdown_radial_force"]
    a4 = registration["agent4_stress_audit_only"]
    for label, obj in (("parent_a5", parent_a5), ("agent3", a3), ("agent4", a4)):
        _require_hex40(obj["head"], f"{label}.head")
        for blob_key in ("source_blob", "test_blob", "workflow_blob"):
            _require_hex40(obj[blob_key], f"{label}.{blob_key}")

    if a3["parent_stress_pr"] != parent.AGENT3_AXIAL_SHUTDOWN_RADIAL_STRESS["pr"]:
        raise ValueError("A3 #1022 parent stress PR mismatch")
    if a3["parent_stress_head"] != parent.AGENT3_AXIAL_SHUTDOWN_RADIAL_STRESS["head"]:
        raise ValueError("A3 #1022 parent stress head mismatch")
    if a3["parent_stress_source_blob"] != parent.AGENT3_AXIAL_SHUTDOWN_RADIAL_STRESS["source_blob"]:
        raise ValueError("A3 #1022 parent stress source mismatch")
    if not a3["radial_force_through_X_2_materialized"]:
        raise ValueError("A3 #1022 X_2 radial force must remain materialized")
    if a3["independent_third_radial_moment_inverse_introduced"]:
        raise ValueError("radial force must remain partial_z sigma_1, not a third inverse")
    if a3["recorded_radial_mean_equated_to_radial_force"]:
        raise ValueError("recorded radial mean cannot be relabeled as radial force")
    if a3["authorized_as_ns_correction_target"] or a3["candidate_residual_evidence"]:
        raise ValueError("scoped radial force cannot be promoted to NS correction/residual evidence")
    if a3["consumes_parent_x4_composite"]:
        raise ValueError("X_2 A3 radial force cannot be laundered into X_4 consumption")

    if a4["audits_agent3_1022_radial_force"] or a4["radial_force_independent_audit_authority"]:
        raise ValueError("A4 #1017 stress audit cannot be promoted into radial-force evidence")
    if a4["audited_agent3_pr"] != 1015:
        raise ValueError("A4 #1017 must remain bound to A3 #1015 stress")

    truth = registration["truth_boundary"]
    required_true = (
        "current_nonlinear_radial_stress_through_axial_shutdown_materialized",
        "current_nonlinear_radial_force_through_axial_shutdown_materialized",
        "agent3_1022_axial_shutdown_radial_force_registered",
        "agent4_1017_stress_audit_does_not_audit_radial_force",
    )
    for key in required_true:
        if not truth[key]:
            raise ValueError(f"required truth boundary regressed: {key}")
    forbidden_true = (
        "agent4_dedicated_axial_shutdown_radial_force_audit_present",
        "agent4_independent_axial_shutdown_radial_force_audit_registered",
        "agent4_independent_axial_shutdown_radial_force_audit_admitted",
        "scoped_axial_shutdown_radial_force_authorized_as_ns_correction_target",
        "current_nonlinear_mean_through_power_law_materialized",
        "current_nonlinear_radial_stress_through_power_law_materialized",
        "current_nonlinear_radial_force_through_power_law_materialized",
        "velocity_after_X_4_materialized",
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "current_nonlinear_mean_authorized_as_correction_target",
        "real_agent3_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "complete_candidate_api_ready",
        "heldout_normalized_ns_residual_assessed",
        "same_protocol_st006_comparison_available_now",
        "residual_reduction_claimed",
        "scientific_admission",
        "pde_validated",
    )
    for key in forbidden_true:
        if key in truth and truth[key]:
            raise ValueError(f"truth boundary illegally promoted: {key}")

    if registration["readiness"] != parent.READINESS:
        raise ValueError("readiness state drifted")
    if registration["final_gate"] != parent.FINAL_GATE:
        raise ValueError("final scientific gate drifted")
    if registration["st006_baseline"] != parent.ST006_BASELINE:
        raise ValueError("ST006 baseline drifted")
    if not registration["frozen_science"]["residual_defined_free_forcing_forbidden"]:
        raise ValueError("residual-defined free forcing must remain forbidden")


def materialize_registration_receipt() -> dict[str, Any]:
    registration = copy.deepcopy(_EXPECTED)
    _enforce_internal_relations(registration)
    return {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "registration": registration,
        "registration_sha256": _sha256(registration),
    }


def enforce_registration_receipt(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA_NAME:
        raise ValueError("registration schema drifted")
    if receipt.get("task_id") != TASK_ID:
        raise ValueError("registration task id drifted")
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
    parser = argparse.ArgumentParser(
        description="Materialize the frozen Kokuno A5 axial-shutdown radial-force registration receipt."
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    path = write_registration_receipt(args.output)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
