"""Kokuno Agent-5 registration for the RF40 axial-shutdown radial-stress seam.

Integration/provenance glue only.  The production velocity lineage remains the
A5 #1013 X_4 leading+oscillatory composite.  This increment advances only the
lagging correction-side lane:

    A3 #1015: current X_2 nonlinear m=0 means -> compact theta/axial stress
      -> A4 #1017: implementation-distinct public-receipt stress audit.

The stress is not an authorized NS correction target.  No radial-force step,
pressure, forcing, complete NS defect, Cartesian correction velocity, full
candidate export, or PDE admission is created here.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from openai_ns_reconstruction import kokuno_a5_current_rf40_power_law_composite_ingest_contract as parent

SCHEMA_NAME = "kokuno-agent5-rf40-axial-shutdown-radial-stress-ingest-v1"
TASK_ID = "KOKUNO-A5-RF40-AXIAL-SHUTDOWN-RADIAL-STRESS-INGEST-101"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 1013,
    "head": "df75e6994ca4ae1fa405d56ef162b5a89d5e6d58",
    "branch": "codex/kokuno-a5-rf40-power-law-composite-ingest-100",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_rf40_power_law_composite_ingest_contract.py",
    "source_blob": "3d7c80f5f8e7f46f3da39402ea163f13b9120a8e",
    "test_path": "tests/test_kokuno_a5_current_rf40_power_law_composite_ingest_contract.py",
    "test_blob": "ad2dae6c862b5170314f95f800fe5c356867ee15",
    "workflow_path": ".github/workflows/kokuno-agent5-current-rf40-power-law-composite-ingest.yml",
    "workflow_blob": "7e4ab337ba958c5981033588c138d7f4a57afa4a",
}

AGENT3_AXIAL_SHUTDOWN_RADIAL_STRESS = {
    "pr": 1015,
    "head": "2e7683f06bea810881afbe39d8f0cf64815ddd56",
    "branch": "codex/kokuno-a3-rf40-axial-shutdown-radial-stress-110",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_axial_shutdown_nonlinear_radial_stress.py",
    "source_blob": "9c6082e6dc5edacf7b6a64e87549cbec142b3d4a",
    "test_path": "tests/test_constrained_kokuno_current_rf40_axial_shutdown_nonlinear_radial_stress.py",
    "test_blob": "46ff23d00f666c6d231f77dcf539e226f53c0658",
    "workflow_path": ".github/workflows/kokuno-agent3-current-rf40-axial-shutdown-nonlinear-radial-stress.yml",
    "workflow_blob": "a32c510dbc354bd9ded8108bd964b3bd19a30968",
    "parent_mean_pr": 1011,
    "parent_mean_head": "6aabf7b6dd28ae683a3776d518e9a977507911e9",
    "agent2_composite_pr": 999,
    "agent2_composite_head": "93b99292fcf141f24b7c6d7e4fbf95775e6a07e5",
    "agent1_leading_pr": 993,
    "agent1_leading_head": "2ac6460b483efb1c07f2fa65e7fed781a32f2718",
    "materialized_domain": "nonlinear m=0 mean -> compact tangential e=2 and axial e=1 radial stress through RF40 axial shutdown X_2",
    "radial_force_through_X_2_materialized": False,
    "authorized_as_ns_correction_target": False,
    "consumes_parent_x4_composite": False,
}

AGENT4_AXIAL_SHUTDOWN_RADIAL_STRESS_AUDIT = {
    "pr": 1017,
    "head": "cc3109bd17d518f3414d972f2303a20f2dabc1b1",
    "branch": "codex/kokuno-a4-rf40-axial-shutdown-radial-stress-audit-103",
    "audited_agent3_pr": 1015,
    "audited_agent3_head": "2e7683f06bea810881afbe39d8f0cf64815ddd56",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_rf40_axial_shutdown_radial_stress_independent_audit.py",
    "source_blob": "c886a7e9d222313125f5e242ca9cc9828babdb23",
    "test_path": "tests/test_constrained_kokuno_a4_current_rf40_axial_shutdown_radial_stress_independent_audit.py",
    "test_blob": "102f5072c97815e808dacf9b9faf03bbb351a6a6",
    "workflow_path": ".github/workflows/kokuno-agent4-current-rf40-axial-shutdown-radial-stress-independent-audit.yml",
    "workflow_blob": "9e9ff1b12e3acd42753be5abd3b8cc97fe9771fe",
    "reference_path": "JSON-round-tripped public A3 receipt -> local piecewise-cubic reconstruction + order-8 Gauss-Legendre cell integration",
    "seed": 9173751,
    "time": 0.39,
    "z": -0.08,
    "radial_interval": [0.005, 0.44],
    "nested_counts": [43, 85, 169],
    "offgrid_count": 169,
    "fine_and_offgrid_relative_rms_gate": 5.0e-2,
    "relative_max_gate": 1.5e-1,
    "weighted_moment_relative_error_gate": 3.0e-2,
    "axis_near_normalized_error_gate": 1.5e-1,
    "outer_edge_normalized_stress_gate": 1.0e-8,
    "aggregate_stress_nontriviality_floor": 1.0e-12,
    "uses_a3_production_radial_inverse_for_reference": False,
    "momentum_or_full_ns_evidence": False,
    "scientific_admission": False,
}

OBSERVED_CI_AT_FREEZE = {
    "parent_a5_1013": {"repository_tests_run": 35616395901, "dedicated_run": 35616395776, "status": "queued", "conclusion": None},
    "agent3_1015_axial_shutdown_radial_stress": {"repository_tests_run": 35620634573, "dedicated_run": 35620634575, "status": "queued", "conclusion": None},
    "agent4_1017_axial_shutdown_radial_stress_audit": {"repository_tests_run": 35622184081, "dedicated_run": 35622184156, "status": "queued", "conclusion": None},
}

FROZEN_SCIENCE = copy.deepcopy(parent.FROZEN_SCIENCE)
FINAL_GATE = copy.deepcopy(parent.FINAL_GATE)
ST006_BASELINE = copy.deepcopy(parent.ST006_BASELINE)
READINESS = copy.deepcopy(parent.READINESS)

TRUTH_BOUNDARY = copy.deepcopy(parent.TRUTH_BOUNDARY)
TRUTH_BOUNDARY.update({
    "current_nonlinear_radial_stress_through_axial_shutdown_materialized": True,
    "agent3_1015_axial_shutdown_radial_stress_registered": True,
    "agent4_1017_independent_axial_shutdown_radial_stress_audit_registered": True,
    "agent4_1017_independent_axial_shutdown_radial_stress_audit_admitted": False,
})

PIPELINE_POSITION = {
    "stage": "X_4 identity-bound velocity seam preserved; correction-side X_2 mean->stress plus matching independent A4 operator audit registered",
    "production_velocity": "A1 #1005 -> A2 #1010 through RF40 power-law X_4; A4 #1012 scoped divergence audit remains the current velocity validator",
    "correction_side": "A3 #1015 consumes A2 #999/A1 #993 through X_2 and materializes tangential/axial compact radial stresses; A4 #1017 independently audits that public receipt",
    "new_output": "checksum-bound X_2 radial-stress registration + implementation-distinct A4 stress-audit registration",
    "not_output": "X_2 radial force, X_4 A3 stress, pressure, forcing, complete NS defect, correction authorization, Cartesian correction velocity, full candidate, or PDE validation",
    "next_shortest_blocker": "carry A3 mean/radial plumbing from X_2 to the exact X_4 production composite identity while A1/A2 continue cone/I1-I4 global completion; only then bind matched pressure and preregistered restricted forcing",
}

_EXPECTED = {
    "parent_a5": PARENT_A5,
    "agent3_axial_shutdown_radial_stress": AGENT3_AXIAL_SHUTDOWN_RADIAL_STRESS,
    "agent4_axial_shutdown_radial_stress_audit": AGENT4_AXIAL_SHUTDOWN_RADIAL_STRESS_AUDIT,
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
    for label, obj in (
        ("parent_a5", registration["parent_a5"]),
        ("agent3", registration["agent3_axial_shutdown_radial_stress"]),
        ("agent4", registration["agent4_axial_shutdown_radial_stress_audit"]),
    ):
        _require_hex40(obj["head"], f"{label}.head")
        for blob_key in ("source_blob", "test_blob", "workflow_blob"):
            _require_hex40(obj[blob_key], f"{label}.{blob_key}")

    a3 = registration["agent3_axial_shutdown_radial_stress"]
    a4 = registration["agent4_axial_shutdown_radial_stress_audit"]
    if a4["audited_agent3_pr"] != a3["pr"] or a4["audited_agent3_head"] != a3["head"]:
        raise ValueError("A4/A3 axial-shutdown stress lineage mismatch")
    if a3["consumes_parent_x4_composite"]:
        raise ValueError("X_2 A3 stress cannot be laundered into X_4 consumption")
    if a3["radial_force_through_X_2_materialized"]:
        raise ValueError("A3 #1015 does not materialize the X_2 radial-force derivative")
    if a3["authorized_as_ns_correction_target"]:
        raise ValueError("scoped stress cannot be authorized before complete NS defect")
    if a4["uses_a3_production_radial_inverse_for_reference"]:
        raise ValueError("A4 reference must remain implementation-distinct")
    if a4["momentum_or_full_ns_evidence"] or a4["scientific_admission"]:
        raise ValueError("scoped stress audit cannot be promoted to full-NS evidence")

    truth = registration["truth_boundary"]
    if not truth["current_nonlinear_radial_stress_through_axial_shutdown_materialized"]:
        raise ValueError("registered X_2 radial stress must remain materialized")
    if not truth["agent4_1017_independent_axial_shutdown_radial_stress_audit_registered"]:
        raise ValueError("matching A4 #1017 audit registration cannot be dropped")
    if truth["agent4_1017_independent_axial_shutdown_radial_stress_audit_admitted"]:
        raise ValueError("unresolved exact-head CI cannot be scientific admission")
    forbidden_true = (
        "current_nonlinear_radial_force_through_axial_shutdown_materialized",
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
    return {"schema": SCHEMA_NAME, "task_id": TASK_ID, "registration": registration, "registration_sha256": _sha256(registration)}


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
    parser = argparse.ArgumentParser(description="Materialize the frozen Kokuno A5 axial-shutdown radial-stress registration receipt.")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    path = write_registration_receipt(args.output)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
