"""Kokuno Agent-5 registration of the independent X2 radial-force audit.

Integration/provenance glue only.  Parent A5 #1024 already registers the
current correction-side chain through RF40 axial shutdown X_2:

    nonlinear m=0 means -> compact theta/axial stresses -> partial_z sigma_1.

Agent-4 PR #1023 was opened after that freeze and independently audits the
public A3 #1022 radial-force receipt with a centered five-point FD4 derivative,
not A3's production centered-FD2 helper.  This module registers that later
validator identity without scientific-admitting unresolved CI and without
promoting the scoped X_2 force into a complete-NS correction target.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from openai_ns_reconstruction import kokuno_a5_rf40_axial_shutdown_radial_force_ingest_contract as parent

SCHEMA_NAME = "kokuno-agent5-rf40-axial-shutdown-radial-force-audit-registration-v1"
TASK_ID = "KOKUNO-A5-RF40-AXIAL-SHUTDOWN-RADIAL-FORCE-AUDIT-REG-103"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 1024,
    "head": "931bf8ae11c7ed0f7118c2f868a657d2ee16863b",
    "branch": "codex/kokuno-a5-rf40-axial-shutdown-radial-force-ingest-102",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_rf40_axial_shutdown_radial_force_ingest_contract.py",
    "source_blob": "3bcdb6f3b3d1ee16f6e4c9e2bd53590bf9bbdb51",
    "test_path": "tests/test_kokuno_a5_rf40_axial_shutdown_radial_force_ingest_contract.py",
    "test_blob": "5181287fe5db7244a8c9146292c024c69f995b7a",
    "workflow_path": ".github/workflows/kokuno-agent5-rf40-axial-shutdown-radial-force-ingest.yml",
    "workflow_blob": "fa140eff1ed63b0fa4cde1a8e5fae5442982d7bb",
}

AGENT4_RADIAL_FORCE_AUDIT = {
    "pr": 1023,
    "head": "6865031eb5084d2bcbfe6832c23c73b22b2f8c21",
    "branch": "codex/kokuno-a4-rf40-axial-shutdown-radial-force-independent-audit-104",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_rf40_axial_shutdown_radial_force_independent_audit.py",
    "source_blob": "9dc41bcc1bacffb714acbcaada0ab64423c8d446",
    "test_path": "tests/test_constrained_kokuno_a4_current_rf40_axial_shutdown_radial_force_independent_audit.py",
    "test_blob": "9f1aa8ae19fc711845bb4662207bc2cc3e15b154",
    "workflow_path": ".github/workflows/kokuno-agent4-current-rf40-axial-shutdown-radial-force-independent-audit.yml",
    "workflow_blob": "45e3c611a1da5e05d7a37f25f82204c3f0929c5f",
    "audited_agent3_pr": 1022,
    "audited_agent3_head": "93200f28deae4372c2338e0620b785543af892d4",
    "audited_agent3_source_blob": "222e75ab5ecb1f0625ed3b25d4e2defe6f780b8f",
    "independent_operator": "centered_five_point_fd4_on_public_axial_stress_receipts",
    "production_operator": "centered_fd2_z_derivative_ladder",
    "production_z_steps": [0.02, 0.01, 0.005],
    "frozen_seed": 9173761,
    "frozen_time": 0.39,
    "frozen_z_centers": [-0.137, -0.073, 0.041],
    "radial_interval": [0.005, 0.44],
    "radial_count": 169,
    "fine_relative_rms_gate": 5.0e-2,
    "fine_relative_max_gate": 1.5e-1,
    "fine_relative_integral_l2_gate": 5.0e-2,
    "axis_near_normalized_error_gate": 1.5e-1,
    "independent_medium_to_fine_relative_rms_gate": 5.0e-2,
    "piece_closure_relative_gate": 5.0e-10,
    "nontrivial_force_rms_floor": 1.0e-12,
    "implementation_distinct": True,
    "complete_ns_residual_evidence": False,
    "authorized_as_ns_correction_target": False,
}

OBSERVED_CI_AT_FREEZE = {
    "parent_a5_1024": {
        "repository_tests_run": 35628986390,
        "dedicated_run": 35628986367,
        "status": "queued",
        "conclusion": None,
    },
    "agent3_1022_radial_force": {
        "repository_tests_run": 35627421903,
        "dedicated_run": 35627422067,
        "status": "queued",
        "conclusion": None,
    },
    "agent4_1023_radial_force_audit": {
        "repository_tests_run": 35628615272,
        "dedicated_run": 35628615181,
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
    "agent4_dedicated_axial_shutdown_radial_force_audit_present": True,
    "agent4_independent_axial_shutdown_radial_force_audit_registered": True,
    "agent4_independent_axial_shutdown_radial_force_audit_admitted": False,
    "scoped_axial_shutdown_radial_force_authorized_as_ns_correction_target": False,
})

PIPELINE_POSITION = {
    "stage": "production velocity remains identity-bound through X_4; correction-side X_2 radial force now has a registered implementation-distinct A4 audit",
    "production_velocity": "A1 #1005 -> A2 #1010 through RF40 power-law X_4; A4 #1012 remains the scoped velocity-divergence validator",
    "correction_side": "A3 #1022 X_2 mean->stress->partial_z sigma_1 plus A4 #1023 FD4 audit of the public radial-force/stress receipts",
    "new_output": "checksum-bound registration of A4 #1023 as the dedicated independent audit for the exact A3 #1022 radial-force seam",
    "not_output": "scientific admission, X_4 stress/force coverage, pressure, restricted forcing, complete NS defect, correction authorization, Cartesian correction velocity, finite cycle, full candidate, or PDE validation",
    "next_shortest_blocker": "carry A3 mean->stress->force plumbing onto the exact X_4 production composite while A1/A2 continue outer/global completion; only a global candidate plus matched pressure and preregistered restricted forcing can define the complete NS defect",
}

_EXPECTED = {
    "parent_a5": PARENT_A5,
    "agent3_radial_force": copy.deepcopy(parent.AGENT3_AXIAL_SHUTDOWN_RADIAL_FORCE),
    "agent4_radial_force_audit": AGENT4_RADIAL_FORCE_AUDIT,
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
    p = registration["parent_a5"]
    a3 = registration["agent3_radial_force"]
    a4 = registration["agent4_radial_force_audit"]
    for label, obj in (("parent_a5", p), ("agent3", a3), ("agent4", a4)):
        _require_hex40(obj["head"], f"{label}.head")
        for blob_key in ("source_blob", "test_blob", "workflow_blob"):
            _require_hex40(obj[blob_key], f"{label}.{blob_key}")

    if a4["audited_agent3_pr"] != a3["pr"]:
        raise ValueError("A4 audit must target exact registered A3 radial-force PR")
    if a4["audited_agent3_head"] != a3["head"]:
        raise ValueError("A4 audit must target exact registered A3 radial-force head")
    if a4["audited_agent3_source_blob"] != a3["source_blob"]:
        raise ValueError("A4 audit must target exact registered A3 radial-force source")
    if not a4["implementation_distinct"]:
        raise ValueError("A4 radial-force audit must remain implementation-distinct")
    if a4["independent_operator"] == a4["production_operator"]:
        raise ValueError("independent radial-force operator cannot collapse to production operator")
    if a4["complete_ns_residual_evidence"] or a4["authorized_as_ns_correction_target"]:
        raise ValueError("scoped radial-force audit cannot be promoted to complete-NS evidence")

    truth = registration["truth_boundary"]
    required_true = (
        "current_nonlinear_radial_force_through_axial_shutdown_materialized",
        "agent3_1022_axial_shutdown_radial_force_registered",
        "agent4_dedicated_axial_shutdown_radial_force_audit_present",
        "agent4_independent_axial_shutdown_radial_force_audit_registered",
    )
    for key in required_true:
        if truth.get(key) is not True:
            raise ValueError(f"required truth boundary regressed: {key}")
    forbidden_true = (
        "agent4_independent_axial_shutdown_radial_force_audit_admitted",
        "scoped_axial_shutdown_radial_force_authorized_as_ns_correction_target",
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
        description="Materialize the frozen Kokuno A5 axial-shutdown radial-force audit registration receipt."
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    path = write_registration_receipt(args.output)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
