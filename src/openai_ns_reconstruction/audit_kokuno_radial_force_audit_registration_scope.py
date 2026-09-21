from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

CONTRACT_REL = Path("configs/kokuno_radial_force_audit_registration_scope.json")
CONSTRAINTS_REL = Path("configs/constraints.json")
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"


def _git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def audit_contract(contract: dict[str, Any]) -> list[str]:
    errors: list[str] = []

    if contract.get("contract_id") != "CR002-KOKUNO-RADIAL-FORCE-AUDIT-REGISTRATION-SCOPE-108":
        errors.append("unexpected contract_id")

    post = contract.get("post_freeze_evidence", {})
    expected_post = {
        "agent3_radial_force_pr": 1022,
        "agent3_radial_force_head": "93200f28deae4372c2338e0620b785543af892d4",
        "agent4_radial_force_audit_pr": 1023,
        "agent4_radial_force_audit_head": "6865031eb5084d2bcbfe6832c23c73b22b2f8c21",
        "agent4_audit_exists": True,
        "agent5_1024_consumes_agent4_1023": False,
        "agent5_1024_registered_radial_force_audit": False,
    }
    for key, value in expected_post.items():
        if post.get(key) != value:
            errors.append(f"post_freeze_evidence.{key} must equal {value!r}")

    allowed = contract.get("allowed_true", {})
    for key in (
        "a3_1022_radial_force_materialized_through_x2",
        "a4_1023_independent_radial_force_audit_artifact_exists",
        "a4_1023_uses_implementation_distinct_fd4_reference",
        "a5_1024_registers_a3_1022_radial_force",
    ):
        if allowed.get(key) is not True:
            errors.append(f"allowed_true.{key} must remain true")

    required_false = contract.get("required_false", {})
    for key, value in required_false.items():
        if value is not False:
            errors.append(f"required_false.{key} must remain false")

    semantics = contract.get("evidence_semantics", {})
    if "later repository evidence" not in semantics.get("snapshot_truth_rule", ""):
        errors.append("snapshot truth rule must protect later repository evidence")
    if "distinct" not in semantics.get("registration_rule", ""):
        errors.append("registration rule must distinguish artifact existence from consumption")
    if "queued/running" not in semantics.get("ci_rule", ""):
        errors.append("CI rule must preserve queued/running truth")

    cr001 = contract.get("cr001_snapshot", {})
    expected_cr001 = {
        "constraints_blob": EXPECTED_CONSTRAINTS_BLOB,
        "nu": 0.01,
        "physical_domain": "R^3",
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25, 0.75],
        "forcing_mode": "restricted_two_parameter_family",
        "validation_seed": 914027,
        "held_out_points": 4096,
        "derivative_steps": [0.02, 0.01, 0.005],
        "quadrature_orders_per_axis": [24, 48, 96],
        "momentum_max_l2_threshold": 0.001,
        "divergence_max_l2_threshold": 0.00001,
        "residual_defined_free_force_forbidden": True,
        "candidate_collapse_forbidden": True,
        "post_hoc_threshold_relaxation_forbidden": True,
    }
    for key, value in expected_cr001.items():
        if cr001.get(key) != value:
            errors.append(f"cr001_snapshot.{key} drifted")
    if cr001.get("evaluation_box") != [[-2, 2], [-2, 2], [-2, 2]]:
        errors.append("CR001 evaluation box drifted")

    eq45 = contract.get("canonical_delivery_independence", {})
    expected_eq45 = {
        "candidate": "eq45_supported_velocity_candidate_v1",
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    for key, value in expected_eq45.items():
        if eq45.get(key) != value:
            errors.append(f"canonical_delivery_independence.{key} drifted")

    return errors


def audit_repository(root: Path) -> list[str]:
    contract = _load_json(root / CONTRACT_REL)
    errors = audit_contract(contract)
    constraints_path = root / CONSTRAINTS_REL
    if not constraints_path.exists():
        errors.append("canonical configs/constraints.json missing")
        return errors
    actual_blob = _git_blob_sha(constraints_path.read_bytes())
    if actual_blob != EXPECTED_CONSTRAINTS_BLOB:
        errors.append(
            f"canonical constraints blob drifted: expected {EXPECTED_CONSTRAINTS_BLOB}, got {actual_blob}"
        )
    return errors


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    errors = audit_repository(root)
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    print("PASS: CR002 radial-force audit registration scope is fail-closed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
