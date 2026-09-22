"""Fail-closed CR002 audit for exact #1154 relative-swirl root-branch scope."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping
import numpy as np

EXPECTED_SCHEMA = "cr002-kokuno-relative-swirl-root-branch-scope-v1"
EXPECTED_TASK = "CR002-KOKUNO-RELATIVE-SWIRL-ROOT-BRANCH-SCOPE-141"
EXPECTED_PARENT_HEAD = "d0e855e37c8f5d6d58f1a85b1dcb874651d20838"
EXPECTED_PARENT_BLOB = "7715dbf9fadb5d9ed37213677703b45f21d29724"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"

def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload

def _git_blob_sha(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()

def _eq(actual: Any, expected: Any, label: str, errors: list[str]) -> None:
    if actual != expected:
        errors.append(f"{label}: expected {expected!r}, got {actual!r}")

def _true(mapping: Mapping[str, Any], key: str, errors: list[str]) -> None:
    if mapping.get(key) is not True:
        errors.append(f"{key} must remain true")

def _false(mapping: Mapping[str, Any], key: str, errors: list[str]) -> None:
    if mapping.get(key) is not False:
        errors.append(f"{key} must remain false")

def multi_root_witness(angular_row: Any, pressure_linear_row: Any, pressure_quadratic_matrix: Any, target: float, root_a: Any, root_b: Any) -> dict[str, float | bool]:
    """Autonomous algebra witness only; not Kokuno/OpenAI/PDE evidence."""
    a = np.asarray(angular_row, dtype=float)
    b = np.asarray(pressure_linear_row, dtype=float)
    q = np.asarray(pressure_quadratic_matrix, dtype=float)
    ca = np.asarray(root_a, dtype=float)
    cb = np.asarray(root_b, dtype=float)
    if a.shape != (2,) or b.shape != (2,) or q.shape != (2, 2) or ca.shape != (2,) or cb.shape != (2,):
        raise ValueError("witness rows/matrix/roots have wrong shape")
    vals = np.concatenate((a, b, q.reshape(-1), ca, cb, np.array([target], dtype=float)))
    if not np.all(np.isfinite(vals)):
        raise ValueError("witness inputs must be finite")
    def residual(c: np.ndarray) -> tuple[float, float]:
        return float(a @ c - target), float(2.0 * (b @ c) + c @ q @ c)
    ar_a, pr_a = residual(ca)
    ar_b, pr_b = residual(cb)
    return {
        "root_distance": float(np.linalg.norm(ca - cb)),
        "root_a_angular_abs": abs(ar_a),
        "root_a_pressure_abs": abs(pr_a),
        "root_b_angular_abs": abs(ar_b),
        "root_b_pressure_abs": abs(pr_b),
        "two_distinct_exact_roots": float(np.linalg.norm(ca - cb)) > 0.0 and max(abs(ar_a), abs(pr_a), abs(ar_b), abs(pr_b)) <= 1e-12,
    }

def _audit_parent(path: Path, errors: list[str]) -> None:
    if _git_blob_sha(path) != EXPECTED_PARENT_BLOB:
        errors.append("exact #1154 candidate module blob changed")
        return
    text = path.read_text(encoding="utf-8")
    required = (
        '"quadratic_root": "real pressure root closest to the zero-target linearized branch"',
        'candidates.append',
        'linearized = np.linalg.solve',
        'coeff = min(candidates, key=lambda c: float(np.linalg.norm(c - linearized)))',
        'if float(np.sum(np.abs(coeff))) >= 1.0:',
        '"current_lineage_angular_entry_I_materialized": False',
        '"current_cartesian_relative_swirl_composed": False',
        '"restricted_forcing_materialized": False',
        '"heldout_ns_residual_assessed": False',
        '"pde_validated": False',
    )
    for fragment in required:
        if fragment not in text:
            errors.append(f"#1154 root/truth fragment missing: {fragment}")

def audit_contract(contract_path: Path, constraints_path: Path, candidate_module_path: Path, project_status_path: Path) -> list[str]:
    errors: list[str] = []
    contract = _load_json(contract_path)
    _eq(contract.get("schema"), EXPECTED_SCHEMA, "schema", errors)
    _eq(contract.get("task_id"), EXPECTED_TASK, "task_id", errors)
    parent = contract.get("exact_parent", {})
    _eq(parent.get("pr"), 1154, "exact_parent.pr", errors)
    _eq(parent.get("head_sha"), EXPECTED_PARENT_HEAD, "exact_parent.head_sha", errors)
    _eq(parent.get("candidate_module_blob_sha"), EXPECTED_PARENT_BLOB, "parent blob", errors)

    provenance = contract.get("four_way_provenance", {})
    classes = {"user_requirements", "public_source_facts", "autonomous_design", "pending_unknown"}
    _eq(set(provenance), classes, "four-way provenance classes", errors)
    for key in classes:
        if not isinstance(provenance.get(key), list) or not provenance.get(key):
            errors.append(f"four_way_provenance.{key} must be nonempty")

    scope = contract.get("root_branch_scope", {})
    for key in ("public_relative_swirl_two_row_algebra_materialized", "repository_autonomous_pointwise_bump_materialized", "repository_autonomous_quadrature_materialized", "current_root_selection_is_repository_autonomous", "current_l1_guard_is_repository_autonomous"):
        _true(scope, key, errors)
    _eq(scope.get("current_l1_guard_strict_upper_bound"), 1.0, "current L1 guard", errors)
    for key in ("source_exact_numeric_coefficients_recovered", "current_lineage_angular_entry_I_materialized", "current_cartesian_relative_swirl_composed", "complete_terminal_global_leading_materialized", "unified_global_kokuno_velocity_export_ready", "matched_global_pressure_materialized", "restricted_forcing_materialized", "heldout_complete_ns_residual_assessed", "pde_validated", "paper_exact", "openai_field_identified"):
        _false(scope, key, errors)

    promotion = contract.get("branch_rule_change_promotion", {})
    for key in ("change_root_selection_requires_new_semantic_candidate_identity", "change_numeric_guard_requires_new_semantic_candidate_identity", "change_pointwise_bump_requires_new_semantic_candidate_identity", "fresh_save_load_replay_required", "fresh_cartesian_composition_audit_required_if_composed", "fresh_downstream_numerical_audit_required"):
        _true(promotion, key, errors)
    _false(promotion, "reuse_old_coefficient_receipts_after_branch_rule_change", errors)
    _false(promotion, "reuse_old_cartesian_or_divergence_receipts_after_branch_rule_change", errors)

    firewall = contract.get("evidence_transfer_firewall", {}).get("forbidden_implications", [])
    joined = "\n".join(map(str, firewall))
    for fragment in ("unique numerical c1,c2 branch", "source-exact branch recovery", "public/source-exact numerical coefficient bound", "entering I or r_I materialized", "composed into Cartesian velocity", "terminal/exterior/global velocity completion", "held-out complete NS residual or PDE validation", "paper exactness or exact OpenAI-field identity"):
        if fragment not in joined:
            errors.append(f"missing forbidden implication fragment: {fragment}")

    witness = contract.get("mechanics_witness", {})
    _eq(witness.get("classification"), "autonomous_mechanics_only", "witness classification", errors)
    try:
        receipt = multi_root_witness(witness.get("angular_row"), witness.get("pressure_linear_row"), witness.get("pressure_quadratic_matrix"), witness.get("angular_target"), witness.get("root_a"), witness.get("root_b"))
    except (TypeError, ValueError) as exc:
        errors.append(f"invalid mechanics witness: {exc}")
    else:
        if receipt["two_distinct_exact_roots"] is not True:
            errors.append("mechanics witness must exhibit two distinct exact roots")

    cr001 = contract.get("cr001_invariants", {})
    expected = {
        "constraints_blob_sha": EXPECTED_CONSTRAINTS_BLOB,
        "nu": 0.01,
        "physical_domain": "R^3",
        "evaluation_box": [[-2,2],[-2,2],[-2,2]],
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25,0.75],
        "forcing_mode": "restricted_two_parameter_family",
        "reference_energy": 1.0,
        "reference_energy_abs_tolerance": 0.001,
        "validation_seed": 914027,
        "held_out_points": 4096,
        "derivative_steps": [0.02,0.01,0.005],
        "quadrature_orders_per_axis": [24,48,96],
        "thresholds": {"divergence_max":1e-5,"divergence_L2":1e-5,"pde_residual_max":1e-3,"pde_residual_L2":1e-3},
    }
    for key, value in expected.items():
        _eq(cr001.get(key), value, f"cr001.{key}", errors)
    for key in ("forcing_family_change_this_increment", "threshold_relaxation_this_increment", "candidate_velocity_change_this_increment"):
        _false(cr001, key, errors)

    if _git_blob_sha(constraints_path) != EXPECTED_CONSTRAINTS_BLOB:
        errors.append("canonical configs/constraints.json blob changed")
    constraints = _load_json(constraints_path)
    _eq(constraints.get("nu"), 0.01, "live nu", errors)
    _eq(constraints.get("domain", {}).get("physical"), "R^3", "live physical domain", errors)
    _eq(constraints.get("domain", {}).get("support"), "r < 2 and abs(z) < 2", "live support", errors)
    restriction = str(constraints.get("forcing", {}).get("restriction", ""))
    if "No residual-dependent basis or pointwise free force" not in restriction:
        errors.append("live forcing restriction no longer forbids residual-dependent free force")

    _audit_parent(candidate_module_path, errors)

    states = contract.get("independent_delivery_states", {})
    _true(states, "canonical_eq45_velocity_export_ready", errors)
    _true(states, "kokuno_relative_swirl_algebra_callable_ready", errors)
    for key in ("canonical_eq45_visual_correspondence_verified", "canonical_eq45_pde_validated", "kokuno_global_velocity_export_ready", "kokuno_visual_correspondence_verified", "kokuno_pde_validated", "kokuno_paper_exact", "kokuno_openai_field_identified"):
        _false(states, key, errors)
    live = _load_json(project_status_path).get("states", {})
    _true(live, "velocity_export_ready", errors)
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _false(live, key, errors)
    return errors

def assert_contract(contract_path: Path, constraints_path: Path, candidate_module_path: Path, project_status_path: Path) -> None:
    errors = audit_contract(contract_path, constraints_path, candidate_module_path, project_status_path)
    if errors:
        raise ValueError("CR002 relative-swirl root-branch scope audit failed:\n- " + "\n- ".join(errors))

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=Path("configs/kokuno_relative_swirl_root_branch_scope.json"))
    parser.add_argument("--constraints", type=Path, default=Path("configs/constraints.json"))
    parser.add_argument("--candidate-module", type=Path, default=Path("src/openai_ns_reconstruction/kokuno_public_relative_swirl_compensator.py"))
    parser.add_argument("--project-status", type=Path, default=Path("project_status.json"))
    args = parser.parse_args()
    errors = audit_contract(args.contract, args.constraints, args.candidate_module, args.project_status)
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    contract = _load_json(args.contract)
    w = contract["mechanics_witness"]
    receipt = multi_root_witness(w["angular_row"], w["pressure_linear_row"], w["pressure_quadratic_matrix"], w["angular_target"], w["root_a"], w["root_b"])
    print(json.dumps({"status":"pass","mechanics_witness":receipt}, sort_keys=True))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
