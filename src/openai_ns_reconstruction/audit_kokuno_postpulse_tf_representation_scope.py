"""Fail-closed CR002 audit for exact #1140 post-pulse T_f provenance/identity scope.

No velocity, pressure, forcing, or scientific threshold is changed here.  The
audit separates the public qualitative requirement (a fixed sufficiently large
T_f), the repository-autonomous numeric realization T_f=100, and any
source-exact numerical claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

EXPECTED_SCHEMA = "cr002-kokuno-postpulse-tf-representation-scope-v1"
EXPECTED_TASK = "CR002-KOKUNO-POSTPULSE-TF-SCOPE-140"
EXPECTED_PARENT_HEAD = "c5442c11165d7f0976889b826bf58dce38a3d3dd"
EXPECTED_PARENT_BLOB = "32cbfc00406aa8d47904e2492ad8916730541d28"
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


def _smooth_step(s: float) -> float:
    s = float(s)
    if not math.isfinite(s):
        raise ValueError("s must be finite")
    if s <= 0.0:
        return 0.0
    if s >= 1.0:
        return 1.0
    logit = -1.0 / (s * s) + 1.0 / ((1.0 - s) * (1.0 - s))
    if logit >= 0.0:
        return 1.0 / (1.0 + math.exp(-logit))
    value = math.exp(logit)
    return value / (1.0 + value)


def tf_sensitivity_witness(T_f_a: float, T_f_b: float, y: float, eta: float) -> dict[str, Any]:
    """Autonomous mechanics witness; not Kokuno/OpenAI or PDE evidence."""
    T_f_a, T_f_b, y, eta = map(float, (T_f_a, T_f_b, y, eta))
    if not all(math.isfinite(v) for v in (T_f_a, T_f_b, y, eta)):
        raise ValueError("witness inputs must be finite")
    if T_f_a <= 0.0 or T_f_b <= 0.0 or T_f_a == T_f_b:
        raise ValueError("T_f values must be distinct positive numbers")
    if not (0.0 < y < min(T_f_a, T_f_b)):
        raise ValueError("y must be inside both flattening intervals")
    sigma_a = _smooth_step(y / T_f_a)
    sigma_b = _smooth_step(y / T_f_b)
    coefficient = math.log1p(eta * eta) - math.log(2.0)
    if abs(coefficient) <= 1e-15:
        raise ValueError("eta probe must make J0-log(2) nonzero")
    correction_a = sigma_a * coefficient
    correction_b = sigma_b * coefficient
    return {
        "sigma_a": sigma_a,
        "sigma_b": sigma_b,
        "correction_a": correction_a,
        "correction_b": correction_b,
        "different_flattening_weight": sigma_a != sigma_b,
        "different_logE_correction": correction_a != correction_b,
        "absolute_correction_delta": abs(correction_a - correction_b),
    }


def _audit_parent(candidate_module_path: Path, errors: list[str]) -> None:
    if _git_blob_sha(candidate_module_path) != EXPECTED_PARENT_BLOB:
        errors.append("exact #1140 candidate module blob changed")
        return
    text = candidate_module_path.read_text(encoding="utf-8")
    fragments = (
        "T_F_AUTONOMOUS = 100.0",
        "fixed sufficiently large T_f",
        "it does not supply a unique numerical value",
        '"repository_autonomous_T_f_materialized": True',
        '"source_exact_T_f_recovered": False',
        '"current_cartesian_postpulse_eta_flattening_composed": True',
        '"source_terminal_hold_after_eta_flattening_materialized": False',
        '"outer_global_leading_velocity_materialized": False',
        '"unified_global_cartesian_velocity_export_ready": False',
        '"restricted_forcing_materialized": False',
        '"heldout_ns_residual_assessed": False',
        '"pde_validated": False',
        '"T_f_role": "repository_autonomous_not_source_exact"',
    )
    for fragment in fragments:
        if fragment not in text:
            errors.append(f"#1140 source/truth fragment missing: {fragment}")


def audit_contract(contract_path: Path, constraints_path: Path, candidate_module_path: Path, project_status_path: Path) -> list[str]:
    errors: list[str] = []
    contract = _load_json(contract_path)
    _eq(contract.get("schema"), EXPECTED_SCHEMA, "schema", errors)
    _eq(contract.get("task_id"), EXPECTED_TASK, "task_id", errors)

    parent = contract.get("exact_parent", {})
    if not isinstance(parent, dict):
        errors.append("exact_parent must be an object")
        parent = {}
    _eq(parent.get("pr"), 1140, "exact_parent.pr", errors)
    _eq(parent.get("head_sha"), EXPECTED_PARENT_HEAD, "exact_parent.head_sha", errors)
    _eq(parent.get("candidate_module_blob_sha"), EXPECTED_PARENT_BLOB, "exact parent blob", errors)

    provenance = contract.get("four_way_provenance", {})
    classes = {"user_requirements", "public_source_facts", "autonomous_design", "pending_unknown"}
    if not isinstance(provenance, dict):
        errors.append("four_way_provenance must be an object")
        provenance = {}
    _eq(set(provenance), classes, "four_way provenance classes", errors)
    for key in classes:
        if not isinstance(provenance.get(key), list) or not provenance.get(key):
            errors.append(f"four_way_provenance.{key} must be a nonempty list")

    scope = contract.get("tf_representation_scope", {})
    if not isinstance(scope, dict):
        errors.append("tf_representation_scope must be an object")
        scope = {}
    _true(scope, "source_requires_fixed_sufficiently_large_T_f", errors)
    _false(scope, "source_supplies_unique_numeric_T_f", errors)
    _eq(scope.get("current_autonomous_T_f"), 100.0, "current autonomous T_f", errors)
    for key in (
        "current_autonomous_T_f_materialized",
        "postpulse_stage_velocity_callable_xyzt",
        "postpulse_stage_save_load_identity_materialized",
        "postpulse_stage_hard_fail_closed_after_T_f",
        "deterministic_ell_window_receipt_materialized",
    ):
        _true(scope, key, errors)
    for key in (
        "source_exact_T_f_recovered",
        "terminal_global_leading_materialized",
        "unified_global_kokuno_velocity_export_ready",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_complete_ns_residual_assessed",
    ):
        _false(scope, key, errors)

    forbidden = contract.get("evidence_transfer_firewall", {}).get("forbidden_implications")
    if not isinstance(forbidden, list) or len(forbidden) < 8:
        errors.append("evidence_transfer_firewall.forbidden_implications is incomplete")
    else:
        joined = "\n".join(map(str, forbidden))
        for fragment in (
            "T_f=100 is a public/source-exact value",
            "unique source T_f recovered",
            "terminal/exterior/global leading completion",
            "unqualified Kokuno velocity_export_ready",
            "different numerical T_f identity",
            "held-out complete NS residual or PDE validation",
            "visual correspondence",
            "paper exactness or exact OpenAI-field identity",
        ):
            if fragment not in joined:
                errors.append(f"missing forbidden implication fragment: {fragment}")

    promotion = contract.get("numeric_tf_change_promotion_rule", {})
    if not isinstance(promotion, dict):
        errors.append("numeric_tf_change_promotion_rule must be an object")
        promotion = {}
    for key in (
        "numeric_T_f_change_requires_new_semantic_candidate_identity",
        "fresh_save_load_audit_required",
        "fresh_finite_cartesian_domain_audit_required",
        "fresh_downstream_numerical_audit_required",
    ):
        _true(promotion, key, errors)
    _false(promotion, "reuse_old_ell_or_divergence_receipts_for_new_T_f", errors)

    cr001 = contract.get("cr001_invariants", {})
    if not isinstance(cr001, dict):
        errors.append("cr001_invariants must be an object")
        cr001 = {}
    expected_pairs = {
        "constraints_blob_sha": EXPECTED_CONSTRAINTS_BLOB,
        "nu": 0.01,
        "physical_domain": "R^3",
        "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25, 0.75],
        "forcing_mode": "restricted_two_parameter_family",
        "reference_energy": 1.0,
        "reference_energy_abs_tolerance": 0.001,
        "validation_seed": 914027,
        "held_out_points": 4096,
        "derivative_steps": [0.02, 0.01, 0.005],
        "quadrature_orders_per_axis": [24, 48, 96],
        "thresholds": {"divergence_max": 1e-5, "divergence_L2": 1e-5, "pde_residual_max": 1e-3, "pde_residual_L2": 1e-3},
    }
    for key, expected in expected_pairs.items():
        _eq(cr001.get(key), expected, f"cr001.{key}", errors)
    _false(cr001, "forcing_family_change_this_increment", errors)
    _false(cr001, "threshold_relaxation_this_increment", errors)

    if _git_blob_sha(constraints_path) != EXPECTED_CONSTRAINTS_BLOB:
        errors.append("canonical configs/constraints.json blob changed")
    constraints = _load_json(constraints_path)
    _eq(constraints.get("nu"), 0.01, "live constraints.nu", errors)
    domain = constraints.get("domain", {})
    _eq(domain.get("physical"), "R^3", "live physical domain", errors)
    _eq(domain.get("evaluation_box"), [[-2, 2], [-2, 2], [-2, 2]], "live evaluation box", errors)
    _eq(domain.get("support"), "r < 2 and abs(z) < 2", "live support", errors)
    _eq(domain.get("time_interval"), [0.25, 0.75], "live time interval", errors)
    forcing = constraints.get("forcing", {})
    _eq(forcing.get("mode"), "restricted_two_parameter_family", "live forcing mode", errors)
    if "No residual-dependent basis or pointwise free force" not in str(forcing.get("restriction", "")):
        errors.append("live forcing restriction no longer forbids residual-dependent free force")
    validation = constraints.get("validation", {})
    _eq(validation.get("seed"), 914027, "live validation seed", errors)
    _eq(validation.get("held_out_points"), 4096, "live held-out points", errors)
    _eq(validation.get("derivative_steps"), [0.02, 0.01, 0.005], "live derivative steps", errors)
    _eq(validation.get("quadrature_orders_per_axis"), [24, 48, 96], "live quadrature", errors)
    thresholds = validation.get("thresholds", {})
    for key, expected in (("divergence_max", 1e-5), ("divergence_L2", 1e-5), ("pde_residual_max", 1e-3), ("pde_residual_L2", 1e-3)):
        _eq(thresholds.get(key), expected, f"live {key}", errors)

    _audit_parent(candidate_module_path, errors)

    states = contract.get("independent_delivery_states", {})
    _true(states, "canonical_eq45_velocity_export_ready", errors)
    _true(states, "kokuno_postpulse_stage_callable_ready", errors)
    for key in (
        "canonical_eq45_visual_correspondence_verified",
        "canonical_eq45_pde_validated",
        "kokuno_global_velocity_export_ready",
        "kokuno_visual_correspondence_verified",
        "kokuno_pde_validated",
        "kokuno_paper_exact",
        "kokuno_openai_field_identified",
    ):
        _false(states, key, errors)

    live = _load_json(project_status_path).get("states", {})
    _true(live, "velocity_export_ready", errors)
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _false(live, key, errors)

    witness = contract.get("mechanics_witness", {})
    _eq(witness.get("classification"), "autonomous_mechanics_only", "witness classification", errors)
    try:
        receipt = tf_sensitivity_witness(witness.get("T_f_a"), witness.get("T_f_b"), witness.get("y_probe"), witness.get("eta_probe"))
    except (TypeError, ValueError) as exc:
        errors.append(f"invalid mechanics witness: {exc}")
    else:
        if receipt["different_flattening_weight"] is not True or receipt["different_logE_correction"] is not True or not receipt["absolute_correction_delta"] > 0.0:
            errors.append("mechanics witness must demonstrate numerical T_f sensitivity")
    return errors


def assert_contract(contract_path: Path, constraints_path: Path, candidate_module_path: Path, project_status_path: Path) -> None:
    errors = audit_contract(contract_path, constraints_path, candidate_module_path, project_status_path)
    if errors:
        raise ValueError("CR002 post-pulse T_f scope audit failed:\n- " + "\n- ".join(errors))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=Path("configs/kokuno_postpulse_tf_representation_scope.json"))
    parser.add_argument("--constraints", type=Path, default=Path("configs/constraints.json"))
    parser.add_argument("--candidate-module", type=Path, default=Path("src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_postpulse_eta_flattening.py"))
    parser.add_argument("--project-status", type=Path, default=Path("project_status.json"))
    args = parser.parse_args()
    errors = audit_contract(args.contract, args.constraints, args.candidate_module, args.project_status)
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    print(json.dumps({"status": "pass", "mechanics_witness": tf_sensitivity_witness(100.0, 120.0, 50.0, 0.5)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
