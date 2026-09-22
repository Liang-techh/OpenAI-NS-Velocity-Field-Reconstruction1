"""Fail-closed CR002 audit for the exact #1133 xi=13 delivery-domain boundary.

This module does not modify or numerically retune a velocity field. It guards
the distinction between a stage-scoped Cartesian callable through the public
source-coordinate pulse endpoint and a globally deliverable project-domain
Kokuno velocity.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

EXPECTED_SCHEMA = "cr002-kokuno-xi13-global-delivery-domain-scope-v1"
EXPECTED_TASK = "CR002-KOKUNO-XI13-GLOBAL-DELIVERY-DOMAIN-139"
EXPECTED_PARENT_HEAD = "3b4af71c4ab547bc11f9f3e320afa4b62c25b930"
EXPECTED_PARENT_BLOB = "4b8fad0781d74c934b6cc6fe7805a11d16807a92"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _git_blob_sha(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode()
    return hashlib.sha1(header + data).hexdigest()


def _eq(actual: Any, expected: Any, label: str, errors: list[str]) -> None:
    if actual != expected:
        errors.append(f"{label}: expected {expected!r}, got {actual!r}")


def _false(mapping: Mapping[str, Any], key: str, errors: list[str]) -> None:
    if mapping.get(key) is not False:
        errors.append(f"{key} must remain false")


def _true(mapping: Mapping[str, Any], key: str, errors: list[str]) -> None:
    if mapping.get(key) is not True:
        errors.append(f"{key} must remain true")


def prefix_nonuniqueness_witness(endpoint: float, outer_probe: float) -> dict[str, Any]:
    """Autonomous mechanics witness: a finite prefix cannot identify its extension."""
    endpoint = float(endpoint)
    outer_probe = float(outer_probe)
    if not (math.isfinite(endpoint) and math.isfinite(outer_probe)):
        raise ValueError("mechanics witness inputs must be finite")
    if not outer_probe > endpoint:
        raise ValueError("outer_probe must lie strictly beyond endpoint")
    prefix_probe = endpoint - 0.5
    prefix_a = math.sin(prefix_probe)
    prefix_b = math.sin(prefix_probe)
    endpoint_a = math.sin(endpoint)
    endpoint_b = math.sin(endpoint)
    outer_a = math.sin(endpoint)
    outer_b = math.sin(endpoint) + (outer_probe - endpoint) ** 2
    return {
        "same_prefix_value": prefix_a == prefix_b,
        "same_endpoint_value": endpoint_a == endpoint_b,
        "different_outer_value": outer_a != outer_b,
        "outer_delta": outer_b - outer_a,
    }


def _audit_source_shape(candidate_module_path: Path, errors: list[str]) -> None:
    if _git_blob_sha(candidate_module_path) != EXPECTED_PARENT_BLOB:
        errors.append("exact #1133 candidate module blob changed")
        return
    text = candidate_module_path.read_text()
    required_fragments = (
        '"current_cartesian_end_compensation_composed": True',
        '"current_pulse_endpoint_xi13_materialized": True',
        '"current_endpoint_MJ_closure_numerically_assessed": True',
        '"source_terminal_tail_schedule_bound_into_current_velocity": False',
        '"source_exterior_heat_replacement_materialized": False',
        '"outer_global_leading_velocity_materialized": False',
        '"unified_global_cartesian_velocity_export_ready": False',
        '"matched_global_pressure_materialized": False',
        '"restricted_forcing_materialized": False',
        '"heldout_ns_residual_assessed": False',
        '"pde_validated": False',
        'beyond the public xi=13 pulse endpoint',
    )
    for fragment in required_fragments:
        if fragment not in text:
            errors.append(f"#1133 source/truth fragment missing: {fragment}")
    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        errors.append(f"#1133 candidate module does not parse: {exc}")
        return
    if not any(
        isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == "velocity"
        for node in ast.walk(tree)
    ):
        errors.append("#1133 candidate module must expose a velocity method")


def audit_contract(
    contract_path: Path,
    constraints_path: Path,
    candidate_module_path: Path,
    project_status_path: Path,
) -> list[str]:
    errors: list[str] = []
    contract = _load_json(contract_path)
    _eq(contract.get("schema"), EXPECTED_SCHEMA, "schema", errors)
    _eq(contract.get("task_id"), EXPECTED_TASK, "task_id", errors)

    parent = contract.get("exact_parent", {})
    if not isinstance(parent, dict):
        errors.append("exact_parent must be an object")
        parent = {}
    _eq(parent.get("pr"), 1133, "exact_parent.pr", errors)
    _eq(parent.get("head_sha"), EXPECTED_PARENT_HEAD, "exact_parent.head_sha", errors)
    _eq(parent.get("candidate_module_blob_sha"), EXPECTED_PARENT_BLOB, "exact parent blob", errors)
    _eq(
        parent.get("candidate_module"),
        "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_pulse_end_compensated.py",
        "exact_parent.candidate_module",
        errors,
    )

    provenance = contract.get("four_way_provenance", {})
    expected_classes = {"user_requirements", "public_source_facts", "autonomous_design", "pending_unknown"}
    if not isinstance(provenance, dict):
        errors.append("four_way_provenance must be an object")
        provenance = {}
    _eq(set(provenance), expected_classes, "four_way_provenance classes", errors)
    for key in expected_classes:
        value = provenance.get(key)
        if not isinstance(value, list) or not value:
            errors.append(f"four_way_provenance.{key} must be a nonempty list")

    scope = contract.get("stage_delivery_scope", {})
    if not isinstance(scope, dict):
        errors.append("stage_delivery_scope must be an object")
        scope = {}
    _eq(scope.get("public_source_pulse_endpoint_xi"), 13.0, "public pulse endpoint", errors)
    for key in (
        "current_cartesian_pulse_endpoint_xi13_materialized",
        "current_endpoint_MJ_closure_numerically_assessed",
        "current_stage_velocity_callable_xyzt",
        "current_stage_save_load_identity_materialized",
        "hard_fail_closed_beyond_xi13",
    ):
        _true(scope, key, errors)
    for key in (
        "xi13_is_canonical_physical_support_endpoint",
        "xi13_is_canonical_evaluation_box_boundary",
        "terminal_exterior_tail_materialized",
        "outer_global_leading_velocity_materialized",
        "global_project_domain_coverage_verified",
        "unified_global_kokuno_velocity_export_ready",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_complete_ns_residual_assessed",
    ):
        _false(scope, key, errors)

    firewall = contract.get("evidence_transfer_firewall", {})
    if not isinstance(firewall, dict):
        errors.append("evidence_transfer_firewall must be an object")
        firewall = {}
    forbidden = firewall.get("forbidden_implications")
    if not isinstance(forbidden, list) or len(forbidden) < 8:
        errors.append("evidence_transfer_firewall.forbidden_implications is incomplete")
    else:
        joined = "\n".join(str(v) for v in forbidden)
        for fragment in (
            "terminal/exterior/global leading completion",
            "canonical physical support endpoint",
            "canonical evaluation-box boundary",
            "held-out complete NS residual or PDE validation",
            "global project-domain coverage",
            "unqualified Kokuno velocity_export_ready",
            "visual correspondence",
            "paper exactness or exact OpenAI-field identity",
        ):
            if fragment not in joined:
                errors.append(f"missing forbidden implication fragment: {fragment}")

    promotion = contract.get("future_extension_promotion_rule", {})
    if not isinstance(promotion, dict):
        errors.append("future_extension_promotion_rule must be an object")
        promotion = {}
    _false(promotion, "reuse_xi13_prefix_receipts_as_outer_domain_evidence", errors)
    for key in (
        "new_identity_required_for_terminal_or_exterior_extension",
        "fresh_save_load_audit_required",
        "fresh_finite_cartesian_domain_audit_required",
        "fresh_downstream_numerical_audit_required",
    ):
        _true(promotion, key, errors)

    cr001 = contract.get("cr001_invariants", {})
    if not isinstance(cr001, dict):
        errors.append("cr001_invariants must be an object")
        cr001 = {}
    _eq(cr001.get("constraints_blob_sha"), EXPECTED_CONSTRAINTS_BLOB, "constraints blob", errors)
    _eq(cr001.get("nu"), 0.01, "nu", errors)
    _eq(cr001.get("physical_domain"), "R^3", "physical domain", errors)
    _eq(cr001.get("evaluation_box"), [[-2, 2], [-2, 2], [-2, 2]], "evaluation box", errors)
    _eq(cr001.get("support"), "r < 2 and abs(z) < 2", "support", errors)
    _eq(cr001.get("time_interval"), [0.25, 0.75], "time interval", errors)
    _eq(cr001.get("forcing_mode"), "restricted_two_parameter_family", "forcing mode", errors)
    _eq(cr001.get("validation_seed"), 914027, "validation seed", errors)
    _eq(cr001.get("held_out_points"), 4096, "held-out points", errors)
    _eq(cr001.get("derivative_steps"), [0.02, 0.01, 0.005], "derivative steps", errors)
    _eq(cr001.get("quadrature_orders_per_axis"), [24, 48, 96], "quadrature orders", errors)
    _false(cr001, "forcing_family_change_this_increment", errors)
    _false(cr001, "threshold_relaxation_this_increment", errors)
    _eq(cr001.get("thresholds"), {
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
        "pde_residual_max": 1e-3,
        "pde_residual_L2": 1e-3,
    }, "CR001 thresholds", errors)

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
    nontriviality = constraints.get("nontriviality", {})
    _eq(nontriviality.get("reference_energy"), 1.0, "live reference energy", errors)
    _eq(nontriviality.get("reference_energy_abs_tolerance"), 0.001, "live reference energy tolerance", errors)
    validation = constraints.get("validation", {})
    _eq(validation.get("seed"), 914027, "live validation seed", errors)
    _eq(validation.get("held_out_points"), 4096, "live held-out points", errors)
    _eq(validation.get("derivative_steps"), [0.02, 0.01, 0.005], "live derivative steps", errors)
    _eq(validation.get("quadrature_orders_per_axis"), [24, 48, 96], "live quadrature orders", errors)
    thresholds = validation.get("thresholds", {})
    _eq(thresholds.get("pde_residual_max"), 0.001, "live momentum max gate", errors)
    _eq(thresholds.get("pde_residual_L2"), 0.001, "live momentum L2 gate", errors)
    _eq(thresholds.get("divergence_max"), 1e-5, "live divergence max gate", errors)
    _eq(thresholds.get("divergence_L2"), 1e-5, "live divergence L2 gate", errors)

    _audit_source_shape(candidate_module_path, errors)

    states = contract.get("independent_delivery_states", {})
    if not isinstance(states, dict):
        errors.append("independent_delivery_states must be an object")
        states = {}
    _true(states, "canonical_eq45_velocity_export_ready", errors)
    _true(states, "kokuno_stage_callable_ready", errors)
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

    project_status = _load_json(project_status_path)
    live_states = project_status.get("states", {})
    _true(live_states, "velocity_export_ready", errors)
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _false(live_states, key, errors)

    witness = contract.get("mechanics_witness", {})
    if not isinstance(witness, dict):
        errors.append("mechanics_witness must be an object")
        witness = {}
    _eq(witness.get("classification"), "autonomous_mechanics_only", "witness classification", errors)
    try:
        receipt = prefix_nonuniqueness_witness(witness.get("endpoint", float("nan")), witness.get("outer_probe", float("nan")))
    except (TypeError, ValueError) as exc:
        errors.append(f"invalid mechanics witness: {exc}")
    else:
        if receipt["same_prefix_value"] is not True:
            errors.append("mechanics witness must preserve prefix data")
        if receipt["same_endpoint_value"] is not True:
            errors.append("mechanics witness must preserve endpoint data")
        if receipt["different_outer_value"] is not True or not receipt["outer_delta"] > 0.0:
            errors.append("mechanics witness must exhibit non-unique outer continuation")
    return errors


def assert_contract(contract_path: Path, constraints_path: Path, candidate_module_path: Path, project_status_path: Path) -> None:
    errors = audit_contract(contract_path, constraints_path, candidate_module_path, project_status_path)
    if errors:
        raise ValueError("CR002 xi13 delivery-domain audit failed:\n- " + "\n- ".join(errors))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=Path("configs/kokuno_xi13_global_delivery_domain_scope.json"))
    parser.add_argument("--constraints", type=Path, default=Path("configs/constraints.json"))
    parser.add_argument("--candidate-module", type=Path, default=Path("src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_pulse_end_compensated.py"))
    parser.add_argument("--project-status", type=Path, default=Path("project_status.json"))
    args = parser.parse_args()
    errors = audit_contract(args.contract, args.constraints, args.candidate_module, args.project_status)
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    print("PASS: CR002 xi=13 pulse stage remains distinct from global delivery")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
