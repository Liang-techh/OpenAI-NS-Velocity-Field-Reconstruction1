"""Fail-closed CR002 audit for the Kokuno relative-swirl delta-target seam.

This module changes no candidate mathematics.  It verifies that exact A1 #1159
continues to expose only a current-minus-imported angular-target correction and
that its missing absolute base target cannot be silently promoted to zero or to
a current total target.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


CONTRACT_PATH = Path("configs/kokuno_relative_swirl_delta_target_scope.json")
UPSTREAM_PATH = Path(
    "src/openai_ns_reconstruction/kokuno_current_relative_swirl_angular_correction.py"
)
CONSTRAINTS_PATH = Path("configs/constraints.json")
EXPECTED_UPSTREAM_BLOB = "37e6b0814fbad14084a8c8c4815122111252ac9b"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
EXPECTED_A1_HEAD = "6ddf35b9dc726c720033ee44161e29efc77d755e"
EXPECTED_INTEGRATION_HEAD = "6a293b3c870d70b8d8aece1b12d68685c9b8b010"


class ScopeAuditError(RuntimeError):
    """Raised when the registered CR002 scope has drifted or been promoted."""


def _git_blob_sha1(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ScopeAuditError(f"{path} must decode to a JSON object")
    return payload


def validate_contract_payload(payload: Mapping[str, Any]) -> None:
    if payload.get("schema_version") != 1:
        raise ScopeAuditError("unexpected delta-target scope schema version")
    if payload.get("scope_id") != "cr002-kokuno-relative-swirl-delta-target-scope-142":
        raise ScopeAuditError("unexpected delta-target scope id")

    integration = payload.get("active_integration_reference")
    if not isinstance(integration, Mapping):
        raise ScopeAuditError("active integration reference is missing")
    if integration.get("branch") != "codex/cr001-constraints":
        raise ScopeAuditError("active integration branch drifted")
    if integration.get("head") != EXPECTED_INTEGRATION_HEAD:
        raise ScopeAuditError("active integration reference head drifted")

    dependency = payload.get("exact_dependency")
    if not isinstance(dependency, Mapping):
        raise ScopeAuditError("exact A1 dependency is missing")
    expected_dependency = {
        "pr": 1159,
        "head": EXPECTED_A1_HEAD,
        "source_path": str(UPSTREAM_PATH),
        "source_blob": EXPECTED_UPSTREAM_BLOB,
    }
    if dict(dependency) != expected_dependency:
        raise ScopeAuditError("exact A1 #1159 dependency identity drifted")

    cr001 = payload.get("canonical_cr001")
    if not isinstance(cr001, Mapping):
        raise ScopeAuditError("canonical CR001 contract is missing")
    expected_cr001 = {
        "path": str(CONSTRAINTS_PATH),
        "blob": EXPECTED_CONSTRAINTS_BLOB,
        "nu": 0.01,
        "physical_domain": "R^3",
        "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25, 0.75],
        "forcing_mode": "restricted_two_parameter_family",
        "free_residual_force_forbidden": True,
        "reference_energy": 1.0,
        "reference_energy_abs_tolerance": 0.001,
        "validation_seed": 914027,
        "held_out_points": 4096,
        "derivative_steps": [0.02, 0.01, 0.005],
        "quadrature_orders_per_axis": [24, 48, 96],
        "momentum_max": 0.001,
        "momentum_L2": 0.001,
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
    }
    if dict(cr001) != expected_cr001:
        raise ScopeAuditError("canonical CR001 values or gates drifted")

    classes = payload.get("provenance_classes")
    if not isinstance(classes, Mapping):
        raise ScopeAuditError("provenance classes are missing")
    if set(classes) != {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }:
        raise ScopeAuditError("four-way provenance taxonomy drifted")
    if any(not isinstance(classes[key], list) or not classes[key] for key in classes):
        raise ScopeAuditError("every provenance class must remain nonempty")

    state = payload.get("materialization_state")
    if not isinstance(state, Mapping):
        raise ScopeAuditError("materialization state is missing")
    expected_state = {
        "current_minus_imported_delta_target_materialized": True,
        "imported_base_absolute_target_materialized": False,
        "current_absolute_entering_I_materialized": False,
        "current_total_relative_swirl_target_materialized": False,
        "current_cartesian_relative_swirl_composed": False,
        "terminal_exterior_global_leading_materialized": False,
        "kokuno_velocity_export_ready": False,
        "visual_correspondence_verified": False,
        "heldout_complete_ns_residual_assessed": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    if dict(state) != expected_state:
        raise ScopeAuditError("delta-only materialization boundary was promoted or drifted")

    rules = payload.get("evidence_transfer_rules")
    if not isinstance(rules, Mapping) or not rules:
        raise ScopeAuditError("evidence-transfer rules are missing")
    if set(rules.values()) != {True}:
        raise ScopeAuditError("every evidence-transfer firewall must remain enabled")

    witness = payload.get("autonomous_mechanics_witness")
    if not isinstance(witness, Mapping):
        raise ScopeAuditError("autonomous mechanics witness metadata is missing")
    if witness.get("role") != "autonomous_mechanics_only":
        raise ScopeAuditError("mechanics witness was laundered into source/candidate evidence")
    for key in ("not_source_evidence", "not_openai_numerical_evidence", "not_pde_evidence"):
        if witness.get(key) is not True:
            raise ScopeAuditError(f"mechanics witness firewall {key} was relaxed")

    forbidden = payload.get("forbidden_promotions")
    if not isinstance(forbidden, list) or len(forbidden) < 5:
        raise ScopeAuditError("forbidden-promotion list was weakened")


def _literal_assignment(module: ast.Module, name: str) -> Any:
    for node in module.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    return ast.literal_eval(node.value)
    raise ScopeAuditError(f"required upstream assignment {name} not found")


def _find_method(module: ast.Module, class_name: str, method_name: str) -> ast.FunctionDef:
    for node in module.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for child in node.body:
                if isinstance(child, ast.FunctionDef) and child.name == method_name:
                    return child
    raise ScopeAuditError(f"required upstream method {class_name}.{method_name} not found")


def audit_upstream_source(source_text: str) -> dict[str, Any]:
    module = ast.parse(source_text)
    truth = _literal_assignment(module, "_TRUTH_BOUNDARY")
    if not isinstance(truth, dict):
        raise ScopeAuditError("upstream _TRUTH_BOUNDARY must remain a literal dict")

    expected_truth = {
        "current_lineage_angular_target_correction_materialized": True,
        "imported_base_absolute_angular_target_materialized": False,
        "current_lineage_angular_entry_I_materialized": False,
        "current_lineage_absolute_r_I_materialized": False,
        "current_relative_swirl_total_target_materialized": False,
        "current_cartesian_relative_swirl_composed": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "heldout_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    for key, expected in expected_truth.items():
        if truth.get(key) is not expected:
            raise ScopeAuditError(f"upstream truth boundary promoted/drifted at {key}")

    method = _find_method(
        module,
        "KokunoCurrentRelativeSwirlAngularCorrection",
        "apply_to_materialized_base_target",
    )
    kw_names = [arg.arg for arg in method.args.kwonlyargs]
    required = {"base_target", "base_target_eta"}
    if not required.issubset(kw_names):
        raise ScopeAuditError("upstream base-target helper no longer requires both explicit base values")
    defaults = dict(zip(kw_names, method.args.kw_defaults, strict=True))
    for name in required:
        if defaults[name] is not None:
            raise ScopeAuditError(f"{name} acquired a default; missing absolute base may be hidden")

    doc = ast.get_docstring(method) or ""
    if "no default" not in doc.lower() or "does not authenticate" not in doc.lower():
        raise ScopeAuditError("upstream explicit-base provenance warning drifted")

    return {
        "delta_target_materialized": True,
        "absolute_base_materialized": False,
        "total_target_materialized": False,
        "explicit_base_target_required": True,
        "explicit_base_target_eta_required": True,
    }


def autonomous_delta_nonuniqueness_witness(
    *, delta: float = 0.0125, base_a: float = -0.03125, base_b: float = 0.046875
) -> dict[str, float | bool | str]:
    """Demonstrate, with autonomous numbers only, that a delta cannot identify a total."""

    values = (float(delta), float(base_a), float(base_b))
    if not all(value == value and abs(value) != float("inf") for value in values):
        raise ValueError("witness values must be finite")
    if delta == 0.0 or base_a == base_b:
        raise ValueError("witness requires a nonzero delta and two distinct bases")
    total_a = base_a + delta
    total_b = base_b + delta
    if total_a == total_b:
        raise AssertionError("distinct bases unexpectedly produced the same total")
    return {
        "role": "autonomous_mechanics_only",
        "delta": delta,
        "base_a": base_a,
        "base_b": base_b,
        "total_a": total_a,
        "total_b": total_b,
        "same_delta_distinct_totals": True,
        "source_evidence": False,
        "openai_numerical_evidence": False,
        "pde_evidence": False,
    }


def audit_repository(root: str | Path) -> dict[str, Any]:
    root_path = Path(root)
    contract_path = root_path / CONTRACT_PATH
    upstream_path = root_path / UPSTREAM_PATH
    constraints_path = root_path / CONSTRAINTS_PATH

    contract = _read_json(contract_path)
    validate_contract_payload(contract)

    upstream_bytes = upstream_path.read_bytes()
    upstream_blob = _git_blob_sha1(upstream_bytes)
    if upstream_blob != EXPECTED_UPSTREAM_BLOB:
        raise ScopeAuditError(
            "exact A1 #1159 source blob drifted; a changed source requires a new governance identity"
        )
    upstream_receipt = audit_upstream_source(upstream_bytes.decode("utf-8"))

    constraints_bytes = constraints_path.read_bytes()
    constraints_blob = _git_blob_sha1(constraints_bytes)
    if constraints_blob != EXPECTED_CONSTRAINTS_BLOB:
        raise ScopeAuditError("canonical CR001 constraints blob drifted")
    constraints = json.loads(constraints_bytes.decode("utf-8"))
    if constraints["forcing"]["mode"] != "restricted_two_parameter_family":
        raise ScopeAuditError("CR001 restricted forcing family drifted")
    if "No residual-dependent basis or pointwise free force" not in constraints["forcing"]["restriction"]:
        raise ScopeAuditError("CR001 no-free-residual-force restriction drifted")
    if constraints["validation"]["thresholds"]["pde_residual_max"] != 0.001:
        raise ScopeAuditError("CR001 momentum max threshold drifted")
    if constraints["validation"]["thresholds"]["pde_residual_L2"] != 0.001:
        raise ScopeAuditError("CR001 momentum L2 threshold drifted")
    if constraints["validation"]["thresholds"]["divergence_max"] != 1e-5:
        raise ScopeAuditError("CR001 divergence max threshold drifted")
    if constraints["validation"]["thresholds"]["divergence_L2"] != 1e-5:
        raise ScopeAuditError("CR001 divergence L2 threshold drifted")

    witness = autonomous_delta_nonuniqueness_witness()
    return {
        "scope_id": contract["scope_id"],
        "exact_dependency_head": EXPECTED_A1_HEAD,
        "exact_dependency_source_blob": upstream_blob,
        "canonical_constraints_blob": constraints_blob,
        "upstream_scope": upstream_receipt,
        "autonomous_mechanics_witness": witness,
        "scientific_state_promoted": False,
        "candidate_bytes_changed_by_this_audit": False,
        "cr001_threshold_changed_by_this_audit": False,
    }


__all__ = [
    "ScopeAuditError",
    "audit_repository",
    "audit_upstream_source",
    "autonomous_delta_nonuniqueness_witness",
    "validate_contract_payload",
]
