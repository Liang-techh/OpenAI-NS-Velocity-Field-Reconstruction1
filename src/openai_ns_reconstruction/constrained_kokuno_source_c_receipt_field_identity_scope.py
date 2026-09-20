"""Fail-closed CR002 audit of #796 receipt identity versus field identity.

Agent-1 PR #796 adds a source-compatible PA.10 physical-center configuration
with C=1000 and a complex-domain self-certificate.  Its public ``sha256``
property is deliberately a report/receipt digest: the hashed report includes
certificate execution settings such as ``certificate_intervals`` and
``chunk_size``.  Those settings do not enter ``values`` or ``derivatives``.

This audit therefore prevents that receipt digest from being promoted to a
stable candidate/velocity identity.  It changes no field, C, source formula,
forcing, threshold, or scientific state.
"""

from __future__ import annotations

import ast
import copy
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_source_c_normalized_physical_center import (
    DEFAULT_CERTIFICATE_INTERVALS,
    KokunoPA10SourceCNormalizedPhysicalCenter,
    SCHEMA as SOURCE_SCHEMA,
    SELECTED_SOURCE_C,
)

SCHEMA = "cr002-kokuno-source-c-receipt-field-identity-scope-v1"
TASK_ID = "CR002-KOKUNO-SOURCE-C-RECEIPT-FIELD-IDENTITY-071"
UPSTREAM_PR = 796
UPSTREAM_HEAD = "393abd5976d4c621ce37dad8544a1b645fdbe83b"
UPSTREAM_SOURCE_BLOB = "722640b0e6e43a2bba521776398d2f9cc8d2a052"

_ROOT = Path(__file__).resolve().parents[2]
_CONTRACT_PATH = _ROOT / "configs" / "kokuno_source_c_receipt_field_identity_scope.json"
_CONSTRAINTS_PATH = _ROOT / "configs" / "constraints.json"
_SOURCE_PATH = (
    _ROOT
    / "src"
    / "openai_ns_reconstruction"
    / "kokuno_pa10_source_c_normalized_physical_center.py"
)


class GovernanceError(ValueError):
    """Raised when the frozen representation/claim boundary drifts."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GovernanceError(message)


def load_contract() -> dict[str, Any]:
    return json.loads(_CONTRACT_PATH.read_text())


def load_constraints() -> dict[str, Any]:
    return json.loads(_CONSTRAINTS_PATH.read_text())


def _class_method(tree: ast.AST, class_name: str, method_name: str) -> ast.FunctionDef:
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == method_name:
                    return item
    raise GovernanceError(f"missing {class_name}.{method_name}")


def _attribute_chain(node: ast.AST) -> tuple[str, ...]:
    parts: list[str] = []
    current = node
    while isinstance(current, ast.Attribute):
        parts.append(current.attr)
        current = current.value
    if isinstance(current, ast.Name):
        parts.append(current.id)
    return tuple(reversed(parts))


def _single_return_call_chain(method: ast.FunctionDef) -> tuple[str, ...]:
    returns = [node for node in ast.walk(method) if isinstance(node, ast.Return)]
    _require(len(returns) == 1, f"{method.name} must have one return")
    value = returns[0].value
    _require(isinstance(value, ast.Call), f"{method.name} must directly return a call")
    return _attribute_chain(value.func)


def validate_source_structure(source_text: str) -> dict[str, Any]:
    """Bind the identity distinction to the exact #796 representation shape."""
    tree = ast.parse(source_text)

    selected_c = None
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if any(isinstance(t, ast.Name) and t.id == "SELECTED_SOURCE_C" for t in node.targets):
                selected_c = ast.literal_eval(node.value)
                break
    _require(selected_c == 1000.0, "selected source C drifted")

    values_method = _class_method(tree, "KokunoPA10SourceCNormalizedPhysicalCenter", "values")
    derivatives_method = _class_method(
        tree, "KokunoPA10SourceCNormalizedPhysicalCenter", "derivatives"
    )
    _require(
        _single_return_call_chain(values_method)
        == ("self", "physical_profiles", "values"),
        "values() no longer delegates only to physical_profiles",
    )
    _require(
        _single_return_call_chain(derivatives_method)
        == ("self", "physical_profiles", "derivatives"),
        "derivatives() no longer delegates only to physical_profiles",
    )

    configuration_method = _class_method(
        tree, "KokunoPA10SourceCNormalizedPhysicalCenter", "configuration"
    )
    configuration_source = ast.get_source_segment(source_text, configuration_method) or ""
    for token in (
        '"physical_profile_configuration"',
        '"certificate_intervals"',
        '"chunk_size"',
    ):
        _require(token in configuration_source, f"configuration missing {token}")

    report_method = _class_method(tree, "KokunoPA10SourceCNormalizedPhysicalCenter", "report")
    report_source = ast.get_source_segment(source_text, report_method) or ""
    for token in (
        "self.source_C_certificate",
        '"configuration": self.configuration()',
        '"receipt_sha256"',
        "_canonical_json(payload)",
    ):
        _require(token in report_source, f"report receipt binding missing {token}")

    sha_method = _class_method(tree, "KokunoPA10SourceCNormalizedPhysicalCenter", "sha256")
    sha_source = ast.get_source_segment(source_text, sha_method) or ""
    _require(
        'self.report()["receipt_sha256"]' in sha_source,
        "sha256 property no longer returns report receipt_sha256",
    )

    return {
        "selected_source_C": float(selected_c),
        "values_delegate_only_to_physical_profiles": True,
        "derivatives_delegate_only_to_physical_profiles": True,
        "configuration_binds_certificate_execution_knobs": True,
        "sha256_is_report_receipt_sha256": True,
    }


def validate_contract(contract: Mapping[str, Any]) -> None:
    _require(contract.get("schema_version") == 1, "contract schema_version drift")
    _require(contract.get("task_id") == TASK_ID, "task id drift")

    upstream = contract.get("upstream_binding", {})
    _require(upstream.get("pr") == UPSTREAM_PR, "upstream PR drift")
    _require(upstream.get("head") == UPSTREAM_HEAD, "upstream head drift")
    _require(upstream.get("source_blob_sha") == UPSTREAM_SOURCE_BLOB, "source blob drift")
    _require(upstream.get("schema") == SOURCE_SCHEMA, "source schema drift")
    _require(upstream.get("selected_source_C") == 1000.0, "source C drift")

    classes = contract.get("source_classification", {})
    _require(
        set(classes)
        == {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"},
        "source classes must remain exactly the canonical four",
    )
    _require(all(classes[key] for key in classes), "every source class must be nonempty")

    facts = contract.get("current_representation_facts", {})
    required_true_facts = (
        "values_delegate_only_to_physical_profiles",
        "derivatives_delegate_only_to_physical_profiles",
        "configuration_contains_physical_profile_configuration",
        "configuration_contains_certificate_intervals",
        "configuration_contains_chunk_size",
        "report_contains_source_C_certificate",
        "report_contains_full_configuration",
        "report_receipt_sha256_hashes_report_payload",
        "sha256_property_returns_report_receipt_sha256",
        "same_physical_profile_can_have_different_evidence_execution_configuration",
    )
    for key in required_true_facts:
        _require(facts.get(key) is True, f"representation fact {key} must stay true")
    _require(
        facts.get("certificate_intervals_enter_field_evaluation") is False,
        "certificate_intervals cannot be promoted to field semantics",
    )
    _require(
        facts.get("chunk_size_enters_field_evaluation") is False,
        "chunk_size cannot be promoted to field semantics",
    )

    roles = contract.get("identity_roles", {})
    _require(roles.get("current_sha256_role") == "evidence_receipt_identity", "sha role drift")
    for key in (
        "current_sha256_is_stable_field_semantic_identity",
        "current_sha256_is_candidate_sha256",
        "certificate_execution_parameters_are_field_parameters",
        "source_C_self_certificate_is_candidate_identity",
        "stable_field_semantic_digest_exposed_by_796",
    ):
        _require(roles.get(key) is False, f"identity laundering: {key}")
    for key in (
        "future_field_identity_must_bind_selected_C",
        "future_field_identity_must_bind_physical_profile_configuration",
        "future_field_identity_must_not_change_only_because_certificate_intervals_change",
        "future_field_identity_must_not_change_only_because_chunk_size_changes",
        "future_evidence_receipt_may_bind_certificate_execution_parameters",
    ):
        _require(roles.get(key) is True, f"missing future identity rule: {key}")

    boundaries = contract.get("claim_boundaries", {})
    for key, value in boundaries.items():
        _require(value is False, f"claim boundary promoted: {key}")

    truth = contract.get("truth_boundary", {})
    for key, value in truth.items():
        _require(value is False, f"truth-boundary promotion: {key}")


def validate_cr001(constraints: Mapping[str, Any], contract: Mapping[str, Any]) -> None:
    lock = contract["CR001_lock"]
    _require(constraints["nu"] == lock["nu"] == 0.01, "nu drift")
    domain = constraints["domain"]
    _require(domain["physical"] == lock["physical_domain"] == "R^3", "physical domain drift")
    _require(domain["evaluation_box"] == lock["evaluation_box"], "evaluation box drift")
    _require(domain["support"] == lock["support"], "support drift")
    _require(domain["time_interval"] == lock["time_interval"], "time interval drift")

    forcing = constraints["forcing"]
    _require(forcing["mode"] == lock["forcing_mode"], "forcing mode drift")
    _require(forcing["parameters"] == lock["forcing_parameter_bounds"], "forcing bounds drift")
    restriction = forcing["restriction"].lower()
    _require("no residual-dependent" in restriction, "free residual-defined force allowed")
    _require("pointwise free force" in restriction, "pointwise free force allowed")

    nontriviality = constraints["nontriviality"]
    _require(nontriviality["reference_energy"] == lock["reference_energy"], "energy target drift")
    _require(
        nontriviality["reference_energy_abs_tolerance"]
        == lock["reference_energy_abs_tolerance"],
        "energy tolerance drift",
    )
    _require(
        "reject collapsed candidates" in nontriviality["enforcement"],
        "amplitude-collapse rejection drift",
    )

    validation = constraints["validation"]
    _require(validation["seed"] == lock["validation_seed"], "validation seed drift")
    _require(validation["held_out_points"] == lock["held_out_points"], "held-out count drift")
    _require(validation["times"] == lock["validation_times"], "validation times drift")
    _require(validation["derivative_steps"] == lock["derivative_steps"], "derivative ladder drift")
    _require(
        validation["quadrature_orders_per_axis"] == lock["quadrature_orders_per_axis"],
        "quadrature ladder drift",
    )
    thresholds = validation["thresholds"]
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds[key] == lock[key], f"threshold drift: {key}")
    _require(
        "changing thresholds requires a new experiment version"
        in validation["failure_policy"].lower(),
        "post-hoc threshold relaxation guard drift",
    )


def behavioral_identity_witness() -> dict[str, Any]:
    """Show cheaply that evidence chunking changes config, not field semantics.

    This intentionally does not execute the expensive source-C certificate.
    """
    left = KokunoPA10SourceCNormalizedPhysicalCenter(
        certificate_intervals=DEFAULT_CERTIFICATE_INTERVALS,
        chunk_size=1024,
    )
    right = KokunoPA10SourceCNormalizedPhysicalCenter(
        certificate_intervals=DEFAULT_CERTIFICATE_INTERVALS,
        chunk_size=2048,
    )
    left_config = left.configuration()
    right_config = right.configuration()
    _require(
        left_config["physical_profile_configuration"]
        == right_config["physical_profile_configuration"],
        "physical profile changed with evidence chunking",
    )
    _require(left_config != right_config, "evidence execution configurations should differ")

    x0, x1 = left.physical_profiles.source_X_interval
    X = np.asarray([x0, 0.5 * (x0 + x1), x1], dtype=float)
    eta = np.asarray([-0.5, 0.0, 0.5], dtype=float)
    left_values = left.values(X, eta)
    right_values = right.values(X, eta)
    left_derivatives = left.derivatives(X, eta)
    right_derivatives = right.derivatives(X, eta)
    for key in left_values:
        _require(
            np.array_equal(left_values[key], right_values[key]),
            f"field value {key} changed with chunk_size",
        )
    for key in left_derivatives:
        _require(
            np.array_equal(left_derivatives[key], right_derivatives[key]),
            f"field derivative {key} changed with chunk_size",
        )

    return {
        "same_physical_profile_configuration": True,
        "different_full_evidence_configuration": True,
        "same_values_on_frozen_probe": True,
        "same_derivatives_on_frozen_probe": True,
        "left_chunk_size": 1024,
        "right_chunk_size": 2048,
        "certificate_executed": False,
    }


def audit(
    *,
    contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    source_text: str | None = None,
) -> dict[str, Any]:
    contract_payload = copy.deepcopy(dict(contract)) if contract is not None else load_contract()
    constraints_payload = (
        copy.deepcopy(dict(constraints)) if constraints is not None else load_constraints()
    )
    source_payload = source_text if source_text is not None else _SOURCE_PATH.read_text()

    validate_contract(contract_payload)
    validate_cr001(constraints_payload, contract_payload)
    source_receipt = validate_source_structure(source_payload)
    witness = behavioral_identity_witness()

    _require(SELECTED_SOURCE_C == 1000.0, "runtime selected source C drift")
    _require(SOURCE_SCHEMA == contract_payload["upstream_binding"]["schema"], "runtime schema drift")

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "upstream_pr": UPSTREAM_PR,
        "upstream_head": UPSTREAM_HEAD,
        "source_structure": source_receipt,
        "behavioral_identity_witness": witness,
        "identity_state": {
            "field_semantics_callable": True,
            "current_sha256_is_evidence_receipt_identity": True,
            "stable_field_semantic_digest_exposed_by_796": False,
            "receipt_sha256_is_candidate_sha256": False,
            "velocity_export_ready_promoted": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False
        },
        "CR001_unchanged": True
    }


def _main() -> None:
    result = audit()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
