from __future__ import annotations

import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "eq45_supported_energy_failure_scope_v1"
ENERGY_AUDIT_SCHEMA = "eq45_supported_energy_audit_record_v1"

_REQUIRED_REVALIDATION = {
    "child_candidate_identity",
    "velocity_export_ready",
    "physical_support_validated",
    "cr001_energy_normalization",
    "divergence_validation",
    "pde_validated",
    "visualization_ready",
    "visual_resolution_stability",
    "visual_correspondence_verified",
}

_EXPECTED_CLASSIFICATION = {
    "velocity_delivery_requirement": "user_requirement",
    "cr001_energy_normalization": "autonomous_design",
    "support_connected_child_representation": "autonomous_design",
    "exact_openai_field_identity": "pending_unknown",
    "paper_exact_identity": "pending_unknown",
}

_EXPECTED_POLICY = {
    "reference_energy_failure_is_real_failure": True,
    "quadrature_pass_overrides_reference_failure": False,
    "per_time_energy_range_pass_overrides_reference_failure": False,
    "reference_energy_failure_blocks_velocity_export": False,
    "reference_energy_failure_blocks_truth_bounded_visualization_or_research_use": False,
    "reference_energy_failure_allows_pde_validated": False,
    "relax_reference_tolerance_to_reclassify_current_child": False,
    "implicit_export_time_amplitude_renormalization_allowed": False,
    "amplitude_rescale_preserves_candidate_identity": False,
    "amplitude_rescale_is_metadata_only": False,
    "renormalized_child_may_inherit_parent_validation_without_recheck": False,
    "renormalization_factor_must_be_serialized": True,
    "renormalization_factor_origin": "autonomous_design",
}

_EXPECTED_FALSE_TRUTH = {
    "current_supported_child_energy_normalized",
    "physical_support_validated",
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
}


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def _require_finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label} must be finite")
    return number


def _same_number(left: Any, right: Any, label: str) -> None:
    a = _require_finite_number(left, f"contract {label}")
    b = _require_finite_number(right, f"source {label}")
    if a != b:
        raise ValueError(f"energy contract drifted: {label}")


def _find_forbidden_inference(delivery_contract: Mapping[str, Any]) -> bool:
    rows = delivery_contract.get("forbidden_inferences")
    if not isinstance(rows, list):
        raise ValueError("delivery_state_contract.forbidden_inferences must be a list")
    for row in rows:
        if not isinstance(row, dict):
            continue
        source = row.get("from")
        if source == ["pde_validation_failed"] and row.get("to") == "velocity_export_not_allowed":
            return True
    return False


def audit_documents(
    contract: Mapping[str, Any],
    constraints: Mapping[str, Any],
    energy_audit: Mapping[str, Any],
    delivery_contract: Mapping[str, Any],
) -> dict[str, Any]:
    if contract.get("schema") != SCHEMA:
        raise ValueError("unexpected supported-child energy failure contract schema")
    if energy_audit.get("schema") != ENERGY_AUDIT_SCHEMA:
        raise ValueError("unexpected supported-child energy audit schema")

    classification = contract.get("classification")
    if classification != _EXPECTED_CLASSIFICATION:
        raise ValueError("supported-child energy source classification drifted")
    vocabulary = delivery_contract.get("classification_vocabulary")
    if not isinstance(vocabulary, dict):
        raise ValueError("delivery-state contract must define classification_vocabulary")
    if not set(classification.values()).issubset(vocabulary):
        raise ValueError("supported-child energy contract uses noncanonical classification")

    dependency = contract.get("dependency")
    if not isinstance(dependency, dict):
        raise ValueError("dependency must be an object")
    if dependency.get("energy_revalidation_pr") != 122:
        raise ValueError("energy revalidation dependency drifted")
    if dependency.get("energy_revalidation_head") != "a327c0a0ccee7bbe6127e0b6b66c6d25bfb3d30f":
        raise ValueError("energy revalidation exact head drifted")
    if dependency.get("energy_audit_artifact") != "artifacts/constrained/eq45_supported_energy_audit.json":
        raise ValueError("energy audit artifact path drifted")

    identity = contract.get("candidate_identity")
    if not isinstance(identity, dict):
        raise ValueError("candidate_identity must be an object")
    if identity.get("supported_child_sha256") != energy_audit.get("candidate_sha256"):
        raise ValueError("supported child identity drifted from independent energy audit")
    if identity.get("parent_sha256") != energy_audit.get("parent_sha256"):
        raise ValueError("parent identity drifted from independent energy audit")
    if identity.get("supported_child_sha256") == identity.get("parent_sha256"):
        raise ValueError("support-transformed child must retain a distinct candidate identity")

    nontriviality = constraints.get("nontriviality")
    validation = constraints.get("validation")
    if not isinstance(nontriviality, dict) or not isinstance(validation, dict):
        raise ValueError("constraints must define nontriviality and validation")
    thresholds = validation.get("thresholds")
    if not isinstance(thresholds, dict):
        raise ValueError("constraints.validation.thresholds must be an object")

    registered = contract.get("registered_energy_contract")
    if not isinstance(registered, dict):
        raise ValueError("registered_energy_contract must be an object")
    expected_registered = {
        "definition": nontriviality.get("kinetic_energy_definition"),
        "reference_time": nontriviality.get("reference_time"),
        "reference_energy": nontriviality.get("reference_energy"),
        "reference_energy_abs_tolerance": nontriviality.get("reference_energy_abs_tolerance"),
        "minimum_energy_each_validation_time": nontriviality.get("minimum_energy_each_validation_time"),
        "maximum_energy_each_validation_time": nontriviality.get("maximum_energy_each_validation_time"),
        "quadrature_orders_per_axis": validation.get("quadrature_orders_per_axis"),
        "quadrature_relative_change_threshold": thresholds.get("energy_quadrature_relative_change"),
    }
    if registered != expected_registered:
        raise ValueError("supported-child energy contract drifted from CR001 preregistration")

    for key in (
        "kinetic_energy_definition",
        "reference_time",
        "reference_energy",
        "reference_energy_abs_tolerance",
        "quadrature_orders_per_axis",
        "quadrature_relative_change_threshold",
    ):
        audit_key = {
            "kinetic_energy_definition": "kinetic_energy_definition",
            "reference_time": "reference_time",
            "reference_energy": "reference_energy",
            "reference_energy_abs_tolerance": "reference_energy_abs_tolerance",
            "quadrature_orders_per_axis": "quadrature_orders_per_axis",
            "quadrature_relative_change_threshold": "quadrature_relative_change_threshold",
        }[key]
        contract_value = {
            "kinetic_energy_definition": registered["definition"],
            "reference_time": registered["reference_time"],
            "reference_energy": registered["reference_energy"],
            "reference_energy_abs_tolerance": registered["reference_energy_abs_tolerance"],
            "quadrature_orders_per_axis": registered["quadrature_orders_per_axis"],
            "quadrature_relative_change_threshold": registered["quadrature_relative_change_threshold"],
        }[key]
        if isinstance(contract_value, (int, float)) and not isinstance(contract_value, bool):
            _same_number(contract_value, energy_audit.get(audit_key), key)
        elif contract_value != energy_audit.get(audit_key):
            raise ValueError(f"energy audit drifted from registered {key}")

    measured = contract.get("measured_supported_child_state")
    if not isinstance(measured, dict):
        raise ValueError("measured_supported_child_state must be an object")
    for key in ("reference_observed_finest_energy", "reference_energy_abs_error"):
        _same_number(measured.get(key), energy_audit.get(key), key)
    for key in ("reference_energy_gate_pass", "all_quadrature_gates_pass", "all_energy_range_gates_pass"):
        if measured.get(key) is not energy_audit.get(key):
            raise ValueError(f"measured supported-child state drifted: {key}")
    if measured.get("reference_energy_gate_pass") is not False:
        raise ValueError("the checked supported child must retain its CR001 reference-energy failure")
    if measured.get("all_quadrature_gates_pass") is not True:
        raise ValueError("this governance receipt expects the checked quadrature gate to pass")
    if measured.get("all_energy_range_gates_pass") is not True:
        raise ValueError("this governance receipt expects the checked per-time energy range to pass")
    if measured.get("cr001_energy_normalization_status") != "fail":
        raise ValueError("CR001 energy normalization status must be explicit fail")

    observed = _require_finite_number(measured.get("reference_observed_finest_energy"), "observed reference energy")
    target = _require_finite_number(registered.get("reference_energy"), "reference energy")
    tolerance = _require_finite_number(registered.get("reference_energy_abs_tolerance"), "reference tolerance")
    error = abs(observed - target)
    if error <= tolerance:
        raise ValueError("checked reference energy no longer violates the preregistered tolerance")
    if error != _require_finite_number(measured.get("reference_energy_abs_error"), "recorded reference error"):
        raise ValueError("recorded reference energy error is inconsistent")

    policy = contract.get("policy")
    if policy != _EXPECTED_POLICY:
        raise ValueError("supported-child energy failure policy drifted")
    if policy["reference_energy_failure_blocks_velocity_export"]:
        raise ValueError("CR001 energy failure cannot become a velocity-export blocker")
    if policy["implicit_export_time_amplitude_renormalization_allowed"]:
        raise ValueError("implicit export-time renormalization is forbidden")
    if policy["amplitude_rescale_preserves_candidate_identity"]:
        raise ValueError("amplitude rescaling must create a new candidate identity")
    if policy["amplitude_rescale_is_metadata_only"]:
        raise ValueError("amplitude rescaling changes [u,v,w] and cannot be metadata-only")
    if policy["renormalized_child_may_inherit_parent_validation_without_recheck"]:
        raise ValueError("renormalized child must re-establish validation evidence")
    if policy["renormalization_factor_origin"] != "autonomous_design":
        raise ValueError("renormalization factor cannot be reclassified as a public-source fact")

    revalidation = contract.get("renormalized_child_requires_revalidation")
    if not isinstance(revalidation, list) or set(revalidation) != _REQUIRED_REVALIDATION:
        raise ValueError("renormalized-child revalidation set drifted")

    truth = contract.get("truth_boundary")
    if not isinstance(truth, dict):
        raise ValueError("truth_boundary must be an object")
    if truth.get("velocity_changed_by_this_governance_increment") is not False:
        raise ValueError("governance increment must not claim a velocity change")
    for key in _EXPECTED_FALSE_TRUTH:
        if truth.get(key) is not False:
            raise ValueError(f"unsupported scientific promotion: {key}")

    audit_truth = energy_audit.get("truth_boundary")
    if not isinstance(audit_truth, dict):
        raise ValueError("energy audit must define truth_boundary")
    for key in (
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if audit_truth.get(key) is not False:
            raise ValueError(f"dependency energy audit unexpectedly promoted {key}")

    if not _find_forbidden_inference(delivery_contract):
        raise ValueError("delivery-state contract lost PDE-failure/export separation")

    return {
        "schema": SCHEMA,
        "candidate_sha256": identity["supported_child_sha256"],
        "reference_energy": observed,
        "reference_energy_abs_error": error,
        "reference_energy_gate": "fail",
        "quadrature_gate": "pass",
        "per_time_energy_range_gate": "pass",
        "velocity_export_blocked_by_energy_failure": False,
        "implicit_renormalization_allowed": False,
        "renormalization_requires_new_candidate_identity": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
    }


def audit_repository(root: str | Path = ".") -> dict[str, Any]:
    root_path = Path(root)
    contract = _load_json(root_path / "configs/eq45_supported_energy_failure_scope.json")
    constraints = _load_json(root_path / "configs/constraints.json")
    energy_audit = _load_json(root_path / "artifacts/constrained/eq45_supported_energy_audit.json")
    delivery_contract = _load_json(root_path / "configs/delivery_state_contract.json")
    return audit_documents(contract, constraints, energy_audit, delivery_contract)


__all__ = ["SCHEMA", "audit_documents", "audit_repository"]
