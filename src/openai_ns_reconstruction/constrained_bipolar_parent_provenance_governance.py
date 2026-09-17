"""Fail-closed governance for the bipolar artifact's parent-provenance token.

The experimental bipolar candidate predates the canonical four-class information
source vocabulary in one narrow place: ``classification.eq45_parent`` stores the
relation token ``parent_provenance``.  That token describes artifact ancestry; it
is not a fifth source class and supplies no scientific acceptance evidence.

This module audits metadata only.  It does not change or validate numerical
velocity values, pressure, forcing, optimization, visual correspondence, or the
Navier--Stokes equations.
"""
from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

_CANONICAL_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_LEGACY_RELATION_PATH = "classification.eq45_parent"
_LEGACY_RELATION_TOKEN = "parent_provenance"
_EXPECTED_INFORMATION_CLASSIFICATION = {
    "callable_velocity_delivery": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "bipolar_phi01_profile_choice": "autonomous_design",
    "hidden_openai_numerical_profiles": "pending_unknown",
}
_EXPECTED_REQUIRED_CLASSIFICATIONS = {
    "classification.openai_hidden_profiles": "pending_unknown",
    "classification.physical_support_transform": "autonomous_design",
    "classification.taper_parameters": "autonomous_design",
    "classification.velocity_interface": "user_requirement",
    "parent_candidate.classification.eq45_backbone": "public_source_fact",
    "parent_candidate.classification.h": "autonomous_design",
    "parent_candidate.classification.openai_hidden_profiles": "pending_unknown",
    "parent_candidate.classification.profile_coefficients": "autonomous_design",
    "parent_candidate.classification.profile_family": "autonomous_design",
    "parent_candidate.classification.time_interval": "autonomous_design",
    "parent_candidate.classification.velocity_interface": "user_requirement",
    "parent_candidate.profile_basis.classification.basis_shape": "autonomous_design",
    "parent_candidate.profile_basis.classification.coefficient_values": "autonomous_design",
    "physical_taper.classification": "autonomous_design",
}
_REQUIRED_FALSE_TRUTH = (
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _get_path(value: Mapping[str, Any], dotted: str) -> Any:
    current: Any = value
    for part in dotted.split("."):
        if not isinstance(current, Mapping) or part not in current:
            raise ValueError(f"missing governed path: {dotted}")
        current = current[part]
    return current


def _classification_strings(value: Any, prefix: tuple[str, ...] = ()) -> list[tuple[str, str]]:
    found: list[tuple[str, str]] = []
    if not isinstance(value, Mapping):
        return found
    for key, child in value.items():
        path = prefix + (str(key),)
        if key == "classification":
            if isinstance(child, str):
                found.append((".".join(path), child))
            elif isinstance(child, Mapping):
                for subkey, subvalue in child.items():
                    if isinstance(subvalue, str):
                        found.append((".".join(path + (str(subkey),)), subvalue))
        if isinstance(child, Mapping):
            found.extend(_classification_strings(child, path))
    return found


def audit_bipolar_parent_provenance_scope(
    scope: Mapping[str, Any],
    candidate: Mapping[str, Any],
    delivery_contract: Mapping[str, Any],
    constraints: Mapping[str, Any],
    project_status: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit the current bipolar parent relation and surrounding claim boundary."""
    _require(scope.get("schema") == "bipolar_parent_provenance_scope_v1", "scope schema drift")
    _require(scope.get("artifact_mutation") is False, "governance scope must not mutate artifact")

    vocabulary = delivery_contract.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "missing canonical classification vocabulary")
    _require(set(vocabulary) == _CANONICAL_CLASSES, "canonical source vocabulary drift")
    _require(
        set(scope.get("canonical_source_vocabulary", [])) == _CANONICAL_CLASSES,
        "scope source vocabulary drift",
    )
    _require(
        scope.get("information_classification") == _EXPECTED_INFORMATION_CLASSIFICATION,
        "information classification drift",
    )

    relation = scope.get("legacy_parent_relation")
    _require(isinstance(relation, Mapping), "missing parent relation governance")
    _require(relation.get("path") == _LEGACY_RELATION_PATH, "parent relation path drift")
    _require(relation.get("token") == _LEGACY_RELATION_TOKEN, "parent relation token drift")
    _require(relation.get("semantic_kind") == "artifact_parent_relation", "parent relation semantic drift")
    _require(relation.get("is_source_classification") is False, "parent relation cannot be a source class")
    for key in (
        "may_imply_public_source",
        "may_imply_paper_exact",
        "may_imply_openai_field_identity",
        "may_imply_pde_validation",
    ):
        _require(relation.get(key) is False, f"unsupported parent-relation inference: {key}")
    _require(_LEGACY_RELATION_TOKEN not in _CANONICAL_CLASSES, "parent relation became a fifth source class")
    _require(
        _get_path(candidate, _LEGACY_RELATION_PATH) == _LEGACY_RELATION_TOKEN,
        "bipolar artifact parent relation changed",
    )

    required = scope.get("required_artifact_classifications")
    _require(required == _EXPECTED_REQUIRED_CLASSIFICATIONS, "required classification map drift")
    for path, expected in _EXPECTED_REQUIRED_CLASSIFICATIONS.items():
        _require(_get_path(candidate, path) == expected, f"artifact source classification drift: {path}")

    observed_classes: set[str] = set()
    for path, classification in _classification_strings(candidate):
        if path == _LEGACY_RELATION_PATH:
            _require(classification == _LEGACY_RELATION_TOKEN, "unexpected parent relation token")
            continue
        _require(
            classification in _CANONICAL_CLASSES,
            f"noncanonical artifact information class at {path}: {classification!r}",
        )
        observed_classes.add(classification)
    _require(observed_classes == _CANONICAL_CLASSES, "bipolar artifact does not exercise all four source classes")

    binding = scope.get("candidate_binding")
    _require(isinstance(binding, Mapping), "missing candidate binding")
    _require(candidate.get("schema") == binding.get("schema"), "candidate schema drift")
    _require(candidate.get("parent_sha256") == binding.get("parent_sha256"), "candidate parent identity drift")
    phi = _get_path(candidate, "parent_candidate.profile_basis.phi_coefficients")
    _require(phi == binding.get("bipolar_phi_coefficients"), "bipolar Phi(0,1) profile identity drift")

    experimental = project_status.get(binding.get("project_status_key"))
    _require(isinstance(experimental, Mapping), "missing project-status bipolar candidate entry")
    inputs = scope.get("inputs")
    _require(isinstance(inputs, Mapping), "missing scope inputs")
    _require(experimental.get("path") == inputs.get("candidate"), "project-status bipolar path drift")
    _require(experimental.get("sha256") == binding.get("candidate_sha256"), "project-status bipolar SHA drift")
    _require(experimental.get("source_correspondence_verified") is False, "unsupported source-correspondence promotion")

    expected_cr001 = scope.get("cr001_invariants")
    _require(isinstance(expected_cr001, Mapping), "missing CR001 invariant snapshot")
    _require(constraints.get("nu") == expected_cr001.get("nu") == 0.01, "nu drift")
    domain = constraints.get("domain")
    _require(isinstance(domain, Mapping), "missing domain contract")
    _require(domain.get("physical") == expected_cr001.get("physical_domain") == "R^3", "physical domain drift")
    _require(domain.get("evaluation_box") == expected_cr001.get("evaluation_box"), "evaluation box drift")
    _require(domain.get("support") == expected_cr001.get("support"), "support drift")
    _require(domain.get("time_interval") == expected_cr001.get("time_interval"), "time interval drift")

    forcing = constraints.get("forcing")
    _require(isinstance(forcing, Mapping), "missing forcing contract")
    _require(forcing.get("mode") == expected_cr001.get("forcing_mode"), "forcing family drift")
    _require(forcing.get("parameters") == expected_cr001.get("forcing_bounds"), "forcing bounds drift")
    restriction = forcing.get("restriction")
    _require(
        isinstance(restriction, str)
        and "Only a,c may be fitted" in restriction
        and "No residual-dependent basis or pointwise free force" in restriction,
        "free-force prohibition drift",
    )

    nontriviality = constraints.get("nontriviality")
    _require(isinstance(nontriviality, Mapping), "missing nontriviality contract")
    _require(nontriviality.get("reference_energy") == expected_cr001.get("reference_energy") == 1.0, "reference energy drift")
    _require(
        nontriviality.get("reference_energy_abs_tolerance")
        == expected_cr001.get("reference_energy_abs_tolerance")
        == 0.001,
        "energy tolerance drift",
    )

    validation = constraints.get("validation")
    _require(isinstance(validation, Mapping), "missing validation contract")
    _require(validation.get("held_out_points") == expected_cr001.get("held_out_points") == 4096, "held-out count drift")
    _require(validation.get("derivative_steps") == expected_cr001.get("derivative_steps"), "derivative ladder drift")
    thresholds = validation.get("thresholds")
    _require(isinstance(thresholds, Mapping), "missing validation thresholds")
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds.get(key) == expected_cr001.get(key), f"registered threshold drift: {key}")

    boundary = scope.get("delivery_truth_boundary")
    _require(isinstance(boundary, Mapping), "missing delivery truth boundary")
    _require(boundary.get("velocity_export_ready") is True, "velocity export must remain independently allowed")
    for key in _REQUIRED_FALSE_TRUTH:
        _require(boundary.get(key) is False, f"scope truth-state promotion: {key}")

    candidate_truth = candidate.get("truth_boundary")
    _require(isinstance(candidate_truth, Mapping), "missing candidate truth boundary")
    _require(candidate_truth.get("velocity_export_ready") is True, "candidate export state drift")
    for key in _REQUIRED_FALSE_TRUTH:
        _require(candidate_truth.get(key) is False, f"candidate truth-state promotion: {key}")

    status_states = project_status.get("states")
    _require(isinstance(status_states, Mapping), "missing project delivery states")
    _require(status_states.get("velocity_export_ready") is True, "project export state drift")
    for key in _REQUIRED_FALSE_TRUTH:
        _require(status_states.get(key) is False, f"project truth-state promotion: {key}")

    return {
        "scope_pass": True,
        "legacy_parent_relation": _LEGACY_RELATION_TOKEN,
        "legacy_parent_relation_is_source_classification": False,
        "canonical_source_classes": sorted(observed_classes),
        "candidate_sha256": experimental.get("sha256"),
        "velocity_export_ready": True,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def audit_bipolar_parent_provenance_scope_file(
    scope_path: str | Path, *, repo_root: str | Path
) -> dict[str, Any]:
    """Load the governed repository inputs and audit their current semantics."""
    root = Path(repo_root)
    scope = _read_json(Path(scope_path))
    inputs = scope.get("inputs")
    _require(isinstance(inputs, Mapping), "missing scope inputs")
    candidate = _read_json(root / str(inputs.get("candidate")))
    delivery = _read_json(root / str(inputs.get("delivery_state_contract")))
    constraints = _read_json(root / str(inputs.get("constraints")))
    status = _read_json(root / str(inputs.get("project_status")))
    return audit_bipolar_parent_provenance_scope(scope, candidate, delivery, constraints, status)
