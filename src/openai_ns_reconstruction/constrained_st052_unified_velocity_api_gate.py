"""Fail-closed audit for the ST052-M unified-velocity API promotion seam.

This module governs representation and delivery claims only.  It does not alter
ST052-M, Eq45, CR001 thresholds, pressure, forcing, or any numerical field.
"""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any


CONTRACT_PATH = Path("configs/st052m_unified_velocity_api_gate.json")


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _expect(errors: list[str], condition: bool, message: str) -> None:
    if not condition:
        errors.append(message)


def _capsule_api_facts(source: str) -> dict[str, bool]:
    tree = ast.parse(source)
    top_functions = {
        node.name: node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}

    loader = top_functions.get("load_bundle_runtime")
    exact_source_root_required = False
    if loader is not None:
        kw_names = [arg.arg for arg in loader.args.kwonlyargs]
        if "exact_source_root" in kw_names:
            idx = kw_names.index("exact_source_root")
            exact_source_root_required = loader.args.kw_defaults[idx] is None

    candidate_class = classes.get("St052LinearTemporalWholeCandidate")
    candidate_method = False
    if candidate_class is not None:
        candidate_method = any(
            isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == "velocity"
            for node in candidate_class.body
        )

    return {
        "top_level_loader_present": loader is not None,
        "exact_source_root_required": exact_source_root_required,
        "candidate_velocity_method_present": candidate_method,
        "top_level_velocity_function_present": "velocity" in top_functions,
    }


def audit_contract(
    contract: dict[str, Any],
    *,
    acceptance: dict[str, Any],
    project_status: dict[str, Any],
    readiness_contract: dict[str, Any],
    velocity_delivery_contract: dict[str, Any],
    capsule_source: str,
) -> list[str]:
    """Return semantic violations; an empty list is the only passing result."""
    errors: list[str] = []

    _expect(
        errors,
        contract.get("schema") == "st052m-unified-velocity-api-gate/v1",
        "unexpected ST052 unified-API gate schema",
    )
    classes = contract.get("source_classification")
    _expect(errors, isinstance(classes, dict), "source classification missing")
    if isinstance(classes, dict):
        _expect(
            errors,
            set(classes) == {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"},
            "four-way source classification must be exact",
        )
        _expect(
            errors,
            classes.get("public_source_fact") == [],
            "this repository-only increment must not invent a new public-source fact",
        )

    current = contract.get("current_st052_delivery_scope", {})
    _expect(errors, current.get("candidate_id") == "ST052-M-linear-temporal-child-v1", "candidate identity drift")
    for key in (
        "versioned_whole_child_bundle_materialized",
        "identity_preserving_save_load_with_authenticated_exact_source_runtime",
        "external_runtime_dependency_accepted_for_visualization_delivery",
        "sampled_grid_export_materialized",
        "candidate_local_velocity_callable",
        "caller_must_supply_exact_source_root",
    ):
        _expect(errors, current.get(key) is True, f"current ST052 evidence must keep {key}=true")
    for key in (
        "repository_unified_st052_velocity_entrypoint_registered",
        "standalone_package_parent_runtime_ready",
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _expect(errors, current.get(key) is False, f"current ST052 state must keep {key}=false")

    semantics = contract.get("readiness_semantics", {})
    expected_required = [
        "versioned_candidate_artifact",
        "unified_velocity_api",
        "finite_velocity_output",
        "reproducible_export",
    ]
    _expect(
        errors,
        semantics.get("unqualified_velocity_export_ready_required_evidence") == expected_required,
        "unqualified velocity readiness evidence drift",
    )
    for key in (
        "sampled_grid_alone_is_sufficient",
        "external_source_checkout_reconstruction_alone_is_sufficient",
        "accepted_external_runtime_plus_candidate_local_loader_is_unified_api",
        "accepted_external_runtime_dependency_implies_standalone_package_runtime",
        "accepted_external_runtime_dependency_implies_velocity_export_ready",
        "candidate_local_callable_implies_velocity_export_ready",
        "a8_delivery_smoke_pass_alone_implies_unified_api_registration",
        "pde_failure_or_pending_blocks_velocity_delivery",
        "visual_correspondence_pending_blocks_velocity_delivery",
    ):
        _expect(errors, semantics.get(key) is False, f"forbidden readiness inference enabled: {key}")

    promotion = contract.get("promotion_gate", {})
    _expect(
        errors,
        promotion.get("state") == "blocked_on_repository_unified_st052_velocity_api_registration",
        "ST052 promotion blocker must remain the unified API seam",
    )
    _expect(errors, promotion.get("upgrade_requires_explicit_contract_update") is True, "promotion must fail closed")
    not_required = promotion.get("not_required_for_velocity_export_ready", [])
    for item in ("PDE acceptance", "visual correspondence verification", "paper exactness", "OpenAI-field identity"):
        _expect(errors, item in not_required, f"velocity delivery must not be blocked on {item}")

    eq45 = contract.get("canonical_eq45_independence", {})
    _expect(errors, eq45.get("candidate_family") == "eq45_supported_velocity_candidate_v1", "Eq45 family drift")
    _expect(errors, eq45.get("velocity_api") == "openai_ns_reconstruction.eq45_supported_delivery:velocity", "Eq45 API drift")
    _expect(errors, eq45.get("velocity_export_ready") is True, "ST052 governance must not downgrade canonical Eq45 delivery")
    _expect(errors, eq45.get("visualization_ready") is False, "Eq45 visualization state promotion is forbidden here")
    _expect(errors, eq45.get("pde_validated") is False, "Eq45 PDE promotion is forbidden here")

    cr001 = contract.get("cr001_nonmutation", {})
    _expect(errors, cr001.get("constraints_git_blob_sha1") == "6c559e42895a606e2ef025ade4cb448966d75814", "canonical CR001 blob drift")
    _expect(errors, cr001.get("nu") == 0.01, "CR001 viscosity drift")
    _expect(errors, cr001.get("time_interval") == [0.25, 0.75], "CR001 time interval drift")
    _expect(errors, cr001.get("derivative_steps") == [0.02, 0.01, 0.005], "CR001 FD ladder drift")
    _expect(errors, cr001.get("quadrature_orders") == [24, 48, 96], "CR001 quadrature ladder drift")
    for key in ("momentum_max", "momentum_L2"):
        _expect(errors, cr001.get(key) == 0.001, f"CR001 {key} drift")
    for key in ("divergence_max", "divergence_L2"):
        _expect(errors, cr001.get(key) == 0.00001, f"CR001 {key} drift")
    for key in (
        "residual_defined_free_forcing_forbidden",
        "candidate_amplitude_collapse_forbidden",
        "post_hoc_threshold_relaxation_forbidden",
    ):
        _expect(errors, cr001.get(key) is True, f"CR001 prohibition weakened: {key}")
    for key in (
        "scientific_parameters_changed",
        "thresholds_changed",
        "forcing_contract_changed",
        "validation_sample_changed",
    ):
        _expect(errors, cr001.get(key) is False, f"nonmutation contract violated: {key}")

    acceptance_truth = acceptance.get("truth_boundary", {})
    _expect(
        errors,
        acceptance_truth.get("external_runtime_dependency_accepted_for_visualization_delivery") is True,
        "merged A8 runtime acceptance must remain explicit",
    )
    _expect(errors, acceptance_truth.get("standalone_package_parent_runtime_ready") is False, "external runtime acceptance is not standalone-package readiness")
    for key in ("velocity_export_ready", "visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _expect(errors, acceptance_truth.get(key) is False, f"runtime acceptance cannot promote {key}")
    _expect(errors, "A8-DELIVERY-02" in acceptance.get("unblocks", []), "runtime acceptance must unblock the A8 delivery smoke lane")

    project_states = project_status.get("states", {})
    _expect(errors, project_status.get("candidate_family") == "eq45_supported_velocity_candidate_v1", "canonical project candidate drift")
    _expect(errors, project_status.get("velocity_api") == "openai_ns_reconstruction.eq45_supported_delivery:velocity", "canonical project velocity API drift")
    _expect(errors, project_states.get("velocity_export_ready") is True, "canonical Eq45 velocity readiness must remain true")
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _expect(errors, project_states.get(key) is False, f"canonical project state unexpectedly promoted: {key}")
    st052 = project_status.get("latest_st052m_visual_candidate_evidence", {})
    for key in ("source_runtime_backed_callable_save_load_ready", "grid_export_materialized"):
        _expect(errors, st052.get(key) is True, f"project ST052 evidence missing: {key}")
    for key in ("standalone_package_parent_runtime_ready", "velocity_export_ready", "visualization_ready", "visual_correspondence_verified", "pde_validated"):
        _expect(errors, st052.get(key) is False, f"project ST052 state unexpectedly promoted: {key}")

    readiness = readiness_contract.get("readiness_semantics", {}).get("unqualified_velocity_export_ready", {})
    _expect(errors, readiness.get("required_evidence") == expected_required, "repository readiness contract evidence mismatch")
    _expect(errors, readiness.get("external_source_checkout_reconstruction_alone_is_sufficient") is False, "external checkout alone cannot qualify as unqualified velocity delivery")

    delivery_primary = velocity_delivery_contract.get("primary_deliverable", {})
    delivery_claim = velocity_delivery_contract.get("claim_status", {})
    delivery_gates = velocity_delivery_contract.get("claim_gates", {})
    _expect(errors, delivery_primary.get("api") == "openai_ns_reconstruction.eq45_supported_delivery:velocity", "registered repository delivery API drift")
    _expect(errors, delivery_claim.get("velocity_export_ready") is True, "canonical delivery contract must remain callable")
    _expect(errors, delivery_gates.get("pde_validation", {}).get("blocking_for_velocity_delivery") is False, "PDE gate must not block callable delivery")
    _expect(errors, delivery_gates.get("visual_correspondence", {}).get("blocking_for_velocity_delivery") is False, "visual correspondence gate must not block callable delivery")

    api = _capsule_api_facts(capsule_source)
    _expect(errors, api["top_level_loader_present"], "ST052 bundle loader disappeared")
    _expect(errors, api["exact_source_root_required"], "current ST052 loader must still require an explicit exact_source_root")
    _expect(errors, api["candidate_velocity_method_present"], "ST052 candidate-local velocity method disappeared")
    _expect(errors, not api["top_level_velocity_function_present"], "a top-level ST052 velocity API appeared without updating this promotion contract")

    mutation = contract.get("mutation_scope", {})
    for key in (
        "velocity_changed",
        "candidate_changed",
        "pressure_changed",
        "forcing_changed",
        "optimizer_changed",
        "validation_sample_changed",
        "norm_definition_changed",
        "threshold_changed",
        "canonical_eq45_delivery_changed",
        "st052_scientific_state_promoted",
    ):
        _expect(errors, mutation.get(key) is False, f"governance-only scope violated: {key}")

    return errors


def audit_repository(root: str | Path = ".") -> dict[str, Any]:
    """Audit the exact repository snapshot and return a machine-readable receipt."""
    root = Path(root)
    contract = _read_json(root / CONTRACT_PATH)
    pinned = contract["pinned_repository_evidence"]
    blob_errors: list[str] = []
    for path_key, sha_key in (
        ("runtime_acceptance_contract", "runtime_acceptance_contract_git_blob_sha1"),
        ("whole_child_module", "whole_child_module_git_blob_sha1"),
        ("delivery_readiness_contract", "delivery_readiness_contract_git_blob_sha1"),
        ("velocity_delivery_contract", "velocity_delivery_contract_git_blob_sha1"),
        ("project_status", "project_status_git_blob_sha1"),
        ("canonical_constraints", "canonical_constraints_git_blob_sha1"),
    ):
        path = root / pinned[path_key]
        observed = _git_blob_sha1(path)
        if observed != pinned[sha_key]:
            blob_errors.append(f"pinned blob mismatch for {pinned[path_key]}: {observed}")

    errors = audit_contract(
        contract,
        acceptance=_read_json(root / pinned["runtime_acceptance_contract"]),
        project_status=_read_json(root / pinned["project_status"]),
        readiness_contract=_read_json(root / pinned["delivery_readiness_contract"]),
        velocity_delivery_contract=_read_json(root / pinned["velocity_delivery_contract"]),
        capsule_source=(root / pinned["whole_child_module"]).read_text(encoding="utf-8"),
    )
    errors = blob_errors + errors
    return {
        "schema": "st052m-unified-velocity-api-gate-audit/v1",
        "task_id": contract["task_id"],
        "passed": not errors,
        "violations": errors,
        "truth_boundary": {
            "st052_velocity_export_ready": False,
            "canonical_eq45_velocity_export_ready": True,
            "pde_pending_blocks_callable_delivery": False,
            "visual_correspondence_pending_blocks_callable_delivery": False,
        },
    }


def main() -> int:
    receipt = audit_repository(Path.cwd())
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0 if receipt["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
