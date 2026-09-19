"""Fail-closed CR002 governance for ST052-M whole-candidate runtime identity scope.

PR #632 materializes a checksum-bound whole-child bundle and can reload it when
an exact historical ST052 source runtime is supplied.  This module keeps that
useful bundle/save-load receipt separate from a stronger standalone executable
function-identity claim: the external runtime code itself is not part of the
bundle identity and is not authenticated by the loader at this snapshot.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST052M-WHOLE-CANDIDATE-RUNTIME-BINDING-072"
ACTIVE_HEAD = "c76402986533346bc200234e8e27dd2f05ba8ccd"
SOURCE_PR = 632
SOURCE_PR_HEAD = "46bbc10bb6b88ae637ac1559fec58c961013d987"
SOURCE_PARENT_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
WHOLE_ID = "333615d89423a22995d9356b29b604d1c241c1b40a631bd66892c2701c536064"
ALLOWED_SOURCE_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


class GovernanceError(ValueError):
    """Raised when a contract mutation widens the frozen evidence scope."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GovernanceError(message)


def load_contract(path: str | Path | None = None) -> dict[str, Any]:
    if path is None:
        path = (
            Path(__file__).resolve().parents[2]
            / "configs"
            / "st052m_whole_candidate_runtime_binding_contract.json"
        )
    obj = json.loads(Path(path).read_text(encoding="utf-8"))
    _require(isinstance(obj, dict), "contract must be a JSON object")
    return obj


def audit_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the frozen truth boundary and return a detached copy."""
    _require(contract.get("task_id") == TASK_ID, "task identity drift")

    snapshot = contract.get("snapshot", {})
    _require(snapshot.get("active_integration_head") == ACTIVE_HEAD, "active head drift")
    _require(snapshot.get("source_pr") == SOURCE_PR, "source PR drift")
    _require(snapshot.get("source_pr_head") == SOURCE_PR_HEAD, "source PR head drift")
    _require(snapshot.get("source_parent_head") == SOURCE_PARENT_HEAD, "source parent head drift")
    _require(snapshot.get("dedicated_workflow_run") == 35442084024, "dedicated CI identity drift")
    _require(snapshot.get("standard_workflow_run") == 35442084021, "standard CI identity drift")
    _require(snapshot.get("dedicated_workflow_success") is True, "dedicated CI must be recorded success")
    _require(snapshot.get("standard_workflow_success") is True, "standard CI must be recorded success")
    _require(snapshot.get("whole_candidate_identity_sha256") == WHOLE_ID, "bundle identity drift")
    _require(snapshot.get("dedicated_parity_max_component_error") == 0.0, "parity receipt drift")
    _require(snapshot.get("dedicated_parity_threshold") == 5e-12, "parity threshold drift")

    allowed = contract.get("allowed_source_classes")
    _require(isinstance(allowed, list) and set(allowed) == ALLOWED_SOURCE_CLASSES, "source-class vocabulary drift")
    classes = contract.get("source_classification", {})
    _require(set(classes.values()) <= ALLOWED_SOURCE_CLASSES, "unknown source class")
    _require(classes.get("callable_save_load_velocity_delivery") == "user_requirement", "delivery source class drift")
    _require(classes.get("public_openai_visualization_structure") == "public_source_fact", "public visualization source class drift")
    _require(classes.get("st052_parent_and_temporal_capsule_architecture") == "autonomous_design", "autonomous ST052 source class drift")
    _require(classes.get("openai_hidden_velocity_coefficients") == "pending_unknown", "hidden-field source class laundering")

    bundle = contract.get("observed_whole_child_bundle", {})
    for key in (
        "bundle_file_checksums_verified",
        "composite_bundle_identity_exists",
        "parent_candidate_bytes_bound",
        "parent_validation_bytes_bound",
        "parent_replay_manifest_bytes_bound",
        "static_transform_spec_bound",
        "temporal_transform_spec_bound",
        "whole_child_bundle_materialized",
        "whole_child_save_load_ready_with_exact_source_runtime",
    ):
        _require(bundle.get(key) is True, f"expected positive bundle fact missing: {key}")

    runtime = contract.get("external_runtime_binding", {})
    _require(runtime.get("external_exact_source_runtime_required") is True, "external runtime requirement erased")
    _require(runtime.get("source_recipe_blob_sha1_verified_on_load") is True, "recipe verification fact erased")
    _require(runtime.get("dedicated_ci_checks_out_exact_source_head") is True, "exact CI checkout fact erased")
    _require(runtime.get("source_native_parent_probe_binding_applied") is True, "probe binding fact erased")
    for key in (
        "source_git_head_verified_on_load",
        "source_git_tree_sha_bound_by_bundle",
        "replay_module_sha_bound_by_bundle",
        "ancestor_runtime_code_shas_bound_by_bundle",
        "source_runtime_code_embedded_in_bundle",
        "parent_probe_binding_is_whole_domain_function_identity",
        "runtime_identity_closed",
    ):
        _require(runtime.get(key) is False, f"runtime identity laundering: {key}")

    identity = contract.get("identity_scope", {})
    _require(
        identity.get("whole_candidate_identity_sha256_scope")
        == "bundled parent evidence bytes plus static and temporal transform identities",
        "bundle identity scope drift",
    )
    for key in (
        "whole_candidate_identity_sha256_is_complete_executable_function_identity",
        "same_bundle_identity_guarantees_same_external_runtime_code",
        "same_bundle_identity_alone_guarantees_same_velocity_function",
        "dedicated_ci_parity_generalizes_to_arbitrary_exact_source_root_argument",
        "green_ci_implies_standalone_runtime_identity_closed",
    ):
        _require(identity.get(key) is False, f"identity overclaim: {key}")

    states = contract.get("delivery_states", {})
    _require(states.get("canonical_eq45_velocity_export_ready") is True, "canonical Eq45 readiness drift")
    for key in (
        "standalone_package_parent_runtime_ready",
        "standalone_reproducible_velocity_identity_ready",
        "experimental_st052_velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "production_candidate_selected",
        "held_out_temporal_child_pde_residual_evaluated",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(states.get(key) is False, f"premature truth-state promotion: {key}")

    promotion = contract.get("promotion_requirements", {})
    for key in (
        "runtime_identity_must_be_bound_before_standalone_export",
        "loader_must_fail_on_runtime_identity_mismatch",
        "complete_identity_must_include_or_immutably_reference_runtime_identity",
        "reload_velocity_parity_must_be_replayed_after_runtime_closure",
    ):
        _require(promotion.get(key) is True, f"runtime-closure prerequisite removed: {key}")
    _require(promotion.get("pde_validation_required_for_velocity_export_ready") is False, "PDE must remain independent of callable delivery")
    _require(promotion.get("visual_correspondence_required_for_velocity_export_ready") is False, "visual correspondence must remain independent of callable delivery")
    routes = promotion.get("allowed_runtime_closure_routes")
    _require(isinstance(routes, list) and len(routes) == 2, "runtime closure routes drift")

    cr001 = contract.get("cr001_frozen", {})
    expected = {
        "nu": 0.01,
        "physical_domain": "R^3",
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25, 0.75],
        "force_mode": "restricted_two_parameter_family",
        "reference_energy": 1.0,
        "reference_energy_abs_tolerance": 0.001,
        "validation_seed": 914027,
        "held_out_points": 4096,
        "validation_times": [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75],
        "derivative_steps": [0.02, 0.01, 0.005],
        "energy_quadrature_orders": [24, 48, 96],
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
        "momentum_max": 0.001,
        "momentum_L2": 0.001,
        "free_residual_defined_forcing_allowed": False,
        "amplitude_collapse_success_allowed": False,
    }
    for key, value in expected.items():
        _require(cr001.get(key) == value, f"CR001 drift: {key}")
    _require(cr001.get("force_bounds") == {"a": [0.0, 10.0], "c": [0.0, 10.0]}, "force bounds drift")
    _require(
        cr001.get("evaluation_box") == [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
        "evaluation box drift",
    )

    mutation = contract.get("mutation_scope", {})
    for key, value in mutation.items():
        _require(value is False, f"governance-only increment mutated {key}")

    return copy.deepcopy(dict(contract))


def audit_default_contract() -> dict[str, Any]:
    return audit_contract(load_contract())


__all__ = [
    "GovernanceError",
    "TASK_ID",
    "audit_contract",
    "audit_default_contract",
    "load_contract",
]
