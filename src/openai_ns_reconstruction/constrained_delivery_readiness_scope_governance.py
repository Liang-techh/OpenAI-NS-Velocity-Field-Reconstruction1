"""Fail-closed CR002 governance for delivery-readiness scope.

Repository-wide ``velocity_export_ready`` names callable candidate delivery.  A
sampled visualization grid may have its own export readiness, but cannot borrow
that unqualified state merely because it is serializable and reloadable.  This
module changes no velocity, candidate, route, or scientific threshold.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-DELIVERY-READINESS-SCOPE-079"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_CANONICAL_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
_OBSERVED_PR = 476
_OBSERVED_HEAD = "d068e0d4c2cdce23a08aec791af8d2d59e459e65"
_EXACT_SOURCE_WORKFLOW_RUN = 35389278384


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "delivery_readiness_scope_contract.json")


def load_delivery_state_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "delivery_state_contract.json")


def load_grid_identity_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "velocity_grid_delivery_identity_contract.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def audit_delivery_readiness_scope(
    contract: Mapping[str, Any] | None = None,
    delivery_state_contract: Mapping[str, Any] | None = None,
    grid_identity_contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit scoped readiness names against current delivery and CR001 contracts."""

    contract = deepcopy(dict(load_contract() if contract is None else contract))
    delivery_state = deepcopy(
        dict(load_delivery_state_contract() if delivery_state_contract is None else delivery_state_contract)
    )
    grid_identity = deepcopy(
        dict(load_grid_identity_contract() if grid_identity_contract is None else grid_identity_contract)
    )
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(contract.get("schema_version"), 1, "contract schema")
    _require_equal(contract.get("task_id"), TASK_ID, "contract task id")

    snapshot = contract.get("snapshot", {})
    _require_equal(snapshot.get("active_integration_branch"), "codex/cr001-constraints", "integration branch")
    _require_equal(snapshot.get("observed_pr"), _OBSERVED_PR, "observed PR")
    _require_equal(snapshot.get("observed_pr_state"), "open_unintegrated", "observed PR state")
    _require_equal(snapshot.get("observed_pr_base"), "main", "observed PR base")
    _require_equal(snapshot.get("observed_pr_head"), _OBSERVED_HEAD, "observed PR head")

    _require_equal(set(contract.get("allowed_source_classes", ())), _ALLOWED_CLASSES, "source-class vocabulary")
    source = contract.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    for key in ("callable_velocity_x_y_z_t", "save_load_delivery"):
        _require_equal(source.get(key), "user_requirement", f"source class {key}")
    _require_equal(
        source.get("public_openai_visualization_observables"),
        "public_source_fact",
        "public visualization source class",
    )
    for key in (
        "readiness_state_names_and_scope_rules",
        "netcdf_format_grid_resolution_times_and_checksum",
        "st051b_reconstruction_and_swirl_redistribution",
        "cr001_validation_protocol",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    for key in (
        "openai_hidden_numerical_velocity",
        "openai_hidden_parameters",
        "openai_hidden_frame_time_camera_mapping",
    ):
        _require_equal(source.get(key), "pending_unknown", f"source class {key}")

    readiness = contract.get("readiness_semantics", {})
    unqualified = readiness.get("unqualified_velocity_export_ready", {})
    sampled = readiness.get("sampled_grid_export_ready", {})
    local_callable = readiness.get("candidate_local_callable_reconstruction_ready", {})

    canonical_state = delivery_state.get("states", {}).get("velocity_export_ready", {})
    _require_equal(
        unqualified.get("required_evidence"),
        canonical_state.get("positive_evidence"),
        "unqualified readiness evidence",
    )
    _require_equal(unqualified.get("required_api_shape"), "velocity(x,y,z,t)->[u,v,w]", "unified API shape")
    for key in (
        "sampled_grid_alone_is_sufficient",
        "auxiliary_points_time_callable_alone_is_sufficient",
        "external_source_checkout_reconstruction_alone_is_sufficient",
    ):
        _require_equal(unqualified.get(key), False, key)
    for key in (
        "does_not_imply_candidate_velocity_export_ready",
        "does_not_imply_visualization_ready",
        "does_not_imply_visual_correspondence_verified",
        "does_not_imply_pde_validated",
        "does_not_imply_paper_exact",
    ):
        _require_equal(sampled.get(key), True, f"sampled-grid scope {key}")
    for key in (
        "does_not_imply_unified_velocity_api",
        "does_not_imply_save_load_candidate_artifact",
        "does_not_imply_unqualified_velocity_export_ready",
    ):
        _require_equal(local_callable.get(key), True, f"candidate-local scope {key}")

    existing_grid = grid_identity.get("identity_layers", {}).get("sampled_grid", {})
    _require_equal(
        existing_grid.get("sampled_grid_is_continuum_velocity_api"),
        False,
        "existing grid/continuum separation",
    )
    _require_equal(
        existing_grid.get("off_grid_or_unsampled_time_interpolation_requires_explicit_new_numerical_contract"),
        True,
        "existing interpolation contract",
    )

    observed = contract.get("observed_pr476", {})
    _require_equal(observed.get("artifact_role"), "candidate_specific_sampled_visualization_grid", "PR476 artifact role")
    _require_equal(observed.get("netcdf_default_resolution"), 33, "PR476 grid resolution")
    _require_equal(observed.get("netcdf_times"), [0.25, 0.5, 0.75], "PR476 grid times")
    _require_equal(observed.get("auxiliary_callable_shape"), "velocity(points,time)->(n,3)", "PR476 callable shape")
    _require_equal(observed.get("requires_external_candidate_root_checkout"), True, "PR476 external checkout")
    _require_equal(observed.get("creates_versioned_candidate_artifact"), False, "PR476 candidate artifact")
    _require_equal(observed.get("registers_repository_unified_velocity_api"), False, "PR476 unified API registration")
    _require_equal(observed.get("creates_checksum_guarded_sampled_grid"), True, "PR476 grid artifact")
    _require_equal(observed.get("sampled_grid_export_ready_status"), "verified_exact_source_ci", "PR476 grid readiness")
    _require_equal(observed.get("exact_source_workflow_run"), _EXACT_SOURCE_WORKFLOW_RUN, "PR476 exact-source workflow run")
    _require_equal(observed.get("exact_source_workflow_conclusion"), "success", "PR476 exact-source workflow conclusion")
    _require_equal(observed.get("pr_body_claims_unqualified_velocity_export_ready_true"), True, "PR476 body claim")
    _require_equal(
        observed.get("netcdf_metadata_sets_unqualified_velocity_export_ready_one"),
        True,
        "PR476 metadata claim",
    )
    _require_equal(
        observed.get("unqualified_velocity_export_ready_claim_is_contract_compatible"),
        False,
        "PR476 readiness compatibility",
    )
    _require("sampled_grid_export_ready" in observed.get("required_semantic_remediation", ""), "missing scoped readiness remediation")
    for key in (
        "grid_success_may_promote_candidate_velocity_export_ready",
        "grid_success_may_replace_canonical_velocity_delivery",
        "grid_success_may_promote_visual_correspondence_verified",
        "grid_success_may_promote_pde_validated",
    ):
        _require_equal(observed.get(key), False, key)

    canonical = contract.get("canonical_callable_delivery", {})
    _require_equal(canonical.get("candidate_family"), project_status.get("candidate_family"), "canonical family")
    _require_equal(canonical.get("candidate_sha256"), project_status.get("candidate_sha256"), "canonical SHA")
    _require_equal(canonical.get("candidate_sha256"), _CANONICAL_SHA, "known canonical SHA")
    _require_equal(canonical.get("velocity_api"), project_status.get("velocity_api"), "canonical velocity API")
    _require_equal(canonical.get("velocity_api"), _CANONICAL_API, "known canonical velocity API")
    _require_equal(canonical.get("velocity_export_ready"), True, "canonical export readiness")
    _require_equal(canonical.get("open_pr476_changes_canonical_candidate"), False, "PR476 canonical candidate separation")
    _require_equal(canonical.get("open_pr476_changes_canonical_velocity_api"), False, "PR476 canonical API separation")
    _require_equal(canonical.get("pde_failure_blocks_callable_velocity_delivery"), False, "delivery/PDE separation")

    domain = constraints.get("domain", {})
    forcing = constraints.get("forcing", {})
    nontriviality = constraints.get("nontriviality", {})
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
    cr001 = contract.get("canonical_cr001", {})
    expected_pairs = (
        (cr001.get("nu"), constraints.get("nu"), "nu"),
        (cr001.get("physical_domain"), domain.get("physical"), "physical domain"),
        (cr001.get("evaluation_box"), domain.get("evaluation_box"), "evaluation box"),
        (cr001.get("support"), domain.get("support"), "support"),
        (cr001.get("time_interval"), domain.get("time_interval"), "time interval"),
        (cr001.get("force_mode"), forcing.get("mode"), "force mode"),
        (cr001.get("force_bounds"), forcing.get("parameters"), "force bounds"),
        (cr001.get("reference_energy"), nontriviality.get("reference_energy"), "reference energy"),
        (cr001.get("reference_energy_abs_tolerance"), nontriviality.get("reference_energy_abs_tolerance"), "energy tolerance"),
        (cr001.get("validation_seed"), validation.get("seed"), "validation seed"),
        (cr001.get("held_out_points"), validation.get("held_out_points"), "held-out points"),
        (cr001.get("validation_times"), validation.get("times"), "validation times"),
        (cr001.get("derivative_steps"), validation.get("derivative_steps"), "derivative ladder"),
        (cr001.get("divergence_max"), thresholds.get("divergence_max"), "divergence max"),
        (cr001.get("divergence_L2"), thresholds.get("divergence_L2"), "divergence L2"),
        (cr001.get("pde_residual_max"), thresholds.get("pde_residual_max"), "PDE max"),
        (cr001.get("pde_residual_L2"), thresholds.get("pde_residual_L2"), "PDE L2"),
    )
    for expected, actual, label in expected_pairs:
        _require_equal(actual, expected, label)
    _require("No residual-dependent basis" in forcing.get("restriction", ""), "free-force restriction weakened")

    truth = contract.get("truth_boundary", {})
    states = project_status.get("states", {})
    _require_equal(truth.get("canonical_velocity_export_ready"), True, "canonical export truth")
    _require_equal(states.get("velocity_export_ready"), True, "project export truth")
    _require_equal(truth.get("pr476_candidate_velocity_export_ready"), False, "PR476 candidate export truth")
    _require_equal(truth.get("pr476_sampled_grid_export_ready"), True, "PR476 sampled-grid truth")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(truth.get(key), False, f"contract truth {key}")
        _require_equal(states.get(key), False, f"project truth {key}")

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
        "active_route_changed",
        "canonical_velocity_api_changed",
    ):
        _require_equal(mutation.get(key), False, key)

    return {
        "task_id": TASK_ID,
        "status": "governance_pass",
        "observed_pr": _OBSERVED_PR,
        "unqualified_velocity_export_ready_requires_unified_candidate_delivery": True,
        "sampled_grid_has_separate_readiness_scope": True,
        "pr476_unqualified_readiness_claim_contract_compatible": False,
        "pr476_candidate_velocity_export_ready": False,
        "pr476_sampled_grid_export_ready": True,
        "pr476_exact_source_workflow_run": _EXACT_SOURCE_WORKFLOW_RUN,
        "canonical_velocity_export_ready": True,
        "canonical_velocity_api": project_status.get("velocity_api"),
        "canonical_thresholds_unchanged": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
    }


def main() -> None:
    print(json.dumps(audit_delivery_readiness_scope(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
