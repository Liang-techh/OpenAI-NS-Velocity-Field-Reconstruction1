"""Fail-closed CR002 governance for derived-challenger lineage.

A child experiment may be stacked on an open challenger, but branch ancestry is
not scientific identity and parent evidence does not transfer through a changed
velocity representation.  This module changes no field or scientific gate.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-DERIVED-CHALLENGER-LINEAGE-063"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_PARENT_SHA = "6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09"
_CANONICAL_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "derived_challenger_lineage_governance.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def validate_parent_lineage_record(record: Mapping[str, Any]) -> None:
    """Reject branch-only provenance for any child proposed for promotion."""
    commit = record.get("parent_git_commit")
    representation = record.get("parent_representation_identity")
    candidate_sha = record.get("parent_candidate_sha256")
    _require(isinstance(commit, str) and len(commit) == 40, "immutable parent git commit required")
    _require(isinstance(representation, str) and bool(representation.strip()), "parent representation identity required")
    _require(isinstance(candidate_sha, str) and len(candidate_sha) == 64, "parent candidate sha256 required")
    _require(record.get("parent_branch") != commit, "branch name cannot substitute for immutable parent commit")


def audit_derived_challenger_lineage(
    contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit stacked-child lineage, delivery separation, and unchanged CR001 gates."""
    contract = deepcopy(dict(load_contract() if contract is None else contract))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(contract.get("schema_version"), 1, "contract schema")
    _require_equal(contract.get("task_id"), TASK_ID, "contract task id")

    snapshot = contract.get("snapshot", {})
    _require_equal(snapshot.get("active_integration_branch"), "codex/cr001-constraints", "integration branch")
    _require_equal(snapshot.get("active_integration_head"), "a2b215790fc483766be712233a134cd48ad7c149", "integration snapshot head")
    _require_equal(snapshot.get("observed_child_pr"), 407, "observed child PR")
    _require_equal(snapshot.get("observed_child_state"), "open_unintegrated", "observed child state")
    _require_equal(snapshot.get("observed_child_head"), "2d7cb2fe05696c6303e09ad20b91c61aad5f977e", "observed child head")
    _require_equal(snapshot.get("observed_parent_pr"), 398, "observed parent PR")
    _require_equal(snapshot.get("observed_parent_state"), "open_unintegrated", "observed parent state")
    _require_equal(snapshot.get("observed_parent_head"), "107fb9aaa75edcf07225d664a38e0b8a779b2dcf", "observed parent head")
    _require_equal(snapshot.get("parent_analytic_candidate_id"), "ST048-S", "parent candidate id")
    _require_equal(snapshot.get("parent_analytic_candidate_sha256"), _PARENT_SHA, "parent candidate sha")

    validate_parent_lineage_record({
        "parent_branch": snapshot.get("observed_parent_branch"),
        "parent_git_commit": snapshot.get("observed_parent_head"),
        "parent_representation_identity": snapshot.get("parent_analytic_candidate_id"),
        "parent_candidate_sha256": snapshot.get("parent_analytic_candidate_sha256"),
    })

    _require_equal(set(contract.get("allowed_source_classes", ())), _ALLOWED_CLASSES, "source-class vocabulary")
    source = contract.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    _require_equal(source.get("callable_velocity_x_y_z_t"), "user_requirement", "callable velocity class")
    _require_equal(source.get("public_openai_visualization_observables"), "public_source_fact", "public visualization class")
    for key in (
        "st048_candidate_recipe_and_repository_measurements",
        "piola_warp_formula_beta_grid_and_capacity_metrics",
        "derived_challenger_lineage_policy",
        "cr001_validation_protocol",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    for key in ("openai_hidden_coordinate_map", "openai_hidden_numerical_velocity", "openai_hidden_parameters"):
        _require_equal(source.get(key), "pending_unknown", f"source class {key}")

    labels = contract.get("evidence_relation_labels", {})
    for label in ("immutable_parent_commit", "parent_candidate_identity", "derived_child_identity", "protocol_scoped_challenger"):
        _require(label in labels, f"missing evidence relation label {label}")
    _require_equal(contract.get("relation_labels_are_source_classes"), False, "relation/source-class separation")

    lineage = contract.get("lineage_contract", {})
    for key in (
        "mutable_branch_name_is_sufficient_parent_identity",
        "child_may_reuse_parent_candidate_sha256",
        "parent_validation_receipt_transfers_automatically",
        "parent_pressure_or_forcing_transfers_automatically",
        "parent_retained_baseline_status_transfers_automatically",
        "parent_canonical_delivery_status_transfers_automatically",
        "parent_pde_status_transfers_automatically",
        "parent_visual_status_transfers_automatically",
        "parent_paper_exact_status_transfers_automatically",
        "parent_openai_identity_status_transfers_automatically",
    ):
        _require_equal(lineage.get(key), False, key)
    for key in (
        "promotion_requires_parent_git_commit",
        "promotion_requires_parent_representation_identity",
        "promotion_requires_parent_candidate_sha256",
        "materialized_child_requires_new_representation_identity",
        "materialized_child_requires_new_candidate_sha256",
    ):
        _require_equal(lineage.get(key), True, key)

    dependency = contract.get("unintegrated_parent_dependency", {})
    for key in (
        "downstream_green_ci_promotes_parent",
        "downstream_merge_promotes_parent",
        "downstream_morphology_improvement_promotes_parent",
        "downstream_divergence_preflight_promotes_parent",
        "branch_tip_only_dependency_is_acceptable_for_promotion",
    ):
        _require_equal(dependency.get(key), False, key)
    _require_equal(dependency.get("promotion_requires_explicit_dependency_resolution"), True, "dependency resolution requirement")
    _require_equal(
        set(dependency.get("allowed_dependency_resolution_modes", ())),
        {
            "integrate_or_vendor_exact_parent_representation_with_immutable_identity",
            "materialize_self_contained_child_from_immutable_parent_identity",
        },
        "dependency resolution modes",
    )
    _require_equal(
        dependency.get("unresolved_stacked_dependency_role"),
        "protocol_scoped_challenger_evidence_only",
        "unresolved dependency role",
    )

    promotion = contract.get("child_promotion_boundary", {})
    for key in (
        "open_or_stacked_child_replaces_retained_baseline",
        "open_or_stacked_child_replaces_canonical_velocity",
        "capacity_crossing_selects_production_parameter",
        "divergence_preservation_implies_pde_validated",
    ):
        _require_equal(promotion.get(key), False, key)
    for key in (
        "fresh_child_validation_required_after_materialization",
        "explicit_repository_promotion_required_for_retained_baseline",
        "explicit_status_update_required_for_canonical_velocity_replacement",
        "formal_scientific_acceptance_requires_unchanged_cr001_gates",
    ):
        _require_equal(promotion.get(key), True, key)

    child = contract.get("current_child_snapshot", {})
    _require_equal(child.get("diagnostic_beta_crossing"), 0.075, "diagnostic beta crossing")
    for key in (
        "production_beta_selected",
        "child_candidate_sha_created",
        "child_representation_identity_created",
        "fresh_child_full_momentum_evaluated",
        "parent_full_momentum_may_be_inherited",
        "canonical_velocity_changed",
        "retained_baseline_changed",
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
    ):
        _require_equal(child.get(key), False, key)

    delivery = contract.get("canonical_delivery_guard", {})
    _require_equal(delivery.get("candidate_family"), project_status.get("candidate_family"), "canonical delivery family")
    _require_equal(delivery.get("candidate_sha256"), project_status.get("candidate_sha256"), "canonical delivery sha")
    _require_equal(delivery.get("candidate_sha256"), _CANONICAL_SHA, "known canonical delivery sha")
    _require_equal(delivery.get("velocity_api"), project_status.get("velocity_api"), "canonical velocity API")
    _require_equal(delivery.get("velocity_api"), _CANONICAL_API, "known canonical velocity API")
    _require_equal(delivery.get("velocity_export_ready"), True, "velocity export state")
    _require_equal(delivery.get("pde_failure_blocks_callable_velocity_delivery"), False, "delivery/PDE separation")
    _require_equal(delivery.get("stacked_challenger_changes_delivery_automatically"), False, "stacked/delivery separation")

    domain = constraints.get("domain", {})
    forcing = constraints.get("forcing", {})
    nontriviality = constraints.get("nontriviality", {})
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
    cr001 = contract.get("canonical_cr001", {})
    for expected, actual, label in (
        (cr001.get("nu"), constraints.get("nu"), "nu"),
        (cr001.get("physical_domain"), domain.get("physical"), "physical domain"),
        (cr001.get("evaluation_box"), domain.get("evaluation_box"), "evaluation box"),
        (cr001.get("support"), domain.get("support"), "support"),
        (cr001.get("time_interval"), domain.get("time_interval"), "time interval"),
        (cr001.get("force_mode"), forcing.get("mode"), "force mode"),
        (cr001.get("force_bounds"), forcing.get("parameters"), "force bounds"),
        (cr001.get("reference_energy"), nontriviality.get("reference_energy"), "reference energy"),
        (cr001.get("reference_energy_abs_tolerance"), nontriviality.get("reference_energy_abs_tolerance"), "reference energy tolerance"),
        (cr001.get("validation_seed"), validation.get("seed"), "validation seed"),
        (cr001.get("held_out_points"), validation.get("held_out_points"), "held-out points"),
        (cr001.get("validation_times"), validation.get("times"), "validation times"),
        (cr001.get("derivative_steps"), validation.get("derivative_steps"), "derivative ladder"),
        (cr001.get("divergence_max"), thresholds.get("divergence_max"), "divergence max"),
        (cr001.get("divergence_L2"), thresholds.get("divergence_L2"), "divergence L2"),
        (cr001.get("pde_residual_max"), thresholds.get("pde_residual_max"), "PDE max"),
        (cr001.get("pde_residual_L2"), thresholds.get("pde_residual_L2"), "PDE L2"),
    ):
        _require_equal(actual, expected, label)
    _require("No residual-dependent basis" in forcing.get("restriction", ""), "free-force restriction weakened")

    truth = contract.get("truth_boundary", {})
    states = project_status.get("states", {})
    _require_equal(truth.get("velocity_export_ready"), True, "contract export truth")
    _require_equal(states.get("velocity_export_ready"), True, "project export truth")
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

    _require_equal(
        project_status.get("active_scientific_route"),
        "materialize_integrated_axial_cap_poloidal_child_then_fresh_validate",
        "active scientific route",
    )

    return {
        "task_id": TASK_ID,
        "status": "governance_pass",
        "observed_child_pr": 407,
        "observed_parent_pr": 398,
        "parent_identity_is_immutable": True,
        "branch_tip_alone_is_insufficient": True,
        "stacked_child_is_protocol_scoped_only": True,
        "parent_evidence_transfers_automatically": False,
        "canonical_velocity_unchanged": True,
        "canonical_thresholds_unchanged": True,
        "velocity_export_ready": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
    }


def main() -> None:
    print(json.dumps(audit_derived_challenger_lineage(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
