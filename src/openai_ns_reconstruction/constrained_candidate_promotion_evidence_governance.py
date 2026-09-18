"""Fail-closed CR002 governance for candidate-promotion evidence lifecycle.

This module changes no velocity, candidate, pressure, forcing, optimizer, sample,
norm, threshold, or active route. It keeps parent-local improvement,
protocol-scoped benchmark challenge, retained-baseline replacement, canonical
velocity delivery, and scientific acceptance as separate claims.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-CANDIDATE-PROMOTION-LIFECYCLE-061"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_ALLOWED_LIFECYCLE_STATES = {
    "protocol_scoped_challenger",
    "repository_retained_baseline",
    "scientifically_accepted_candidate",
}
_REQUIRED_PROTOCOL_FIELDS = {
    "candidate_freeze_and_identity",
    "sample_points_and_seed",
    "evaluation_times",
    "spatial_and_time_derivative_operator_and_ladders",
    "momentum_residual_definition",
    "reported_norm_definitions",
    "forcing_family_and_bounds",
    "pressure_convention",
    "nu",
    "physical_domain_and_support",
}
_ST006_SHA = "6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3"
_DELIVERY_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "candidate_promotion_evidence_contract.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def audit_candidate_promotion_evidence(
    contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit evidence comparability, lifecycle separation, and truth boundaries."""

    contract = deepcopy(dict(load_contract() if contract is None else contract))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(contract.get("schema_version"), 2, "contract schema")
    _require_equal(contract.get("task_id"), TASK_ID, "contract task id")
    snapshot = contract.get("snapshot", {})
    _require_equal(snapshot.get("active_integration_branch"), "codex/cr001-constraints", "integration branch")
    _require_equal(snapshot.get("stale_superseded_governance_pr"), 308, "superseded governance PR")

    _require_equal(set(contract.get("allowed_source_classes", ())), _ALLOWED_CLASSES, "source-class vocabulary")
    source = contract.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    _require_equal(source.get("callable_velocity_delivery"), "user_requirement", "velocity delivery class")
    _require_equal(source.get("public_openai_visualization_structure"), "public_source_fact", "public visualization class")
    _require_equal(source.get("cr001_validation_protocol"), "autonomous_design", "CR001 protocol class")
    _require_equal(source.get("st006_validation_protocol"), "autonomous_design", "ST006 protocol class")
    _require_equal(source.get("st006_measured_residual_record"), "autonomous_design", "ST006 measurement class")
    _require_equal(source.get("challenger_branch_measurements"), "autonomous_design", "challenger measurement class")
    _require_equal(source.get("candidate_promotion_lifecycle_policy"), "autonomous_design", "promotion policy class")
    _require_equal(source.get("openai_hidden_numerical_profile"), "pending_unknown", "hidden profile class")

    labels = contract.get("evidence_relation_labels", {})
    for label in (
        "project_derived_measurement",
        "frozen_parent_paired_comparison",
        "repository_retained_benchmark",
        "protocol_scoped_challenger",
    ):
        _require(label in labels, f"missing evidence relation label {label}")
    _require_equal(contract.get("relation_labels_are_source_classes"), False, "relation/source-class separation")

    paired = contract.get("paired_parent_claim_scope", {})
    _require_equal(
        paired.get("allowed_claim"),
        "improved over this frozen parent under this exact declared protocol",
        "paired allowed claim",
    )
    for key in (
        "allows_repository_st006_superiority_claim",
        "allows_retained_baseline_replacement",
        "allows_canonical_velocity_replacement",
        "allows_pde_validation_claim",
        "allows_visual_correspondence_claim",
        "allows_paper_exact_claim",
        "allows_openai_field_identification_claim",
    ):
        _require_equal(paired.get(key), False, key)

    benchmark = contract.get("repository_baseline", {})
    _require_equal(benchmark.get("benchmark_name"), "ST006", "benchmark name")
    _require_equal(benchmark.get("lifecycle_state"), "repository_retained_baseline", "benchmark lifecycle")
    _require_equal(benchmark.get("candidate_sha256"), _ST006_SHA, "ST006 candidate sha")
    _require_equal(benchmark.get("retained_seed"), 9172801, "ST006 seed")
    _require_equal(benchmark.get("retained_held_out_points"), 4096, "ST006 held-out count")
    _require_equal(
        benchmark.get("retained_times"),
        [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75],
        "ST006 times",
    )
    _require_equal(benchmark.get("retained_finest_spatial_step"), 0.005, "ST006 finest spatial step")
    _require_equal(benchmark.get("retained_momentum_sampled_max"), 0.1082289305112118, "ST006 momentum max")
    _require_equal(benchmark.get("retained_momentum_volume_L2"), 0.10758432876230622, "ST006 momentum L2")
    _require_equal(benchmark.get("retained_pde_validated"), False, "ST006 PDE state")
    _require_equal(benchmark.get("formal_momentum_target"), 0.001, "formal momentum target")
    _require_equal(
        benchmark.get("direct_superiority_claim_requires_identical_protocol"),
        True,
        "benchmark comparability gate",
    )
    _require_equal(
        set(benchmark.get("required_protocol_identity_fields", ())),
        _REQUIRED_PROTOCOL_FIELDS,
        "protocol identity fields",
    )
    _require_equal(benchmark.get("historical_or_cross_route_mismatch_role"), "indicative_only", "mismatch evidence role")
    _require_equal(benchmark.get("unknown_protocol_field_fails_closed"), True, "unknown protocol field policy")

    validation = constraints.get("validation", {})
    forcing = constraints.get("forcing", {})
    domain = constraints.get("domain", {})
    nontriviality = constraints.get("nontriviality", {})
    thresholds = validation.get("thresholds", {})
    canonical = contract.get("canonical_cr001_comparability", {})
    _require_equal(canonical.get("validation_seed"), validation.get("seed"), "CR001 validation seed")
    _require_equal(canonical.get("held_out_points"), validation.get("held_out_points"), "CR001 held-out count")
    _require_equal(canonical.get("times"), validation.get("times"), "CR001 validation times")
    _require_equal(canonical.get("derivative_steps"), validation.get("derivative_steps"), "CR001 derivative ladder")
    _require(validation.get("seed") != benchmark.get("retained_seed"), "ST006 and CR001 seeds unexpectedly match")
    _require_equal(canonical.get("seed_matches_st006_retained_protocol"), False, "seed mismatch declaration")
    _require_equal(canonical.get("directly_comparable_to_st006_without_replay"), False, "direct ST006 comparability")
    _require_equal(
        canonical.get("required_action_for_repository_level_superiority_claim"),
        "reevaluate both candidates under one identical frozen protocol before comparing",
        "repository-level comparison action",
    )

    lifecycle = contract.get("lifecycle_separation", {})
    _require_equal(set(lifecycle.get("allowed_states", ())), _ALLOWED_LIFECYCLE_STATES, "lifecycle state vocabulary")
    for key in (
        "open_or_draft_pr_changes_retained_baseline",
        "protocol_scoped_superiority_changes_retained_baseline",
        "protocol_scoped_superiority_changes_canonical_velocity",
        "retained_baseline_replacement_implies_pde_validated",
    ):
        _require_equal(lifecycle.get(key), False, key)
    for key in (
        "retained_baseline_replacement_requires_explicit_integration_or_promotion",
        "retained_baseline_replacement_requires_immutable_candidate_identity",
        "retained_baseline_replacement_requires_governed_validation_receipt",
        "retained_baseline_replacement_requires_status_document_update",
        "canonical_velocity_replacement_requires_explicit_api_or_candidate_status_update",
        "scientific_acceptance_requires_all_unchanged_formal_gates",
    ):
        _require_equal(lifecycle.get(key), True, key)

    challengers = contract.get("current_challenger_snapshot", {})
    for name in ("st047_e", "st048"):
        challenger = challengers.get(name, {})
        _require_equal(challenger.get("observed_state"), "draft_open_unintegrated", f"{name} observed state")
        _require_equal(challenger.get("below_st006_numbers_on_own_fresh_holdouts"), True, f"{name} numerical relation")
        _require_equal(challenger.get("same_protocol_st006_superiority_established_here"), False, f"{name} direct ST006 comparison")
        _require_equal(challenger.get("retained_baseline_replaced"), False, f"{name} retained promotion")
        _require_equal(challenger.get("canonical_velocity_replaced"), False, f"{name} delivery promotion")
        _require_equal(challenger.get("pde_validated"), False, f"{name} PDE state")
        _require(
            0.001 < float(challenger.get("best_reported_fresh_sampled_max")) < float(benchmark["retained_momentum_sampled_max"]),
            f"{name} sampled-max snapshot no longer has the governed challenger ordering",
        )
        _require(
            0.001 < float(challenger.get("best_reported_fresh_volume_L2")) < float(benchmark["retained_momentum_volume_L2"]),
            f"{name} L2 snapshot no longer has the governed challenger ordering",
        )

    snap = contract.get("canonical_thresholds", {})
    expected_pairs = (
        (snap.get("nu"), constraints.get("nu"), "nu"),
        (snap.get("physical_domain"), domain.get("physical"), "physical domain"),
        (snap.get("support"), domain.get("support"), "support"),
        (snap.get("time_interval"), domain.get("time_interval"), "time interval"),
        (snap.get("force_mode"), forcing.get("mode"), "force mode"),
        (snap.get("force_bounds"), forcing.get("parameters"), "force bounds"),
        (snap.get("reference_energy"), nontriviality.get("reference_energy"), "reference energy"),
        (
            snap.get("reference_energy_abs_tolerance"),
            nontriviality.get("reference_energy_abs_tolerance"),
            "reference energy tolerance",
        ),
        (snap.get("divergence_max"), thresholds.get("divergence_max"), "divergence max"),
        (snap.get("divergence_L2"), thresholds.get("divergence_L2"), "divergence L2"),
        (snap.get("pde_residual_max"), thresholds.get("pde_residual_max"), "PDE max"),
        (snap.get("pde_residual_L2"), thresholds.get("pde_residual_L2"), "PDE L2"),
    )
    for expected, actual, label in expected_pairs:
        _require_equal(actual, expected, label)
    _require("No residual-dependent basis" in forcing.get("restriction", ""), "free-force restriction weakened")

    delivery = contract.get("delivery_separation", {})
    _require_equal(delivery.get("canonical_delivery_family"), project_status.get("candidate_family"), "delivery family")
    _require_equal(delivery.get("canonical_delivery_sha256"), project_status.get("candidate_sha256"), "delivery candidate sha")
    _require_equal(delivery.get("canonical_delivery_sha256"), _DELIVERY_SHA, "known delivery candidate sha")
    _require_equal(delivery.get("velocity_export_ready"), True, "delivery export state")
    _require_equal(delivery.get("delivery_candidate_is_same_object_as_retained_pde_benchmark"), False, "delivery/benchmark identity")
    _require_equal(delivery.get("benchmark_challenge_changes_velocity_api_automatically"), False, "benchmark/delivery separation")
    _require_equal(delivery.get("pde_failure_blocks_callable_velocity_delivery"), False, "delivery/PDE separation")

    truth = contract.get("scientific_truth_boundary", {})
    _require_equal(truth.get("velocity_export_ready_is_independent_of_pde_validation"), True, "export/PDE independence")
    for key in (
        "paired_improvement_implies_pde_validated",
        "repository_baseline_improvement_implies_pde_validated",
        "retained_baseline_replacement_implies_pde_validated",
        "green_ci_implies_pde_validated",
        "visual_similarity_implies_pde_validated",
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(truth.get(key), False, key)

    states = project_status.get("states", {})
    _require_equal(states.get("velocity_export_ready"), True, "current velocity export state")
    for key in (
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(states.get(key), False, f"project state {key}")

    mutation = contract.get("mutation_scope", {})
    for key in (
        "velocity_changed",
        "candidate_changed",
        "pressure_changed",
        "forcing_changed",
        "optimizer_changed",
        "sampling_changed",
        "norm_definition_changed",
        "threshold_changed",
        "active_route_changed",
    ):
        _require_equal(mutation.get(key), False, key)

    return {
        "task_id": TASK_ID,
        "status": "governance_pass",
        "retained_baseline": "ST006",
        "retained_baseline_pde_validated": False,
        "current_challengers_are_unintegrated": True,
        "direct_st006_comparison_requires_replay": True,
        "parent_paired_claim_is_parent_local_only": True,
        "canonical_delivery_is_separate_from_pde_benchmark": True,
        "canonical_thresholds_unchanged": True,
        "velocity_export_ready": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
    }


def main() -> None:
    print(json.dumps(audit_candidate_promotion_evidence(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
