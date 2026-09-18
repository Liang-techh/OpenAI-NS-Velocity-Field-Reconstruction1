"""Fail-closed CR002 governance for candidate-promotion evidence semantics.

This module changes no velocity, candidate, pressure, forcing, optimizer, sample,
norm, or threshold.  It keeps parent-local improvement, repository-baseline
comparison, and scientific-state promotion as separate claims.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-CANDIDATE-PROMOTION-EVIDENCE-050"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
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
    """Audit evidence comparability and claim-promotion boundaries."""

    contract = deepcopy(dict(load_contract() if contract is None else contract))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(contract.get("task_id"), TASK_ID, "contract task id")
    _require_equal(set(contract.get("allowed_source_classes", ())), _ALLOWED_CLASSES, "source-class vocabulary")
    source = contract.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    _require_equal(source.get("callable_velocity_delivery"), "user_requirement", "velocity delivery class")
    _require_equal(source.get("public_openai_visualization_structure"), "public_source_fact", "public visualization class")
    _require_equal(source.get("cr001_validation_protocol"), "autonomous_design", "CR001 protocol class")
    _require_equal(source.get("st006_validation_protocol"), "autonomous_design", "ST006 protocol class")
    _require_equal(source.get("st006_measured_residual_record"), "autonomous_design", "ST006 measurement class")
    _require_equal(source.get("paired_parent_child_measurements"), "autonomous_design", "paired measurement class")
    _require_equal(source.get("openai_hidden_numerical_profile"), "pending_unknown", "hidden profile class")

    labels = contract.get("evidence_relation_labels", {})
    _require("project_derived_measurement" in labels, "project measurement provenance label missing")
    _require("frozen_parent_paired_comparison" in labels, "paired-comparison relation label missing")
    _require("repository_retained_benchmark" in labels, "repository benchmark relation label missing")
    _require_equal(contract.get("relation_labels_are_source_classes"), False, "relation/source-class separation")

    paired = contract.get("paired_parent_claim_scope", {})
    _require_equal(
        paired.get("allowed_claim"),
        "improved over this frozen parent under this exact declared protocol",
        "paired allowed claim",
    )
    for key in (
        "allows_repository_st006_superiority_claim",
        "allows_pde_validation_claim",
        "allows_visual_correspondence_claim",
        "allows_paper_exact_claim",
        "allows_openai_field_identification_claim",
    ):
        _require_equal(paired.get(key), False, key)

    benchmark = contract.get("repository_baseline_promotion", {})
    _require_equal(benchmark.get("benchmark_name"), "ST006", "benchmark name")
    _require_equal(
        benchmark.get("candidate_sha256"),
        "6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3",
        "ST006 candidate sha",
    )
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
    _require_equal(benchmark.get("direct_superiority_claim_requires_identical_protocol"), True, "benchmark comparability gate")
    _require_equal(set(benchmark.get("required_protocol_identity_fields", ())), _REQUIRED_PROTOCOL_FIELDS, "protocol identity fields")
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

    truth = contract.get("scientific_truth_boundary", {})
    _require_equal(truth.get("velocity_export_ready_is_independent_of_pde_validation"), True, "export/PDE independence")
    for key in (
        "paired_improvement_implies_pde_validated",
        "repository_baseline_improvement_implies_pde_validated",
        "green_ci_implies_pde_validated",
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
    ):
        _require_equal(mutation.get(key), False, key)

    return {
        "task_id": TASK_ID,
        "status": "governance_pass",
        "st006_seed": benchmark["retained_seed"],
        "cr001_seed": validation["seed"],
        "direct_st006_comparison_requires_replay": True,
        "parent_paired_claim_is_parent_local_only": True,
        "canonical_thresholds_unchanged": True,
        "velocity_export_ready": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
    }


def main() -> None:
    print(json.dumps(audit_candidate_promotion_evidence(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
