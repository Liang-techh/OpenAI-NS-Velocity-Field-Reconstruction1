"""Fail-closed governance for the bipolar F(2,0) radial-swirl capacity screen.

Agent 7's F(2,0) screen measures representation capacity around the frozen capped
bipolar field.  It does not materialize or select a nonzero F(2,0) candidate and
it does not evaluate held-out Navier--Stokes residuals.  In particular, a
swirl-pure *velocity* response is not a statement that only the theta momentum
equation changes: nonlinear momentum and kinetic energy can change when swirl is
changed.

This module audits that evidence/claim boundary and the unchanged CR001 contract.
It changes no velocity, coefficient, pressure, force, optimizer, or threshold.
"""
from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .constrained_eq45_bipolar_f20_capacity import audit_bipolar_f20_capacity

_CANONICAL_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_REQUIRED_FALSE_REPORT_TRUTH = (
    "velocity_changed",
    "candidate_artifact_changed",
    "force_or_pressure_fitted",
    "held_out_pde_residual_evaluated",
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
)
_REQUIRED_FALSE_MATERIALIZATION = (
    "f20_nonzero_candidate_materialized",
    "f20_nonzero_velocity_export_ready",
    "candidate_selection_resolved",
    "pde_validated",
    "visual_correspondence_verified",
    "paper_exact",
    "openai_field_identified",
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


def audit_bipolar_f20_capacity_scope(
    scope: Mapping[str, Any],
    constraints: Mapping[str, Any],
    delivery_contract: Mapping[str, Any],
    source_candidate: Mapping[str, Any],
    capacity_report: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit F20 capacity evidence without promoting a candidate or PDE claim."""
    _require(scope.get("schema") == "bipolar_f20_capacity_scope_v1", "scope schema drift")
    _require(scope.get("task_id") == "CR002-BIPOLAR-F20-CAPACITY-SCOPE-030", "task id drift")

    vocabulary = delivery_contract.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "missing source vocabulary")
    _require(set(vocabulary) == _CANONICAL_CLASSES, "canonical source vocabulary drift")
    information_source = scope.get("information_source")
    _require(isinstance(information_source, Mapping), "missing information-source map")
    _require(set(information_source.values()) <= _CANONICAL_CLASSES, "noncanonical information source class")
    _require(information_source.get("f20_basis_mode") == "autonomous_design", "F20 basis must remain autonomous")
    _require(information_source.get("capacity_probe_layout") == "autonomous_design", "probe layout provenance drift")
    _require(information_source.get("capacity_thresholds") == "autonomous_design", "capacity thresholds provenance drift")
    _require(information_source.get("hidden_openai_f20_coefficient") == "pending_unknown", "hidden coefficient cannot be invented")

    upstream = scope.get("upstream")
    _require(isinstance(upstream, Mapping), "missing upstream binding")
    _require(upstream.get("task_id") == "CR003-BIPOLAR-F20-RADIAL-SWIRL-CAPACITY-035", "upstream task drift")
    _require(upstream.get("pr") == 212, "upstream PR drift")
    _require(upstream.get("head") == "33f710de04e7c0d36db666e3863efcaaa711aa3e", "upstream head drift")
    _require(capacity_report.get("task_id") == upstream.get("task_id"), "capacity report task drift")
    _require(
        capacity_report.get("source_candidate_sha256") == upstream.get("source_candidate_sha256"),
        "capacity source candidate identity drift",
    )

    representation = scope.get("representation")
    _require(isinstance(representation, Mapping), "missing representation contract")
    basis = source_candidate.get("parent_candidate", {}).get("profile_basis", {})
    _require(isinstance(basis, Mapping), "missing source profile basis")
    _require(basis.get("radial_degree") == representation.get("source_radial_degree") == 1, "source radial degree drift")
    _require(representation.get("embedded_radial_degree") == 2, "embedded radial degree drift")
    _require(basis.get("eta_degree") == representation.get("eta_degree") == 3, "eta degree drift")
    _require(basis.get("coefficient_limit") == representation.get("coefficient_limit") == 4.0, "coefficient limit drift")
    modes = [tuple(mode) for mode in basis.get("mode_indices", [])]
    _require(tuple(representation.get("existing_reference_swirl_mode", ())) == (1, 0), "reference mode drift")
    _require(tuple(representation.get("new_swirl_mode", ())) == (2, 0), "F20 mode drift")
    _require((1, 0) in modes, "source candidate lost F10")
    _require((2, 0) not in modes, "source candidate already contains F20; capacity semantics changed")
    _require(representation.get("f20_coefficient_during_capacity_screen") == 0.0, "capacity screen selected F20")
    _require(representation.get("zero_embedding_preserves_source_velocity") is True, "zero embedding contract drift")
    growth = representation.get("parameter_growth_if_materialized")
    _require(growth == {"spatial_modes": 1, "scalar_coefficients": 1}, "parameter growth drift")
    _require(
        capacity_report.get("parameter_growth")
        == {"spatial_modes_added": 1, "scalar_coefficients_added_if_selected": 1},
        "capacity report parameter-growth drift",
    )

    evidence = scope.get("capacity_evidence_contract")
    _require(isinstance(evidence, Mapping), "missing capacity evidence contract")
    _require(evidence.get("kind") == "target_free_representation_capacity_only", "capacity kind drift")
    _require(
        float(capacity_report.get("embedding_max_abs_velocity_error", float("inf")))
        < float(evidence.get("embedding_max_abs_velocity_error_max")),
        "zero-mode embedding no longer replays the source velocity",
    )
    _require(capacity_report.get("response_rank") == evidence.get("response_rank") == 2, "response rank drift")
    _require(
        float(capacity_report.get("response_condition_number", float("inf")))
        < float(evidence.get("response_condition_number_max")),
        "response conditioning drift",
    )
    _require(
        float(capacity_report.get("finite_difference_refinement_relative_change", float("inf")))
        < float(evidence.get("finite_difference_refinement_relative_change_max")),
        "capacity finite-difference refinement drift",
    )
    _require(
        float(capacity_report.get("F20_novelty_outside_F10_span", -1.0))
        > float(evidence.get("f20_novelty_outside_f10_span_min")),
        "F20 novelty evidence drift",
    )
    loc = capacity_report.get("localization_at_t_050")
    _require(isinstance(loc, Mapping), "missing localization evidence")
    for name in ("F10", "F20"):
        entry = loc.get(name)
        _require(isinstance(entry, Mapping), f"missing {name} localization evidence")
        _require(
            float(entry.get("swirl_fraction", -1.0)) > float(evidence.get("swirl_fraction_min")),
            f"{name} velocity response no longer swirl-pure",
        )
    _require(
        float(capacity_report.get("F20_vs_F10_outer_selectivity_ratio", -1.0))
        > float(evidence.get("f20_vs_f10_outer_selectivity_ratio_min")),
        "F20 radial-selectivity evidence drift",
    )
    for key in ("held_out_pde_residual_evaluated", "public_reference_used", "coefficient_selected"):
        _require(evidence.get(key) is False, f"capacity scope cannot promote {key}")

    report_truth = capacity_report.get("truth_boundary")
    _require(isinstance(report_truth, Mapping), "missing capacity truth boundary")
    for key in _REQUIRED_FALSE_REPORT_TRUTH:
        _require(report_truth.get(key) is False, f"capacity report truth promotion: {key}")

    forbidden = scope.get("forbidden_inferences")
    _require(isinstance(forbidden, Mapping), "missing forbidden-inference map")
    _require(forbidden and all(value is False for value in forbidden.values()), "forbidden inference was enabled")
    for required_key in (
        "swirl_pure_velocity_response_implies_theta_only_momentum_change",
        "swirl_pure_velocity_response_implies_energy_unchanged",
        "radial_selectivity_implies_pde_improvement",
        "capacity_rank_implies_nonzero_f20_should_be_selected",
        "zero_embedding_implies_nonzero_f20_candidate_validated",
    ):
        _require(required_key in forbidden, f"missing forbidden inference: {required_key}")

    materialization = scope.get("materialization_state")
    _require(isinstance(materialization, Mapping), "missing F20 materialization state")
    _require(materialization.get("source_candidate_velocity_export_ready") is True, "source export state drift")
    _require(materialization.get("f20_nonzero_candidate_sha256") is None, "capacity screen invented an F20 candidate SHA")
    for key in _REQUIRED_FALSE_MATERIALIZATION:
        _require(materialization.get(key) is False, f"premature F20 materialization/promotion: {key}")

    source_truth = source_candidate.get("truth_boundary")
    _require(isinstance(source_truth, Mapping), "missing source candidate truth boundary")
    _require(source_truth.get("velocity_export_ready") is True, "source candidate exportability drift")
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _require(source_truth.get(key) is False, f"source candidate truth promotion: {key}")

    requirements = scope.get("requirements_before_nonzero_f20_promotion")
    _require(isinstance(requirements, Mapping), "missing nonzero-F20 promotion requirements")
    _require(requirements and all(value is True for value in requirements.values()), "nonzero-F20 promotion requirement weakened")

    expected = scope.get("cr001_invariants")
    _require(isinstance(expected, Mapping), "missing CR001 invariant snapshot")
    _require(constraints.get("nu") == expected.get("nu") == 0.01, "nu drift")
    domain = constraints.get("domain")
    _require(isinstance(domain, Mapping), "missing domain contract")
    _require(domain.get("physical") == expected.get("physical_domain") == "R^3", "physical domain drift")
    _require(domain.get("evaluation_box") == expected.get("evaluation_box"), "evaluation box drift")
    _require(domain.get("support") == expected.get("support"), "support drift")
    _require(domain.get("time_interval") == expected.get("time_interval"), "time interval drift")

    forcing = constraints.get("forcing")
    _require(isinstance(forcing, Mapping), "missing forcing contract")
    _require(forcing.get("mode") == expected.get("forcing_mode") == "restricted_two_parameter_family", "forcing mode drift")
    _require(forcing.get("parameters") == expected.get("forcing_bounds"), "forcing bounds drift")
    restriction = forcing.get("restriction")
    _require(
        isinstance(restriction, str)
        and "Only a,c may be fitted" in restriction
        and "No residual-dependent basis or pointwise free force" in restriction,
        "free-force prohibition drift",
    )

    nontriviality = constraints.get("nontriviality")
    _require(isinstance(nontriviality, Mapping), "missing nontriviality contract")
    _require(nontriviality.get("reference_energy") == expected.get("reference_energy") == 1.0, "reference energy drift")
    _require(
        nontriviality.get("reference_energy_abs_tolerance")
        == expected.get("reference_energy_abs_tolerance")
        == 0.001,
        "reference energy tolerance drift",
    )

    validation = constraints.get("validation")
    _require(isinstance(validation, Mapping), "missing validation contract")
    _require(validation.get("held_out_points") == expected.get("held_out_points") == 4096, "held-out count drift")
    _require(validation.get("derivative_steps") == expected.get("derivative_steps"), "derivative ladder drift")
    thresholds = validation.get("thresholds")
    _require(isinstance(thresholds, Mapping), "missing validation thresholds")
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds.get(key) == expected.get(key), f"registered threshold drift: {key}")

    return {
        "scope_pass": True,
        "source_candidate_sha256": capacity_report.get("source_candidate_sha256"),
        "f20_mode": [2, 0],
        "f20_source_classification": "autonomous_design",
        "f20_nonzero_candidate_materialized": False,
        "candidate_selection_resolved": False,
        "held_out_pde_residual_evaluated": False,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def audit_bipolar_f20_capacity_scope_file(
    scope_path: str | Path = "configs/bipolar_f20_capacity_scope.json",
    *,
    repo_root: str | Path = ".",
) -> dict[str, Any]:
    """Load governed files, rerun the target-free capacity report, and audit."""
    root = Path(repo_root)
    scope = _read_json(root / Path(scope_path))
    inputs = scope.get("inputs")
    _require(isinstance(inputs, Mapping), "missing scope inputs")
    constraints = _read_json(root / str(inputs.get("constraints")))
    delivery = _read_json(root / str(inputs.get("delivery_state_contract")))
    candidate = _read_json(root / str(inputs.get("source_candidate")))
    report = audit_bipolar_f20_capacity()
    return audit_bipolar_f20_capacity_scope(scope, constraints, delivery, candidate, report)
