"""Fail-closed governance for residual-map-guided spatial-swirl routing.

A measured residual/error map may be useful for deciding which *predeclared*
representation family deserves the next capacity screen.  Once that same map or
its samples influence the location, width, support, basis shape, or coefficient
of a correction, however, they are model-selection data and cannot also be
presented as independent PDE acceptance evidence.

This module audits that boundary for the capped bipolar F(2,0) overlap result.
It changes no velocity, coefficient, pressure, force, optimizer, sample, or
threshold.
"""
from __future__ import annotations

import json
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any

_CANONICAL_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_EXPECTED_SOURCE_SHA = "7875214ad65f5e1ba4f7c5e362217f05ef47bc3d893e8641e8e40a85fc4e7609"
_EXPECTED_UPSTREAM_HEAD = "7d742148415a5333eb016496996f281c58186544"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _close(actual: Any, expected: float, *, name: str, atol: float = 1e-12) -> None:
    _require(isinstance(actual, (int, float)) and math.isfinite(float(actual)), f"invalid numeric evidence: {name}")
    _require(math.isclose(float(actual), expected, rel_tol=1e-10, abs_tol=atol), f"routing evidence drift: {name}")


def audit_bipolar_error_map_routing_scope(
    scope: Mapping[str, Any],
    constraints: Mapping[str, Any],
    delivery_contract: Mapping[str, Any],
    source_candidate: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit error-map routing without laundering model-selection data."""
    _require(scope.get("schema") == "bipolar_error_map_routing_scope_v1", "scope schema drift")
    _require(scope.get("task_id") == "CR002-BIPOLAR-ERROR-MAP-ROUTING-SCOPE-031", "task id drift")

    vocabulary = delivery_contract.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "missing canonical source vocabulary")
    _require(set(vocabulary) == _CANONICAL_CLASSES, "canonical source vocabulary drift")
    _require(set(scope.get("classification_vocabulary", ())) == _CANONICAL_CLASSES, "scope source vocabulary drift")

    source = scope.get("source_candidate")
    _require(isinstance(source, Mapping), "missing source-candidate binding")
    _require(source.get("artifact") == "artifacts/bipolar_joint_capped/candidate.json", "source artifact drift")
    _require(source.get("candidate_sha256") == _EXPECTED_SOURCE_SHA, "source candidate SHA drift")
    _require(source.get("classification") == "autonomous_design", "source candidate classification drift")
    _require(source.get("velocity_export_ready") is True, "source candidate export state drift")
    _require(source.get("pde_validated") is False, "source candidate PDE promotion")

    source_truth = source_candidate.get("truth_boundary")
    _require(isinstance(source_truth, Mapping), "missing source candidate truth boundary")
    _require(source_truth.get("velocity_export_ready") is True, "source artifact exportability drift")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(source_truth.get(key) is False, f"source artifact truth promotion: {key}")

    upstream = scope.get("upstream_routing_evidence")
    _require(isinstance(upstream, Mapping), "missing upstream routing evidence")
    _require(upstream.get("pr") == 220, "upstream PR drift")
    _require(upstream.get("head") == _EXPECTED_UPSTREAM_HEAD, "upstream head drift")
    _require(upstream.get("classification") == "autonomous_design", "routing evidence provenance drift")
    _require(upstream.get("evidence_role") == "routing_diagnostic_only", "routing evidence role drift")
    _close(upstream.get("time"), 0.5, name="routing time")

    baseline = upstream.get("baseline_theta_error_rms")
    _require(isinstance(baseline, Mapping), "missing baseline theta-error map receipt")
    _close(baseline.get("interior"), 0.7362262073, name="interior theta RMS")
    _close(baseline.get("radial_support_flank"), 0.0789660930, name="flank theta RMS")
    _close(baseline.get("flank_over_interior"), 0.107257922, name="flank/interior ratio")
    _require(float(baseline.get("interior")) > float(baseline.get("radial_support_flank")), "theta error no longer interior dominated")

    overlap = upstream.get("f20_vs_f10")
    _require(isinstance(overlap, Mapping), "missing F20/F10 overlap receipt")
    _require(overlap.get("both_velocity_responses_swirl_pure") is True, "swirl-response receipt drift")
    _close(overlap.get("f20_over_f10_flank_selectivity"), 2.74087783, name="F20/F10 flank selectivity")
    _close(overlap.get("global_error_enrichment_ratio_f20_over_f10"), 0.843068186, name="F20/F10 error enrichment")
    _close(overlap.get("response_error_cosine_f10"), 0.643200558, name="F10 error cosine")
    _close(overlap.get("response_error_cosine_f20"), 0.497506209, name="F20 error cosine")
    _close(overlap.get("top_quartile_error_response_energy_share_f10"), 0.1216896, name="F10 top-quartile share")
    _close(overlap.get("top_quartile_error_response_energy_share_f20"), 0.0601040, name="F20 top-quartile share")
    _require(float(overlap.get("f20_over_f10_flank_selectivity")) > 1.0, "F20 lost outer-selectivity evidence")
    _require(float(overlap.get("global_error_enrichment_ratio_f20_over_f10")) < 1.0, "F20 no longer has the governed weaker global alignment")
    _require(float(overlap.get("response_error_cosine_f20")) < float(overlap.get("response_error_cosine_f10")), "F20/F10 alignment ordering drift")
    _require(
        float(overlap.get("top_quartile_error_response_energy_share_f20"))
        < float(overlap.get("top_quartile_error_response_energy_share_f10")),
        "F20/F10 top-error overlap ordering drift",
    )

    allowed = upstream.get("allowed_conclusion")
    _require(isinstance(allowed, str) and "route the next predeclared capacity screen" in allowed, "allowed routing conclusion drift")
    forbidden = upstream.get("forbidden_conclusions")
    _require(isinstance(forbidden, list) and len(forbidden) >= 5, "missing forbidden routing conclusions")
    forbidden_text = "\n".join(str(item) for item in forbidden)
    for required in (
        "F20 coefficient selected",
        "F20 finite perturbation reduces PDE residual",
        "interior-localized correction parameters have been identified",
        "OpenAI hidden swirl profile recovered",
    ):
        _require(required in forbidden_text, f"missing forbidden conclusion: {required}")

    boundary = scope.get("model_selection_boundary")
    _require(isinstance(boundary, Mapping), "missing model-selection boundary")
    _require(boundary.get("residual_map_may_route_predeclared_family") is True, "routing use was disabled")
    _require(boundary.get("residual_map_may_fit_spatial_envelope") is False, "residual map may not fit correction envelope")
    _require(boundary.get("predeclared_interior_geometry_source") == "existing plateau/support geometry only", "interior geometry provenance drift")
    for key in (
        "forbid_residual_fitted_center",
        "forbid_residual_fitted_width",
        "forbid_residual_fitted_support",
        "forbid_residual_fitted_basis_shape",
    ):
        _require(boundary.get(key) is True, f"residual-fit leakage enabled: {key}")
    _require(boundary.get("if_samples_influence_mode_location_shape_or_coefficient") == "consumed_for_model_selection", "sample-consumption rule drift")
    _require(boundary.get("consumed_samples_may_count_as_independent_validation") is False, "model-selection samples laundered into validation")
    _require(boundary.get("required_after_candidate_freeze") == "fresh_or_pristine_independent_validation", "fresh validation requirement drift")

    registered = scope.get("registered_contract")
    _require(isinstance(registered, Mapping), "missing registered contract snapshot")
    _require(constraints.get("nu") == registered.get("nu") == 0.01, "nu drift")
    domain = constraints.get("domain")
    _require(isinstance(domain, Mapping), "missing domain contract")
    _require(domain.get("physical") == registered.get("physical_domain") == "R^3", "physical domain drift")
    _require(domain.get("evaluation_box") == registered.get("evaluation_box"), "evaluation box drift")
    _require(domain.get("support") == registered.get("support"), "support drift")
    _require(domain.get("time_interval") == registered.get("time_interval"), "time interval drift")

    forcing = constraints.get("forcing")
    _require(isinstance(forcing, Mapping), "missing forcing contract")
    _require(forcing.get("mode") == registered.get("forcing_mode") == "restricted_two_parameter_family", "forcing mode drift")
    _require(forcing.get("parameters") == registered.get("forcing_parameters"), "forcing bounds drift")
    _require(registered.get("residual_defined_free_force_allowed") is False, "free-force route enabled")
    restriction = forcing.get("restriction")
    _require(isinstance(restriction, str) and "No residual-dependent basis or pointwise free force" in restriction, "free-force prohibition drift")

    nontriviality = constraints.get("nontriviality")
    _require(isinstance(nontriviality, Mapping), "missing nontriviality contract")
    _require(nontriviality.get("reference_energy") == registered.get("reference_energy") == 1.0, "reference energy drift")
    _require(
        nontriviality.get("reference_energy_abs_tolerance")
        == registered.get("reference_energy_abs_tolerance")
        == 0.001,
        "reference energy tolerance drift",
    )

    validation = constraints.get("validation")
    _require(isinstance(validation, Mapping), "missing validation contract")
    _require(validation.get("seed") == registered.get("validation_seed") == 914027, "validation seed drift")
    _require(validation.get("held_out_points") == registered.get("held_out_points") == 4096, "held-out count drift")
    _require(validation.get("derivative_steps") == registered.get("derivative_steps"), "derivative ladder drift")
    thresholds = validation.get("thresholds")
    _require(isinstance(thresholds, Mapping), "missing registered thresholds")
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds.get(key) == registered.get(key), f"registered threshold drift: {key}")

    requirements = scope.get("promotion_requirements_for_any_nonzero_spatial_correction")
    _require(isinstance(requirements, list), "missing correction promotion requirements")
    required_phrases = (
        "explicit bounded coefficient(s)",
        "serialized candidate identity",
        "initial energy revalidation",
        "nonzero center rotation and sign revalidation",
        "bipolar radial/axial direction and parity revalidation",
        "axis/support regularity revalidation",
        "fresh_or_pristine full momentum and divergence validation",
        "training/model-selection samples identified separately from acceptance samples",
    )
    for phrase in required_phrases:
        _require(phrase in requirements, f"missing promotion requirement: {phrase}")

    states = scope.get("states")
    _require(isinstance(states, Mapping), "missing state boundary")
    _require(states.get("velocity_export_ready_for_source_candidate") is True, "source delivery state drift")
    for key in (
        "f20_selected",
        "interior_localized_mode_materialized",
        "candidate_selection_resolved",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(states.get(key) is False, f"premature state promotion: {key}")

    return {
        "scope_pass": True,
        "source_candidate_sha256": _EXPECTED_SOURCE_SHA,
        "upstream_pr": 220,
        "error_map_role": "routing_diagnostic_only",
        "f20_selected": False,
        "interior_localized_mode_materialized": False,
        "model_selection_samples_reusable_as_independent_validation": False,
        "fresh_or_pristine_validation_required_after_freeze": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def audit_bipolar_error_map_routing_scope_file(
    scope_path: str | Path = "configs/bipolar_error_map_routing_scope.json",
    *,
    repo_root: str | Path = ".",
) -> dict[str, Any]:
    """Load the governed repository files and audit the routing boundary."""
    root = Path(repo_root)
    scope = _read_json(root / Path(scope_path))
    constraints = _read_json(root / "configs/constraints.json")
    delivery = _read_json(root / "configs/delivery_state_contract.json")
    source = _read_json(root / str(scope["source_candidate"]["artifact"]))
    return audit_bipolar_error_map_routing_scope(scope, constraints, delivery, source)
