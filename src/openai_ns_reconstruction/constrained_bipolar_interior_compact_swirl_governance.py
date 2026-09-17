"""Fail-closed governance for the interior compact swirl capacity screen.

The integrated INTERIOR_C4_SWIRL direction is an autonomous additive physical-
space correction.  It is not an existing Eq45 profile-basis coefficient, and a
nonzero materialization must therefore receive its own representation identity
and candidate SHA before scientific evidence can be attached to it.

This module changes no velocity value, coefficient, force, pressure, sample, or
acceptance threshold.
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
_EXPECTED_UPSTREAM_HEAD = "009185a05133359a089cb7489b67a32a812883ed"


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
    _require(math.isclose(float(actual), expected, rel_tol=1e-9, abs_tol=atol), f"capacity evidence drift: {name}")


def audit_bipolar_interior_compact_swirl_scope(
    scope: Mapping[str, Any],
    constraints: Mapping[str, Any],
    delivery_contract: Mapping[str, Any],
    source_candidate: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit representation, diagnostic-trial, and promotion semantics."""
    _require(scope.get("schema") == "bipolar_interior_compact_swirl_capacity_scope_v1", "scope schema drift")
    _require(scope.get("task_id") == "CR002-BIPOLAR-INTERIOR-COMPACT-SWIRL-SCOPE-032", "task id drift")

    vocabulary = delivery_contract.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "missing canonical source vocabulary")
    _require(set(vocabulary) == _CANONICAL_CLASSES, "canonical source vocabulary drift")
    _require(set(scope.get("classification_vocabulary", ())) == _CANONICAL_CLASSES, "scope source vocabulary drift")

    source = scope.get("source_candidate")
    _require(isinstance(source, Mapping), "missing source candidate binding")
    _require(source.get("artifact") == "artifacts/bipolar_joint_capped/candidate.json", "source artifact drift")
    _require(source.get("candidate_sha256") == _EXPECTED_SOURCE_SHA, "source candidate SHA drift")
    _require(source.get("classification") == "autonomous_design", "source classification drift")
    _require(source.get("velocity_export_ready") is True, "source export state drift")
    _require(source.get("pde_validated") is False, "source PDE promotion")

    source_truth = source_candidate.get("truth_boundary")
    _require(isinstance(source_truth, Mapping), "missing source truth boundary")
    _require(source_truth.get("velocity_export_ready") is True, "source artifact export state drift")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(source_truth.get(key) is False, f"source truth promotion: {key}")

    taper = source_candidate.get("physical_taper")
    _require(isinstance(taper, Mapping), "missing physical taper")
    radial_plateau = float(taper.get("radial_support")) * math.sqrt(float(taper.get("radial_plateau_q")))
    axial_plateau = float(taper.get("axial_half_height")) * math.sqrt(float(taper.get("axial_plateau_q")))
    _close(radial_plateau, 1.6, name="source radial plateau half width")
    _close(axial_plateau, 1.6, name="source axial plateau half width")
    parent = source_candidate.get("parent_candidate")
    _require(isinstance(parent, Mapping), "missing source parent candidate")
    profile_basis = parent.get("profile_basis")
    _require(isinstance(profile_basis, Mapping), "missing source profile basis")
    _close(profile_basis.get("coefficient_limit"), 4.0, name="inherited profile coefficient limit")

    upstream = scope.get("upstream_capacity_evidence")
    _require(isinstance(upstream, Mapping), "missing upstream capacity evidence")
    _require(upstream.get("pr") == 228, "upstream PR drift")
    _require(upstream.get("head") == _EXPECTED_UPSTREAM_HEAD, "upstream head drift")
    _require(upstream.get("task_id") == "CR003-BIPOLAR-INTERIOR-COMPACT-SWIRL-CAPACITY-037", "upstream task drift")
    _require(upstream.get("classification") == "autonomous_design", "capacity evidence provenance drift")
    _require(upstream.get("evidence_role") == "target_free_representation_capacity_only", "capacity evidence role drift")
    _require(upstream.get("mode") == "INTERIOR_C4_SWIRL", "mode identity drift")
    _require(upstream.get("shape_parameters_fitted") == 0, "capacity shape was fitted")
    _require(upstream.get("residual_map_used_to_fit_basis_shape") is False, "residual-fitted basis shape enabled")
    _close(upstream.get("time"), 0.5, name="capacity time")
    _close(upstream.get("radial_plateau_half_width"), 1.6, name="capacity radial plateau")
    _close(upstream.get("axial_plateau_half_width"), 1.6, name="capacity axial plateau")
    _require(upstream.get("normalized_response_rank") == 3, "capacity rank drift")
    _close(upstream.get("normalized_response_condition_number"), 7.344633235, name="capacity condition")
    _close(upstream.get("compact_novelty_outside_F10_F20_span"), 0.667113975, name="capacity novelty")
    _close(upstream.get("coefficient_step_refinement_relative_change"), 5.77e-14, name="capacity refinement", atol=1e-15)
    _close(upstream.get("swirl_fraction"), 1.0, name="swirl fraction")
    _close(upstream.get("core_response_rms"), 0.141456627, name="core response RMS")
    _close(upstream.get("interior_response_rms"), 0.0647900904, name="interior response RMS")
    _close(upstream.get("radial_flank_response_rms"), 0.0, name="radial flank response RMS")
    _close(upstream.get("axial_flank_response_rms"), 0.0, name="axial flank response RMS")
    _close(upstream.get("compact_collar_response_energy_fraction"), 0.0, name="collar response energy")

    semantics = scope.get("representation_semantics")
    _require(isinstance(semantics, Mapping), "missing representation semantics")
    _require(semantics.get("mode_origin") == "autonomous_design", "mode origin drift")
    _require(
        semantics.get("geometry_source") == "existing AxisymmetricPhysicalTaper identity-plateau geometry",
        "geometry provenance drift",
    )
    _require(semantics.get("is_eq45_public_source_profile_mode") is False, "autonomous correction laundered as public Eq45 mode")
    _require(semantics.get("is_existing_eq45_profile_basis_coefficient") is False, "physical-space correction laundered as profile coefficient")
    _require(semantics.get("is_additive_physical_space_swirl_correction") is True, "physical-space correction semantics drift")
    _require(semantics.get("axis_regular_by_r_factor") is True, "axis-regular construction drift")
    _require(semantics.get("axisymmetric_pure_swirl_is_analytically_divergence_free") is True, "structural divergence identity drift")
    _require(semantics.get("zero_before_physical_taper_collars") is True, "collar-localization semantics drift")
    _require(semantics.get("nonzero_materialization_changes_representation_family_identity") is True, "representation identity change disabled")
    _require(semantics.get("parent_candidate_identity_may_be_reused_for_nonzero_child") is False, "parent identity reuse enabled")

    trial = scope.get("diagnostic_trial_semantics")
    _require(isinstance(trial, Mapping), "missing diagnostic-trial semantics")
    _require(trial.get("symmetric_trial_coefficients") == [-0.25, 0.25], "diagnostic trial grid drift")
    _require(trial.get("trial_role") == "capacity_and_target_free_morphology_probe_only", "diagnostic trial role drift")
    for key in (
        "trial_coefficients_are_selected",
        "trial_coefficients_define_materialization_bound",
        "visual_correspondence_evidence",
        "pde_evidence",
        "energy_acceptance_evidence",
    ):
        _require(trial.get(key) is False, f"diagnostic trial promotion: {key}")

    forbidden = scope.get("forbidden_inferences")
    _require(isinstance(forbidden, list) and len(forbidden) >= 10, "missing forbidden inferences")
    forbidden_text = "\n".join(str(item) for item in forbidden)
    for phrase in (
        "theta-only momentum impact",
        "initial energy is unchanged",
        "registered numerical divergence gate passed",
        "selected",
        "materialization bound",
        "public-source Eq45 profile coefficient",
        "parent candidate SHA or family identity",
        "public OpenAI correspondence",
        "PDE validation",
    ):
        _require(phrase in forbidden_text, f"missing forbidden inference: {phrase}")

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
    _require(isinstance(thresholds, Mapping), "missing validation thresholds")
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds.get(key) == registered.get(key), f"registered threshold drift: {key}")

    requirements = scope.get("promotion_requirements_for_nonzero_materialization")
    _require(isinstance(requirements, list), "missing materialization requirements")
    required_phrases = (
        "explicit coefficient value and explicit coefficient bound",
        "new representation-family label describing Eq45 parent plus autonomous physical-space swirl correction",
        "new serialized candidate SHA distinct from the unchanged parent",
        "initial energy E(0.25)=1 within registered tolerance revalidated",
        "minimum/maximum validation-time energy gates revalidated",
        "nonzero center rotation and required sign revalidated",
        "bipolar radial/axial direction and parity revalidated",
        "axis/support regularity revalidated",
        "full per-component momentum and divergence recomputed",
        "registered 4096-point fresh/pristine validation remains separate from model-selection evidence",
        "restricted forcing family and bounds unchanged unless a new preregistered experiment is declared",
    )
    for phrase in required_phrases:
        _require(phrase in requirements, f"missing materialization requirement: {phrase}")

    states = scope.get("states")
    _require(isinstance(states, Mapping), "missing state boundary")
    _require(states.get("source_candidate_velocity_export_ready") is True, "source export state drift")
    _require(states.get("interior_compact_mode_capacity_screened") is True, "capacity-screen state drift")
    for key in (
        "interior_compact_nonzero_candidate_materialized",
        "interior_compact_coefficient_selected",
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
        "upstream_pr": 228,
        "mode": "INTERIOR_C4_SWIRL",
        "mode_origin": "autonomous_design",
        "is_existing_eq45_profile_basis_coefficient": False,
        "diagnostic_trial_selected": False,
        "nonzero_materialization_requires_new_identity": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def audit_bipolar_interior_compact_swirl_scope_file(
    scope_path: str | Path = "configs/bipolar_interior_compact_swirl_capacity_scope.json",
    *,
    repo_root: str | Path = ".",
) -> dict[str, Any]:
    """Load governed repository files and audit the compact-swirl scope."""
    root = Path(repo_root)
    scope = _read_json(root / Path(scope_path))
    constraints = _read_json(root / "configs/constraints.json")
    delivery = _read_json(root / "configs/delivery_state_contract.json")
    source = _read_json(root / str(scope["source_candidate"]["artifact"]))
    return audit_bipolar_interior_compact_swirl_scope(scope, constraints, delivery, source)
