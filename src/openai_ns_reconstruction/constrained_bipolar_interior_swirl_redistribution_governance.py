"""Fail-closed governance for the interior swirl-redistribution capacity screen.

The Agent-7 redistribution direction is an autonomous, sign-changing physical-
space swirl basis.  Its analytic cylindrical L2 orthogonality is a property of
one raw basis pair; it is not a full-candidate energy or momentum decoupling
statement.  A finite nonzero child must therefore receive an explicit governed
bound and a new candidate identity before scientific evidence can attach to it.

This module changes no velocity value, coefficient, force, pressure, sample,
residual norm, or acceptance threshold.
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
_EXPECTED_UPSTREAM_HEAD = "78e4eea0d08116e4138e81db58edebba25706621"
_EXPECTED_ORTHOGONALIZER = 13.0 / 2.0


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
    _require(
        isinstance(actual, (int, float)) and math.isfinite(float(actual)),
        f"invalid numeric evidence: {name}",
    )
    _require(
        math.isclose(float(actual), expected, rel_tol=1e-9, abs_tol=atol),
        f"capacity evidence drift: {name}",
    )


def audit_bipolar_interior_swirl_redistribution_scope(
    scope: Mapping[str, Any],
    constraints: Mapping[str, Any],
    delivery_contract: Mapping[str, Any],
    source_candidate: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit provenance, orthogonality, coefficient, and promotion semantics."""
    _require(
        scope.get("schema") == "bipolar_interior_swirl_redistribution_capacity_scope_v1",
        "scope schema drift",
    )
    _require(
        scope.get("task_id") == "CR002-BIPOLAR-INTERIOR-SWIRL-REDISTRIBUTION-SCOPE-033",
        "task id drift",
    )

    vocabulary = delivery_contract.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "missing canonical source vocabulary")
    _require(set(vocabulary) == _CANONICAL_CLASSES, "canonical source vocabulary drift")
    _require(
        set(scope.get("classification_vocabulary", ())) == _CANONICAL_CLASSES,
        "scope source vocabulary drift",
    )

    source = scope.get("source_candidate")
    _require(isinstance(source, Mapping), "missing source candidate binding")
    _require(
        source.get("artifact") == "artifacts/bipolar_joint_capped/candidate.json",
        "source artifact drift",
    )
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
    _require(upstream.get("pr") == 236, "upstream PR drift")
    _require(upstream.get("head") == _EXPECTED_UPSTREAM_HEAD, "upstream head drift")
    _require(
        upstream.get("task_id") == "CR003-BIPOLAR-INTERIOR-SWIRL-REDISTRIBUTION-CAPACITY-038",
        "upstream task drift",
    )
    _require(upstream.get("classification") == "autonomous_design", "capacity evidence provenance drift")
    _require(
        upstream.get("evidence_role") == "target_free_representation_capacity_only",
        "capacity evidence role drift",
    )
    _require(upstream.get("mode") == "INTERIOR_C4_SWIRL_REDISTRIBUTION", "mode identity drift")
    _require(upstream.get("parent_compact_mode") == "INTERIOR_C4_SWIRL", "parent compact mode drift")
    _close(upstream.get("time"), 0.5, name="capacity time")
    _close(upstream.get("orthogonalizer"), _EXPECTED_ORTHOGONALIZER, name="orthogonalizer")
    _require(
        upstream.get("orthogonality_measure")
        == "raw compact-vs-redistribution basis response under cylindrical radial L2 inner product with weight r dr",
        "orthogonality measure drift",
    )
    _require(
        upstream.get("orthogonality_target") == "INTERIOR_C4_SWIRL raw basis response only",
        "orthogonality target drift",
    )
    _require(upstream.get("shape_parameters_fitted") == 0, "capacity shape was fitted")
    _require(upstream.get("residual_map_used_to_fit_basis_shape") is False, "residual-fitted basis shape enabled")
    _require(upstream.get("public_image_used_to_fit_basis_shape") is False, "public-image-fitted basis shape enabled")
    _require(upstream.get("normalized_response_rank") == 4, "capacity rank drift")
    _close(upstream.get("normalized_response_condition_number"), 9.84211, name="capacity condition", atol=1e-8)
    _close(
        upstream.get("redistribution_novelty_outside_existing_three_mode_span"),
        0.578515,
        name="capacity novelty",
        atol=1e-8,
    )
    _close(
        upstream.get("coefficient_step_refinement_relative_change"),
        5.86e-14,
        name="capacity refinement",
        atol=1e-15,
    )
    _require(upstream.get("pure_swirl_response") is True, "swirl-purity evidence drift")
    _close(upstream.get("radial_flank_response_rms"), 0.0, name="radial flank response")
    _close(upstream.get("axial_flank_response_rms"), 0.0, name="axial flank response")
    _close(upstream.get("collar_response_energy_fraction"), 0.0, name="collar response energy")
    _close(
        upstream.get("analytic_sign_change_r_over_Rp"),
        math.sqrt(2.0 / 13.0),
        name="analytic sign-change fraction",
    )
    _close(
        upstream.get("analytic_sign_change_radius"),
        1.6 * math.sqrt(2.0 / 13.0),
        name="analytic sign-change radius",
    )
    _require(float(upstream.get("inner_signed_swirl_mean")) > 0.0, "inner redistribution sign drift")
    _require(float(upstream.get("shoulder_signed_swirl_mean")) < 0.0, "shoulder redistribution sign drift")

    orthogonality = scope.get("orthogonality_semantics")
    _require(isinstance(orthogonality, Mapping), "missing orthogonality semantics")
    _require(orthogonality.get("orthogonality_is_analytic_representation_fact") is True, "analytic orthogonality fact disabled")
    _require(orthogonality.get("orthogonality_scope_is_raw_basis_pair_only") is True, "orthogonality scope widened")
    for key in (
        "orthogonality_is_to_full_source_candidate_velocity",
        "orthogonality_implies_full_candidate_energy_neutrality",
        "orthogonality_implies_finite_coefficient_energy_neutrality",
        "orthogonality_implies_radial_or_axial_momentum_neutrality",
        "pure_swirl_velocity_response_implies_theta_only_nonlinear_momentum_impact",
        "zero_collar_response_implies_full_energy_unchanged",
        "zero_collar_response_implies_physical_support_validated",
    ):
        _require(orthogonality.get(key) is False, f"orthogonality/collar inference promoted: {key}")
    _require(
        orthogonality.get("finite_nonzero_coefficient_has_positive_quadratic_energy_term") is True,
        "finite-coefficient quadratic energy term hidden",
    )

    semantics = scope.get("representation_semantics")
    _require(isinstance(semantics, Mapping), "missing representation semantics")
    _require(semantics.get("mode_origin") == "autonomous_design", "mode origin drift")
    _require(
        semantics.get("geometry_source")
        == "existing INTERIOR_C4_SWIRL and AxisymmetricPhysicalTaper identity-plateau geometry",
        "geometry provenance drift",
    )
    _require(semantics.get("is_eq45_public_source_profile_mode") is False, "autonomous correction laundered as public Eq45 mode")
    _require(semantics.get("is_existing_eq45_profile_basis_coefficient") is False, "redistribution laundered as profile coefficient")
    _require(semantics.get("is_additive_physical_space_swirl_correction") is True, "physical-space correction semantics drift")
    _require(semantics.get("axis_regular_by_r_factor") is True, "axis-regular construction drift")
    _require(semantics.get("sign_change_is_fixed_analytic_design_not_fit") is True, "analytic sign change relabeled as fit")
    _require(semantics.get("zero_before_physical_taper_collars") is True, "collar-localization semantics drift")
    _require(semantics.get("nonzero_materialization_changes_representation_family_identity") is True, "representation identity change disabled")
    _require(semantics.get("parent_candidate_identity_may_be_reused_for_nonzero_child") is False, "parent identity reuse enabled")

    coefficient = scope.get("coefficient_semantics")
    _require(isinstance(coefficient, Mapping), "missing coefficient semantics")
    _close(coefficient.get("capacity_evaluator_technical_guard"), 4.0, name="capacity technical guard")
    _require(
        coefficient.get("technical_guard_origin")
        == "borrowed implementation guard from parent profile coefficient_limit",
        "technical guard provenance drift",
    )
    _require(
        coefficient.get("technical_guard_is_preregistered_materialization_bound") is False,
        "technical guard laundered into materialization bound",
    )
    _require(coefficient.get("symmetric_trial_coefficients") == [-0.1, 0.1], "diagnostic trial grid drift")
    _require(coefficient.get("trial_role") == "target_free_morphology_probe_only", "diagnostic trial role drift")
    for key in (
        "trial_coefficients_are_selected",
        "trial_coefficients_define_materialization_bound",
        "trial_collar_energy_invariance_is_global_energy_acceptance",
    ):
        _require(coefficient.get(key) is False, f"diagnostic coefficient promotion: {key}")

    forbidden = scope.get("forbidden_inferences")
    _require(isinstance(forbidden, list) and len(forbidden) >= 12, "missing forbidden inferences")
    forbidden_text = "\n".join(str(item) for item in forbidden)
    for phrase in (
        "full source candidate is energy-orthogonal",
        "finite nonzero redistribution coefficient preserves E(0.25)",
        "theta-only nonlinear momentum impact",
        "full kinetic energy is unchanged",
        "physical-support validation",
        "should be selected",
        "preregistered materialization bound",
        "select a sign, value, or bound",
        "public-source Eq45 profile coefficient",
        "unchanged parent candidate SHA or family identity",
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
    _require(
        forcing.get("mode") == registered.get("forcing_mode") == "restricted_two_parameter_family",
        "forcing mode drift",
    )
    _require(forcing.get("parameters") == registered.get("forcing_parameters"), "forcing bounds drift")
    _require(registered.get("residual_defined_free_force_allowed") is False, "free-force route enabled")
    restriction = forcing.get("restriction")
    _require(
        isinstance(restriction, str) and "No residual-dependent basis or pointwise free force" in restriction,
        "free-force prohibition drift",
    )

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
        "explicit redistribution coefficient value and explicitly preregistered autonomous materialization bound",
        "new representation-family label describing Eq45 parent plus autonomous sign-changing physical-space swirl redistribution correction",
        "new serialized candidate SHA distinct from the unchanged parent",
        "initial energy E(0.25)=1 within registered tolerance revalidated on the full child",
        "minimum/maximum validation-time energy gates revalidated on the full child",
        "nonzero center rotation and required rotation sign revalidated",
        "bipolar radial/axial direction and parity revalidated",
        "axis/support regularity revalidated",
        "full per-component momentum and divergence recomputed because nonlinear coupling is not theta-only",
        "registered 4096-point fresh/pristine validation remains separate from model-selection evidence",
        "restricted forcing family and bounds unchanged unless a new preregistered experiment is declared",
    )
    for phrase in required_phrases:
        _require(phrase in requirements, f"missing materialization requirement: {phrase}")

    states = scope.get("states")
    _require(isinstance(states, Mapping), "missing state boundary")
    _require(states.get("source_candidate_velocity_export_ready") is True, "source export state drift")
    _require(states.get("interior_redistribution_mode_capacity_screened") is True, "capacity-screen state drift")
    for key in (
        "interior_redistribution_nonzero_candidate_materialized",
        "interior_redistribution_coefficient_selected",
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
        "upstream_pr": 236,
        "mode": "INTERIOR_C4_SWIRL_REDISTRIBUTION",
        "mode_origin": "autonomous_design",
        "raw_basis_pair_orthogonality_only": True,
        "full_candidate_energy_neutrality_claimed": False,
        "technical_guard_is_materialization_bound": False,
        "nonzero_materialization_requires_new_identity": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def audit_bipolar_interior_swirl_redistribution_scope_file(
    scope_path: str | Path = "configs/bipolar_interior_swirl_redistribution_capacity_scope.json",
    *,
    repo_root: str | Path = ".",
) -> dict[str, Any]:
    """Load governed repository files and audit the redistribution scope."""
    root = Path(repo_root)
    scope = _read_json(root / Path(scope_path))
    constraints = _read_json(root / "configs/constraints.json")
    delivery = _read_json(root / "configs/delivery_state_contract.json")
    source = _read_json(root / str(scope["source_candidate"]["artifact"]))
    return audit_bipolar_interior_swirl_redistribution_scope(scope, constraints, delivery, source)
