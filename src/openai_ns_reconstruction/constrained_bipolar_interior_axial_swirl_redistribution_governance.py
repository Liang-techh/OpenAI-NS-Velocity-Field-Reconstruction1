"""Fail-closed governance for the axial interior swirl-redistribution screen.

The screened INTERIOR_C4_SWIRL_AXIAL_REDISTRIBUTION direction is an autonomous
additive physical-space correction.  Its raw response is constructed to be
axially L2-orthogonal to the one-sign compact response, but that fact does not
make a finite corrected *candidate* energy-neutral and does not select a
coefficient.  A nonzero materialization must get its own representation
identity and must be revalidated under the unchanged CR001 contract.

This module changes no velocity value, coefficient, force, pressure, sample,
or acceptance threshold.
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
_EXPECTED_UPSTREAM_HEAD = "c0c50c72ead554c0e041a0d2476a838d13fdccc7"
_EXPECTED_MODE = "INTERIOR_C4_SWIRL_AXIAL_REDISTRIBUTION"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _close(actual: Any, expected: float, *, name: str, atol: float = 1e-12, rtol: float = 1e-9) -> None:
    _require(isinstance(actual, (int, float)) and math.isfinite(float(actual)), f"invalid numeric evidence: {name}")
    _require(math.isclose(float(actual), expected, rel_tol=rtol, abs_tol=atol), f"capacity evidence drift: {name}")


def audit_bipolar_interior_axial_swirl_redistribution_scope(
    scope: Mapping[str, Any],
    constraints: Mapping[str, Any],
    delivery_contract: Mapping[str, Any],
    source_candidate: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit provenance, orthogonality, sign, coefficient, and promotion semantics."""
    _require(scope.get("schema") == "bipolar_interior_axial_swirl_redistribution_scope_v1", "scope schema drift")
    _require(
        scope.get("task_id") == "CR002-BIPOLAR-INTERIOR-AXIAL-SWIRL-REDISTRIBUTION-SCOPE-034",
        "task id drift",
    )

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
    _require(upstream.get("pr") == 248, "upstream PR drift")
    _require(upstream.get("head") == _EXPECTED_UPSTREAM_HEAD, "upstream head drift")
    _require(
        upstream.get("task_id") == "CR003-BIPOLAR-INTERIOR-AXIAL-SWIRL-REDISTRIBUTION-CAPACITY-039",
        "upstream task drift",
    )
    _require(upstream.get("classification") == "autonomous_design", "capacity evidence provenance drift")
    _require(upstream.get("evidence_role") == "target_free_representation_capacity_only", "capacity evidence role drift")
    _require(upstream.get("mode") == _EXPECTED_MODE, "mode identity drift")
    _require(upstream.get("shape_parameters_fitted") == 0, "capacity shape was fitted")
    _require(upstream.get("residual_map_used_to_fit_basis_shape") is False, "residual-fitted basis shape enabled")
    _require(upstream.get("public_image_used_to_fit_basis_shape") is False, "public-image-fitted basis shape enabled")
    _close(upstream.get("orthogonalizer"), 23.0, name="axial orthogonalizer")
    _close(upstream.get("radial_plateau_half_width"), 1.6, name="capacity radial plateau")
    _close(upstream.get("axial_plateau_half_width"), 1.6, name="capacity axial plateau")
    _close(upstream.get("analytic_sign_change_abs_z_over_Zp"), 1.0 / math.sqrt(23.0), name="sign-change fraction")
    _close(upstream.get("analytic_sign_change_abs_z"), 1.6 / math.sqrt(23.0), name="sign-change height")
    _require(upstream.get("normalized_response_rank") == 4, "capacity rank drift")
    _close(upstream.get("normalized_response_condition_number"), 10.0479, name="capacity condition", atol=1e-4, rtol=1e-5)
    _close(upstream.get("novelty_outside_F10_F20_compact_span"), 0.844934, name="capacity novelty", atol=1e-6, rtol=1e-6)
    _close(upstream.get("coefficient_step_refinement_relative_change"), 4.08e-14, name="capacity refinement", atol=1e-15, rtol=1e-2)
    _close(upstream.get("swirl_fraction"), 1.0, name="swirl fraction")
    _require(float(upstream.get("center_signed_swirl_mean")) > 0.0, "center response sign drift")
    _require(float(upstream.get("axial_shoulder_signed_swirl_mean")) < 0.0, "axial-shoulder response sign drift")
    _close(upstream.get("radial_flank_response_rms"), 0.0, name="radial flank response RMS")
    _close(upstream.get("axial_flank_response_rms"), 0.0, name="axial flank response RMS")
    _close(upstream.get("collar_response_energy_fraction"), 0.0, name="collar response energy")
    _close(upstream.get("diagnostic_coefficient"), 0.10, name="diagnostic coefficient")
    _require(upstream.get("diagnostic_coefficient_is_candidate_bound") is False, "diagnostic coefficient laundered into bound")
    _require(upstream.get("visual_correspondence_evidence") is False, "capacity evidence promoted to visual correspondence")
    _require(upstream.get("perturbed_candidate_pde_residual_evaluated") is False, "capacity screen laundered into PDE evidence")

    orthogonality = scope.get("orthogonality_semantics")
    _require(isinstance(orthogonality, Mapping), "missing orthogonality semantics")
    _require(orthogonality.get("raw_compact_vs_axial_redistribution_axial_L2_orthogonal") is True, "raw axial orthogonality drift")
    _require(orthogonality.get("orthogonality_is_to_complete_candidate_velocity") is False, "raw orthogonality laundered to complete candidate")
    _require(orthogonality.get("orthogonality_implies_finite_child_energy_unchanged") is False, "orthogonality laundered to energy neutrality")
    weighting = orthogonality.get("weighting")
    _require(isinstance(weighting, str) and "(1-q^2)^10" in weighting, "orthogonality weighting drift")

    representation = scope.get("representation_semantics")
    _require(isinstance(representation, Mapping), "missing representation semantics")
    _require(representation.get("mode_origin") == "autonomous_design", "mode origin drift")
    _require(
        representation.get("geometry_source") == "existing AxisymmetricPhysicalTaper identity-plateau geometry",
        "geometry provenance drift",
    )
    _require(representation.get("is_eq45_public_source_profile_mode") is False, "autonomous correction laundered as public Eq45 mode")
    _require(representation.get("is_existing_eq45_profile_basis_coefficient") is False, "physical-space correction laundered as profile coefficient")
    _require(representation.get("is_additive_physical_space_swirl_correction") is True, "physical-space correction semantics drift")
    _require(representation.get("axis_regular_by_parent_r_factor") is True, "axis-regular construction drift")
    _require(representation.get("axisymmetric_pure_swirl_is_analytically_divergence_free") is True, "structural divergence identity drift")
    _require(representation.get("zero_before_physical_taper_collars") is True, "collar-localization semantics drift")
    _require(representation.get("nonzero_materialization_changes_representation_family_identity") is True, "representation identity change disabled")
    _require(representation.get("parent_candidate_identity_may_be_reused_for_nonzero_child") is False, "parent identity reuse enabled")

    signs = scope.get("local_rotation_sign_semantics")
    _require(isinstance(signs, Mapping), "missing local rotation-sign semantics")
    _require(signs.get("basis_response_changes_sign_between_axial_center_and_interior_shoulders") is True, "basis sign-change semantics drift")
    _require(signs.get("capacity_sign_change_is_not_a_required_global_counter_rotation") is True, "capacity sign change promoted to requirement")
    _require(signs.get("canonical_core_requirement") == "u_theta > 0 at registered core_probe", "core-sign requirement drift")
    _require(signs.get("global_u_theta_positivity_is_preregistered") is False, "invented global swirl-positivity gate")
    _require(signs.get("future_nonzero_child_must_recheck_registered_core_sign") is True, "core-sign recheck disabled")
    _require(signs.get("future_nonzero_child_should_report_predeclared_center_and_shoulder_signed_swirl") is True, "local sign reporting disabled")
    _require(signs.get("capacity_probe_alone_may_not_select_coefficient_sign") is True, "capacity sign probe allowed to select coefficient")
    structure = constraints.get("structure")
    _require(isinstance(structure, Mapping), "missing canonical structure contract")
    core_signs = structure.get("core_sign_requirements")
    _require(isinstance(core_signs, list) and "u_theta > 0" in core_signs, "registered core swirl-sign requirement drift")

    coefficient = scope.get("coefficient_semantics")
    _require(isinstance(coefficient, Mapping), "missing coefficient semantics")
    _close(coefficient.get("capacity_evaluator_inherited_abs_guard"), 4.0, name="capacity evaluator guard")
    _require(coefficient.get("inherited_guard_is_preregistered_bound_for_this_autonomous_mode") is False, "inherited guard laundered into mode bound")
    _require(coefficient.get("symmetric_diagnostic_coefficients") == [-0.10, 0.10], "diagnostic coefficient grid drift")
    _require(coefficient.get("diagnostic_coefficients_are_selected") is False, "diagnostic coefficients promoted to selection")
    _require(coefficient.get("diagnostic_coefficients_define_materialization_bound") is False, "diagnostic coefficients promoted to materialization bound")
    _require(coefficient.get("future_nonzero_materialization_requires_explicit_autonomous_bound") is True, "explicit autonomous bound requirement disabled")

    forbidden = scope.get("forbidden_inferences")
    _require(isinstance(forbidden, list) and len(forbidden) >= 10, "missing forbidden inferences")
    forbidden_text = "\n".join(str(item) for item in forbidden)
    for phrase in (
        "full candidate energy is unchanged",
        "theta-only nonlinear momentum impact",
        "full physical-support validation",
        "registered numerical divergence gate passed",
        "selects a nonzero coefficient",
        "materialization bound",
        "counter-rotation",
        "global u_theta positivity",
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
    _require(nontriviality.get("reference_energy_abs_tolerance") == registered.get("reference_energy_abs_tolerance") == 0.001, "reference energy tolerance drift")
    _require(nontriviality.get("minimum_energy_each_validation_time") == registered.get("minimum_energy_each_validation_time") == 0.1, "minimum validation energy drift")
    _require(nontriviality.get("maximum_energy_each_validation_time") == registered.get("maximum_energy_each_validation_time") == 10.0, "maximum validation energy drift")

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
    for phrase in (
        "explicit autonomous coefficient value and explicit autonomous coefficient bound",
        "new representation-family label describing Eq45 parent plus autonomous axial swirl redistribution correction",
        "new serialized candidate SHA distinct from the unchanged parent",
        "initial energy E(0.25)=1 within registered tolerance revalidated",
        "minimum/maximum validation-time energy gates revalidated",
        "registered core u_theta positive and nonzero center rotation revalidated",
        "bipolar radial/axial direction and parity revalidated",
        "predeclared center and axial-shoulder signed-swirl diagnostics reported without inventing a global positivity gate",
        "axis/support regularity revalidated",
        "full per-component momentum and divergence recomputed",
        "registered 4096-point fresh/pristine validation remains separate from model-selection evidence",
        "restricted forcing family and bounds unchanged unless a new preregistered experiment is declared",
    ):
        _require(phrase in requirements, f"missing materialization requirement: {phrase}")

    states = scope.get("states")
    _require(isinstance(states, Mapping), "missing state boundary")
    _require(states.get("source_candidate_velocity_export_ready") is True, "source export state drift")
    _require(states.get("axial_redistribution_capacity_screened") is True, "capacity-screen state drift")
    for key in (
        "axial_redistribution_nonzero_candidate_materialized",
        "axial_redistribution_coefficient_selected",
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
        "upstream_pr": 248,
        "mode": _EXPECTED_MODE,
        "mode_origin": "autonomous_design",
        "raw_axial_orthogonality_only": True,
        "global_swirl_positivity_preregistered": False,
        "diagnostic_coefficient_selected": False,
        "nonzero_materialization_requires_new_identity": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def audit_bipolar_interior_axial_swirl_redistribution_scope_file(
    scope_path: str | Path = "configs/bipolar_interior_axial_swirl_redistribution_scope.json",
    *,
    repo_root: str | Path = ".",
) -> dict[str, Any]:
    """Load governed repository files and audit the axial redistribution scope."""
    root = Path(repo_root)
    scope = _read_json(root / Path(scope_path))
    constraints = _read_json(root / "configs/constraints.json")
    delivery = _read_json(root / "configs/delivery_state_contract.json")
    source = _read_json(root / str(scope["source_candidate"]["artifact"]))
    return audit_bipolar_interior_axial_swirl_redistribution_scope(scope, constraints, delivery, source)
