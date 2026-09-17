"""Fail-closed governance for the joint radial/axial interior-swirl stopping screen.

The upstream screen is target-free representation capacity evidence.  Public
velocity rank, response orthogonality, and morphology rank are not residual-
Jacobian completeness, coefficient selection, visual correspondence, or PDE
validation.  The current no-mixed-mode routing is therefore a local stopping
rule only, not a global basis-completeness theorem.
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
_EXPECTED_HEAD = "04829f63cbbcd093610684ed92f52b411a90356a"
_EXPECTED_CANDIDATE_SHA = "7875214ad65f5e1ba4f7c5e362217f05ef47bc3d893e8641e8e40a85fc4e7609"
_ROOT = Path(__file__).resolve().parents[2]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _close(actual: Any, expected: float, name: str, *, atol: float = 1e-12, rtol: float = 1e-8) -> None:
    _require(isinstance(actual, (int, float)) and math.isfinite(float(actual)), f"invalid numeric evidence: {name}")
    _require(math.isclose(float(actual), expected, rel_tol=rtol, abs_tol=atol), f"joint capacity evidence drift: {name}")


def audit_bipolar_interior_swirl_2d_stopping_scope(
    scope: Mapping[str, Any],
    constraints: Mapping[str, Any],
    delivery_contract: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate the joint-capacity stopping rule without promoting its evidence."""
    _require(scope.get("schema") == "bipolar_interior_swirl_2d_stopping_scope_v1", "scope schema drift")
    _require(scope.get("task_id") == "CR002-BIPOLAR-INTERIOR-SWIRL-2D-STOPPING-SCOPE-035", "task id drift")

    vocabulary = delivery_contract.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "missing canonical source vocabulary")
    _require(set(vocabulary) == _CANONICAL_CLASSES, "canonical source vocabulary drift")
    _require(set(scope.get("classification_vocabulary", ())) == _CANONICAL_CLASSES, "scope source vocabulary drift")

    upstream = scope.get("upstream_capacity_evidence")
    _require(isinstance(upstream, Mapping), "missing upstream joint-capacity evidence")
    _require(upstream.get("pr") == 259, "upstream PR drift")
    _require(upstream.get("head") == _EXPECTED_HEAD, "upstream head drift")
    _require(upstream.get("task_id") == "CR003-BIPOLAR-INTERIOR-SWIRL-2D-JACOBIAN-040", "upstream task drift")
    _require(upstream.get("classification") == "autonomous_design", "upstream source classification drift")
    _require(upstream.get("evidence_role") == "target_free_joint_representation_capacity_and_routing_only", "upstream evidence role drift")
    _require(upstream.get("source_candidate_sha256") == _EXPECTED_CANDIDATE_SHA, "source candidate identity drift")
    _require(upstream.get("modes") == ["INTERIOR_C4_SWIRL_REDISTRIBUTION", "INTERIOR_C4_SWIRL_AXIAL_REDISTRIBUTION"], "joint mode identity drift")
    for key in ("new_basis_mode_added", "coefficients_selected", "residual_map_used_to_fit_basis_shape", "public_image_used_to_fit_basis_shape", "perturbed_candidate_pde_residual_evaluated"):
        _require(upstream.get(key) is False, f"upstream capacity promoted: {key}")

    _require(upstream.get("normalized_public_velocity_rank") == 5, "public velocity rank drift")
    _close(upstream.get("normalized_public_velocity_condition_number"), 13.7220, "public velocity condition", atol=1e-4, rtol=1e-5)
    _close(upstream.get("response_step_refinement_relative_change"), 4.67e-14, "response refinement", atol=1e-15, rtol=1e-2)
    _close(upstream.get("radial_novelty"), 0.66427, "radial novelty", atol=1e-5, rtol=1e-5)
    _close(upstream.get("axial_novelty"), 0.66331, "axial novelty", atol=1e-5, rtol=1e-5)
    _require(upstream.get("radial_axial_pair_rank") == 2, "pair rank drift")
    _close(upstream.get("radial_axial_pair_condition_number"), 1.07811, "pair condition", atol=1e-5, rtol=1e-5)
    _close(upstream.get("radial_axial_pair_cosine"), 0.07506, "pair cosine", atol=1e-5, rtol=1e-5)
    _require(upstream.get("morphology_rank") == 2, "morphology rank drift")
    _close(upstream.get("morphology_normalized_condition_number"), 22.3982, "morphology condition", atol=1e-4, rtol=1e-5)
    _close(upstream.get("centroid_only_condition_number"), 24.5752, "centroid condition", atol=1e-4, rtol=1e-5)
    _close(upstream.get("morphology_step_refinement_relative_change"), 7.47e-5, "morphology refinement", atol=1e-7, rtol=1e-3)
    _close(upstream.get("radial_primary_centroid_fraction"), 0.9841, "radial primary-centroid fraction", atol=1e-4, rtol=1e-4)
    _close(upstream.get("axial_primary_centroid_fraction"), 0.1055, "axial primary-centroid fraction", atol=1e-4, rtol=1e-4)
    _close(upstream.get("diagnostic_corner_coefficient_abs"), 0.10, "diagnostic coefficient")
    _require(upstream.get("diagnostic_coefficient_is_materialization_bound") is False, "diagnostic coefficient laundered into bound")
    _close(upstream.get("radial_centroid_span"), 0.0087360, "radial centroid span", atol=1e-7, rtol=1e-5)
    _close(upstream.get("abs_z_centroid_span"), 0.00034105, "axial centroid span", atol=1e-8, rtol=1e-5)
    _close(upstream.get("physical_plateau_collar_swirl_energy_change"), 0.0, "collar energy change")

    semantics = scope.get("joint_capacity_semantics")
    _require(isinstance(semantics, Mapping), "missing joint-capacity semantics")
    required_false = (
        "public_velocity_rank_is_residual_jacobian_rank",
        "public_velocity_rank_implies_pde_leverage",
        "morphology_rank_implies_visual_correspondence",
        "raw_response_orthogonality_implies_finite_child_energy_neutrality",
        "weak_global_axial_centroid_leverage_proves_axial_swirl_is_useless",
        "weak_global_axial_centroid_leverage_proves_poloidal_support_fix_is_correct",
        "zero_audited_collar_response_implies_physical_support_validated",
    )
    _require(semantics.get("raw_radial_and_axial_responses_are_independent_in_public_velocity_space") is True, "joint public-velocity independence drift")
    for key in required_false:
        _require(semantics.get(key) is False, f"forbidden joint inference enabled: {key}")

    stop = scope.get("stopping_rule_semantics")
    _require(isinstance(stop, Mapping), "missing stopping-rule semantics")
    _require(stop.get("current_rule") == "do not add a mixed r-z swirl basis on public-velocity capacity-rank grounds alone", "stopping rule drift")
    _require(stop.get("rule_scope") == "current target-free two-mode capacity screen only", "stopping-rule scope drift")
    for key in ("global_no_mixed_mode_theorem", "global_basis_completeness_claim", "candidate_selection_resolved", "optimizer_stagnation_alone_is_sufficient_to_reopen_basis_growth", "green_ci_alone_is_sufficient_to_reopen_or_close_basis_growth", "routing_to_poloidal_support_lane_is_a_scientific_acceptance_result"):
        _require(stop.get(key) is False, f"stopping rule over-promoted: {key}")
    _require(stop.get("future_mixed_mode_may_be_reopened_only_by_new_independent_blocker_evidence") is True, "independent evidence requirement disabled")

    materialization = scope.get("materialization_semantics")
    _require(isinstance(materialization, Mapping), "missing materialization semantics")
    for key in ("nonzero_radial_child_materialized", "nonzero_axial_child_materialized", "joint_nonzero_child_materialized", "radial_coefficient_selected", "axial_coefficient_selected", "diagnostic_plus_minus_0p10_selects_sign_or_value", "inherited_abs_4_implementation_guard_is_preregistered_bound_for_autonomous_modes"):
        _require(materialization.get(key) is False, f"capacity evidence laundered into materialization: {key}")
    for key in ("future_nonzero_child_requires_explicit_autonomous_bounds", "future_nonzero_child_requires_new_family_identity_and_sha", "future_nonzero_child_requires_initial_and_validation_time_energy_recheck", "future_nonzero_child_requires_core_rotation_and_bipolar_flow_sign_recheck", "future_nonzero_child_requires_axis_support_regularity_recheck", "future_nonzero_child_requires_full_per_component_momentum_and_divergence_recheck", "future_acceptance_data_must_be_pristine_from_model_selection"):
        _require(materialization.get(key) is True, f"future child guard disabled: {key}")

    registered = scope.get("registered_contract")
    _require(isinstance(registered, Mapping), "missing registered CR001 snapshot")
    _require(constraints.get("nu") == registered.get("nu") == 0.01, "nu drift")
    domain = constraints.get("domain")
    _require(isinstance(domain, Mapping), "missing canonical domain")
    _require(domain.get("physical") == registered.get("physical_domain") == "R^3", "physical domain drift")
    _require(domain.get("evaluation_box") == registered.get("evaluation_box"), "evaluation box drift")
    _require(domain.get("support") == registered.get("support"), "support drift")
    _require(domain.get("time_interval") == registered.get("time_interval"), "time interval drift")
    forcing = constraints.get("forcing")
    _require(isinstance(forcing, Mapping), "missing canonical forcing")
    _require(forcing.get("mode") == registered.get("forcing_mode") == "restricted_two_parameter_family", "forcing mode drift")
    _require(forcing.get("parameters") == registered.get("forcing_parameters"), "forcing bounds drift")
    _require(registered.get("residual_defined_free_force_allowed") is False, "free-force route enabled")
    _require("No residual-dependent basis or pointwise free force" in str(forcing.get("restriction")), "canonical free-force prohibition drift")
    nontriviality = constraints.get("nontriviality")
    _require(isinstance(nontriviality, Mapping), "missing canonical nontriviality")
    for key in ("reference_energy", "reference_energy_abs_tolerance", "minimum_energy_each_validation_time", "maximum_energy_each_validation_time"):
        _require(nontriviality.get(key) == registered.get(key), f"nontriviality drift: {key}")
    validation = constraints.get("validation")
    _require(isinstance(validation, Mapping), "missing canonical validation")
    _require(validation.get("seed") == registered.get("validation_seed") == 914027, "validation seed drift")
    _require(validation.get("held_out_points") == registered.get("held_out_points") == 4096, "held-out count drift")
    _require(validation.get("derivative_steps") == registered.get("derivative_steps"), "derivative ladder drift")
    thresholds = validation.get("thresholds")
    _require(isinstance(thresholds, Mapping), "missing canonical thresholds")
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds.get(key) == registered.get(key), f"threshold drift: {key}")

    truth = scope.get("truth_boundary")
    _require(isinstance(truth, Mapping), "missing truth boundary")
    _require(truth.get("velocity_export_ready") is True, "delivery state regressed")
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified", "blowup_proved"):
        _require(truth.get(key) is False, f"unsupported truth promotion: {key}")

    forbidden = scope.get("forbidden_inferences")
    _require(isinstance(forbidden, list) and len(forbidden) >= 10, "missing forbidden inferences")
    text = "\n".join(str(item) for item in forbidden)
    for phrase in ("residual-Jacobian", "energy-neutral", "visual correspondence", "can never be useful", "poloidal/support correction", "physical-support validation", "select coefficient", "absolute 4", "no mixed r-z swirl mode", "PDE validation", "paper-exact"):
        _require(phrase in text, f"missing forbidden inference: {phrase}")

    return {
        "task_id": scope["task_id"],
        "upstream_head": upstream["head"],
        "stopping_rule_scope": stop["rule_scope"],
        "registered_validation_points": registered["held_out_points"],
        "pde_validated": False,
    }


def audit_repository() -> dict[str, Any]:
    return audit_bipolar_interior_swirl_2d_stopping_scope(
        _read_json(_ROOT / "configs/bipolar_interior_swirl_2d_stopping_scope.json"),
        _read_json(_ROOT / "configs/constraints.json"),
        _read_json(_ROOT / "configs/delivery_state_contract.json"),
    )


def main() -> None:
    print(json.dumps(audit_repository(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
