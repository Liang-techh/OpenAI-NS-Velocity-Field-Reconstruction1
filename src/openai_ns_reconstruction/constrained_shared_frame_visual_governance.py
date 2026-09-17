"""Govern evidentiary use of the shared-frame multi-candidate visual comparator.

This is a governance-only layer.  It changes no velocity values, candidate
parameters, support transform, forcing, thresholds, optimizer, or canonical
routing.  Its purpose is to keep fair visual comparability distinct from public
visual correspondence, true 3-D streamline evidence, PDE acceptance, and exact
OpenAI-field identity.
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path
from typing import Any

from .constrained_candidate_visual_compare import (
    SharedFrameCandidateComparison,
    sample_shared_frame_candidates,
)

EXPECTED_SCHEMA = "shared_frame_visual_evidence_scope_v1"
EXPECTED_TASK = "CR002-SHARED-FRAME-VISUAL-EVIDENCE-SCOPE-028"
EXPECTED_SOURCE_TASK = "CR-A9-023"
EXPECTED_SOURCE_PR = 192
EXPECTED_SOURCE_HEAD = "47c3fb66ebb19d46565d1f1dce159bc9597c1a7b"
EXPECTED_INTEGRATION_HEAD = "7dd05bb0f1f68e593e37e5f3a4eee926a89fbdf5"
EXPECTED_SOURCE_VOCABULARY = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


def _load(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _same_float(left: Any, right: Any, *, atol: float = 1.0e-12) -> bool:
    try:
        return abs(float(left) - float(right)) <= atol
    except (TypeError, ValueError):
        return False


def audit_shared_frame_visual_evidence_scope(
    *,
    contract_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
    delivery_state_path: str | Path | None = None,
    project_status_path: str | Path | None = None,
) -> dict[str, Any]:
    """Fail closed if fair-render evidence is promoted beyond its declared scope."""
    root = Path(__file__).resolve().parents[2]
    contract = _load(contract_path or root / "configs/shared_frame_visual_evidence_scope.json")
    constraints = _load(constraints_path or root / "configs/constraints.json")
    delivery = _load(delivery_state_path or root / "configs/delivery_state_contract.json")
    status = _load(project_status_path or root / "project_status.json")

    _require(contract.get("schema") == EXPECTED_SCHEMA, "shared-frame governance schema drifted")
    _require(contract.get("task_id") == EXPECTED_TASK, "shared-frame governance task drifted")

    canonical_vocab = set(delivery["classification_vocabulary"])
    local_vocab = set(contract.get("source_vocabulary", []))
    _require(canonical_vocab == EXPECTED_SOURCE_VOCABULARY, "delivery-state source vocabulary drifted")
    _require(local_vocab == EXPECTED_SOURCE_VOCABULARY, "shared-frame source vocabulary drifted")
    classifications = contract.get("source_classification", {})
    _require(set(classifications.values()) <= EXPECTED_SOURCE_VOCABULARY, "unknown source classification")
    expected_classifications = {
        "callable_saveable_loadable_3d_time_varying_velocity": "user_requirement",
        "independent_velocity_visual_pde_and_exactness_states": "user_requirement",
        "same_frame_same_time_candidate_comparison": "autonomous_design",
        "single_global_speed_normalization": "autonomous_design",
        "projected_meridional_streamline_overlay": "autonomous_design",
        "renderer_activity_floor": "autonomous_design",
        "public_openai_reference_correspondence": "pending_unknown",
        "paper_exact_or_hidden_openai_field_identity": "pending_unknown",
    }
    _require(classifications == expected_classifications, "shared-frame source classification drifted")

    binding = contract["integration_binding"]
    _require(binding["source_task"] == EXPECTED_SOURCE_TASK, "source task drifted")
    _require(int(binding["source_pr"]) == EXPECTED_SOURCE_PR, "source PR drifted")
    _require(binding["source_head"] == EXPECTED_SOURCE_HEAD, "source head drifted")
    _require(binding["integration_head"] == EXPECTED_INTEGRATION_HEAD, "integration head drifted")
    _require(
        binding["module"] == "openai_ns_reconstruction.constrained_candidate_visual_compare",
        "shared-frame module binding drifted",
    )
    _require(binding["sampler"] == "sample_shared_frame_candidates", "sampler binding drifted")
    _require(binding["component_order"] == ["u", "v", "w"], "velocity component order drifted")

    semantics = contract["comparison_semantics"]
    required_true = (
        "same_physical_frame_for_all_candidates",
        "same_time_values_for_all_candidates",
        "full_cartesian_velocity_retained",
        "single_global_speed_scale_across_all_candidates_and_times",
        "swirl_component_retained_in_saved_velocity",
    )
    for key in required_true:
        _require(semantics[key] is True, f"{key} must remain true")
    required_false = (
        "per_panel_autoscaling",
        "projected_streamlines_are_true_3d_streamlines",
        "camera_fit_performed",
        "registration_fit_performed",
        "hidden_time_alignment_performed",
        "candidate_selection_performed",
        "public_reference_comparison_performed",
    )
    for key in required_false:
        _require(semantics[key] is False, f"{key} must remain false")
    _require(semantics["projected_streamline_components"] == ["u", "w"], "projected streamline components drifted")

    fields = SharedFrameCandidateComparison.__dataclass_fields__
    _require(fields["per_panel_autoscaling"].default is False, "live comparator enabled per-panel autoscaling")
    _require(fields["camera_fitted"].default is False, "live comparator claims camera fitting")
    _require(fields["registration_fitted"].default is False, "live comparator claims registration fitting")
    _require(fields["visualization_ready"].default is False, "live comparator promoted visualization readiness")
    _require(
        fields["visual_correspondence_verified"].default is False,
        "live comparator promoted visual correspondence",
    )
    _require(fields["pde_validated"].default is False, "live comparator promoted PDE validity")
    _require(fields["paper_exact"].default is False, "live comparator promoted paper exactness")
    _require(fields["openai_field_identified"].default is False, "live comparator identified hidden OpenAI field")
    _require(
        fields["projected_streamline_semantics"].default
        == "meridional_(u,w)_projection_only_not_true_3d_streamline",
        "live projected-streamline semantics drifted",
    )
    signature = inspect.signature(sample_shared_frame_candidates)
    _require(
        _same_float(signature.parameters["min_panel_activity"].default, 1.0e-12),
        "live renderer activity-floor default drifted",
    )

    guard = contract["activity_guard"]
    _require(_same_float(guard["default_min_panel_activity"], 1.0e-12), "activity guard value drifted")
    _require(
        guard["role"] == "fail_closed_rendering_guard_against_zero_or_numerically_inactive_panels",
        "activity guard role drifted",
    )
    for key in (
        "scientific_acceptance_threshold",
        "may_replace_cr001_energy_nontriviality",
        "may_select_candidate",
        "may_promote_visual_correspondence_verified",
    ):
        _require(guard[key] is False, f"activity guard incorrectly promoted: {key}")

    cr001 = contract["cr001_binding"]
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]
    _require(_same_float(cr001["nu"], constraints["nu"]), "nu drifted")
    _require(cr001["physical_domain"] == domain["physical"], "physical domain drifted")
    _require(cr001["evaluation_box"] == domain["evaluation_box"], "evaluation box drifted")
    _require(cr001["support"] == domain["support"], "support contract drifted")
    _require(cr001["time_interval"] == domain["time_interval"], "time interval drifted")
    _require(cr001["forcing_mode"] == forcing["mode"], "forcing mode drifted")
    _require(cr001["forcing_parameters"] == sorted(forcing["parameters"]), "forcing parameter set drifted")
    _require(
        "No residual-dependent basis or pointwise free force" in forcing["restriction"],
        "restricted forcing no-free-force guard drifted",
    )
    _require(_same_float(cr001["reference_energy_time"], nontriviality["reference_time"]), "energy time drifted")
    _require(_same_float(cr001["reference_energy"], nontriviality["reference_energy"]), "energy target drifted")
    _require(
        _same_float(cr001["reference_energy_abs_tolerance"], nontriviality["reference_energy_abs_tolerance"]),
        "energy tolerance drifted",
    )
    _require(int(cr001["held_out_points"]) == int(validation["held_out_points"]), "held-out count drifted")
    _require(
        [float(x) for x in cr001["derivative_steps"]]
        == [float(x) for x in validation["derivative_steps"]],
        "derivative ladder drifted",
    )
    _require(cr001["registered_norms"] == validation["norms"], "registered residual norms drifted")
    for local_key, registered_key in (
        ("divergence_max_threshold", "divergence_max"),
        ("divergence_L2_threshold", "divergence_L2"),
        ("pde_residual_max_threshold", "pde_residual_max"),
        ("pde_residual_L2_threshold", "pde_residual_L2"),
    ):
        _require(_same_float(cr001[local_key], thresholds[registered_key]), f"{registered_key} threshold drifted")

    promotion = contract["promotion_contract"]
    for key in (
        "shared_frame_comparability_is_visual_correspondence",
        "global_normalization_is_visual_correspondence",
        "projected_2d_streamlines_are_sufficient_for_3d_streamline_acceptance",
        "green_ci_or_render_success_may_select_candidate",
        "serialization_or_exportability_may_select_candidate",
        "formal_pde_pass_required_for_candidate_local_save_load_export",
        "formal_pde_pass_required_for_truth_bounded_visualization_candidate",
        "candidate_selection_resolved",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _require(promotion[key] is False, f"promotion shortcut enabled: {key}")
    for key in (
        "public_reference_comparison_required_for_visual_correspondence",
        "independent_visual_resolution_stability_required_for_visual_correspondence",
        "true_3d_streamline_or_equivalent_full_vector_diagnostic_required_when_streamline_structure_is_claimed",
        "formal_pde_gate_remains_independent",
    ):
        _require(promotion[key] is True, f"required promotion evidence disabled: {key}")

    states = contract["truth_states"]
    expected_states = {
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    _require(states == expected_states, "shared-frame truth-state boundary drifted")
    for state_name, state_value in expected_states.items():
        if state_name in status["states"]:
            _require(status["states"][state_name] is state_value, f"project status disagrees on {state_name}")

    forbidden = {
        "same_frame_same_time_comparison=>candidate_selection_resolved",
        "single_global_speed_scale=>visual_correspondence_verified",
        "projected_meridional_streamlines=>true_3d_streamline_acceptance",
        "renderer_activity_guard_pass=>cr001_nontriviality_pass",
        "green_ci_or_render_success=>visual_correspondence_verified",
        "serialization_or_exportability=>candidate_selection_resolved",
        "visual_preference_without_public_reference=>openai_field_identified",
        "visual_similarity=>pde_validated_or_paper_exact",
        "pde_failed_or_pending=>velocity_export_not_allowed",
    }
    _require(set(contract["forbidden_inferences"]) == forbidden, "forbidden inference set drifted")

    return {
        "integration_head": EXPECTED_INTEGRATION_HEAD,
        "source_pr": EXPECTED_SOURCE_PR,
        "velocity_export_ready": states["velocity_export_ready"],
        "candidate_selection_resolved": promotion["candidate_selection_resolved"],
        "visual_correspondence_verified": states["visual_correspondence_verified"],
        "pde_validated": states["pde_validated"],
        "per_panel_autoscaling": semantics["per_panel_autoscaling"],
        "projected_streamlines_are_true_3d_streamlines": semantics[
            "projected_streamlines_are_true_3d_streamlines"
        ],
    }


def main() -> int:
    print(json.dumps(audit_shared_frame_visual_evidence_scope(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
