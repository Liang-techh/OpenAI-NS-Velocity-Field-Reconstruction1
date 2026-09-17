"""Govern evidence use for compact-quartic Phi(1,0) blend selection.

The materialized family may be evaluated, serialized, and exported for any caller-
declared blend weight in [0,1].  That capability is intentionally separate from a
repository-level decision that one weight is the visual candidate, from formal PDE
acceptance, and from any paper-exact/OpenAI-field identity claim.

This module is governance only.  It changes no velocity value, temporal schedule,
support transform, force, pressure, optimizer, threshold, or canonical binding.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

EXPECTED_SCHEMA = "eq45_phi10_blend_selection_evidence_scope_v1"
EXPECTED_RECIPE_TASK = "CR011-EQ45-PHI10-COMPACT-QUARTIC-BLEND-DELIVERY-CAPSULE-018"
EXPECTED_RECIPE_HEAD = "42c2a43e3a60fecec3302d4dfa7d6f61c787dabb"
EXPECTED_RECIPE_BASE = "37be0d6ee963cd9b98b8853df00be731d8e28bf4"
EXPECTED_SUPPORTED_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
EXPECTED_SOURCE_VOCABULARY = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
EXPECTED_EVIDENCE = {
    "delivery_capsule": (190, EXPECTED_RECIPE_HEAD, "ancestry_consumed"),
    "temporal_derivative": (187, EXPECTED_RECIPE_BASE, "ancestry_consumed"),
    "target_free_morphology": (
        184,
        "fec36f21f5ee8945385b7070f19ed0127885c2d4",
        "open_unconsumed_sibling",
    ),
    "restricted_force_pde_side": (
        188,
        "fc6ad7b2f0c0a892c1197718db7b225066a60d8e",
        "open_unconsumed_sibling",
    ),
    "spectral_fingerprint": (
        189,
        "0286a04ecb541ba02b247c44545044b3c7298691",
        "open_unconsumed_sibling",
    ),
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


def audit_phi10_blend_selection_evidence_scope(
    *,
    contract_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
    delivery_state_path: str | Path | None = None,
    recipe_path: str | Path | None = None,
) -> dict[str, Any]:
    """Fail closed if blend evidence is promoted beyond its registered scope."""
    root = Path(__file__).resolve().parents[2]
    contract = _load(
        contract_path
        or root / "configs/eq45_phi10_blend_selection_evidence_scope.json"
    )
    constraints = _load(constraints_path or root / "configs/constraints.json")
    delivery = _load(delivery_state_path or root / "configs/delivery_state_contract.json")
    recipe = _load(
        recipe_path
        or root / "artifacts/constrained/eq45_supported_phi10_compact_quartic_blend_recipe.json"
    )

    _require(contract.get("schema") == EXPECTED_SCHEMA, "blend-selection governance schema drifted")

    canonical_vocab = set(delivery["classification_vocabulary"])
    local_vocab = set(contract.get("source_vocabulary", []))
    _require(canonical_vocab == EXPECTED_SOURCE_VOCABULARY, "delivery-state source vocabulary drifted")
    _require(local_vocab == EXPECTED_SOURCE_VOCABULARY, "blend-selection source vocabulary drifted")
    classifications = contract.get("source_classification", {})
    _require(set(classifications.values()) <= EXPECTED_SOURCE_VOCABULARY, "unknown source classification")
    expected_classifications = {
        "callable_saveable_loadable_3d_time_varying_velocity": "user_requirement",
        "independent_velocity_visual_pde_and_exactness_states": "user_requirement",
        "compact_quartic_phi10_temporal_interpolation": "autonomous_design",
        "caller_declared_blend_weight": "autonomous_design",
        "target_free_morphology_screen": "autonomous_design",
        "restricted_force_pressure_free_pde_side_screen": "autonomous_design",
        "public_velocity_temporal_derivative_audit": "autonomous_design",
        "finite_box_spectral_fingerprint": "autonomous_design",
        "public_openai_observable_correspondence_for_any_blend_member": "pending_unknown",
        "paper_exact_or_hidden_openai_field_identity": "pending_unknown",
    }
    _require(classifications == expected_classifications, "blend-selection source classification drifted")

    binding = contract["cr001_binding"]
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]
    _require(_same_float(binding["nu"], constraints["nu"]), "nu drifted")
    _require(binding["physical_domain"] == domain["physical"], "physical domain drifted")
    _require(binding["evaluation_box"] == domain["evaluation_box"], "evaluation box drifted")
    _require(binding["support"] == domain["support"], "support contract drifted")
    _require(binding["time_interval"] == domain["time_interval"], "time interval drifted")
    _require(binding["forcing_mode"] == forcing["mode"], "forcing mode drifted")
    _require(binding["forcing_parameters"] == sorted(forcing["parameters"]), "forcing parameter set drifted")
    _require(
        "No residual-dependent basis or pointwise free force" in forcing["restriction"],
        "restricted forcing no-free-force guard drifted",
    )
    _require(_same_float(binding["reference_energy_time"], nontriviality["reference_time"]), "energy time drifted")
    _require(_same_float(binding["reference_energy"], nontriviality["reference_energy"]), "energy target drifted")
    _require(
        _same_float(binding["reference_energy_abs_tolerance"], nontriviality["reference_energy_abs_tolerance"]),
        "energy tolerance drifted",
    )
    _require(int(binding["held_out_points"]) == int(validation["held_out_points"]), "held-out count drifted")
    _require(
        [float(x) for x in binding["derivative_steps"]]
        == [float(x) for x in validation["derivative_steps"]],
        "derivative ladder drifted",
    )
    _require(binding["registered_norms"] == validation["norms"], "registered residual norms drifted")
    threshold_pairs = (
        ("divergence_max_threshold", "divergence_max"),
        ("divergence_L2_threshold", "divergence_L2"),
        ("pde_residual_max_threshold", "pde_residual_max"),
        ("pde_residual_L2_threshold", "pde_residual_L2"),
    )
    for local_key, registered_key in threshold_pairs:
        _require(
            _same_float(binding[local_key], thresholds[registered_key]),
            f"{registered_key} threshold drifted",
        )

    family = contract["family_binding"]
    _require(family["delivery_recipe_task_id"] == EXPECTED_RECIPE_TASK, "delivery recipe task drifted")
    _require(family["delivery_recipe_base_head"] == EXPECTED_RECIPE_BASE, "delivery recipe base drifted")
    _require(family["delivery_recipe_implementation_head"] == EXPECTED_RECIPE_HEAD, "delivery recipe head drifted")
    _require(family["base_supported_sha256"] == EXPECTED_SUPPORTED_SHA, "supported parent identity drifted")
    _require(recipe["task_id"] == EXPECTED_RECIPE_TASK, "live recipe task drifted")
    _require(recipe["base_dependency_head"] == EXPECTED_RECIPE_BASE, "live recipe dependency drifted")
    _require(recipe["base_supported_sha256"] == EXPECTED_SUPPORTED_SHA, "live recipe supported parent drifted")
    representation = recipe["representation"]
    _require(family["family"] == representation["family"] == "phi", "blend family drifted")
    _require(family["index"] == representation["index"] == [1, 0], "blend mode index drifted")
    _require(family["blend_formula"] == representation["blend_formula"], "blend formula drifted")
    _require(family["blend_weight_interval"] == representation["blend_weight_interval"] == [0.0, 1.0], "blend interval drifted")
    _require(
        family["representative_candidate_sha256"] == representation["representative_candidate_sha256"],
        "representative candidate identities drifted",
    )
    _require(family["component_order"] == recipe["delivery"]["component_order"] == ["u", "v", "w"], "component order drifted")
    _require(family["grid_layout"] == recipe["delivery"]["grid_layout"], "grid layout drifted")
    _require(family["blend_weight_selected"] is False, "governed family must remain unselected")
    _require(representation["blend_weight_selected"] is False, "delivery recipe selected a blend weight")
    _require(family["candidate_selection_resolved"] is False, "candidate selection must remain unresolved")
    _require(recipe["status"]["candidate_selection_resolved"] is False, "delivery recipe resolved candidate selection")
    _require(family["caller_declared_blend_weight"] is True, "caller-declared blend semantics drifted")
    _require(recipe["status"]["caller_declared_blend_weight"] is True, "delivery recipe caller-declared semantics drifted")

    evidence = contract["evidence_ledger"]
    _require(set(evidence) == set(EXPECTED_EVIDENCE), "evidence ledger membership drifted")
    for name, (pr, head, status) in EXPECTED_EVIDENCE.items():
        item = evidence[name]
        _require(int(item["pr"]) == pr, f"{name} PR drifted")
        _require(item["head"] == head, f"{name} head drifted")
        _require(item["status"] == status, f"{name} consumption state drifted")
        _require(item["may_select_blend_weight"] is False, f"{name} may not select blend weight")
    _require(evidence["restricted_force_pde_side"]["formal_cr001_pde_gate"] is False, "PDE-side routing diagnostic became formal gate")
    for name in ("delivery_capsule", "target_free_morphology", "spectral_fingerprint"):
        _require(
            evidence[name]["may_promote_visual_correspondence_verified"] is False,
            f"{name} may not promote visual correspondence",
        )
    for name in ("delivery_capsule", "temporal_derivative", "restricted_force_pde_side", "spectral_fingerprint"):
        _require(evidence[name]["may_promote_pde_validated"] is False, f"{name} may not promote PDE validity")
    _require(
        recipe["sibling_evidence"]["restricted_force_status"] == "open_unconsumed"
        and recipe["sibling_evidence"]["spectral_fingerprint_status"] == "open_unconsumed",
        "delivery recipe sibling-consumption boundary drifted",
    )
    _require(recipe["sibling_evidence"]["may_select_blend_weight"] is False, "delivery recipe sibling evidence selected lambda")

    selection = contract["selection_contract"]
    required_true = (
        "caller_declared_weights_remain_research_candidates",
        "joint_tradeoff_review_may_occur_after_evidence_replay",
        "public_observable_comparison_required_for_visual_correspondence_selection",
        "independent_visual_resolution_stability_required_for_visual_correspondence_selection",
        "selected_member_must_serialize_blend_weight_and_candidate_sha256",
        "selected_member_requires_candidate_specific_revalidation",
        "formal_pde_gate_remains_independent",
    )
    for key in required_true:
        _require(selection[key] is True, f"{key} must remain true")
    required_false = (
        "blend_weight_selected",
        "candidate_selection_resolved",
        "formal_pde_pass_required_for_candidate_local_save_load_export",
        "formal_pde_pass_required_for_truth_bounded_visualization_candidate",
        "target_free_morphology_alone_may_select_blend_weight",
        "restricted_force_pde_side_diagnostic_alone_may_select_blend_weight",
        "temporal_derivative_audit_alone_may_select_blend_weight",
        "spectral_fingerprint_alone_may_select_blend_weight",
        "green_ci_or_serialization_alone_may_select_blend_weight",
    )
    for key in required_false:
        _require(selection[key] is False, f"{key} must remain false")

    states = contract["truth_states"]
    expected_states = {
        "velocity_export_ready": True,
        "visualization_candidate_only": True,
        "blend_weight_selected": False,
        "candidate_selection_resolved": False,
        "physical_support_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    _require(states == expected_states, "blend-family truth-state boundary drifted")

    required_forbidden = {
        "caller_declared_blend_weight=>repository_selected_blend_weight",
        "target_free_morphology_improvement=>visual_correspondence_verified",
        "restricted_force_pde_side_ordering=>pde_validated",
        "temporal_derivative_convergence=>pde_validated",
        "spectral_softness_or_resolution_stability=>visual_correspondence_verified",
        "green_ci_or_serialization=>candidate_selection_resolved",
        "single_evidence_lane=>blend_weight_selected",
        "visual_similarity=>pde_validated_or_paper_exact_or_openai_field_identified",
        "pde_failed_or_pending=>velocity_export_not_allowed",
    }
    _require(set(contract["forbidden_inferences"]) == required_forbidden, "forbidden inference set drifted")

    return {
        "base_supported_sha256": EXPECTED_SUPPORTED_SHA,
        "representative_candidate_sha256": family["representative_candidate_sha256"],
        "velocity_export_ready": states["velocity_export_ready"],
        "blend_weight_selected": states["blend_weight_selected"],
        "candidate_selection_resolved": states["candidate_selection_resolved"],
        "visual_correspondence_verified": states["visual_correspondence_verified"],
        "pde_validated": states["pde_validated"],
        "unconsumed_sibling_prs": sorted(
            int(item["pr"])
            for item in evidence.values()
            if item["status"] == "open_unconsumed_sibling"
        ),
    }


def main() -> int:
    print(json.dumps(audit_phi10_blend_selection_evidence_scope(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
