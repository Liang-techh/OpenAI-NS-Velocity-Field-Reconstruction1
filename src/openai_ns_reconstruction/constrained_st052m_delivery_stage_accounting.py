"""Fail-closed governance audit for ST052-M delivery-stage accounting.

This module changes no scientific object.  It checks that a fixed render replay of
an exact source PR cannot be credited as materialization/save-load/API delivery on
the live constrained ancestry, and that canonical CR001/delivery truth remains
unchanged.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = ROOT / "configs" / "st052m_delivery_stage_accounting_contract.json"


class GovernanceError(RuntimeError):
    """Raised when the governed delivery-stage boundary is violated."""


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GovernanceError(message)


def audit(
    root: Path = ROOT,
    *,
    contract_path: Path | None = None,
) -> dict[str, Any]:
    """Audit current repository state against the frozen stage-accounting contract."""

    contract = _load(contract_path or (root / "configs" / CONTRACT_PATH.name))
    status = _load(root / "project_status.json")
    constraints = _load(root / "configs" / "constraints.json")
    delivery = _load(root / "configs" / "delivery_state_contract.json")
    render_scope = _load(root / "configs" / "render_readiness_scope_contract.json")

    expected_classes = {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    _require(set(contract["allowed_source_classes"]) == expected_classes, "source vocabulary drift")
    _require(set(contract["source_classification"].values()) <= expected_classes, "invalid source class")

    snap = contract["snapshot"]
    _require(snap["active_integration_branch"] == "codex/cr001-constraints", "integration branch drift")
    _require(
        snap["active_integration_head"] == "2d97cfaece49d62b163637348da56a545ab66fa7",
        "snapshot base head drift",
    )
    _require(snap["source_child_pr"] == 559, "source child PR drift")
    _require(
        snap["source_child_head"] == "39b106ad8cb8df2064cabead3a12682089575e74",
        "source child identity drift",
    )
    _require(snap["render_pr"] == 583, "render PR drift")
    _require(
        snap["render_head"] == "ee80de11775bd45157d114000378b6ee3422ff46",
        "render head drift",
    )
    _require(snap["render_workflow_conclusion"] == "success", "render receipt is not green")
    _require(snap["standard_tests_conclusion"] == "success", "render standard tests are not green")

    _require(
        status["delivery_shortest_route"] == snap["project_status_delivery_shortest_route"],
        "live delivery route drift",
    )
    next_delivery = status["next_delivery_task"]
    _require(next_delivery.startswith("replay and materialize the exact ST052-M child"), "materialization is no longer first")
    _require("assign a new versioned candidate identity/SHA plus save-load and unified velocity(x,y,z,t) adapter" in next_delivery,
             "live route lost versioned identity/save-load/API requirement")
    _require("then run fixed-camera/fixed-seed 3-D render" in next_delivery,
             "live route lost post-materialization render ordering")

    latest = status["latest_st052m_visual_candidate_evidence"]
    _require(latest["exact_child_source_pr"] == 559, "project status source PR drift")
    _require(latest["exact_child_source_head"] == snap["source_child_head"], "project status source head drift")
    _require(latest["redistribution_kappa"] == 0.05, "redistribution kappa drift")
    _require(latest["taper_tau"] == 0.05, "taper tau drift")
    _require(latest["shoulder_swirl_beta"] == 0.08837490297155456, "energy-root beta drift")
    _require(latest["post_transform_common_scale"] == 1.0, "post-transform scale drift")
    _require(latest["production_candidate_selected"] is False, "experimental child promoted")
    _require(latest["velocity_export_ready"] is False, "experimental child export-readiness promoted")
    _require(latest["visualization_ready"] is False, "experimental child visualization-readiness promoted")
    _require(latest["visual_correspondence_verified"] is False, "visual correspondence promoted")
    _require(latest["pde_validated"] is False, "experimental child PDE state promoted")

    child = contract["exact_source_child"]
    for key in (
        "materialized_on_live_constrained_ancestry",
        "versioned_candidate_id_assigned",
        "versioned_candidate_sha_assigned",
        "save_load_contract_registered",
        "unified_velocity_api_registered",
        "velocity_export_ready",
        "production_candidate_selected",
    ):
        _require(child[key] is False, f"source child stage promoted: {key}")

    render = contract["render_evidence"]
    _require(render["source_checkout_is_exact_pr559_head"] is True, "render is not bound to exact source")
    _require(render["source_checkout_is_separate_from_render_branch_ancestry"] is True,
             "cross-ancestry source checkout semantics lost")
    _require(render["render_artifact_exists"] is True and render["render_ci_success"] is True,
             "positive render receipt lost")
    _require(render["openai_numeric_image_target_used"] is False, "hidden/public image numeric target introduced")
    _require(render["visual_acceptance_threshold_used"] is False, "visual threshold introduced")
    for key in (
        "counts_as_live_ancestry_materialization",
        "counts_as_child_velocity_export_ready",
        "counts_as_canonical_candidate_replacement",
        "counts_as_visual_correspondence_verified",
        "counts_as_pde_validation",
    ):
        _require(render[key] is False, f"render evidence over-promoted: {key}")

    stages = contract["delivery_stage_accounting"]
    _require(stages["source_child_defined"] is True, "source child definition lost")
    _require(stages["target_free_exact_source_render_available"] is True, "render evidence lost")
    _require(stages["render_before_materialization_is_allowed_as_diagnostic"] is True,
             "diagnostic render incorrectly forbidden")
    for key in (
        "materialization_stage_complete",
        "identity_and_save_load_stage_complete",
        "post_materialization_render_stage_complete",
        "fresh_pde_validation_stage_complete",
        "render_before_materialization_satisfies_materialization_prerequisite",
        "render_before_materialization_satisfies_post_materialization_render_stage",
        "delivery_route_complete",
    ):
        _require(stages[key] is False, f"delivery stage laundering: {key}")

    transfer = contract["later_materialized_child_evidence_transfer"]
    _require(transfer["automatic_transfer_allowed"] is False, "render evidence made automatically inheritable")
    _require(transfer["same_parameters_alone_are_sufficient"] is False, "parameters alone made identity proof")
    _require(transfer["same_source_commit_alone_is_sufficient"] is False, "source commit alone made delivery proof")
    _require(transfer["same_visual_appearance_alone_is_sufficient"] is False, "appearance made identity proof")
    required_transfer = set(transfer["required_before_transfer"])
    _require(
        {
            "explicit_functional_representation_binding_to_pr559_exact_source_identity",
            "new_versioned_candidate_id_and_sha",
            "deterministic_save_load_replay",
            "unified_velocity_api_replay",
            "render_or_equivalence_replay_on_the_materialized_identity",
        } <= required_transfer,
        "later-child evidence-transfer prerequisites weakened",
    )

    canonical = contract["canonical_delivery"]
    _require(status["candidate_family"] == canonical["candidate_family"], "canonical family drift")
    _require(status["candidate_sha256"] == canonical["candidate_sha256"], "canonical candidate SHA drift")
    _require(status["velocity_api"] == canonical["velocity_api"], "canonical velocity API drift")
    _require(status["states"]["velocity_export_ready"] is True, "canonical export readiness lost")
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _require(status["states"][key] is canonical[key] is False, f"canonical truth state promoted: {key}")

    _require(
        delivery["states"]["velocity_export_ready"]["positive_evidence"]
        == ["versioned_candidate_artifact", "unified_velocity_api", "finite_velocity_output", "reproducible_export"],
        "velocity_export_ready evidence contract drift",
    )
    _require("visualization_ready" in delivery["states"]["velocity_export_ready"]["does_not_imply"],
             "velocity/export and visualization states coupled")
    _require(render_scope["readiness_semantics"]["render_artifact_ready"]["does_not_imply_visualization_ready"] is True,
             "render artifact improperly implies visualization readiness")
    _require(render_scope["readiness_semantics"]["render_artifact_ready"]["does_not_imply_pde_validated"] is True,
             "render artifact improperly implies PDE validation")

    cr = contract["canonical_cr001"]
    _require(constraints["nu"] == cr["nu"], "nu drift")
    _require(constraints["domain"]["physical"] == cr["physical_domain"], "physical domain drift")
    _require(constraints["domain"]["evaluation_box"] == cr["evaluation_box"], "evaluation box drift")
    _require(constraints["domain"]["support"] == cr["support"], "support drift")
    _require(constraints["domain"]["time_interval"] == cr["time_interval"], "time interval drift")
    _require(constraints["forcing"]["mode"] == cr["force_mode"], "forcing family drift")
    _require(constraints["forcing"]["parameters"] == cr["force_bounds"], "forcing bounds drift")
    _require(constraints["nontriviality"]["reference_energy"] == cr["reference_energy"], "reference energy drift")
    _require(constraints["nontriviality"]["reference_energy_abs_tolerance"] == cr["reference_energy_abs_tolerance"],
             "reference-energy tolerance drift")
    validation = constraints["validation"]
    _require(validation["seed"] == cr["validation_seed"], "validation seed drift")
    _require(validation["held_out_points"] == cr["held_out_points"], "held-out point count drift")
    _require(validation["times"] == cr["validation_times"], "validation-time drift")
    _require(validation["derivative_steps"] == cr["derivative_steps"], "derivative ladder drift")
    _require(validation["quadrature_orders_per_axis"] == cr["quadrature_orders_per_axis"], "quadrature ladder drift")
    thresholds = validation["thresholds"]
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds[key] == cr[key], f"CR001 threshold drift: {key}")

    truth = contract["truth_boundary"]
    _require(truth["canonical_velocity_export_ready"] is True, "canonical export state corrupted")
    for key in (
        "experimental_child_velocity_export_ready",
        "production_candidate_selected",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth[key] is False, f"truth-state promotion in contract: {key}")

    _require(not any(contract["mutation_scope"].values()), "governance increment mutates scientific state")

    return {
        "task_id": contract["task_id"],
        "status": "pass",
        "source_child_head": snap["source_child_head"],
        "render_head": snap["render_head"],
        "render_receipt_available": True,
        "live_materialization_complete": False,
        "experimental_child_velocity_export_ready": False,
        "delivery_route_complete": False,
        "canonical_candidate_unchanged": True,
        "cr001_unchanged": True,
    }


def main() -> None:
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
