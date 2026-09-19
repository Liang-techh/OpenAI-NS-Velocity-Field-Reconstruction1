"""Fail-closed CR002 audit for ST052-M parent/candidate identity binding.

The live branch now contains the exact PR559 post-processing transform kernel, but
that kernel accepts a caller-supplied parent velocity callable.  This audit keeps
three identities separate:

1. the frozen transform specification/checksum;
2. the intended ST052-M parent metadata;
3. a complete materialized child candidate identity.

It changes no velocity, pressure, forcing, validation protocol, or threshold.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path
from typing import Any

import numpy as np

from openai_ns_reconstruction.st052_local_swirl_transform import (
    DEFAULT_SPEC,
    TRUTH_BOUNDARY as TRANSFORM_TRUTH_BOUNDARY,
    St052LocalSwirlAdapter,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = ROOT / "configs" / "st052m_parent_identity_binding_contract.json"
PRIOR_STAGE_PATH = ROOT / "configs" / "st052m_delivery_stage_accounting_contract.json"


class GovernanceError(RuntimeError):
    """Raised when identity/provenance boundaries are promoted without evidence."""


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GovernanceError(message)


def _diagnostic_unbound_parent_acceptance() -> bool:
    """Confirm that current construction authenticates metadata, not callable bytes.

    This deliberately uses a finite non-ST052 diagnostic callable.  Acceptance is
    evidence about the current adapter interface only; it is never a candidate or
    a scientific field.
    """

    def diagnostic_parent(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        values = np.array([0.123, -0.456, 0.789], dtype=float)
        return np.broadcast_to(values, points.shape).copy()

    try:
        adapter = St052LocalSwirlAdapter(
            diagnostic_parent,
            parent_candidate_id=DEFAULT_SPEC.parent_candidate_id,
            parent_source_head=DEFAULT_SPEC.parent_source_head,
        )
        out = adapter.velocity(0.2, -0.1, 0.3, 0.5)
    except (TypeError, ValueError, RuntimeError):
        return False
    return bool(np.asarray(out).shape == (3,) and np.isfinite(out).all())


def audit(
    root: Path = ROOT,
    *,
    contract_path: Path | None = None,
) -> dict[str, Any]:
    """Audit the live repository against the frozen parent-identity contract."""

    contract = _load(contract_path or (root / "configs" / CONTRACT_PATH.name))
    status = _load(root / "project_status.json")
    constraints = _load(root / "configs" / "constraints.json")
    prior_stage = _load(root / "configs" / PRIOR_STAGE_PATH.name)

    expected_classes = {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    _require(set(contract["allowed_source_classes"]) == expected_classes, "source vocabulary drift")
    _require(
        set(contract["source_classification"].values()) <= expected_classes,
        "invalid source classification",
    )

    snap = contract["snapshot"]
    _require(snap["active_integration_branch"] == "codex/cr001-constraints", "integration branch drift")
    _require(
        snap["active_integration_head"] == "ea2b99903db40775f9397837eed3f26622dd09d3",
        "integration snapshot head drift",
    )
    _require(snap["source_migration_pr"] == 591, "migration PR identity drift")
    _require(
        snap["source_migration_head"] == "fd3929b6950e2821ac41a1930cd0620e43f7f2d6",
        "migration source head drift",
    )
    _require(snap["source_migration_dedicated_conclusion"] == "success", "migration parity CI not green")
    _require(snap["exact_child_source_pr"] == 559, "exact child PR drift")
    _require(
        snap["exact_child_source_head"] == "39b106ad8cb8df2064cabead3a12682089575e74",
        "exact child head drift",
    )
    _require(
        snap["parent_source_head"] == "b3b8bfdbe1077f9ec967d158602951997d81e17d",
        "parent source head drift",
    )
    _require(
        snap["transform_spec_sha256"] == "470103b68fd5ccece5d43d054711c894e40e4cb2f7ce890d52d04c4623bf7c25",
        "registered transform checksum drift",
    )
    _require(DEFAULT_SPEC.sha256() == snap["transform_spec_sha256"], "live transform spec checksum mismatch")
    _require(DEFAULT_SPEC.parent_candidate_id == "ST052-M", "live parent candidate metadata drift")
    _require(DEFAULT_SPEC.parent_source_head == snap["parent_source_head"], "live parent source metadata drift")
    _require(DEFAULT_SPEC.source_child_pr == 559, "live source child PR drift")
    _require(DEFAULT_SPEC.source_child_head == snap["exact_child_source_head"], "live source child head drift")

    migrated = contract["migrated_transform_state"]
    for key in (
        "package_transform_kernel_materialized_on_live_ancestry",
        "package_vectorized_adapter_materialized_on_live_ancestry",
        "transform_spec_save_load_ready",
        "transform_spec_checksum_verified",
        "exact_source_parity_receipt_passed",
    ):
        _require(migrated[key] is True, f"migrated transform state weakened: {key}")
    for key in (
        "full_st052m_parent_evaluator_materialized_by_this_kernel",
        "complete_child_candidate_materialized",
        "complete_child_save_load_ready",
        "experimental_child_velocity_export_ready",
    ):
        _require(migrated[key] is False, f"incomplete child promoted: {key}")

    binding = contract["parent_binding_current_state"]
    _require(binding["intended_parent_candidate_id"] == DEFAULT_SPEC.parent_candidate_id, "parent id drift")
    _require(binding["intended_parent_source_head"] == DEFAULT_SPEC.parent_source_head, "parent head drift")
    for key in (
        "adapter_accepts_external_base_velocity_callable",
        "adapter_checks_parent_candidate_id_string",
        "adapter_checks_parent_source_head_string",
        "arbitrary_finite_callable_with_correct_metadata_can_construct_current_adapter",
    ):
        _require(binding[key] is True, f"current adapter-interface fact drift: {key}")
    for key in (
        "adapter_constructor_has_parent_artifact_sha256_parameter",
        "adapter_constructor_has_parent_callable_digest_parameter",
        "adapter_verifies_parent_artifact_bytes",
        "adapter_verifies_parent_callable_behavior",
        "correct_metadata_strings_are_cryptographic_parent_proof",
        "transform_spec_sha256_is_parent_artifact_sha256",
        "transform_spec_sha256_is_complete_child_candidate_sha256",
        "exact_source_parity_of_migration_implies_all_future_injected_parents_are_equivalent",
    ):
        _require(binding[key] is False, f"unsupported parent-binding claim promoted: {key}")

    init_params = inspect.signature(St052LocalSwirlAdapter.__init__).parameters
    _require("base_velocity" in init_params, "adapter lost external parent callable")
    _require("parent_candidate_id" in init_params, "adapter lost parent metadata check")
    _require("parent_source_head" in init_params, "adapter lost parent source metadata check")
    _require("parent_artifact_sha256" not in init_params, "contract stale: parent artifact binding was added")
    _require("parent_callable_digest" not in init_params, "contract stale: parent callable binding was added")
    _require(
        _diagnostic_unbound_parent_acceptance(),
        "contract stale: arbitrary finite callable no longer passes metadata-only construction",
    )

    semantics = contract["identity_semantics"]
    _require(semantics["callability_without_complete_identity_allowed"] is True, "callability incorrectly forbidden")
    _require(
        semantics["callability_without_complete_identity_may_be_labelled_complete_delivery"] is False,
        "callability laundered into complete delivery",
    )
    _require(semantics["pde_failure_blocks_honest_callable_delivery"] is False, "PDE pending incorrectly blocks delivery")

    required = contract["required_before_complete_child_promotion"]
    for key, value in required.items():
        _require(value is True, f"promotion prerequisite disabled: {key}")

    _require(
        prior_stage["exact_source_child"]["versioned_candidate_id_assigned"] is False,
        "prior stage contract already claims child id",
    )
    _require(
        prior_stage["exact_source_child"]["versioned_candidate_sha_assigned"] is False,
        "prior stage contract already claims child sha",
    )
    _require(
        prior_stage["exact_source_child"]["save_load_contract_registered"] is False,
        "prior stage contract already claims full save/load",
    )
    _require(
        prior_stage["exact_source_child"]["unified_velocity_api_registered"] is False,
        "prior stage contract already claims unified child API",
    )
    _require(
        prior_stage["truth_boundary"]["experimental_child_velocity_export_ready"] is False,
        "prior stage contract promoted experimental export readiness",
    )

    latest = status["latest_st052m_visual_candidate_evidence"]
    _require(latest["exact_child_source_pr"] == 559, "project status child PR drift")
    _require(latest["exact_child_source_head"] == snap["exact_child_source_head"], "project status child head drift")
    _require(latest["production_candidate_selected"] is False, "experimental child production-selected")
    _require(latest["velocity_export_ready"] is False, "experimental child export-readiness promoted")
    _require(latest["visualization_ready"] is False, "experimental child visualization-readiness promoted")
    _require(latest["visual_correspondence_verified"] is False, "visual correspondence promoted")
    _require(latest["pde_validated"] is False, "experimental child PDE-promoted")
    _require(
        status["delivery_shortest_route"] == "materialize_exact_st052m_local_swirl_energy_child_then_render_and_fresh_validate",
        "delivery route drift",
    )
    _require(
        "assign a new versioned candidate identity/SHA plus save-load and unified velocity(x,y,z,t) adapter"
        in status["next_delivery_task"],
        "project status lost full-child identity/save-load requirement",
    )

    canonical = contract["canonical_delivery"]
    _require(status["candidate_family"] == canonical["candidate_family"], "canonical candidate family drift")
    _require(status["candidate_sha256"] == canonical["candidate_sha256"], "canonical candidate sha drift")
    _require(status["velocity_api"] == canonical["velocity_api"], "canonical velocity API drift")
    _require(status["states"]["velocity_export_ready"] is True, "canonical export readiness lost")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _require(status["states"][key] is False, f"canonical truth state promoted: {key}")

    cr001 = contract["canonical_cr001"]
    _require(constraints["nu"] == cr001["nu"], "nu drift")
    _require(constraints["domain"]["physical"] == cr001["physical_domain"], "physical domain drift")
    _require(constraints["domain"]["evaluation_box"] == cr001["evaluation_box"], "evaluation box drift")
    _require(constraints["domain"]["support"] == cr001["support"], "support drift")
    _require(constraints["domain"]["time_interval"] == cr001["time_interval"], "time interval drift")
    _require(constraints["forcing"]["mode"] == cr001["force_mode"], "force mode drift")
    _require(constraints["forcing"]["parameters"] == cr001["force_bounds"], "force bounds drift")
    _require(constraints["nontriviality"]["reference_energy"] == cr001["reference_energy"], "energy target drift")
    _require(
        constraints["nontriviality"]["reference_energy_abs_tolerance"]
        == cr001["reference_energy_abs_tolerance"],
        "energy tolerance drift",
    )
    validation = constraints["validation"]
    _require(validation["seed"] == cr001["validation_seed"], "validation seed drift")
    _require(validation["held_out_points"] == cr001["held_out_points"], "held-out count drift")
    _require(validation["times"] == cr001["validation_times"], "validation times drift")
    _require(validation["derivative_steps"] == cr001["derivative_steps"], "derivative ladder drift")
    _require(
        validation["quadrature_orders_per_axis"] == cr001["quadrature_orders_per_axis"],
        "energy quadrature ladder drift",
    )
    thresholds = validation["thresholds"]
    _require(thresholds["divergence_max"] == cr001["divergence_max"], "divergence max threshold drift")
    _require(thresholds["divergence_L2"] == cr001["divergence_L2"], "divergence L2 threshold drift")
    _require(thresholds["pde_residual_max"] == cr001["pde_residual_max"], "momentum max threshold drift")
    _require(thresholds["pde_residual_L2"] == cr001["pde_residual_L2"], "momentum L2 threshold drift")

    transform_truth = TRANSFORM_TRUTH_BOUNDARY
    for key in (
        "full_st052_parent_materialized_here",
        "complete_candidate_save_load_ready",
        "production_candidate_selected",
        "held_out_pde_residual_evaluated",
        "parent_pde_receipt_transferred",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(transform_truth[key] is False, f"transform truth boundary promoted: {key}")

    truth = contract["truth_boundary"]
    _require(truth["transform_kernel_materialized"] is True, "integrated transform kernel erased")
    for key in (
        "full_parent_materialized_and_bound",
        "complete_child_candidate_materialized",
        "experimental_child_velocity_export_ready",
        "production_candidate_selected",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth[key] is False, f"unsupported truth promotion: {key}")

    mutations = contract["mutation_scope"]
    _require(all(value is False for value in mutations.values()), "governance increment mutates scientific state")

    return {
        "status": "PASS",
        "task_id": contract["task_id"],
        "transform_spec_sha256": DEFAULT_SPEC.sha256(),
        "metadata_only_parent_construction_observed": True,
        "transform_kernel_materialized": True,
        "full_parent_bound": False,
        "complete_child_velocity_export_ready": False,
        "canonical_velocity_export_ready": True,
        "pde_validated": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
