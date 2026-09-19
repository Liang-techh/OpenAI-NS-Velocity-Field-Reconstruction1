"""Fail-closed CR002 audit for the fresh #714 ST052-M tip-path protocol.

This module governs *what kind* of evidence a future exact-head #714 result may
support.  It does not execute the expensive ODE replay and it does not change a
candidate.  The key distinction is between data-disjoint post-freeze
trajectory generalization and a stronger blind/global visual-validation claim.
"""
from __future__ import annotations

import ast
import json
from pathlib import Path
from typing import Any


CONTRACT_REL = Path("configs/st052m_tip_holdout_protocol_scope.json")
HOLDOUT_REL = Path("experiments/root_st052/agent7_st052m_tip_odd_poloidal_holdout.py")
SCREEN_REL = Path("experiments/root_st052/agent7_st052m_tip_odd_poloidal_screen.py")
LOCALIZATION_REL = Path("experiments/root_st052/agent7_st052m_radial_split_localization.py")
CR001_REL = Path("configs/constraints.json")

EXPECTED_HEAD = "a352336d4b18771620fe0f5f554c3f3a7fb9d157"
EXPECTED_701_HEAD = "e3ecc8654f2c9a2324265bbd01df596b4e503404"
ALLOWED_SOURCE_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


class GovernanceError(AssertionError):
    """Raised when a governed truth boundary drifts."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GovernanceError(message)


def repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def load_contract(root: Path | None = None) -> dict[str, Any]:
    root = repository_root() if root is None else Path(root)
    return json.loads((root / CONTRACT_REL).read_text(encoding="utf-8"))


def _literal_assignments(source: str) -> dict[str, Any]:
    tree = ast.parse(source)
    values: dict[str, Any] = {}
    for node in tree.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        target = None
        value = None
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            target = node.targets[0].id
            value = node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            target = node.target.id
            value = node.value
        if target is None or value is None:
            continue
        try:
            values[target] = ast.literal_eval(value)
        except (ValueError, TypeError):
            pass
    return values


def audit_contract(contract: dict[str, Any]) -> None:
    upstream = contract["audited_upstream"]
    _require(upstream["pull_request"] == 714, "wrong audited PR")
    _require(upstream["head_sha"] == EXPECTED_HEAD, "#714 head identity drift")
    _require(upstream["frozen_candidate_source_pr"] == 701, "wrong frozen candidate source")
    _require(upstream["frozen_candidate_source_head"] == EXPECTED_701_HEAD, "#701 head identity drift")
    _require(upstream["selection_data_governance_pr"] == 712, "selection-governance linkage drift")
    _require(upstream["prior_band_localization_pr"] == 692, "band-localization provenance drift")
    _require(upstream["prior_band_localization_classification"] == "band_local", "band diagnosis drift")

    source = contract["source_classification"]
    _require(set(source["allowed_classes"]) == ALLOWED_SOURCE_CLASSES, "source vocabulary drift")
    for item in source["items"]:
        _require(item["class"] in ALLOWED_SOURCE_CLASSES, f"invalid source class: {item}")
    by_item = {item["item"]: item["class"] for item in source["items"]}
    _require(
        by_item["fresh_radii_angles_z_bands_solver_tolerances_and_directional_gates"]
        == "autonomous_design",
        "fresh holdout protocol must remain autonomous design",
    )
    _require(
        by_item["hidden_openai_path_seeds_camera_frame_time_mapping_coefficients_and_numeric_field"]
        == "pending_unknown",
        "hidden OpenAI numerical data cannot be promoted",
    )

    sep = contract["selection_evaluation_separation"]
    for key in (
        "candidate_parameter_frozen_before_fresh_path_evaluation",
        "alpha_is_derived_only_from_historical_development_protocol",
        "fresh_protocol_is_data_disjoint_from_alpha_calibration",
        "fresh_protocol_is_post_freeze_directional_generalization_evidence_if_frozen_gate_passes",
    ):
        _require(sep[key] is True, f"selection/evaluation separation lost: {key}")
    _require(sep["fresh_path_count"] == 64, "fresh path count drift")
    _require(sep["development_path_count"] == 48, "development path count drift")
    _require(sep["exact_seed_overlap_with_development_required"] == 0, "seed-overlap guard relaxed")
    for key in (
        "fresh_paths_used_to_fit_alpha",
        "fresh_results_feed_back_into_alpha",
        "parameter_scan_or_retuning_on_fresh_paths",
    ):
        _require(sep[key] is False, f"fresh-data leakage allowed: {key}")

    scope = contract["protocol_independence_scope"]
    _require(scope["protocol_designed_after_prior_band_local_diagnosis"] is True, "prior-diagnosis timing erased")
    _require(scope["protocol_reuses_prior_shoulder_tip_semantic_partition"] is True, "band semantics provenance erased")
    for key in (
        "protocol_is_blind_to_prior_band_diagnosis",
        "protocol_is_random_statistical_holdout",
        "protocol_is_global_visual_validation",
        "protocol_is_cr001_pde_acceptance_sample",
        "fresh_seed_disjointness_alone_establishes_visual_correspondence",
        "fresh_seed_disjointness_alone_establishes_pde_validation",
    ):
        _require(scope[key] is False, f"holdout scope over-promoted: {key}")

    identity = contract["candidate_identity_scope"]
    _require(identity["frozen_candidate_is_procedurally_reconstructed_from_exact_701_source"] is True, "procedural freeze erased")
    _require(identity["alpha_is_recomputed_deterministically_from_frozen_historical_development_data"] is True, "alpha reconstruction drift")
    for key in (
        "standalone_serialized_714_candidate_identity_created",
        "velocity_export_ready_promoted_by_714",
        "canonical_velocity_changed_by_714",
        "saved_velocity_changed_by_714",
    ):
        _require(identity[key] is False, f"candidate identity/delivery over-promotion: {key}")

    claims = contract["claim_states"]
    # The upstream scientific result is unresolved at this governance commit.
    for key, value in claims.items():
        _require(value is False, f"premature claim promotion: {key}")

    lock = contract["cr001_lock"]
    _require(lock["nu"] == 0.01, "nu drift")
    _require(lock["physical_domain"] == "R^3", "physical domain drift")
    _require(lock["evaluation_box"] == [[-2, 2], [-2, 2], [-2, 2]], "box drift")
    _require(lock["support"] == "r < 2 and abs(z) < 2", "support drift")
    _require(lock["time_interval"] == [0.25, 0.75], "time interval drift")
    _require(lock["forcing_mode"] == "restricted_two_parameter_family", "forcing family drift")
    _require(lock["forcing_bounds"] == {"a": [0.0, 10.0], "c": [0.0, 10.0]}, "forcing bounds drift")
    _require(lock["reference_energy"] == 1.0 and lock["reference_energy_abs_tolerance"] == 0.001, "energy normalization drift")
    _require(lock["validation_seed"] == 914027 and lock["held_out_points"] == 4096, "formal validation sample drift")
    _require(lock["derivative_steps"] == [0.02, 0.01, 0.005], "derivative ladder drift")
    _require(lock["quadrature_orders_per_axis"] == [24, 48, 96], "quadrature ladder drift")
    _require(lock["divergence_max"] == 1e-5 and lock["divergence_L2"] == 1e-5, "divergence gate drift")
    _require(lock["pde_residual_max"] == 1e-3 and lock["pde_residual_L2"] == 1e-3, "momentum gate drift")
    _require(lock["residual_defined_or_pointwise_free_force_allowed"] is False, "free-force shortcut enabled")
    _require(lock["amplitude_collapse_success_allowed"] is False, "amplitude-collapse shortcut enabled")
    _require(lock["threshold_relaxation_allowed"] is False, "threshold relaxation enabled")


def audit_live(root: Path | None = None, contract: dict[str, Any] | None = None) -> dict[str, Any]:
    root = repository_root() if root is None else Path(root)
    contract = load_contract(root) if contract is None else contract
    audit_contract(contract)

    holdout = (root / HOLDOUT_REL).read_text(encoding="utf-8")
    screen = (root / SCREEN_REL).read_text(encoding="utf-8")
    localization = (root / LOCALIZATION_REL).read_text(encoding="utf-8")
    cr001 = json.loads((root / CR001_REL).read_text(encoding="utf-8"))

    h = _literal_assignments(holdout)
    s = _literal_assignments(screen)

    _require(h["PREREG_ISSUE"] == 713, "#714 preregistration identity drift")
    _require(h["SOURCE_PARENT_PR"] == 701 and h["SOURCE_PARENT_HEAD"] == EXPECTED_701_HEAD, "#714 source parent drift")
    _require(h["SELECTION_GOVERNANCE_PR"] == 712, "#714 no longer binds #712 selection governance")
    _require(h["HOLDOUT_RADII"] == (0.50, 0.80, 1.10, 1.40), "fresh radii drift")
    _require(h["HOLDOUT_BANDS"] == (("shoulder", 0.75), ("tip", 1.25)), "fresh band geometry drift")
    _require(h["HOLDOUT_ANGLES"] == 4 and h["OUTPUT_SAMPLES"] == 37, "fresh path discretization drift")
    _require(h["CHECKPOINT_TIMES"] == (0.25, 0.375, 0.50, 0.625, 0.75), "fresh checkpoint ladder drift")
    _require(h["SOLVER_METHOD"] == "DOP853", "fresh integrator drift")
    _require(h["SOLVER_RTOL"] == 5e-10 and h["SOLVER_ATOL"] == 5e-12 and h["SOLVER_MAX_STEP"] == 0.008, "fresh solver tolerance drift")

    truth = h["TRUTH"]
    for key in (
        "source_701_candidate_retuned",
        "heldout_evaluation_feedback_into_alpha",
        "parameter_grid_scan_performed",
        "optimization_performed",
        "pressure_or_force_changed",
        "held_out_pde_residual_evaluated",
        "parent_pde_receipt_transferred",
        "public_image_numeric_target_used",
        "pixel_similarity_objective_used",
        "visual_acceptance_threshold_defined",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth[key] is False, f"#714 truth boundary promoted: {key}")

    _require("if _development_seed_overlap_count(seeds) != 0:" in holdout, "fresh exact-overlap guard removed")
    _require("dev_times, dev_positions, dev_metadata = screen._integrate(child_fn)" in holdout, "historical alpha-calibration path changed")
    _require("screen.derive_alpha_from_frozen_midpoint(" in holdout, "#701 frozen alpha rule not reused")
    _require("candidate_fn = lambda p, t: screen.corrected_velocity(child_fn, p, t, alpha)" in holdout, "frozen candidate construction drift")
    _require("np.random" not in holdout and "random." not in holdout, "fresh protocol became random without governance update")

    # Bind the temporal ordering/protocol semantics to the stacked source itself.
    _require(s["SOURCE_PARENT_PR"] == 692, "#701 no longer descends from band-localization PR #692")
    _require("import agent7_st052m_radial_split_localization as loc" in screen, "#701 band-localization dependency removed")
    _require("tip-band-local odd-poloidal" in screen, "#701 no longer identifies the correction as tip-band-local")
    _require("classification = \"band_local\"" in localization, "#692 band-local classification route missing")
    _require("shoulder/tip window" in localization, "#692 band-local routing semantics missing")

    # Bind the governance snapshot back to canonical CR001 rather than creating a new threshold set.
    lock = contract["cr001_lock"]
    _require(cr001["nu"] == lock["nu"], "contract/live nu mismatch")
    _require(cr001["domain"]["physical"] == lock["physical_domain"], "contract/live physical domain mismatch")
    _require(cr001["domain"]["evaluation_box"] == lock["evaluation_box"], "contract/live evaluation box mismatch")
    _require(cr001["domain"]["support"] == lock["support"], "contract/live support mismatch")
    _require(cr001["domain"]["time_interval"] == lock["time_interval"], "contract/live time mismatch")
    _require(cr001["forcing"]["mode"] == lock["forcing_mode"], "contract/live forcing mode mismatch")
    _require(cr001["forcing"]["parameters"] == lock["forcing_bounds"], "contract/live forcing bounds mismatch")
    _require(cr001["validation"]["seed"] == lock["validation_seed"], "contract/live validation seed mismatch")
    _require(cr001["validation"]["held_out_points"] == lock["held_out_points"], "contract/live held-out count mismatch")
    _require(cr001["validation"]["derivative_steps"] == lock["derivative_steps"], "contract/live derivative ladder mismatch")
    _require(cr001["validation"]["quadrature_orders_per_axis"] == lock["quadrature_orders_per_axis"], "contract/live quadrature ladder mismatch")
    thresholds = cr001["validation"]["thresholds"]
    _require(thresholds["divergence_max"] == lock["divergence_max"] and thresholds["divergence_L2"] == lock["divergence_L2"], "contract/live divergence gate mismatch")
    _require(thresholds["pde_residual_max"] == lock["pde_residual_max"] and thresholds["pde_residual_L2"] == lock["pde_residual_L2"], "contract/live momentum gate mismatch")
    _require("No residual-dependent basis or pointwise free force" in cr001["forcing"]["restriction"], "canonical no-free-force rule disappeared")

    return {
        "status": "PASS",
        "audited_pr": 714,
        "fresh_path_count": contract["selection_evaluation_separation"]["fresh_path_count"],
        "exact_seed_overlap_required": 0,
        "post_freeze_directional_generalization_eligible_if_upstream_gate_passes": True,
        "blind_global_visual_validation": False,
        "cr001_pde_acceptance_sample": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit_live(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
