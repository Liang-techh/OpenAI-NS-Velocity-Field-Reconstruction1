"""Fail-closed CR002 audit for #740 ST052-M reshaped-child identity scope.

PR #740 makes a new p=9,m=4 tip-envelope velocity formula procedurally
evaluable and deterministically rederives its scalar coefficient from historical
development paths.  That is not yet the same delivery claim as a standalone,
immutable, saved/loadable candidate identity.  This audit keeps those two facts
separate before a future path/render result is called post-freeze evidence.

No velocity, pressure, forcing, scientific threshold, or public API is changed.
"""
from __future__ import annotations

import ast
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

CONTRACT_REL = Path("configs/st052m_tip_reshaped_child_identity_scope.json")
REPLAY_REL = Path("experiments/root_st052/agent7_st052m_tip_reshaped_nonlinear_replay.py")
CR001_REL = Path("configs/constraints.json")

EXPECTED_LIVE_HEAD = "7b32d9c0fa46fe369d0fc636eb5a7530aadebe04"
EXPECTED_740_HEAD = "21019df839557e64b8d8689f68cdf8fc685a3a36"
EXPECTED_738_HEAD = "8072b099b8886e7c8716347ef2a60066d5985206"
ALLOWED_SOURCE_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


class GovernanceError(AssertionError):
    """Raised when the governed identity/truth boundary drifts."""


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


def _audit_cr001(root: Path, lock: dict[str, Any]) -> None:
    canonical = json.loads((root / CR001_REL).read_text(encoding="utf-8"))
    _require(canonical["nu"] == lock["nu"] == 0.01, "CR001 viscosity drift")
    _require(canonical["domain"]["physical"] == lock["physical_domain"] == "R^3", "CR001 physical domain drift")
    _require(canonical["domain"]["evaluation_box"] == lock["evaluation_box"], "CR001 evaluation box drift")
    _require(canonical["domain"]["support"] == lock["support"], "CR001 support drift")
    _require(canonical["domain"]["time_interval"] == lock["time_interval"], "CR001 time interval drift")
    _require(canonical["forcing"]["mode"] == lock["forcing_mode"], "CR001 forcing mode drift")
    _require(canonical["forcing"]["parameters"] == lock["forcing_bounds"], "CR001 forcing bounds drift")
    _require(
        "No residual-dependent basis or pointwise free force" in canonical["forcing"]["restriction"],
        "CR001 free-force prohibition drift",
    )
    nontrivial = canonical["nontriviality"]
    _require(nontrivial["reference_energy"] == lock["reference_energy"], "CR001 reference energy drift")
    _require(
        nontrivial["reference_energy_abs_tolerance"] == lock["reference_energy_abs_tolerance"],
        "CR001 energy tolerance drift",
    )
    _require("reject collapsed candidates" in nontrivial["enforcement"], "CR001 collapse guard drift")
    validation = canonical["validation"]
    _require(validation["seed"] == lock["validation_seed"], "CR001 validation seed drift")
    _require(validation["held_out_points"] == lock["held_out_points"], "CR001 held-out count drift")
    _require(validation["times"] == lock["validation_times"], "CR001 validation times drift")
    _require(validation["derivative_steps"] == lock["derivative_steps"], "CR001 derivative ladder drift")
    _require(
        validation["quadrature_orders_per_axis"] == lock["quadrature_orders_per_axis"],
        "CR001 quadrature ladder drift",
    )
    thresholds = validation["thresholds"]
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds[key] == lock[key], f"CR001 {key} drift")
    _require(
        "changing thresholds requires a new experiment version" in validation["failure_policy"],
        "CR001 threshold-relaxation policy drift",
    )
    _require(lock["residual_defined_or_pointwise_free_force_allowed"] is False, "contract allows free force")
    _require(lock["amplitude_collapse_success_allowed"] is False, "contract allows amplitude collapse")
    _require(lock["threshold_relaxation_allowed"] is False, "contract allows threshold relaxation")


def audit(
    root: Path | None = None,
    contract: dict[str, Any] | None = None,
    replay_source: str | None = None,
) -> dict[str, Any]:
    root = repository_root() if root is None else Path(root)
    contract = load_contract(root) if contract is None else deepcopy(contract)

    _require(
        contract["contract_id"] == "cr002-st052m-tip-reshaped-child-identity-scope-v1",
        "unexpected contract id",
    )
    _require(contract["status"] == "governance_scope_for_exact_740_head", "unexpected contract status")
    _require(
        contract["live_integration_authority"]
        == {"branch": "codex/cr001-constraints", "head_sha": EXPECTED_LIVE_HEAD},
        "live integration authority drift",
    )

    upstream = contract["audited_upstream"]
    _require(upstream["pull_request"] == 740, "wrong audited PR")
    _require(upstream["head_sha"] == EXPECTED_740_HEAD, "wrong exact #740 head")
    _require(upstream["coordinate_governance_pull_request"] == 738, "missing #738 coordinate governance")
    _require(upstream["reshape_preflight_pull_request"] == 732, "wrong reshape preflight PR")
    _require(upstream["selection_data_governance_pull_request"] == 712, "missing #712 selection governance")
    _require(upstream["fresh_path_scope_governance_pull_request"] == 721, "missing #721 fresh-path governance")
    _require(upstream["source_visual_child_pull_request"] == 652, "wrong source child")
    _require(upstream["old_tip_child_pull_request"] == 701, "wrong old tip child")

    classes = contract["source_classification"]
    _require(set(classes["allowed_classes"]) == ALLOWED_SOURCE_CLASSES, "source-class vocabulary drift")
    used_classes = {item["class"] for item in classes["items"]}
    _require(used_classes == ALLOWED_SOURCE_CLASSES, "four-way source classification incomplete")

    source = replay_source if replay_source is not None else (root / REPLAY_REL).read_text(encoding="utf-8")
    values = _literal_assignments(source)
    _require(values["PREREG_ISSUE"] == 739, "#740 preregistration issue drift")
    _require(values["SOURCE_PARENT_PR"] == 738, "#740 source-parent PR drift")
    _require(values["SOURCE_PARENT_HEAD"] == EXPECTED_738_HEAD, "#740 source-parent exact head drift")
    _require(values["SELECTED_P"] == 9 and values["SELECTED_M"] == 4, "#740 p/m identity drift")
    _require(values["RELATED_FRESH_HOLDOUT_PR"] == 714, "#740 fresh holdout identity drift")
    _require(values["TARGET_PROBE_ABS_Z"] == 1.55, "#740 targeted probe z drift")

    truth = values["TRUTH"]
    _require(truth["canonical_velocity_changed"] is False, "#740 promoted canonical velocity")
    _require(truth["saved_velocity_changed"] is False, "#740 promoted saved velocity")
    _require(truth["production_candidate_selected"] is False, "#740 selected production candidate")
    _require(truth["historical_development_paths_used_for_coefficient"] is True, "#740 alpha provenance drift")
    _require(truth["fresh_714_path_data_used"] is False, "#714 data leaked into #740")
    _require(truth["fresh_714_path_data_used_for_retuning"] is False, "#714 data retuned #740")
    _require(truth["pressure_or_force_changed"] is False, "#740 changed pressure/forcing")
    _require(truth["held_out_pde_residual_evaluated"] is False, "#740 laundered PDE validation")
    for key in (
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth[key] is False, f"forbidden #740 truth promotion: {key}")

    for snippet in (
        "alpha = float(-mean_base / mean_unit)",
        '"historical_development_protocol_used": True',
        '"fresh_714_data_used": False',
        "candidate_fn = lambda p, t: reshaped_velocity(child_fn, p, t, alpha_94)",
        "out.write_text(json.dumps(report",
    ):
        _require(snippet in source, f"#740 procedural-child implementation drift: {snippet}")

    for forbidden in ("save_candidate(", "load_candidate(", "candidate_identity_sha256"):
        _require(forbidden not in source, f"#740 unexpectedly gained standalone identity mechanism: {forbidden}")

    scope = contract["procedural_child_scope"]
    for key in (
        "velocity_formula_evaluable_in_740",
        "new_alpha_rederived_from_historical_development_paths",
        "candidate_lambda_constructed_inside_run",
        "alpha_recorded_in_report_when_run_completes",
    ):
        _require(scope[key] is True, f"expected procedural fact drifted: {key}")
    for key in (
        "fresh_714_data_used_for_alpha",
        "old_701_alpha_reused_as_same_coordinate",
        "basis_dimension_changed",
        "pressure_or_force_rebuilt",
        "held_out_pde_residual_evaluated",
        "standalone_candidate_payload_emitted",
        "standalone_candidate_load_api_provided",
        "standalone_candidate_identity_digest_bound",
        "canonical_or_saved_velocity_promoted",
    ):
        _require(scope[key] is False, f"forbidden procedural/identity promotion: {key}")

    rules = contract["identity_rule"]
    for key in (
        "procedural_reconstruction_is_not_the_same_claim_as_saved_loadable_candidate_identity",
        "a_future_fresh_post_freeze_evaluation_must_not_rederive_or_retune_candidate_parameters_from_the_fresh_data",
        "before_calling_a_future_740_child_check_post_freeze_evidence_freeze_the_child_identity",
        "740_targeted_and_historical_development_evidence_can_select_or_reject_the_child",
        "pde_validation_requires_cr001_independent_acceptance",
    ):
        _require(rules[key] is True, f"identity rule weakened: {key}")
    for key in (
        "740_targeted_and_historical_development_evidence_is_independent_visual_validation",
        "fresh_path_or_render_evidence_may_promote_visual_correspondence_without_public_target_protocol",
    ):
        _require(rules[key] is False, f"evidence scope laundered: {key}")
    minimum = set(rules["minimum_frozen_identity_fields"])
    _require(
        minimum
        == {
            "parent/source candidate identity",
            "p=9,m=4 envelope identity",
            "numeric alpha_94",
            "time activation g(t)=2*(t-.25)",
            "support r<1.6 and 1<|z|<1.8",
            "exact source/runtime or equivalent reproducibility binding",
        },
        "minimum frozen child identity fields drift",
    )
    _require(len(rules["acceptable_freeze_forms"]) == 2, "freeze-form contract drift")

    states = contract["claim_states"]
    _require(states["procedural_velocity_evaluable"] is True, "procedural evaluability was lost")
    for key in (
        "standalone_saved_loadable_child_ready",
        "velocity_export_ready_promoted_by_740",
        "visualization_ready_promoted_by_740",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(states[key] is False, f"forbidden claim-state promotion: {key}")

    _audit_cr001(root, contract["cr001_lock"])
    return {
        "contract_id": contract["contract_id"],
        "audited_pr": upstream["pull_request"],
        "audited_head": upstream["head_sha"],
        "procedural_velocity_evaluable": True,
        "standalone_saved_loadable_child_ready": False,
        "fresh_evidence_requires_frozen_child_identity": True,
        "pde_validated": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
