"""CR002 fail-closed audit for ST052-M tip-basis calibration reuse.

The upstream CR003 screen is allowed to use a deterministic coefficient solved
from its frozen development trajectory protocol.  That makes the screen useful
expression-capacity/routing evidence, but it does not make the same path replay
an independent held-out validation receipt.  This module keeps that distinction
machine-readable without changing the screened velocity field.
"""
from __future__ import annotations

import ast
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = ROOT / "configs" / "st052m_tip_basis_selection_data_contract.json"
SCREEN_PATH = ROOT / "experiments" / "root_st052" / "agent7_st052m_tip_odd_poloidal_screen.py"
CONSTRAINTS_PATH = ROOT / "configs" / "constraints.json"

EXPECTED_UPSTREAM_HEAD = "e3ecc8654f2c9a2324265bbd01df596b4e503404"
EXPECTED_SCREEN_TASK = "CR003-ST052M-TIP-ODD-POLOIDAL-SCREEN-098"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    _require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def _literal_assignment(tree: ast.Module, name: str) -> Any:
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    return ast.literal_eval(node.value)
        if isinstance(node, ast.AnnAssign):
            if isinstance(node.target, ast.Name) and node.target.id == name:
                return ast.literal_eval(node.value)
    raise ValueError(f"missing literal assignment {name}")


def _function(tree: ast.Module, name: str) -> ast.FunctionDef:
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise ValueError(f"missing function {name}")


def _call_name(call: ast.Call) -> str | None:
    if isinstance(call.func, ast.Name):
        return call.func.id
    if isinstance(call.func, ast.Attribute):
        return call.func.attr
    return None


def _calls(function: ast.FunctionDef, name: str) -> list[ast.Call]:
    return [
        node
        for node in ast.walk(function)
        if isinstance(node, ast.Call) and _call_name(node) == name
    ]


def _arg_names(call: ast.Call) -> list[str]:
    return [ast.unparse(arg) for arg in call.args]


def _audit_upstream_source(screen_path: Path) -> dict[str, Any]:
    text = screen_path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(screen_path))

    _require(_literal_assignment(tree, "TASK_ID") == EXPECTED_SCREEN_TASK, "upstream task drifted")
    _require(_literal_assignment(tree, "PREREG_ISSUE") == 700, "upstream preregistration drifted")
    _require(_literal_assignment(tree, "MID_TIME") == 0.50, "midpoint calibration time drifted")
    _require(_literal_assignment(tree, "MID_INDEX") == 16, "midpoint calibration index drifted")

    truth = _literal_assignment(tree, "TRUTH")
    _require(isinstance(truth, dict), "upstream TRUTH must stay a literal mapping")
    for key in (
        "new_spatial_basis_promoted",
        "production_candidate_selected",
        "held_out_pde_residual_evaluated",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth.get(key) is False, f"upstream truth promotion requires CR002 re-audit: {key}")
    _require(truth.get("parameter_grid_scan_performed") is False, "unexpected alpha grid scan")
    _require(truth.get("optimization_performed") is False, "unexpected optimizer use")
    _require(truth.get("post_result_damping_or_retuning_performed") is False, "unexpected retuning")

    calibration = _function(tree, "derive_alpha_from_frozen_midpoint")
    calibration_source = ast.unparse(calibration)
    calibration_args = [arg.arg for arg in calibration.args.args]
    _require("child_positions" in calibration_args, "calibration no longer consumes parent path positions")
    _require("metadata" in calibration_args, "calibration no longer consumes frozen path metadata")
    _require('str(meta["band"]) == "tip"' in calibration_source, "tip-band calibration filter drifted")
    _require("np.mean(base_ur)" in calibration_source, "parent radial-mean calibration drifted")
    _require("np.mean(unit_ur)" in calibration_source, "unit radial-mean calibration drifted")
    _require("alpha = float(-mean_base / mean_unit)" in calibration_source, "frozen alpha rule drifted")
    _require('"tip_path_count": int(len(tip_idx))' in calibration_source, "tip path-count receipt drifted")

    run = _function(tree, "run")
    run_source = ast.unparse(run)
    _require(
        "child_times, child_positions, metadata = _integrate(child_fn)" in run_source,
        "source child path integration drifted",
    )

    derive_calls = _calls(run, "derive_alpha_from_frozen_midpoint")
    _require(len(derive_calls) == 1, "expected exactly one frozen alpha calibration call")
    _require(
        _arg_names(derive_calls[0]) == ["child_fn", "child_times", "child_positions", "metadata"],
        "alpha calibration inputs drifted",
    )

    grouped_calls = _calls(run, "_grouped_records_from_positions")
    _require(
        any(
            _arg_names(call)[:4] == ["child_fn", "child_times", "child_positions", "metadata"]
            for call in grouped_calls
        ),
        "calibration parent paths are no longer reused in source-side path comparison",
    )
    path_calls = _calls(run, "_path_records_from_positions")
    _require(
        any(_arg_names(call)[:2] == ["child_times", "child_positions"] for call in path_calls),
        "calibration parent paths are no longer reused in source-side turns comparison",
    )
    _require(
        "candidate_times, candidate_positions, candidate_metadata = _integrate(candidate_fn)" in run_source,
        "candidate nonlinear replay drifted",
    )
    _require(
        "if not np.array_equal(times, child_times) or meta != metadata" in run_source,
        "same frozen trajectory protocol guard drifted",
    )
    _require(
        '"tip_odd_poloidal_basis_growth_justified": justified' in run_source,
        "screen routing gate drifted",
    )

    return {
        "task_id": EXPECTED_SCREEN_TASK,
        "calibration_time": 0.50,
        "parent_positions_consumed": True,
        "same_protocol_replayed": True,
        "upstream_truth_boundary_fail_closed": True,
    }


def _audit_cr001(constraints_path: Path, contract: dict[str, Any]) -> dict[str, Any]:
    constraints = _load_json(constraints_path)
    freeze = contract.get("cr001_freeze")
    _require(isinstance(freeze, dict), "missing CR001 freeze receipt")

    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]

    expected = {
        "nu": constraints["nu"],
        "physical_domain": domain["physical"],
        "evaluation_box": domain["evaluation_box"],
        "support": domain["support"],
        "time_interval": domain["time_interval"],
        "forcing_mode": forcing["mode"],
        "forcing_parameter_bounds": forcing["parameters"],
        "reference_energy": nontriviality["reference_energy"],
        "reference_energy_abs_tolerance": nontriviality["reference_energy_abs_tolerance"],
        "validation_seed": validation["seed"],
        "held_out_points": validation["held_out_points"],
        "validation_times": validation["times"],
        "derivative_steps": validation["derivative_steps"],
        "quadrature_orders_per_axis": validation["quadrature_orders_per_axis"],
        "divergence_max": thresholds["divergence_max"],
        "divergence_L2": thresholds["divergence_L2"],
        "momentum_residual_max": thresholds["pde_residual_max"],
        "momentum_residual_L2": thresholds["pde_residual_L2"],
    }
    for key, value in expected.items():
        _require(freeze.get(key) == value, f"CR001 drift in {key}")

    _require(
        "No residual-dependent basis or pointwise free force" in forcing["restriction"],
        "CR001 free-force prohibition drifted",
    )
    _require(
        freeze.get("residual_defined_or_pointwise_free_force_forbidden") is True,
        "contract must keep free-force route forbidden",
    )
    _require(
        "reject outside [1e-4,100]" in constraints["optimization"]["candidate_parameter_bounds"]["amplitude"],
        "CR001 amplitude-collapse guard drifted",
    )
    _require(
        freeze.get("amplitude_collapse_success_route_forbidden") is True,
        "contract must keep amplitude-collapse success forbidden",
    )
    return expected


def audit(
    contract_path: str | Path = CONTRACT_PATH,
    screen_path: str | Path = SCREEN_PATH,
    constraints_path: str | Path = CONSTRAINTS_PATH,
) -> dict[str, Any]:
    """Audit selection-data accounting and the independent truth boundaries."""
    contract_path = Path(contract_path)
    screen_path = Path(screen_path)
    constraints_path = Path(constraints_path)
    contract = _load_json(contract_path)

    _require(contract.get("schema") == "st052m-tip-basis-selection-data-governance/v1", "schema drift")
    _require(contract.get("task_id") == "CR002", "task ownership drift")
    upstream = contract.get("audited_upstream")
    _require(isinstance(upstream, dict), "missing audited upstream receipt")
    _require(upstream.get("pull_request") == 701, "audited PR drifted")
    _require(upstream.get("head_sha") == EXPECTED_UPSTREAM_HEAD, "audited head drifted")
    _require(upstream.get("task_id") == EXPECTED_SCREEN_TASK, "audited screen task drifted")
    _require(upstream.get("preregistration_issue") == 700, "preregistration issue drifted")

    classification = contract.get("source_classification")
    _require(isinstance(classification, dict), "missing four-way source classification")
    _require(
        set(classification) == {"user_requirement", "public_source_fact", "autonomous_design", "pending_or_unknown"},
        "source classification must remain four-way and explicit",
    )
    for key, values in classification.items():
        _require(isinstance(values, list) and values, f"empty source-classification bucket: {key}")

    selection = contract.get("selection_data_accounting")
    _require(isinstance(selection, dict), "missing selection-data accounting")
    required_selection = {
        "calibration_time": 0.5,
        "calibration_band": "tip",
        "calibration_parent_path_count": 24,
        "calibration_consumes_frozen_parent_tip_path_positions": True,
        "coefficient_rule": "alpha=-mean(parent_tip_u_r)/mean(unit_tip_correction_u_r)",
        "parameter_grid_scan_performed": False,
        "optimizer_used": False,
        "post_result_damping_or_retuning_allowed": False,
        "post_calibration_path_screen_reuses_same_frozen_seed_protocol": True,
        "parent_path_positions_used_for_calibration_are_also_reused_in_source_side_path_comparisons": True,
        "path_screen_is_independent_held_out_validation": False,
        "morphology_grid_used_to_solve_alpha": False,
        "screen_is_development_expression_capacity_evidence": True,
        "screen_can_support_routing_without_production_promotion": True,
    }
    for key, value in required_selection.items():
        _require(selection.get(key) == value, f"selection-data truth boundary drifted: {key}")

    promotion = contract.get("promotion_boundary")
    _require(isinstance(promotion, dict), "missing promotion boundary")
    for key in (
        "new_spatial_basis_promoted",
        "production_candidate_selected",
        "velocity_export_ready_promoted_by_this_screen",
        "independent_held_out_visual_validation_performed",
        "visual_correspondence_verified",
        "held_out_temporal_child_pde_residual_evaluated",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(promotion.get(key) is False, f"unsupported promotion: {key}")
    _require(
        promotion.get("fresh_post_freeze_evidence_required_before_production_or_visual_correspondence_promotion") is True,
        "fresh-evidence requirement was removed",
    )

    fresh = contract.get("fresh_evidence_requirement")
    _require(isinstance(fresh, dict), "missing fresh-evidence contract")
    for key in (
        "coefficient_must_be_frozen_before_new_evaluation",
        "fresh_path_seeds_or_other_explicitly_held_out_path_protocol_required_for_independent_path_claim",
        "fresh_visual_protocol_must_not_be_used_to_refit_alpha_if_later_called_independent",
        "pde_validation_remains_separate_and_must_use_the_registered_cr001_held_out_protocol",
    ):
        _require(fresh.get(key) is True, f"fresh-evidence guard drifted: {key}")

    upstream_result = _audit_upstream_source(screen_path)
    cr001_result = _audit_cr001(constraints_path, contract)
    return {
        "contract_schema": contract["schema"],
        "audited_upstream_head": upstream["head_sha"],
        "upstream": upstream_result,
        "cr001": cr001_result,
        "selection_data_scope": "development/model-selection evidence; not independent held-out validation",
        "passed": True,
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
