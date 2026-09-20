"""Fail-closed CR002 audit for PR #750 axial return-concentration semantics.

PR #750 evaluates only the one-dimensional axial response shape G_p(z) of the
repository-designed compact tip correction channel.  Its unweighted dz metrics
are useful representation diagnostics, but they are not an absolute correction
magnitude, a cylindrical 3-D energy norm, or a total parent-plus-correction
velocity statement.  This audit keeps those claims separate.

No velocity, pressure, forcing, candidate parameter, or scientific threshold is
changed here.
"""
from __future__ import annotations

import ast
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

CONTRACT_REL = Path("configs/st052m_tip_return_concentration_metric_scope.json")
AUDIT_REL = Path("experiments/root_st052/agent7_st052m_tip_return_concentration_audit.py")
CR001_REL = Path("configs/constraints.json")

EXPECTED_LIVE_HEAD = "7b32d9c0fa46fe369d0fc636eb5a7530aadebe04"
EXPECTED_750_HEAD = "376223a3bb2dae0d53756ea48da26937a95b6e58"
EXPECTED_740_HEAD = "21019df839557e64b8d8689f68cdf8fc685a3a36"
ALLOWED_SOURCE_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


class GovernanceError(AssertionError):
    """Raised when the governed metric/truth boundary drifts."""


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


def _function(tree: ast.Module, name: str) -> ast.FunctionDef:
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise GovernanceError(f"missing function {name}")


def _audit_metric_implementation(source: str) -> None:
    tree = ast.parse(source)
    axial = _function(tree, "axial_response")
    metrics = _function(tree, "_metrics")

    _require([arg.arg for arg in axial.args.args] == ["z", "p"], "axial response gained non-axial inputs")
    _require([arg.arg for arg in metrics.args.args] == ["p", "n"], "metric gained physical-field inputs")

    metric_names = {node.id for node in ast.walk(metrics) if isinstance(node, ast.Name)}
    for forbidden in ("alpha", "alpha_94", "velocity", "points", "radius", "radii", "jacobian"):
        _require(forbidden not in metric_names, f"#750 metric unexpectedly uses {forbidden}")

    trapezoids = [
        node
        for node in ast.walk(metrics)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "trapezoid"
    ]
    _require(len(trapezoids) >= 5, "expected frozen dz quadratures missing")
    for call in trapezoids:
        _require(len(call.args) >= 2, "trapezoid call missing coordinate")
        _require(isinstance(call.args[1], ast.Name) and call.args[1].id == "z", "metric is no longer unweighted dz")

    snippets = (
        "g = axial_response(z, p)",
        "energy = float(np.trapezoid(g * g, z))",
        "outward_energy = float(np.trapezoid(np.where(g < 0.0, g * g, 0.0), z))",
        '"normalized_derivative_sharpness"',
        '"normalized_peak_response"',
        '"return_flow_concentrated_by_reshape"',
    )
    for snippet in snippets:
        _require(snippet in source, f"#750 metric implementation drift: {snippet}")


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
    audit_source: str | None = None,
) -> dict[str, Any]:
    root = repository_root() if root is None else Path(root)
    contract = load_contract(root) if contract is None else deepcopy(contract)

    _require(
        contract["contract_id"] == "cr002-st052m-tip-return-concentration-metric-scope-v1",
        "unexpected contract id",
    )
    _require(contract["status"] == "governance_scope_for_exact_750_head", "unexpected contract status")
    _require(
        contract["live_integration_authority"]
        == {"branch": "codex/cr001-constraints", "head_sha": EXPECTED_LIVE_HEAD},
        "live integration authority drift",
    )

    upstream = contract["audited_upstream"]
    _require(upstream["pull_request"] == 750, "wrong audited PR")
    _require(upstream["head_sha"] == EXPECTED_750_HEAD, "wrong exact #750 head")
    _require(upstream["nonlinear_child_pull_request"] == 740, "wrong nonlinear child PR")
    _require(upstream["nonlinear_child_head_sha"] == EXPECTED_740_HEAD, "wrong exact #740 head")
    _require(upstream["reshape_preflight_pull_request"] == 732, "wrong reshape preflight")
    _require(upstream["prior_return_flow_scope_pull_request"] == 730, "missing prior correction-vs-total scope")
    _require(upstream["coefficient_coordinate_scope_pull_request"] == 738, "missing coefficient-coordinate scope")
    _require(upstream["child_identity_scope_pull_request"] == 748, "missing child-identity scope")

    classes = contract["source_classification"]
    _require(set(classes["allowed_classes"]) == ALLOWED_SOURCE_CLASSES, "source-class vocabulary drift")
    _require({item["class"] for item in classes["items"]} == ALLOWED_SOURCE_CLASSES, "four-way source classification incomplete")

    source = audit_source if audit_source is not None else (root / AUDIT_REL).read_text(encoding="utf-8")
    values = _literal_assignments(source)
    _require(values["TASK_ID"] == "CR003-ST052M-TIP-RETURN-CONCENTRATION-AUDIT-104", "#750 task identity drift")
    _require(values["ISSUE"] == 749, "#750 preregistration issue drift")
    _require(values["SOURCE_NONLINEAR_PR"] == 740, "#750 nonlinear source PR drift")
    _require(values["SOURCE_NONLINEAR_HEAD"] == EXPECTED_740_HEAD, "#750 nonlinear source head drift")
    _require(values["BASELINE_P"] == 4 and values["RESHAPED_P"] == 9 and values["M"] == 4, "#750 envelope identity drift")
    _require(tuple(values["GRID_SIZES"]) == (25001, 50001, 100001), "#750 grid ladder drift")

    truth = values["TRUTH"]
    for key in (
        "canonical_velocity_changed",
        "saved_velocity_changed",
        "candidate_velocity_changed",
        "basis_dimension_changed",
        "second_poloidal_basis_added",
        "new_temporal_basis_added",
        "coefficient_fit_performed",
        "trajectory_replay_performed",
        "fresh_714_path_data_used",
        "held_out_pde_residual_evaluated",
        "pressure_or_force_changed",
        "public_image_numeric_target_used",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth[key] is False, f"forbidden #750 truth promotion: {key}")

    _audit_metric_implementation(source)

    semantics = contract["metric_semantics"]
    _require(semantics["diagnostic_object"] == "one_dimensional_axial_response_shape_of_autonomous_tip_correction_channel", "metric object drift")
    _require(semantics["integration_variable"] == "z", "metric integration variable drift")
    _require(semantics["integration_measure"] == "unweighted_dz_on_positive_tip_lobe_[1,1.8]", "metric measure drift")
    _require(semantics["grid_sizes"] == [25001, 50001, 100001], "metric grid contract drift")
    for key in (
        "compares_p4_m4_to_p9_m4",
        "metric_is_shape_level_representation_diagnostic",
        "metric_is_not_a_physical_3d_energy_norm",
    ):
        _require(semantics[key] is True, f"expected metric-scope fact drifted: {key}")
    for key in (
        "coefficient_alpha_used",
        "recalibrated_alpha_94_used",
        "radial_profile_used_in_metric",
        "cylindrical_volume_jacobian_used",
        "cartesian_or_cylindrical_3d_volume_integral_used",
        "nonlinear_child_velocity_evaluated",
        "total_parent_plus_correction_velocity_evaluated",
        "absolute_correction_magnitude_compared",
        "three_dimensional_correction_energy_compared",
    ):
        _require(semantics[key] is False, f"metric scope laundered: {key}")

    rules = contract["interpretation_rule"]
    for key in (
        "a_successful_750_receipt_may_support_axial_shape_return_concentration",
        "a_successful_750_receipt_may_support_thinner_geometric_return_collar",
        "a_successful_750_receipt_may_support_p9_has_higher_shape_normalized_outward_fraction_or_sharpness",
        "later_absolute_or_3d_claim_requires_alpha_bound_and_physical_measure",
        "later_total_field_claim_requires_parent_plus_correction_evaluation",
        "queued_740_ci_must_not_be_laundered_into_pass",
    ):
        _require(rules[key] is True, f"interpretation guard weakened: {key}")
    for key in (
        "axial_shape_concentration_implies_larger_absolute_physical_return_velocity",
        "axial_shape_concentration_implies_three_dimensional_correction_energy_concentration",
        "axial_shape_concentration_implies_total_child_return_flow_concentration",
        "axial_shape_concentration_implies_nonlinear_child_passes_740",
        "axial_shape_concentration_implies_visual_correspondence",
        "axial_shape_concentration_implies_pde_validation",
    ):
        _require(rules[key] is False, f"forbidden inference promoted: {key}")

    states = contract["claim_states"]
    _require(states["axial_response_shape_scope_defined"] is True, "axial metric scope lost")
    for key in (
        "absolute_correction_magnitude_verified_by_750",
        "three_dimensional_correction_energy_concentration_verified",
        "cylindrical_volume_weighted_concentration_verified",
        "nonlinear_child_total_field_return_flow_concentration_verified",
        "nonlinear_child_740_passed_at_contract_freeze",
        "production_candidate_selected",
        "velocity_export_ready_promoted_by_750",
        "visualization_ready_promoted_by_750",
        "visual_correspondence_verified",
        "pde_validated",
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
        "diagnostic_scope": "one_dimensional_axial_response_shape",
        "absolute_or_3d_concentration_verified": False,
        "total_child_return_flow_concentration_verified": False,
        "pde_validated": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
