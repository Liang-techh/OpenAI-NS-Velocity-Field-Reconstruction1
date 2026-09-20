"""Fail-closed CR002 audit for ST052-M reservoir axial-capacity scope.

PR #785 derives exact/numerical geometry for the *correction channel* used by
the outer-reservoir child.  This module prevents that unit-channel geometry
from being promoted to total-child geometry, public OpenAI numerical data,
visual correspondence, or PDE evidence.
"""
from __future__ import annotations

import ast
import copy
import json
import math
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = REPO_ROOT / "configs/st052m_reservoir_axial_capacity_scope_contract.json"
SOURCE_PATH = REPO_ROOT / "experiments/root_st052/agent7_st052m_reservoir_axial_stretch_audit.py"
CONSTRAINTS_PATH = REPO_ROOT / "configs/constraints.json"

EXPECTED_TASK_ID = "CR002-ST052M-RESERVOIR-AXIAL-CAPACITY-SCOPE-112"
EXPECTED_SOURCE_TASK_ID = "CR003-ST052M-RESERVOIR-AXIAL-STRETCH-AUDIT-111"
EXPECTED_SOURCE_HEAD = "55747bb813fedccd3df5635668fce9099a2c13dc"
EXPECTED_PARENT_HEAD = "93887a59729d22113badf2ae4dac2f7d868e2703"
EXPECTED_CLASSES = [
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
]
EXPECTED_Q_ROOT = (13.0 + 5.0 * math.sqrt(17.0)) / 30.0
EXPECTED_R_ROOT = math.sqrt(EXPECTED_Q_ROOT)
EXPECTED_AREA_FRACTION = EXPECTED_Q_ROOT / (64.0 / 25.0)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _literal_assignments(source_text: str) -> dict[str, Any]:
    tree = ast.parse(source_text)
    out: dict[str, Any] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                continue
            name = node.targets[0].id
            value_node = node.value
        elif isinstance(node, ast.AnnAssign):
            if not isinstance(node.target, ast.Name) or node.value is None:
                continue
            name = node.target.id
            value_node = node.value
        else:
            continue
        try:
            out[name] = ast.literal_eval(value_node)
        except (ValueError, TypeError):
            continue
    return out


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def load_default_inputs() -> tuple[dict[str, Any], str, dict[str, Any]]:
    return (
        _load_json(CONTRACT_PATH),
        SOURCE_PATH.read_text(encoding="utf-8"),
        _load_json(CONSTRAINTS_PATH),
    )


def _audit_cr001(contract: dict[str, Any], constraints: dict[str, Any]) -> None:
    snap = contract["cr001_snapshot"]
    _require(constraints["nu"] == snap["nu"] == 0.01, "CR001 nu drift")
    domain = constraints["domain"]
    _require(domain["physical"] == snap["physical_domain"] == "R^3", "CR001 physical domain drift")
    _require(domain["evaluation_box"] == snap["evaluation_box"] == [[-2, 2], [-2, 2], [-2, 2]], "CR001 evaluation box drift")
    _require(domain["support"] == snap["support"] == "r < 2 and abs(z) < 2", "CR001 support drift")
    _require(domain["time_interval"] == snap["time_interval"] == [0.25, 0.75], "CR001 time interval drift")

    force = constraints["forcing"]
    _require(force["mode"] == snap["forcing_mode"] == "restricted_two_parameter_family", "CR001 forcing mode drift")
    _require(force["parameters"] == snap["forcing_parameter_bounds"], "CR001 forcing bounds drift")
    restriction = force["restriction"]
    _require("No residual-dependent basis or pointwise free force" in restriction, "free/residual-defined forcing shortcut detected")

    nontriviality = constraints["nontriviality"]
    _require(nontriviality["reference_energy"] == snap["reference_energy"] == 1.0, "reference energy drift")
    _require(nontriviality["reference_energy_abs_tolerance"] == snap["reference_energy_abs_tolerance"] == 0.001, "energy tolerance drift")
    _require("reject collapsed candidates" in nontriviality["enforcement"], "amplitude-collapse guard drift")

    validation = constraints["validation"]
    _require(validation["seed"] == snap["validation_seed"] == 914027, "validation seed drift")
    _require(validation["held_out_points"] == snap["held_out_points"] == 4096, "held-out point count drift")
    _require(validation["times"] == snap["validation_times"], "validation times drift")
    _require(validation["derivative_steps"] == snap["derivative_steps"] == [0.02, 0.01, 0.005], "derivative ladder drift")
    _require(validation["quadrature_orders_per_axis"] == snap["quadrature_orders_per_axis"] == [24, 48, 96], "quadrature ladder drift")
    thresholds = validation["thresholds"]
    for key, expected in (
        ("divergence_max", 1e-5),
        ("divergence_L2", 1e-5),
        ("pde_residual_max", 1e-3),
        ("pde_residual_L2", 1e-3),
    ):
        _require(thresholds[key] == snap[key] == expected, f"CR001 threshold drift: {key}")
    _require("changing thresholds requires a new experiment version" in validation["failure_policy"], "failure-policy relaxation detected")


def audit(
    contract: dict[str, Any] | None = None,
    source_text: str | None = None,
    constraints: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if contract is None or source_text is None or constraints is None:
        default_contract, default_source, default_constraints = load_default_inputs()
        contract = default_contract if contract is None else contract
        source_text = default_source if source_text is None else source_text
        constraints = default_constraints if constraints is None else constraints

    contract = copy.deepcopy(contract)
    constraints = copy.deepcopy(constraints)
    assignments = _literal_assignments(source_text)

    _require(contract["task_id"] == EXPECTED_TASK_ID, "contract task id drift")
    audited = contract["audited_increment"]
    _require(audited["pr"] == 785, "audited PR drift")
    _require(audited["head"] == EXPECTED_SOURCE_HEAD, "audited exact head drift")
    _require(audited["task_id"] == EXPECTED_SOURCE_TASK_ID, "audited source task drift")
    _require(audited["parent_pr"] == 775 and audited["parent_head"] == EXPECTED_PARENT_HEAD, "audited parent drift")

    _require(contract["classification_vocabulary"] == EXPECTED_CLASSES, "classification vocabulary drift")
    _require(list(contract["classification"].keys()) == EXPECTED_CLASSES, "classification buckets drift")
    for name in EXPECTED_CLASSES:
        _require(bool(contract["classification"][name]), f"empty classification bucket: {name}")

    geometry = contract["frozen_channel_geometry"]
    _require(geometry["quantity_under_audit"] == "outer_reservoir_correction_channel_unit_axial_response", "quantity scope drift")
    _require(geometry["geometry_provenance"] == "autonomous_design", "channel geometry was laundered into source provenance")
    _require(geometry["public_numeric_core_target"] is None, "invented public numerical core target")
    _require(geometry["sign_polynomial_coefficients"] == [-225, 195, 64], "sign polynomial drift")
    _require(abs(geometry["physical_q_root"] - EXPECTED_Q_ROOT) <= 5e-14, "q-root drift")
    _require(abs(geometry["physical_r_root"] - EXPECTED_R_ROOT) <= 5e-14, "r-root drift")
    _require(abs(geometry["core_cross_section_area_fraction"] - EXPECTED_AREA_FRACTION) <= 5e-14, "core-area fraction drift")
    _require(geometry["support_radius"] == 1.6, "correction support radius drift")
    _require(geometry["live_probe_abs_z"] == 1.55, "probe z drift")
    _require(geometry["live_probe_radii"] == [0.6, 0.9, 1.2, 1.5], "probe radii drift")
    _require(geometry["radial_refinement_grids"] == [5001, 20001, 80001], "radial refinement drift")

    scope = contract["channel_vs_total_field"]
    _require(scope["correction_channel_geometry_eligible_if_pr785_passes"] is True, "channel evidence eligibility drift")
    for key in (
        "total_child_axial_velocity_evaluated_by_pr785",
        "parent_velocity_contribution_included_in_axial_shape_factor",
        "frozen_child_coefficient_included_in_axial_shape_factor",
        "channel_core_radius_may_be_called_total_child_core_radius",
        "channel_zero_axial_flux_may_be_called_total_child_zero_axial_flux",
        "channel_core_area_fraction_may_be_called_public_openai_core_fraction",
    ):
        _require(scope[key] is False, f"correction-channel/total-field claim laundering: {key}")

    anti = contract["anti_shortcut"]
    _require(all(value is False for value in anti.values()), "scientific shortcut was enabled")
    truth = contract["truth_boundary"]
    _require(all(value is False for value in truth.values()), "truth-boundary promotion detected")

    _require(assignments.get("TASK_ID") == EXPECTED_SOURCE_TASK_ID, "live #785 task id drift")
    _require(assignments.get("PREREG_ISSUE") == 784, "live #785 preregistration drift")
    _require(assignments.get("SOURCE_PARENT_PR") == 775, "live #785 parent PR drift")
    _require(assignments.get("SOURCE_PARENT_HEAD") == EXPECTED_PARENT_HEAD, "live #785 parent head drift")
    _require(assignments.get("SOURCE_RESERVOIR_PR") == 767, "live #785 reservoir source drift")
    _require(assignments.get("PROBE_ABS_Z") == 1.55, "live #785 probe z drift")
    _require(assignments.get("PROBE_RADII") == (0.6, 0.9, 1.2, 1.5), "live #785 probe radii drift")
    _require(assignments.get("RADIAL_GRIDS") == (5001, 20001, 80001), "live #785 radial grids drift")

    source_truth = assignments.get("TRUTH")
    _require(isinstance(source_truth, dict), "live #785 TRUTH dictionary missing")
    for key in (
        "candidate_velocity_changed",
        "canonical_velocity_changed",
        "saved_velocity_changed",
        "candidate_coefficient_changed",
        "pressure_or_force_changed",
        "held_out_pde_residual_evaluated",
        "public_image_numeric_target_used",
        "core_thickness_source_threshold_defined",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(source_truth.get(key) is False, f"live #785 truth bit promoted: {key}")

    _require("return 2.0 * (radial + q * radial_q)" in source_text, "unit-channel axial response formula drift")
    _require("reservoir.reservoir_correction" in source_text, "live correction-channel probe path missing")
    _require('"reservoir_channel_has_visible_inward_radial_and_core_axial_stretch_capacity": True' in source_text, "channel-capacity conclusion drift")
    _require('"axial_stretch_core_radius_is_fixed_by_current_radial_profile": True' in source_text, "channel core-radius conclusion drift")

    _audit_cr001(contract, constraints)

    return {
        "task_id": EXPECTED_TASK_ID,
        "audited_pr": 785,
        "audited_head": EXPECTED_SOURCE_HEAD,
        "scope": "correction_channel_unit_response_only",
        "physical_r_root": EXPECTED_R_ROOT,
        "core_area_fraction": EXPECTED_AREA_FRACTION,
        "public_numeric_core_target": None,
        "total_child_axial_velocity_evaluated": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "cr001_unchanged": True,
    }


def main() -> int:
    print(json.dumps(audit(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
