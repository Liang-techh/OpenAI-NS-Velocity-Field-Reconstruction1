"""Fail-closed CR002 audit for the #723 ST052-M tip return-flow obstruction.

The upstream analytic observation is deliberately narrow: the autonomous #701
*correction channel* has a compact axial potential, so its radial correction
has zero integral across each connected axial lobe and changes sign.  This does
not, by itself, prove that the total #701 child has return flow everywhere the
correction does, nor does it prove an obstruction for arbitrary divergence-free
fields.  A fresh #714 path PASS at |z|=1.25 is local directional evidence only.

This module changes no candidate, coefficient, force, threshold, or delivery API.
"""
from __future__ import annotations

import ast
import json
import math
from copy import deepcopy
from fractions import Fraction
from pathlib import Path
from typing import Any


CONTRACT_REL = Path("configs/st052m_tip_return_flow_scope.json")
UPSTREAM_REL = Path("experiments/root_st052/agent7_st052m_tip_return_flow_obstruction.py")
CR001_REL = Path("configs/constraints.json")

EXPECTED_723_HEAD = "1b16147b93e7b5279c833e8629e10aff1ef4f5a6"
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


def _independent_sign_data() -> dict[str, float]:
    """Re-derive the #701 axial sign change without calling the #723 helper."""
    q_a = Fraction(1, 1)
    q_b = Fraction(81, 25)
    c2 = Fraction(-17, 1)
    c1 = 9 * (q_a + q_b)
    c0 = -(q_a * q_b)
    _require(c1 == Fraction(954, 25), "independent sign-polynomial linear coefficient drift")
    _require(c0 == Fraction(-81, 25), "independent sign-polynomial constant drift")

    a, b, c = float(c2), float(c1), float(c0)
    disc = b * b - 4.0 * a * c
    _require(disc > 0.0, "sign polynomial lost real roots")
    roots = sorted(((-b - math.sqrt(disc)) / (2.0 * a), (-b + math.sqrt(disc)) / (2.0 * a)))
    inside = [q for q in roots if float(q_a) < q < float(q_b)]
    _require(len(inside) == 1, "expected exactly one physical sign root")
    q_star = float(inside[0])
    z_star = math.sqrt(q_star)
    return {
        "q_small": float(roots[0]),
        "q_star": q_star,
        "z_star": z_star,
        "n_at_1p25": a * (1.25**4) + b * (1.25**2) + c,
        "n_at_1p60": a * (1.60**4) + b * (1.60**2) + c,
    }


def audit_contract(contract: dict[str, Any]) -> None:
    upstream = contract["audited_upstream"]
    _require(upstream["pull_request"] == 723, "wrong audited PR")
    _require(upstream["head_sha"] == EXPECTED_723_HEAD, "#723 exact-head identity drift")
    _require(upstream["source_candidate_pull_request"] == 701, "wrong source candidate PR")
    _require(upstream["source_candidate_head_sha"] == EXPECTED_701_HEAD, "#701 source head drift")
    _require(upstream["fresh_path_pull_request"] == 714, "fresh-path linkage drift")
    _require(upstream["fresh_path_scope_governance_pull_request"] == 721, "fresh-path governance linkage drift")
    _require(upstream["upstream_exact_head_ci_passed_at_contract_freeze"] is False, "queued upstream CI laundered into PASS")

    source = contract["source_classification"]
    _require(set(source["allowed_classes"]) == ALLOWED_SOURCE_CLASSES, "source vocabulary drift")
    for item in source["items"]:
        _require(item["class"] in ALLOWED_SOURCE_CLASSES, f"invalid source class: {item}")
    by_item = {item["item"]: item["class"] for item in source["items"]}
    _require(
        by_item["openai_public_qualitative_inward_spiraling_and_axial_elongation"] == "public_source_fact",
        "public qualitative source classification drift",
    )
    _require(
        by_item["701_tip_local_compact_odd_poloidal_channel_formula_support_orientation_and_alpha_rule"] == "autonomous_design",
        "#701 autonomous basis laundered into public-source fact",
    )
    _require(
        by_item["claim_that_every_point_of_the_compact_tip_lobe_has_inward_radial_velocity"] == "pending_unknown",
        "whole-lobe inwardness prematurely promoted",
    )
    _require(
        by_item["hidden_openai_velocity_coefficients_camera_frame_time_mapping_and_full_numeric_field"] == "pending_unknown",
        "hidden OpenAI numerical information promoted",
    )

    channel = contract["single_channel_scope"]
    _require(channel["vector_potential"] == "A=(-y*f,x*f,0)", "vector-potential identity drift")
    _require(channel["scalar"] == "f=R(r^2)*z*B(z^2)", "scalar channel identity drift")
    _require(channel["radial_component"] == "C_r=-r*f_z", "radial component identity drift")
    _require(channel["axial_q_support"] == [1.0, 3.24], "q support drift")
    _require(channel["axial_abs_z_support"] == [1.0, 1.8], "z support drift")
    _require(channel["fresh_714_tip_abs_z"] == 1.25, "fresh #714 local probe drift")
    _require(channel["outer_return_flow_spot_abs_z"] == 1.6, "outer return-flow probe drift")
    _require(channel["correction_channel_has_inward_and_outward_radial_subbands"] is True, "single-channel sign-change fact erased")
    _require(channel["nontrivial_single_compact_correction_channel_can_be_strictly_inward_on_entire_connected_lobe"] is False, "single-channel obstruction erased")
    _require(channel["fresh_714_tip_seed_lies_in_analytic_inward_subband"] is True, "fresh local probe scope drift")

    sign = _independent_sign_data()
    _require(abs(sign["z_star"] - channel["physical_abs_z_sign_change_expected"]) <= 1e-12, "physical sign-change value drift")
    _require(1.0 < 1.25 < sign["z_star"] < 1.60 < 1.8, "local/return-flow probe ordering drift")
    _require(sign["n_at_1p25"] > 0.0 and sign["n_at_1p60"] < 0.0, "independent sign classification drift")

    scope = contract["scope_limits"]
    _require(scope["obstruction_applies_to_701_correction_channel"] is True, "correction-channel scope lost")
    _require(scope["fresh_714_pass_may_support_local_post_freeze_directional_generalization"] is True, "legitimate local fresh evidence erased")
    for key in (
        "obstruction_proves_total_701_child_radial_velocity_has_return_flow",
        "obstruction_proves_any_divergence_free_compact_field_requires_same_return_flow_pattern",
        "fresh_714_pass_may_support_whole_tip_support_inwardness",
        "fresh_714_pass_may_support_global_visual_correspondence",
        "return_flow_obstruction_is_a_pde_acceptance_result",
        "return_flow_obstruction_is_a_public_openai_field_identity_result",
    ):
        _require(scope[key] is False, f"return-flow scope over-promoted: {key}")

    promotion = contract["promotion_rules"]
    for key, value in promotion.items():
        _require(value is True, f"promotion guard relaxed: {key}")

    claims = contract["claim_states"]
    for key, value in claims.items():
        _require(value is False, f"scientific/delivery claim promoted by representation audit: {key}")

    lock = contract["cr001_lock"]
    _require(lock["nu"] == 0.01, "nu drift")
    _require(lock["physical_domain"] == "R^3", "physical domain drift")
    _require(lock["evaluation_box"] == [[-2, 2], [-2, 2], [-2, 2]], "evaluation box drift")
    _require(lock["support"] == "r < 2 and abs(z) < 2", "support drift")
    _require(lock["time_interval"] == [0.25, 0.75], "time interval drift")
    _require(lock["forcing_mode"] == "restricted_two_parameter_family", "forcing mode drift")
    _require(lock["forcing_bounds"] == {"a": [0.0, 10.0], "c": [0.0, 10.0]}, "forcing bounds drift")
    _require(lock["reference_energy"] == 1.0 and lock["reference_energy_abs_tolerance"] == 0.001, "energy normalization drift")
    _require(lock["validation_seed"] == 914027 and lock["held_out_points"] == 4096, "formal held-out sample drift")
    _require(lock["validation_times"] == [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75], "validation times drift")
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

    upstream_source = (root / UPSTREAM_REL).read_text(encoding="utf-8")
    values = _literal_assignments(upstream_source)
    cr001 = json.loads((root / CR001_REL).read_text(encoding="utf-8"))

    _require(values["PREREG_ISSUE"] == 722, "#723 preregistration identity drift")
    _require(values["SOURCE_PARENT_PR"] == 701, "#723 source parent PR drift")
    _require(values["SOURCE_PARENT_HEAD"] == EXPECTED_701_HEAD, "#723 source parent head drift")
    _require(values["RELATED_FRESH_HOLDOUT_PR"] == 714, "#723 fresh-holdout linkage drift")
    _require(values["RELATED_HOLDOUT_SCOPE_PR"] == 721, "#723 fresh-holdout governance linkage drift")
    _require(values["AXIAL_Z_INNER"] == 1.0 and values["AXIAL_Z_OUTER"] == 1.8, "#723 axial support drift")
    _require(values["FRESH_HOLDOUT_TIP_Z_ABS"] == 1.25, "#723 local holdout z drift")
    _require(values["SPOT_OUTWARD_Z_ABS"] == 1.60, "#723 return-flow spot drift")

    truth = values["TRUTH"]
    for key in (
        "canonical_velocity_changed",
        "saved_velocity_changed",
        "production_candidate_selected",
        "new_spatial_basis_added",
        "new_temporal_basis_added",
        "source_701_candidate_retuned",
        "parameter_grid_scan_performed",
        "optimization_performed",
        "pressure_or_force_changed",
        "trajectory_result_used",
        "held_out_pde_residual_evaluated",
        "parent_pde_receipt_transferred",
        "public_image_numeric_target_used",
        "pixel_similarity_objective_used",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth[key] is False, f"#723 truth boundary promoted: {key}")

    # Bind the governance scope to the executable analytic structure, while
    # independently re-deriving the root rather than trusting #723 output.
    for token in (
        "AXIAL_Q_A = Fraction(1, 1)",
        "AXIAL_Q_B = Fraction(81, 25)",
        "return Fraction(-17, 1), 9 * (a + b), -(a * b)",
        "integral_identity_rhs = -SPOT_RADIUS * (endpoint_f_outer - endpoint_f_inner)",
        '"single_compact_tip_channel_strictly_inward_everywhere_possible": False',
        '"return_flow_required_by_compact_lobe_identity": True',
        '"full_tip_support_inwardness_inferred_from_holdout": False',
    ):
        _require(token in upstream_source, f"#723 analytic/scope binding missing: {token}")

    sign = _independent_sign_data()
    _require(1.0 < sign["z_star"] < 1.8, "independent root outside physical tip support")
    _require(1.25 < sign["z_star"] < 1.60, "fresh/local vs return-flow subband ordering changed")

    # Bind back to canonical CR001; this audit cannot create a new scientific contract.
    lock = contract["cr001_lock"]
    _require(cr001["nu"] == lock["nu"], "contract/live nu mismatch")
    _require(cr001["domain"]["physical"] == lock["physical_domain"], "contract/live domain mismatch")
    _require(cr001["domain"]["evaluation_box"] == lock["evaluation_box"], "contract/live box mismatch")
    _require(cr001["domain"]["support"] == lock["support"], "contract/live support mismatch")
    _require(cr001["domain"]["time_interval"] == lock["time_interval"], "contract/live time mismatch")
    _require(cr001["forcing"]["mode"] == lock["forcing_mode"], "contract/live forcing mode mismatch")
    _require(cr001["forcing"]["parameters"] == lock["forcing_bounds"], "contract/live forcing bounds mismatch")
    _require(cr001["nontriviality"]["reference_energy"] == lock["reference_energy"], "contract/live energy target mismatch")
    _require(cr001["nontriviality"]["reference_energy_abs_tolerance"] == lock["reference_energy_abs_tolerance"], "contract/live energy tolerance mismatch")
    _require(cr001["validation"]["seed"] == lock["validation_seed"], "contract/live validation seed mismatch")
    _require(cr001["validation"]["held_out_points"] == lock["held_out_points"], "contract/live held-out count mismatch")
    _require(cr001["validation"]["times"] == lock["validation_times"], "contract/live validation times mismatch")
    _require(cr001["validation"]["derivative_steps"] == lock["derivative_steps"], "contract/live derivative ladder mismatch")
    _require(cr001["validation"]["quadrature_orders_per_axis"] == lock["quadrature_orders_per_axis"], "contract/live quadrature ladder mismatch")
    thresholds = cr001["validation"]["thresholds"]
    _require(thresholds["divergence_max"] == lock["divergence_max"] and thresholds["divergence_L2"] == lock["divergence_L2"], "contract/live divergence gate mismatch")
    _require(thresholds["pde_residual_max"] == lock["pde_residual_max"] and thresholds["pde_residual_L2"] == lock["pde_residual_L2"], "contract/live momentum gate mismatch")
    _require("No residual-dependent basis or pointwise free force" in cr001["forcing"]["restriction"], "canonical no-free-force rule disappeared")
    _require("reject collapsed candidates" in cr001["nontriviality"]["enforcement"], "canonical anti-collapse rule disappeared")
    _require("changing thresholds requires a new experiment version" in cr001["validation"]["failure_policy"], "canonical threshold-change policy disappeared")

    return {
        "status": "PASS",
        "audited_pr": 723,
        "independent_physical_abs_z_sign_change": sign["z_star"],
        "fresh_714_tip_seed_in_inward_subband": True,
        "single_channel_return_flow_obstruction": True,
        "obstruction_applies_to_total_701_child": False,
        "obstruction_applies_to_arbitrary_divergence_free_fields": False,
        "whole_tip_support_inwardness_verified": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
    }


def mutated_contract(contract: dict[str, Any], path: tuple[str, ...], value: Any) -> dict[str, Any]:
    """Small helper used by mutation regressions without mutating the fixture."""
    copy = deepcopy(contract)
    cursor: Any = copy
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    return copy


def main() -> None:
    print(json.dumps(audit_live(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
