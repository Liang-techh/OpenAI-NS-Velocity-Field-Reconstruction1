"""Fail-closed CR002 audit for #732 ST052-M tip-envelope coefficient semantics.

Agent-7 #732 keeps one poloidal channel and one time law, but changes the axial
basis function from the symmetric p=4,m=4 bump to a peak-normalized p=9,m=4
bump. Keeping the *dimension* of the span at one does not make the old and new
basis vectors identical, and therefore does not make #701's scalar ``alpha`` an
invariant physical coordinate.

This governance module preserves #732's valid sign-topology preflight while
preventing a future nonlinear child from silently treating "same dimension" as
"same unit response". It changes no velocity, coefficient, force, threshold,
or delivery API.
"""
from __future__ import annotations

import ast
import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Any

CONTRACT_REL = Path("configs/st052m_tip_envelope_coordinate_scope.json")
RELOCATION_REL = Path("experiments/root_st052/agent7_st052m_tip_return_flow_relocation.py")
SOURCE_BASIS_REL = Path("experiments/root_st052/agent7_st052m_tip_odd_poloidal_screen.py")
CR001_REL = Path("configs/constraints.json")

EXPECTED_732_HEAD = "1f0391c6e1bad213c2af4b545746d47652ef3aed"
EXPECTED_LIVE_HEAD = "7b32d9c0fa46fe369d0fc636eb5a7530aadebe04"
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


def _peak_normalized_bump(q: float, p: int, m: int, a: float = 1.0, b: float = 3.24) -> tuple[float, float]:
    """Return B and dB/dq for the exact peak-normalized polynomial family."""
    q = float(q)
    p = int(p)
    m = int(m)
    if not (a < q < b):
        return 0.0, 0.0
    peak = (p * b + m * a) / (p + m)
    peak_raw = (peak - a) ** p * (b - peak) ** m
    _require(math.isfinite(peak_raw) and peak_raw > 0.0, "invalid independent bump normalization")
    left = q - a
    right = b - q
    raw = left**p * right**m
    raw_q = p * left ** (p - 1) * right**m - m * left**p * right ** (m - 1)
    return raw / peak_raw, raw_q / peak_raw


def _axial_z_factor(abs_z: float, p: int, m: int = 4) -> float:
    """Factor multiplying the common -r*R(r^2) in the radial correction."""
    q = float(abs_z) ** 2
    bump, bump_q = _peak_normalized_bump(q, p, m)
    return float(bump + 2.0 * q * bump_q)


def _physical_sign_change(p: int, m: int = 4) -> float:
    """Independently solve the #732 sign polynomial and return physical |z|."""
    a = 1.0
    b = 81.0 / 25.0
    c2 = -(1.0 + 2.0 * (p + m))
    c1 = (a + b) + 2.0 * (p * b + m * a)
    c0 = -(a * b)
    disc = c1 * c1 - 4.0 * c2 * c0
    _require(disc > 0.0, "sign polynomial lost real roots")
    roots = sorted(
        (
            (-c1 - math.sqrt(disc)) / (2.0 * c2),
            (-c1 + math.sqrt(disc)) / (2.0 * c2),
        )
    )
    inside = [q for q in roots if a < q < b]
    _require(len(inside) == 1, f"expected one physical sign root for p={p}, m={m}")
    return math.sqrt(float(inside[0]))


def _select_minimal_p() -> tuple[int, float]:
    for p in range(4, 65):
        z_star = _physical_sign_change(p, 4)
        if z_star >= 1.60:
            if p > 4:
                _require(_physical_sign_change(p - 1, 4) < 1.60, "selected p is not minimal")
            return p, z_star
    raise GovernanceError("no discrete p reaches the frozen #732 sign target")


def independent_unit_response_witness() -> dict[str, float]:
    """Show that peak-normalized old/new basis vectors have different unit response."""
    abs_z = 1.25
    baseline = _axial_z_factor(abs_z, 4, 4)
    reshaped = _axial_z_factor(abs_z, 9, 4)
    _require(baseline > 0.0 and reshaped > 0.0, "unit-response witness lost inward orientation")
    ratio = abs(reshaped / baseline)
    _require(abs(reshaped - baseline) > 1.0e-3, "old/new unit responses unexpectedly identical")
    _require(0.0 < ratio < 0.1, "unit-response witness no longer shows material coordinate change")
    return {
        "abs_z": abs_z,
        "q": abs_z**2,
        "baseline_axial_z_factor": baseline,
        "reshaped_axial_z_factor": reshaped,
        "reshaped_to_baseline_magnitude_ratio": ratio,
    }


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
    _require("reject collapsed candidates" in nontrivial["enforcement"], "CR001 amplitude-collapse guard drift")
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
    _require(
        lock["residual_defined_or_pointwise_free_force_allowed"] is False,
        "governance contract allows free forcing",
    )
    _require(lock["amplitude_collapse_success_allowed"] is False, "governance contract allows amplitude collapse")
    _require(lock["threshold_relaxation_allowed"] is False, "governance contract allows threshold relaxation")


def audit(root: Path | None = None, contract: dict[str, Any] | None = None) -> dict[str, Any]:
    root = repository_root() if root is None else Path(root)
    contract = load_contract(root) if contract is None else deepcopy(contract)

    _require(
        contract["contract_id"] == "cr002-st052m-tip-envelope-coordinate-scope-v1",
        "unexpected contract id",
    )
    live = contract["live_integration_authority"]
    _require(
        live == {"branch": "codex/cr001-constraints", "head_sha": EXPECTED_LIVE_HEAD},
        "live integration authority drift",
    )

    upstream = contract["audited_upstream"]
    _require(upstream["pull_request"] == 732, "wrong audited PR")
    _require(upstream["head_sha"] == EXPECTED_732_HEAD, "wrong #732 exact head")
    _require(upstream["parent_obstruction_pull_request"] == 723, "wrong obstruction parent")
    _require(upstream["source_basis_pull_request"] == 701, "wrong source basis PR")
    _require(upstream["selection_data_governance_pull_request"] == 712, "missing #712 calibration governance")
    _require(upstream["fresh_path_pull_request"] == 714, "wrong fresh-path PR")
    _require(upstream["fresh_path_scope_governance_pull_request"] == 721, "missing #721 fresh-path governance")

    classes = contract["source_classification"]
    _require(set(classes["allowed_classes"]) == ALLOWED_SOURCE_CLASSES, "source-class vocabulary drift")
    used_classes = {item["class"] for item in classes["items"]}
    _require(used_classes <= ALLOWED_SOURCE_CLASSES, "unclassified source item")
    _require(ALLOWED_SOURCE_CLASSES == used_classes, "four-way source classification incomplete")

    relocation_source = (root / RELOCATION_REL).read_text(encoding="utf-8")
    source_basis = (root / SOURCE_BASIS_REL).read_text(encoding="utf-8")
    relocation_values = _literal_assignments(relocation_source)
    source_values = _literal_assignments(source_basis)

    _require(
        relocation_values["SOURCE_PARENT_HEAD"] == "1b16147b93e7b5279c833e8629e10aff1ef4f5a6",
        "#732 parent head drift",
    )
    _require(relocation_values["SOURCE_BASIS_PR"] == 701, "#732 source basis drift")
    _require(relocation_values["OUTER_EXPONENT"] == 4, "#732 outer exponent drift")
    _require(relocation_values["MIN_INNER_EXPONENT"] == 4, "#732 minimum p drift")
    _require(relocation_values["TARGET_SIGN_CHANGE_ABS_Z"] == 1.60, "#732 sign target drift")
    _require(
        source_values["Z2_BUMP_A"] == 1.0 and source_values["Z2_BUMP_B"] == 3.24,
        "#701 axial support drift",
    )
    _require(
        source_values["Z_SUPPORT_INNER"] == 1.0 and source_values["Z_SUPPORT_OUTER"] == 1.8,
        "#701 |z| support drift",
    )

    for snippet in (
        "peak = (float(p) * b + float(m) * a) / float(p + m)",
        "out[inside] = raw / peak_raw",
        "deriv[inside] = raw_q / peak_raw",
    ):
        _require(snippet in relocation_source, f"#732 peak-normalized implementation drift: {snippet}")
    _require('"""Normalized ((q-a)(b-q))^4 bump' in source_basis, "#701 baseline normalization doc drift")
    _require("value[mask] = (p / pmax) ** 4" in source_basis, "#701 baseline peak normalization drift")

    p_selected, z_star = _select_minimal_p()
    _require(p_selected == 9, "independent #732 p selection no longer yields p=9")
    _require(
        math.isclose(z_star, 1.602054456052944, rel_tol=0.0, abs_tol=2.0e-15),
        "independent #732 sign root drift",
    )

    scope = contract["representation_scope"]
    for key in (
        "same_basis_dimension",
        "same_vector_potential_form",
        "same_radial_support",
        "same_axial_support",
        "same_time_activation_semantics",
        "both_axial_bumps_peak_normalized_to_one",
    ):
        _require(scope[key] is True, f"expected positive representation fact drifted: {key}")
    for key in (
        "spatial_basis_function_identical",
        "unit_velocity_response_identical",
        "coefficient_coordinate_identity_preserved",
        "old_701_alpha_transfer_validated",
        "new_alpha_calibration_performed_by_732",
        "nonlinear_child_evaluated_by_732",
    ):
        _require(scope[key] is False, f"forbidden coordinate/child promotion: {key}")

    witness = independent_unit_response_witness()
    frozen_witness = contract["unit_response_witness"]
    for key, actual in (
        ("baseline_axial_z_factor_approx", witness["baseline_axial_z_factor"]),
        ("reshaped_axial_z_factor_approx", witness["reshaped_axial_z_factor"]),
        ("reshaped_to_baseline_magnitude_ratio_approx", witness["reshaped_to_baseline_magnitude_ratio"]),
    ):
        _require(
            math.isclose(float(frozen_witness[key]), actual, rel_tol=0.0, abs_tol=2.0e-14),
            f"unit-response witness drift: {key}",
        )

    rules = contract["promotion_rules"]
    _require(
        rules["732_pass_may_justify_one_same_dimension_nonlinear_replay"] is True,
        "#732 valid replay consequence was lost",
    )
    _require(
        rules["732_pass_may_claim_old_alpha_is_invariant_under_reshape"] is False,
        "old alpha was laundered across a changed basis vector",
    )
    _require(
        rules["next_child_must_freeze_alpha_semantics_before_result_evaluation"] is True,
        "next child may silently choose alpha after results",
    )
    _require(
        rules["rederiving_alpha_from_701_historical_tip_paths_remains_development_selection_data"] is True,
        "#712 selection-data boundary lost",
    )
    _require(
        rules["fresh_714_paths_must_not_be_used_to_retune_alpha_if_they_are_to_remain_post_freeze_evidence"] is True,
        "#721 fresh-path independence boundary lost",
    )
    _require(
        rules["same_dimension_must_not_be_relabelled_as_same_basis_vector"] is True,
        "basis dimension was confused with basis identity",
    )
    _require(
        rules["same_support_must_not_be_relabelled_as_same_unit_response"] is True,
        "support identity was confused with response identity",
    )
    _require(
        rules["visual_correspondence_requires_independent_visual_evidence"] is True,
        "visual claim firewall drift",
    )
    _require(
        rules["pde_validation_requires_cr001_held_out_acceptance"] is True,
        "PDE claim firewall drift",
    )

    for key, value in contract["claim_states"].items():
        _require(value is False, f"preflight illegally promotes {key}")

    truth = relocation_values.get("TRUTH")
    _require(isinstance(truth, dict), "#732 TRUTH boundary is no longer a literal dictionary")
    for key in (
        "canonical_velocity_changed",
        "saved_velocity_changed",
        "production_candidate_selected",
        "nonlinear_child_evaluated",
        "coefficient_alpha_fitted",
        "held_out_pde_residual_evaluated",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth.get(key) is False, f"#732 upstream truth promotion drift: {key}")

    _audit_cr001(root, contract["cr001_lock"])
    return {
        "status": "pass",
        "audited_upstream_head": EXPECTED_732_HEAD,
        "selected_inner_exponent_p": p_selected,
        "physical_abs_z_sign_change": z_star,
        "unit_response_witness": witness,
        "old_701_alpha_transfer_validated": False,
        "next_child_requires_explicit_frozen_alpha_semantics": True,
    }
