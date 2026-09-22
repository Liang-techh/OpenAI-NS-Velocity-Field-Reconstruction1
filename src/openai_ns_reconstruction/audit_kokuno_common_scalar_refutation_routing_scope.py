"""CR002 fail-closed routing audit for A1 #1242 common-scalar refutation.

This module governs evidence/representation transfer only.  It does not alter
candidate mathematics, CR001 gates, forcing, or any scientific acceptance flag.
"""
from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

CONTRACT = Path("configs/kokuno_common_scalar_refutation_routing_scope.json")
SOURCE_MODULE = Path(
    "src/openai_ns_reconstruction/kokuno_current_bridge_common_scalar_refutation.py"
)
CONSTRAINTS = Path("configs/constraints.json")
PROJECT_STATUS = Path("project_status.json")

_CLASSES = {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}
_PARENT_FALSE = (
    "full_eta_interval_entry_above_qp_established",
    "full_eta_interval_target_totality_established",
    "full_eta_interval_root_uniqueness_established",
    "full_eta_interval_transversality_established",
    "smooth_eta_target_time_map_established",
    "full_eta_interval_common_scalar_bridge_length_established",
    "current_l_minus_h_matching_bridge_materialized",
    "eta_dependent_cartesian_matching_boundary_materialized",
)
_PARENT_TRUTH_FALSE = (
    "current_cartesian_terminal_multiplier_composed",
    "outer_global_leading_velocity_materialized",
    "unified_global_cartesian_velocity_export_ready",
    "heldout_ns_residual_assessed",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AssertionError(f"{path} must contain a JSON object")
    return value


def _blob(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _literal(path: Path, name: str) -> Any:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == name for target in node.targets
        ):
            return ast.literal_eval(node.value)
    raise AssertionError(f"missing literal assignment {name} in {path}")


def _methods(path: Path, class_name: str) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            return {
                child.name
                for child in node.body
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
            }
    raise AssertionError(f"missing class {class_name}")


def _eq(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected!r}, got {actual!r}")


def _require_parent_false(report: Mapping[str, Any]) -> None:
    for key in _PARENT_FALSE:
        if report.get(key) is not False:
            raise AssertionError(f"dependency report illegally promotes {key}")
    truth = report.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise AssertionError("dependency report must include truth_boundary mapping")
    for key in _PARENT_TRUTH_FALSE:
        if truth.get(key) is not False:
            raise AssertionError(f"dependency truth illegally promotes {key}")


def route_refutation_report(report: Mapping[str, Any]) -> dict[str, Any]:
    """Map either legal A1 #1242 outcome to fail-closed representation routing."""
    if not isinstance(report, Mapping):
        raise TypeError("report must be a mapping")
    if report.get("finite_witness_is_sufficient_only_for_refutation") is not True:
        raise AssertionError("dependency one-way witness semantics drifted")
    if report.get("absence_of_refutation_does_not_establish_common_scalar_bridge") is not True:
        raise AssertionError("dependency negative-outcome semantics drifted")
    _require_parent_false(report)

    refuted = report.get("common_scalar_bridge_refuted_by_disjoint_unique_root_brackets")
    pair = report.get("disjoint_root_pair")
    gap = report.get("disjoint_gap_lower_bound")
    if not isinstance(refuted, bool):
        raise AssertionError("dependency refutation outcome must be boolean")
    if not isinstance(gap, (int, float)) or not math.isfinite(float(gap)):
        raise AssertionError("dependency gap must be finite")

    if refuted:
        if not isinstance(pair, Mapping):
            raise AssertionError("positive refutation requires a witness pair")
        earlier = pair.get("earlier_root_bracket")
        later = pair.get("later_root_bracket")
        if not (
            isinstance(earlier, list)
            and len(earlier) == 2
            and isinstance(later, list)
            and len(later) == 2
        ):
            raise AssertionError("positive refutation requires two root brackets")
        e0, e1 = map(float, earlier)
        l0, l1 = map(float, later)
        if not all(math.isfinite(v) for v in (e0, e1, l0, l1)):
            raise AssertionError("root brackets must be finite")
        if not (e0 <= e1 < l0 <= l1):
            raise AssertionError("positive refutation requires ordered disjoint brackets")
        if float(gap) <= 0.0 or float(pair.get("disjoint_gap_lower_bound", 0.0)) <= 0.0:
            raise AssertionError("positive refutation requires a positive gap")
        status = "refuted_for_exact_current_candidate_identity"
        next_routes = (
            "new_upstream_candidate_identity_then_retest_scalar_compatibility",
            "new_eta_dependent_autonomous_bridge_identity_after_continuum_gates",
        )
    else:
        if pair is not None or float(gap) != 0.0:
            raise AssertionError("negative/inconclusive result must not carry a refutation pair")
        status = "unresolved_not_established"
        next_routes = (
            "continuum_common_scalar_proof_before_promotion",
            "new_explicit_representation_identity_if_scalar_route_is_abandoned",
        )

    return {
        "common_scalar_bridge_status": status,
        "dependency_positive_refutation": refuted,
        "allowed_next_representation_routes": list(next_routes),
        "current_common_scalar_bridge_established": False,
        "current_l_minus_h_matching_bridge_materialized": False,
        "eta_dependent_cartesian_matching_boundary_materialized": False,
        "current_cartesian_terminal_multiplier_composed": False,
        "kokuno_unified_global_velocity_export_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def _mechanics_refutation_witness(contract: Mapping[str, Any]) -> dict[str, Any]:
    witness = contract["mechanics_witness"]
    _eq(witness["role"], "autonomous_representation_mechanics_only", "mechanics role")
    target = float(witness["target"])
    if target <= 0.0:
        raise AssertionError("mechanics target must be positive")
    root_a = math.log(2.0 / target)
    root_b = 0.5 * math.log(2.0 / target)
    if not (root_a > 0.0 and root_b > 0.0 and not math.isclose(root_a, root_b)):
        raise AssertionError("autonomous mechanics witness must have distinct positive roots")
    return {
        "target": target,
        "root_curve_a": root_a,
        "root_curve_b": root_b,
        "absolute_root_gap": abs(root_a - root_b),
        "common_scalar_time_exists_for_both_curves": False,
        "role": witness["role"],
    }


def audit(root: str | Path | None = None) -> dict[str, Any]:
    repo = Path(root) if root is not None else Path(__file__).resolve().parents[2]
    contract = _load_json(repo / CONTRACT)
    constraints = _load_json(repo / CONSTRAINTS)
    status = _load_json(repo / PROJECT_STATUS)
    bindings = contract["bindings"]

    _eq(contract["schema_version"], 1, "contract schema")
    _eq(set(contract["source_classification"]), _CLASSES, "four-class provenance")
    for name, entries in contract["source_classification"].items():
        if not isinstance(entries, list) or not entries or not all(
            isinstance(value, str) and value.strip() for value in entries
        ):
            raise AssertionError(f"source_classification.{name} must be nonempty strings")

    _eq(_blob(repo / SOURCE_MODULE), bindings["dependency_module_blob_sha1"], "#1242 module blob")
    _eq(_blob(repo / CONSTRAINTS), bindings["canonical_constraints_blob_sha1"], "CR001 blob")
    _eq(_blob(repo / PROJECT_STATUS), bindings["project_status_blob_sha1"], "status blob")
    _eq(bindings["dependency_pr"], 1242, "dependency PR")
    _eq(
        bindings["dependency_exact_head"],
        "7e987ff0a3bea743fe3f0fe9d02d5b3c8c593ea3",
        "dependency exact head",
    )

    _eq(_literal(repo / SOURCE_MODULE, "_WITNESS_COUNT"), 65, "#1242 witness count")
    _eq(_literal(repo / SOURCE_MODULE, "_ROOT_EXPANSIONS"), 16, "#1242 root expansions")
    _eq(_literal(repo / SOURCE_MODULE, "_ROOT_BISECTIONS"), 96, "#1242 root bisections")
    _eq(_literal(repo / SOURCE_MODULE, "_MARGIN_REL"), 2.0 ** -40, "#1242 sign margin")
    numerical = _literal(repo / SOURCE_MODULE, "_NUMERICAL_REALIZATION")
    if "absence of such a pair does not establish a common bridge" not in numerical["promotion_policy"]:
        raise AssertionError("#1242 one-way promotion policy drifted")

    source_truth = _literal(repo / SOURCE_MODULE, "_TRUTH_UPDATES")
    _eq(
        source_truth.get("current_common_scalar_bridge_witness_check_materialized"),
        True,
        "witness materialized",
    )
    for key in _PARENT_TRUTH_FALSE:
        _eq(source_truth.get(key), False, f"#1242 truth {key}")

    methods = _methods(repo / SOURCE_MODULE, "KokunoCurrentBridgeCommonScalarRefutation")
    for name in (
        "unique_root_certificate",
        "unique_root_bracket",
        "refutation_report",
        "save_configuration",
        "load_configuration",
    ):
        if name not in methods:
            raise AssertionError(f"#1242 missing interface {name}")
    if "velocity" in methods:
        raise AssertionError("#1242 unexpectedly exposes velocity")

    identity = contract["identity_policy"]
    for key in (
        "positive_refutation_forbids_in_place_mean_min_max_or_representative_eta_repair",
        "upstream_current_state_revision_requires_new_candidate_identity",
        "eta_dependent_bridge_requires_new_representation_identity",
        "save_load_and_downstream_audits_must_bind_the_new_identity",
    ):
        _eq(identity.get(key), True, f"identity policy {key}")
    for key in (
        "prior_refutation_receipt_transfers_to_revised_candidate",
        "source_scalar_stage_label_transfers_to_eta_dependent_autonomous_design",
    ):
        _eq(identity.get(key), False, f"identity policy {key}")

    for key, value in contract["evidence_transfer"].items():
        _eq(value, False, f"forbidden evidence transfer {key}")
    for key, value in contract["truth_state"].items():
        if key == "current_common_scalar_bridge_witness_check_materialized":
            _eq(value, True, f"truth {key}")
        else:
            _eq(value, False, f"truth {key}")

    mechanics = _mechanics_refutation_witness(contract)

    canonical = contract["cr001_guardrails"]
    _eq(constraints["nu"], canonical["nu"], "nu")
    _eq(constraints["domain"]["physical"], canonical["physical_domain"], "physical domain")
    _eq(constraints["domain"]["evaluation_box"], canonical["evaluation_box"], "evaluation box")
    _eq(constraints["domain"]["support"], canonical["support"], "support")
    _eq(constraints["domain"]["time_interval"], canonical["time_interval"], "time interval")
    _eq(constraints["forcing"]["mode"], canonical["forcing_mode"], "forcing mode")
    if "No residual-dependent basis or pointwise free force" not in constraints["forcing"]["restriction"]:
        raise AssertionError("residual-defined/free pointwise force prohibition drifted")
    _eq(constraints["nontriviality"]["reference_energy"], canonical["reference_energy"], "reference energy")
    _eq(
        constraints["nontriviality"]["reference_energy_abs_tolerance"],
        canonical["reference_energy_abs_tolerance"],
        "energy tolerance",
    )
    validation = constraints["validation"]
    _eq(validation["seed"], canonical["validation_seed"], "validation seed")
    _eq(validation["held_out_points"], canonical["held_out_points"], "held-out points")
    _eq(validation["derivative_steps"], canonical["derivative_steps"], "FD ladder")
    _eq(
        validation["quadrature_orders_per_axis"],
        canonical["quadrature_orders_per_axis"],
        "quadrature ladder",
    )
    gates = validation["thresholds"]
    _eq(gates["pde_residual_max"], canonical["momentum_max_threshold"], "momentum max")
    _eq(gates["pde_residual_L2"], canonical["momentum_l2_threshold"], "momentum L2")
    _eq(gates["divergence_max"], canonical["divergence_max_threshold"], "divergence max")
    _eq(gates["divergence_L2"], canonical["divergence_l2_threshold"], "divergence L2")

    states = status["states"]
    independent = contract["canonical_independence"]
    _eq(
        states["velocity_export_ready"],
        independent["canonical_eq45_velocity_export_ready_expected"],
        "Eq45 export",
    )
    _eq(
        states["visual_correspondence_verified"],
        independent["canonical_visual_correspondence_verified_expected"],
        "visual state",
    )
    _eq(states["pde_validated"], independent["canonical_pde_validated_expected"], "PDE state")
    _eq(states["paper_exact"], independent["canonical_paper_exact_expected"], "paper exactness")
    _eq(
        states["openai_field_identified"],
        independent["canonical_openai_field_identified_expected"],
        "OpenAI identity",
    )

    return {
        "contract_id": contract["contract_id"],
        "dependency_pr": bindings["dependency_pr"],
        "dependency_exact_head": bindings["dependency_exact_head"],
        "dependency_module_blob_sha1": bindings["dependency_module_blob_sha1"],
        "mechanics_witness": mechanics,
        "truth_state": contract["truth_state"],
        "canonical_states": {
            "velocity_export_ready": states["velocity_export_ready"],
            "visual_correspondence_verified": states["visual_correspondence_verified"],
            "pde_validated": states["pde_validated"],
            "paper_exact": states["paper_exact"],
            "openai_field_identified": states["openai_field_identified"],
        },
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True, allow_nan=False))
