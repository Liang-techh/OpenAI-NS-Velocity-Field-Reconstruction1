"""CR002 fail-closed audit for current-Qs eta/matching-bridge scope.

This module governs representation and evidence transfer only. It does not alter
Kokuno candidate mathematics, CR001 thresholds, forcing, or scientific states.
"""

from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

CONTRACT = Path("configs/kokuno_current_qs_eta_bridge_scope.json")
SOURCE_MODULE = Path(
    "src/openai_ns_reconstruction/kokuno_current_exterior_qs_release2_state.py"
)
CONSTRAINTS = Path("configs/constraints.json")
PROJECT_STATUS = Path("project_status.json")

_CANONICAL_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}

_REQUIRED_CURRENT_TRUE = (
    "current_q_s_release2_endpoint_materialized",
    "current_q_s_uses_general_radial_identity",
    "current_absolute_I_state_consumed",
    "current_actual_M_M_eta_state_consumed",
    "current_actual_J_J_eta_state_consumed",
)

_REQUIRED_DOWNSTREAM_FALSE = (
    "current_l_minus_h_matching_bridge_materialized",
    "current_cartesian_terminal_multiplier_composed",
    "outer_global_leading_velocity_materialized",
    "unified_global_cartesian_velocity_export_ready",
    "matched_global_pressure_materialized",
    "restricted_forcing_materialized",
    "heldout_ns_residual_assessed",
    "same_protocol_comparable_to_st006",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise AssertionError(f"{path} must contain a JSON object")
    return value


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _literal_assignment(path: Path, name: str) -> Mapping[str, Any]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if any(
                isinstance(target, ast.Name) and target.id == name
                for target in node.targets
            ):
                value = ast.literal_eval(node.value)
                if not isinstance(value, dict):
                    raise AssertionError(f"{name} must be a literal dict")
                return value
    raise AssertionError(f"missing literal assignment {name} in {path}")


def _class_method_names(path: Path, class_name: str) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            return {
                child.name
                for child in node.body
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
            }
    raise AssertionError(f"missing class {class_name} in {path}")


def _string_literals(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return {
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected!r}, got {actual!r}")


def _mechanics_eta_bridge_witness(contract: Mapping[str, Any]) -> dict[str, Any]:
    """Autonomous logic witness; never Kokuno/OpenAI numerical evidence."""
    witness = contract["mechanics_witness"]
    q0 = float(witness["q0"])
    q_p = float(witness["q_p"])
    h = float(witness["h"])
    alpha = float(witness["alpha"])
    eta = [float(v) for v in witness["eta_samples"]]
    if not (q0 > 0.0 and q_p > 0.0 and 0.0 < h < 1.0 and eta):
        raise AssertionError(
            "mechanics witness parameters must define a positive matching problem"
        )
    k = 1.0 - h
    q = [q0 * math.exp(alpha * e) for e in eta]
    t = [math.log(value / q_p) / k for value in q]
    if not all(math.isfinite(value) for value in q + t):
        raise AssertionError("mechanics witness became non-finite")
    spread = max(t) - min(t)
    if spread <= 0.0:
        raise AssertionError("mechanics witness must exhibit eta-dependent matching time")

    representative_index = min(range(len(eta)), key=lambda i: abs(eta[i]))
    scalar_t = t[representative_index]
    terminal = [value * math.exp(-k * scalar_t) for value in q]
    mismatch = [abs(value - q_p) for value in terminal]
    if max(mismatch) <= 0.0:
        raise AssertionError(
            "scalarized mechanics witness unexpectedly matches every eta"
        )

    return {
        "eta": eta,
        "q": q,
        "pointwise_matching_time": t,
        "matching_time_spread": spread,
        "representative_eta": eta[representative_index],
        "representative_scalar_time": scalar_t,
        "terminal_after_representative_scalar_time": terminal,
        "max_abs_terminal_mismatch": max(mismatch),
    }


def audit(root: str | Path | None = None) -> dict[str, Any]:
    """Audit the exact #1225 eta-state / future bridge representation boundary."""
    repo = Path(root) if root is not None else Path(__file__).resolve().parents[2]
    contract = _load_json(repo / CONTRACT)
    constraints = _load_json(repo / CONSTRAINTS)
    status = _load_json(repo / PROJECT_STATUS)
    bindings = contract["bindings"]

    _require_equal(contract["schema_version"], 1, "contract schema")
    _require_equal(
        set(contract["source_classification"]),
        _CANONICAL_CLASSES,
        "four-class provenance vocabulary",
    )
    for class_name, entries in contract["source_classification"].items():
        if (
            not isinstance(entries, list)
            or not entries
            or not all(isinstance(item, str) and item.strip() for item in entries)
        ):
            raise AssertionError(
                f"source_classification.{class_name} must be a nonempty list of strings"
            )

    _require_equal(
        _git_blob_sha1(repo / CONSTRAINTS),
        bindings["canonical_constraints_blob_sha1"],
        "canonical CR001 constraints blob",
    )
    _require_equal(
        _git_blob_sha1(repo / PROJECT_STATUS),
        bindings["project_status_blob_sha1"],
        "project status blob",
    )
    _require_equal(
        _git_blob_sha1(repo / SOURCE_MODULE),
        bindings["dependency_module_blob_sha1"],
        "PR #1225 current-Qs module blob",
    )

    source_truth = _literal_assignment(repo / SOURCE_MODULE, "_TRUTH_UPDATES")
    for key in _REQUIRED_CURRENT_TRUE:
        _require_equal(source_truth.get(key), True, f"#1225 current truth {key}")
    for key in _REQUIRED_DOWNSTREAM_FALSE:
        _require_equal(source_truth.get(key), False, f"#1225 downstream truth {key}")

    methods = _class_method_names(
        repo / SOURCE_MODULE, "KokunoCurrentExteriorQsRelease2State"
    )
    if "state" not in methods or "report" not in methods:
        raise AssertionError("#1225 must expose the eta-resolved state/report surface")
    if "velocity" in methods:
        raise AssertionError(
            "#1225 unexpectedly exposes a velocity method; bridge delivery scope changed"
        )
    literals = _string_literals(repo / SOURCE_MODULE)
    for key in (
        "Q_s_current_release2_end",
        "Q_s_source_ideal_release2_end",
        "Q_s_current_minus_source_ideal",
    ):
        if key not in literals:
            raise AssertionError(f"#1225 eta-state output key missing: {key}")

    truth = contract["truth_state"]
    for key in (
        "current_q_s_release2_endpoint_materialized",
        "current_q_s_exposed_as_eta_resolved_state",
    ):
        _require_equal(truth.get(key), True, f"contract truth {key}")
    for key in (
        "current_q_s_proven_eta_independent",
        "current_scalar_l_minus_h_matching_length_materialized",
        "eta_dependent_cartesian_matching_bridge_materialized",
        "current_cartesian_terminal_multiplier_composed",
        "kokuno_unified_global_velocity_export_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(truth.get(key), False, f"contract truth {key}")

    for key, value in contract["evidence_transfer"].items():
        _require_equal(value, False, f"forbidden evidence transfer {key}")
    if len(contract["promotion_prerequisites"]) < 5:
        raise AssertionError(
            "promotion_prerequisites must retain eta-structure, bridge, Cartesian, replay, and state gates"
        )

    witness = contract["mechanics_witness"]
    _require_equal(
        witness["role"],
        "autonomous_representation_mechanics_only",
        "mechanics witness role",
    )
    mechanics = _mechanics_eta_bridge_witness(contract)

    canonical = contract["cr001_guardrails"]
    _require_equal(constraints["nu"], canonical["nu"], "CR001 nu")
    _require_equal(
        constraints["domain"]["physical"],
        canonical["physical_domain"],
        "CR001 physical domain",
    )
    _require_equal(
        constraints["domain"]["evaluation_box"],
        canonical["evaluation_box"],
        "CR001 evaluation box",
    )
    _require_equal(
        constraints["domain"]["support"], canonical["support"], "CR001 support"
    )
    _require_equal(
        constraints["domain"]["time_interval"],
        canonical["time_interval"],
        "CR001 time interval",
    )
    _require_equal(
        constraints["forcing"]["mode"],
        canonical["forcing_mode"],
        "CR001 forcing mode",
    )
    if (
        "No residual-dependent basis or pointwise free force"
        not in constraints["forcing"]["restriction"]
    ):
        raise AssertionError(
            "CR001 residual-defined/free pointwise force prohibition drifted"
        )
    _require_equal(
        constraints["nontriviality"]["reference_energy"],
        canonical["reference_energy"],
        "CR001 reference energy",
    )
    _require_equal(
        constraints["nontriviality"]["reference_energy_abs_tolerance"],
        canonical["reference_energy_abs_tolerance"],
        "CR001 energy tolerance",
    )
    _require_equal(
        constraints["validation"]["seed"],
        canonical["validation_seed"],
        "CR001 validation seed",
    )
    _require_equal(
        constraints["validation"]["held_out_points"],
        canonical["held_out_points"],
        "CR001 held-out points",
    )
    _require_equal(
        constraints["validation"]["derivative_steps"],
        canonical["derivative_steps"],
        "CR001 derivative steps",
    )
    _require_equal(
        constraints["validation"]["quadrature_orders_per_axis"],
        canonical["quadrature_orders_per_axis"],
        "CR001 quadrature orders",
    )
    thresholds = constraints["validation"]["thresholds"]
    _require_equal(
        thresholds["pde_residual_max"],
        canonical["momentum_max_threshold"],
        "CR001 momentum max threshold",
    )
    _require_equal(
        thresholds["pde_residual_L2"],
        canonical["momentum_l2_threshold"],
        "CR001 momentum L2 threshold",
    )
    _require_equal(
        thresholds["divergence_max"],
        canonical["divergence_max_threshold"],
        "CR001 divergence max threshold",
    )
    _require_equal(
        thresholds["divergence_L2"],
        canonical["divergence_l2_threshold"],
        "CR001 divergence L2 threshold",
    )

    states = status["states"]
    independence = contract["canonical_independence"]
    _require_equal(
        states["velocity_export_ready"],
        independence["canonical_eq45_velocity_export_ready_expected"],
        "canonical Eq45 velocity_export_ready",
    )
    _require_equal(
        states["visual_correspondence_verified"],
        independence["canonical_visual_correspondence_verified_expected"],
        "canonical visual correspondence",
    )
    _require_equal(
        states["pde_validated"],
        independence["canonical_pde_validated_expected"],
        "canonical PDE state",
    )
    _require_equal(
        states["paper_exact"],
        independence["canonical_paper_exact_expected"],
        "canonical paper exactness",
    )
    _require_equal(
        states["openai_field_identified"],
        independence["canonical_openai_field_identified_expected"],
        "canonical OpenAI-field identity",
    )
    _require_equal(
        independence[
            "kokuno_bridge_incompleteness_must_not_revoke_canonical_eq45_callability"
        ],
        True,
        "delivery independence",
    )

    return {
        "contract_id": contract["contract_id"],
        "dependency_pr": bindings["dependency_pr"],
        "dependency_exact_head": bindings["dependency_exact_head"],
        "current_q_s_release2_endpoint_materialized": True,
        "current_q_s_eta_independence_proven": False,
        "current_scalar_matching_length_materialized": False,
        "current_cartesian_matching_bridge_materialized": False,
        "mechanics_only_matching_time_spread": mechanics["matching_time_spread"],
        "mechanics_only_max_abs_terminal_mismatch_after_scalarization": mechanics[
            "max_abs_terminal_mismatch"
        ],
        "canonical_eq45_velocity_export_ready": states["velocity_export_ready"],
        "kokuno_unified_global_velocity_export_ready": truth[
            "kokuno_unified_global_velocity_export_ready"
        ],
    }


def main() -> int:
    print(json.dumps(audit(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
