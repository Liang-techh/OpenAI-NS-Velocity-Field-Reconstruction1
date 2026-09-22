"""CR002 fail-closed audit for source-ideal Q_s evidence transfer.

This module governs *evidence scope*.  It does not alter Kokuno candidate
mathematics, CR001 thresholds, forcing, or any scientific readiness state.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

CONTRACT = Path("configs/kokuno_source_ideal_qs_evidence_scope.json")
SOURCE_MODULE = Path("src/openai_ns_reconstruction/kokuno_public_ideal_exterior_qs_schedule.py")
CONSTRAINTS = Path("configs/constraints.json")
PROJECT_STATUS = Path("project_status.json")

_CANONICAL_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}

_REQUIRED_SOURCE_FALSE = (
    "current_q_s_release2_endpoint_materialized",
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
            if any(isinstance(target, ast.Name) and target.id == name for target in node.targets):
                value = ast.literal_eval(node.value)
                if not isinstance(value, dict):
                    raise AssertionError(f"{name} must be a literal dict")
                return value
    raise AssertionError(f"missing literal assignment {name} in {path}")


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected!r}, got {actual!r}")


def audit(root: str | Path | None = None) -> dict[str, Any]:
    """Audit the exact source/current evidence boundary and CR001 invariants."""
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
        if not isinstance(entries, list) or not entries or not all(isinstance(item, str) and item.strip() for item in entries):
            raise AssertionError(f"source_classification.{class_name} must be a nonempty list of strings")

    _require_equal(
        _git_blob_sha1(repo / CONSTRAINTS),
        bindings["canonical_constraints_blob_sha1"],
        "canonical CR001 constraints blob",
    )
    _require_equal(
        _git_blob_sha1(repo / SOURCE_MODULE),
        bindings["dependency_module_blob_sha1"],
        "PR #1217 source-ideal module blob",
    )

    source_truth = _literal_assignment(repo / SOURCE_MODULE, "_TRUTH_UPDATES")
    for key in (
        "public_source_ideal_exterior_qs_schedule_materialized",
        "public_source_ideal_qs_release2_endpoint_materialized",
        "public_source_ideal_l_minus_h_matching_length_materialized",
    ):
        _require_equal(source_truth.get(key), True, f"#1217 source-ideal truth {key}")
    for key in _REQUIRED_SOURCE_FALSE:
        _require_equal(source_truth.get(key), False, f"#1217 current/global truth {key}")

    truth = contract["truth_state"]
    for key in (
        "source_ideal_qs_schedule_materialized",
        "source_ideal_release2_endpoint_materialized",
        "source_ideal_t_match_materialized",
    ):
        _require_equal(truth.get(key), True, f"contract truth {key}")
    for key in (
        "current_candidate_qs_release2_endpoint_materialized",
        "current_candidate_t_match_materialized",
        "current_cartesian_matching_bridge_composed",
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
        raise AssertionError("promotion_prerequisites must retain the current-state, Cartesian, replay, and domain gates")

    canonical = contract["cr001_guardrails"]
    _require_equal(constraints["nu"], canonical["nu"], "CR001 nu")
    _require_equal(constraints["domain"]["physical"], canonical["physical_domain"], "CR001 physical domain")
    _require_equal(constraints["domain"]["evaluation_box"], canonical["evaluation_box"], "CR001 evaluation box")
    _require_equal(constraints["domain"]["support"], canonical["support"], "CR001 support")
    _require_equal(constraints["domain"]["time_interval"], canonical["time_interval"], "CR001 time interval")
    _require_equal(constraints["forcing"]["mode"], canonical["forcing_mode"], "CR001 forcing mode")
    if "No residual-dependent basis or pointwise free force" not in constraints["forcing"]["restriction"]:
        raise AssertionError("CR001 residual-defined/free pointwise force prohibition drifted")
    _require_equal(constraints["nontriviality"]["reference_energy"], canonical["reference_energy"], "CR001 reference energy")
    _require_equal(
        constraints["nontriviality"]["reference_energy_abs_tolerance"],
        canonical["reference_energy_abs_tolerance"],
        "CR001 energy tolerance",
    )
    _require_equal(constraints["validation"]["seed"], canonical["validation_seed"], "CR001 validation seed")
    _require_equal(constraints["validation"]["held_out_points"], canonical["held_out_points"], "CR001 held-out points")
    _require_equal(constraints["validation"]["derivative_steps"], canonical["derivative_steps"], "CR001 derivative steps")
    _require_equal(
        constraints["validation"]["quadrature_orders_per_axis"],
        canonical["quadrature_orders_per_axis"],
        "CR001 quadrature orders",
    )
    thresholds = constraints["validation"]["thresholds"]
    _require_equal(thresholds["pde_residual_max"], canonical["momentum_max_threshold"], "CR001 momentum max threshold")
    _require_equal(thresholds["pde_residual_L2"], canonical["momentum_l2_threshold"], "CR001 momentum L2 threshold")
    _require_equal(thresholds["divergence_max"], canonical["divergence_max_threshold"], "CR001 divergence max threshold")
    _require_equal(thresholds["divergence_L2"], canonical["divergence_l2_threshold"], "CR001 divergence L2 threshold")

    states = status["states"]
    independence = contract["canonical_independence"]
    _require_equal(states["velocity_export_ready"], independence["canonical_eq45_velocity_export_ready_expected"], "canonical Eq45 velocity_export_ready")
    _require_equal(states["visual_correspondence_verified"], independence["canonical_visual_correspondence_verified_expected"], "canonical visual correspondence")
    _require_equal(states["pde_validated"], independence["canonical_pde_validated_expected"], "canonical PDE state")
    _require_equal(states["paper_exact"], independence["canonical_paper_exact_expected"], "canonical paper exactness")
    _require_equal(states["openai_field_identified"], independence["canonical_openai_field_identified_expected"], "canonical OpenAI-field identity")
    _require_equal(independence["kokuno_incompleteness_must_not_revoke_canonical_eq45_callability"], True, "delivery independence")

    return {
        "contract_id": contract["contract_id"],
        "dependency_pr": bindings["dependency_pr"],
        "dependency_exact_head": bindings["dependency_exact_head"],
        "source_ideal_reference_materialized": True,
        "current_candidate_qs_materialized": False,
        "source_ideal_to_current_transfer_allowed": False,
        "canonical_eq45_velocity_export_ready": states["velocity_export_ready"],
        "kokuno_unified_global_velocity_export_ready": truth["kokuno_unified_global_velocity_export_ready"],
    }


def main() -> int:
    print(json.dumps(audit(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
