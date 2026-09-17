"""Govern residual-norm semantics for the support-connected Eq. (4.5) candidate.

The independent supported-child vorticity diagnostic intentionally evaluates a
small fixed probe set. Its RMS, maximum, and term-normalized RMS are useful
representation/PDE diagnostics, but they are not interchangeable with the CR001
registered validation sampling and volume-weighted spatial L2 acceptance metrics.
This module keeps those evidence scopes separate without making PDE acceptance a
prerequisite for callable/saveable/exportable velocity delivery.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

CANONICAL_SOURCE_VOCABULARY = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _as_float(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not (result == result and abs(result) != float("inf")):
        raise ValueError(f"{name} must be finite")
    return result


def audit_supported_residual_norm_scope(
    contract: Mapping[str, Any],
    constraints: Mapping[str, Any],
    diagnostic: Mapping[str, Any],
    candidate_truth: Mapping[str, Any],
) -> dict[str, Any]:
    """Fail closed on metric-scope drift or scientific-state promotion."""

    _require(contract.get("schema") == "eq45_supported_residual_norm_scope_v1", "contract schema drifted")
    _require(
        contract.get("task_id") == "CR002-EQ45-SUPPORTED-RESIDUAL-NORM-SCOPE-020",
        "task identity drifted",
    )

    vocabulary = set(contract.get("source_vocabulary", ()))
    _require(vocabulary == CANONICAL_SOURCE_VOCABULARY, "source vocabulary drifted")
    classifications = contract.get("classification")
    _require(isinstance(classifications, Mapping), "classification map missing")
    _require(set(classifications.values()) <= vocabulary, "unknown source classification")
    _require(classifications.get("velocity_delivery") == "user_requirement", "velocity delivery misclassified")
    _require(classifications.get("cr001_validation_norms") == "autonomous_design", "CR001 norms misclassified")
    _require(classifications.get("cr001_thresholds") == "autonomous_design", "CR001 thresholds misclassified")
    _require(
        classifications.get("fixed_probe_vorticity_diagnostic") == "autonomous_design",
        "fixed-probe diagnostic misclassified",
    )
    _require(
        classifications.get("openai_visual_correspondence") == "pending_unknown",
        "visual correspondence cannot be promoted by a PDE diagnostic",
    )

    registered = contract.get("registered_cr001_contract")
    _require(isinstance(registered, Mapping), "registered CR001 contract missing")
    _require(constraints.get("experiment_id") == registered.get("experiment_id"), "CR001 experiment identity drifted")
    validation = constraints.get("validation")
    _require(isinstance(validation, Mapping), "validation contract missing")
    _require(int(validation.get("held_out_points", -1)) == int(registered.get("held_out_points", -2)), "held-out sample count drifted")
    _require(list(validation.get("derivative_steps", ())) == list(registered.get("derivative_steps", ())), "derivative ladder drifted")
    _require(list(validation.get("norms", ())) == list(registered.get("norms", ())), "registered residual norms drifted")
    thresholds = validation.get("thresholds")
    _require(isinstance(thresholds, Mapping), "validation thresholds missing")
    _require(
        _as_float(thresholds.get("pde_residual_max"), "pde_residual_max")
        == _as_float(registered.get("pde_residual_max"), "registered pde_residual_max"),
        "PDE max threshold drifted",
    )
    _require(
        _as_float(thresholds.get("pde_residual_L2"), "pde_residual_L2")
        == _as_float(registered.get("pde_residual_L2"), "registered pde_residual_L2"),
        "PDE L2 threshold drifted",
    )

    identity = contract.get("candidate_identity")
    _require(isinstance(identity, Mapping), "candidate identity contract missing")
    _require(identity.get("classification") == "autonomous_design", "candidate identity misclassified")
    _require(diagnostic.get("parent_sha256") == identity.get("parent_sha256"), "diagnostic parent identity drifted")
    _require(
        diagnostic.get("supported_child_sha256") == identity.get("supported_child_sha256"),
        "diagnostic child identity drifted",
    )

    expected_diag = contract.get("supported_child_diagnostic")
    _require(isinstance(expected_diag, Mapping), "supported-child diagnostic contract missing")
    _require(diagnostic.get("schema") == expected_diag.get("schema"), "diagnostic schema drifted")
    _require(diagnostic.get("claim_scope") == expected_diag.get("claim_scope"), "diagnostic claim scope drifted")
    _require(int(diagnostic.get("probe_count", -1)) == int(expected_diag.get("probe_count", -2)), "diagnostic probe count drifted")
    _require(
        int(diagnostic.get("probe_count", 0)) < int(validation.get("held_out_points", 0)),
        "fixed-probe diagnostic must remain distinct from registered held-out sampling",
    )
    _require(list(diagnostic.get("spatial_steps", ())) == list(validation.get("derivative_steps", ())), "diagnostic derivative ladder drifted")
    diag_thresholds = diagnostic.get("pde_thresholds")
    _require(isinstance(diag_thresholds, Mapping), "diagnostic threshold metadata missing")
    _require(_as_float(diag_thresholds.get("max"), "diagnostic max threshold") == _as_float(thresholds.get("pde_residual_max"), "registered max threshold"), "diagnostic threshold metadata changed the registered max")
    _require(_as_float(diag_thresholds.get("L2"), "diagnostic L2 threshold") == _as_float(thresholds.get("pde_residual_L2"), "registered L2 threshold"), "diagnostic threshold metadata changed the registered L2")

    for key in (
        "finest_supported_zero_force_rms",
        "finest_supported_zero_force_max",
        "finest_supported_term_normalized_rms",
    ):
        _as_float(diagnostic.get(key), key)

    _require(diagnostic.get("formal_pde_gate_assessed") is False, "fixed-probe diagnostic cannot assess the formal PDE gate")
    reason = str(diagnostic.get("reason_formal_pde_gate_unassessed", ""))
    _require("probe RMS is not the preregistered volume-weighted spatial L2" in reason, "formal-gate limitation was dropped")
    _require(
        diagnostic.get("supported_force_coefficients_status")
        == expected_diag.get("supported_force_coefficients_status"),
        "supported-child force binding status drifted",
    )
    _require(diagnostic.get("untapered_fitted_force_transferred") is False, "untapered force fit cannot be inherited")
    _require(diagnostic.get("forcing_family") == constraints.get("forcing", {}).get("mode"), "forcing family drifted")
    _require(diagnostic.get("forcing_family") == "restricted_two_parameter_family", "free/residual-defined force is forbidden")

    policy = contract.get("policy")
    _require(isinstance(policy, Mapping), "norm-scope policy missing")
    false_keys = (
        "fixed_probe_rms_may_satisfy_registered_L2",
        "fixed_probe_max_may_satisfy_registered_max",
        "term_normalized_rms_has_registered_acceptance_threshold",
        "zero_force_pressure_free_check_may_set_pde_validated",
        "green_ci_may_set_pde_validated",
        "formal_pde_failure_or_unassessed_blocks_velocity_export",
        "residual_defined_or_pointwise_free_force_allowed",
        "posthoc_threshold_relaxation_allowed",
    )
    for key in false_keys:
        _require(policy.get(key) is False, f"illegal policy promotion: {key}")

    _require(candidate_truth.get("velocity_export_ready") is True, "supported child lost export readiness")
    _require(candidate_truth.get("pde_validated") is False, "PDE state promoted without formal validation")
    _require(candidate_truth.get("visualization_ready") is False, "visualization readiness promoted by norm audit")
    _require(candidate_truth.get("visual_correspondence_verified") is False, "visual correspondence promoted by norm audit")
    for key in ("paper_exact", "openai_field_identified", "blowup_proved"):
        _require(candidate_truth.get(key) is False, f"scientific state illegally promoted: {key}")

    truth = contract.get("truth_boundary")
    _require(isinstance(truth, Mapping), "truth boundary missing")
    for key, expected in truth.items():
        _require(candidate_truth.get(key) is expected, f"candidate truth boundary drifted for {key}")

    return {
        "task_id": contract["task_id"],
        "candidate_sha256": diagnostic["supported_child_sha256"],
        "probe_count": diagnostic["probe_count"],
        "registered_held_out_points": validation["held_out_points"],
        "fixed_probe_metrics_are_formal_pde_acceptance": False,
        "formal_pde_gate_assessed": False,
        "velocity_export_ready": True,
        "pde_validated": False,
    }


def audit_current_repository(root: str | Path | None = None) -> dict[str, Any]:
    """Replay the current supported-child diagnostic and audit its evidence scope."""
    repo = Path(root) if root is not None else Path(__file__).resolve().parents[2]
    contract = json.loads((repo / "configs/eq45_supported_residual_norm_scope.json").read_text(encoding="utf-8"))
    constraints = json.loads((repo / "configs/constraints.json").read_text(encoding="utf-8"))

    from .constrained_eq45_candidate import Eq45VelocityCandidate
    from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
    from .constrained_eq45_supported_vorticity import audit_supported_eq45

    diagnostic = audit_supported_eq45(constraints_path=repo / "configs/constraints.json")
    parent = Eq45VelocityCandidate.load_json(repo / "artifacts/constrained/eq45_velocity_candidate_seed.json")
    child = Eq45SupportedVelocityCandidate(parent=parent)
    candidate_truth = child.to_dict()["truth_boundary"]
    return audit_supported_residual_norm_scope(contract, constraints, diagnostic, candidate_truth)


def main() -> None:
    print(json.dumps(audit_current_repository(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
