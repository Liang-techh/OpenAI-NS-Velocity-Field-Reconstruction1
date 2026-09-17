"""Fail-closed scope governance for Eq45 velocity export readiness.

The integrated Eq45 object is callable and serializable through its own
``at_points`` evaluator. That candidate-local fact is intentionally kept
separate from the repository's package-level unified/default velocity binding.
This module audits only that scope boundary; it does not evaluate PDE residuals
or visual correspondence.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

_REQUIRED_FALSE_EQ45_CLAIMS = {
    "physical_support_validated",
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
}


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def audit_eq45_export_scope_contract(
    contract: Mapping[str, Any], *, repo_root: str | Path
) -> dict[str, Any]:
    """Audit candidate-local versus package-unified velocity readiness.

    A pass means the current repository metadata does not confuse the callable
    Eq45 artifact with the package-level default velocity implementation. It is
    not evidence of visualization correspondence, PDE validity, paper-exact
    identity, OpenAI hidden-field recovery, or blow-up.
    """

    root = Path(repo_root)
    _require(
        contract.get("schema") == "eq45_export_scope_contract_v1",
        "unsupported Eq45 export scope schema",
    )

    package = contract.get("package_unified_velocity")
    eq45 = contract.get("eq45_candidate")
    semantics = contract.get("state_semantics")
    _require(isinstance(package, Mapping), "missing package_unified_velocity")
    _require(isinstance(eq45, Mapping), "missing eq45_candidate")
    _require(isinstance(semantics, Mapping), "missing state_semantics")

    delivery_path = root / str(package.get("delivery_contract", ""))
    state_path = root / str(semantics.get("delivery_state_contract", ""))
    eq45_path = root / str(eq45.get("artifact", ""))
    _require(delivery_path.is_file(), "package delivery contract is missing")
    _require(state_path.is_file(), "delivery state contract is missing")
    _require(eq45_path.is_file(), "Eq45 candidate artifact is missing")

    delivery = _read_json(delivery_path)
    state_contract = _read_json(state_path)
    eq45_artifact = _read_json(eq45_path)

    primary = delivery.get("primary_deliverable")
    _require(isinstance(primary, Mapping), "package delivery contract lacks primary_deliverable")
    _require(
        primary.get("api") == package.get("api"),
        "package unified velocity API drifted from Eq45 scope contract",
    )
    _require(
        primary.get("bundled_candidate_family") == package.get("bundled_candidate_family"),
        "package bundled candidate family drifted from Eq45 scope contract",
    )
    _require(
        package.get("api") == "openai_ns_reconstruction.velocity_components:velocity",
        "unexpected package unified velocity API",
    )
    _require(
        package.get("bundled_candidate_family") == "coupled_velocity_v1",
        "unexpected package unified candidate family",
    )

    _require(eq45_artifact.get("schema") == "eq45_velocity_candidate_v1", "unexpected Eq45 schema")
    _require(eq45.get("scope") == "candidate_local", "Eq45 export scope must remain candidate_local")
    _require(
        eq45.get("evaluator")
        == "openai_ns_reconstruction.constrained_eq45_candidate:Eq45VelocityCandidate.at_points",
        "unexpected Eq45 candidate evaluator",
    )

    truth = eq45_artifact.get("truth_boundary")
    _require(isinstance(truth, Mapping), "Eq45 candidate lacks truth_boundary")
    _require(truth.get("callable_serializable") is True, "Eq45 artifact is not callable/serializable")
    _require(truth.get("velocity_export_ready") is True, "Eq45 artifact is not candidate-locally export ready")
    _require(eq45.get("callable_serializable") is True, "scope contract lost callable/serializable state")
    _require(eq45.get("velocity_export_ready") is True, "scope contract lost candidate-local export readiness")

    _require(
        eq45.get("is_package_unified_default") is False,
        "Eq45 cannot be package unified/default while the package binding points to coupled_velocity_v1",
    )
    _require(
        semantics.get("eq45_local_readiness_is_evidence_for_package_unified_binding") is False,
        "candidate-local Eq45 readiness must not promote the package unified binding",
    )
    _require(
        semantics.get("pde_failure_blocks_candidate_local_export") is False,
        "PDE failure must not block candidate-local velocity export",
    )

    states = state_contract.get("states")
    _require(isinstance(states, Mapping), "delivery state contract lacks states")
    export_state = states.get("velocity_export_ready")
    _require(isinstance(export_state, Mapping), "delivery state contract lacks velocity_export_ready")
    positive_evidence = export_state.get("positive_evidence")
    _require(isinstance(positive_evidence, list), "velocity_export_ready lacks positive_evidence")
    _require(
        "unified_velocity_api" in positive_evidence,
        "canonical velocity_export_ready no longer requires a unified velocity API",
    )
    _require(
        semantics.get("canonical_velocity_export_ready_requires") == "unified_velocity_api",
        "Eq45 scope contract drifted from canonical velocity_export_ready semantics",
    )

    required_false = contract.get("required_false_eq45_claims")
    _require(isinstance(required_false, list), "required_false_eq45_claims must be a list")
    _require(
        set(required_false) == _REQUIRED_FALSE_EQ45_CLAIMS,
        "required false Eq45 claim set drifted",
    )
    for key in _REQUIRED_FALSE_EQ45_CLAIMS:
        _require(truth.get(key) is False, f"unsupported Eq45 claim promotion: {key}")

    package_status = delivery.get("claim_status")
    _require(isinstance(package_status, Mapping), "package delivery contract lacks claim_status")

    return {
        "contract_pass": True,
        "eq45_scope": "candidate_local",
        "eq45_callable_serializable": True,
        "eq45_velocity_export_ready": True,
        "eq45_is_package_unified_default": False,
        "package_unified_api": package.get("api"),
        "package_candidate_family": package.get("bundled_candidate_family"),
        "package_velocity_export_ready": package_status.get("velocity_export_ready"),
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def audit_eq45_export_scope_contract_file(
    path: str | Path, *, repo_root: str | Path
) -> dict[str, Any]:
    return audit_eq45_export_scope_contract(_read_json(Path(path)), repo_root=repo_root)
