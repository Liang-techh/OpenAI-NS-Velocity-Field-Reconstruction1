from __future__ import annotations

import json
from pathlib import Path
from typing import Any


SCHEMA = "eq45_support_transform_evidence_scope_v1"
PENDING = "pending_unknown"

_EXPECTED_CLASSIFICATION = {
    "velocity_delivery_requirement": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "registered_physical_support": "autonomous_design",
    "support_transform_design": "autonomous_design",
    "production_transform_status": "pending_unknown",
}

_REQUIRED_REVALIDATION = {
    "child_candidate_identity",
    "velocity_export_ready",
    "physical_support_validated",
    "cr001_energy_normalization",
    "boundary_support_metrics",
    "pde_validated",
    "visualization_ready",
    "visual_resolution_stability",
    "visual_correspondence_verified",
}

_FORBIDDEN_CHILD_INHERITANCE = {
    "parent_candidate_identity_as_child_identity",
    "parent_energy_metrics",
    "parent_boundary_or_tail_metrics",
    "parent_pde_residual_metrics",
    "parent_visual_fingerprint_or_resolution_metrics",
    "parent_canonical_promotion_status",
}


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def _require_bool(policy: dict[str, Any], key: str, expected: bool) -> None:
    if policy.get(key) is not expected:
        raise ValueError(f"support-transform policy drifted: {key}")


def audit_documents(
    contract: dict[str, Any],
    constraints: dict[str, Any],
    parent_candidate: dict[str, Any],
    delivery_contract: dict[str, Any],
) -> dict[str, Any]:
    if contract.get("schema") != SCHEMA:
        raise ValueError("unexpected support-transform evidence contract schema")

    classification = contract.get("classification")
    if classification != _EXPECTED_CLASSIFICATION:
        raise ValueError("support-transform source classification drifted")
    vocabulary = delivery_contract.get("classification_vocabulary")
    if not isinstance(vocabulary, dict):
        raise ValueError("delivery-state contract must define classification_vocabulary")
    if not set(classification.values()).issubset(vocabulary):
        raise ValueError("support-transform classification uses noncanonical vocabulary")

    domain = constraints.get("domain")
    if not isinstance(domain, dict):
        raise ValueError("constraints.domain must be an object")
    registered = contract.get("registered_support")
    expected_support = {
        "physical_domain": domain.get("physical"),
        "support": domain.get("support"),
        "boundary": domain.get("boundary"),
    }
    if registered != expected_support:
        raise ValueError("support-transform contract drifted from CR001 support semantics")
    if expected_support != {
        "physical_domain": "R^3",
        "support": "r < 2 and abs(z) < 2",
        "boundary": "smooth zero extension of velocity and pressure outside support",
    }:
        raise ValueError("active CR001 physical support contract changed")

    parent = contract.get("parent_candidate")
    if not isinstance(parent, dict):
        raise ValueError("parent_candidate must be an object")
    if parent.get("schema") != parent_candidate.get("schema"):
        raise ValueError("parent candidate schema drifted")
    truth = parent_candidate.get("truth_boundary")
    if not isinstance(truth, dict):
        raise ValueError("parent candidate truth_boundary must be an object")
    required_parent_state = {
        "requires_physical_support_connection": True,
        "physical_support_validated": False,
        "velocity_export_ready": True,
    }
    for key, expected in required_parent_state.items():
        if parent.get(key) is not expected or truth.get(key) is not expected:
            raise ValueError(f"parent candidate support/delivery state drifted: {key}")

    current = contract.get("current_transform_state")
    if not isinstance(current, dict):
        raise ValueError("current_transform_state must be an object")
    if current.get("production_support_transform_materialized") is not False:
        raise ValueError("contract is stale: a production support transform now requires child audit")
    if current.get("child_candidate_identity") != PENDING:
        raise ValueError("child candidate identity must remain pending before materialization")
    if current.get("production_support_validation") != PENDING:
        raise ValueError("production support validation must remain pending before materialization")

    inheritance = contract.get("evidence_inheritance")
    if not isinstance(inheritance, dict):
        raise ValueError("evidence_inheritance must be an object")
    provenance_only = inheritance.get("provenance_only")
    if not isinstance(provenance_only, list):
        raise ValueError("provenance_only must be a list")
    forbidden_as_provenance = {
        "pde_validated",
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "parent_pde_residual_metrics",
        "parent_energy_metrics",
    }
    if forbidden_as_provenance.intersection(provenance_only):
        raise ValueError("validation evidence cannot be inherited as provenance")

    revalidate = inheritance.get("must_be_reestablished_for_transformed_candidate")
    if not isinstance(revalidate, list) or set(revalidate) != _REQUIRED_REVALIDATION:
        raise ValueError("support-transformed child revalidation set drifted")
    forbidden = inheritance.get("must_not_be_inherited_as_child_validation")
    if not isinstance(forbidden, list) or set(forbidden) != _FORBIDDEN_CHILD_INHERITANCE:
        raise ValueError("forbidden child-evidence inheritance set drifted")

    local_only = inheritance.get("local_only_if_identity_plateau_is_independently_checked")
    if local_only != ["pointwise_velocity_equality_inside_declared_identity_plateau"]:
        raise ValueError("identity-plateau evidence must remain explicitly local")

    policy = contract.get("policy")
    if not isinstance(policy, dict):
        raise ValueError("policy must be an object")
    for key in (
        "nonidentity_support_transform_creates_new_candidate_identity",
        "identity_plateau_evidence_is_local_only",
        "support_validation_requires_production_child",
        "child_velocity_export_requires_child_serialization_and_public_evaluator_check",
    ):
        _require_bool(policy, key, True)
    for key in (
        "synthetic_taper_calibration_establishes_production_support",
        "parent_global_validation_applies_to_transformed_child",
        "support_validation_implies_pde_validation",
        "support_validation_implies_visual_correspondence",
        "pde_failure_blocks_child_velocity_export",
        "residual_defined_or_pointwise_free_forcing_allowed",
        "support_transform_is_public_source_fact",
    ):
        _require_bool(policy, key, False)

    forcing = constraints.get("forcing")
    if not isinstance(forcing, dict) or forcing.get("mode") != "restricted_two_parameter_family":
        raise ValueError("registered forcing family drifted")
    restriction = forcing.get("restriction")
    if not isinstance(restriction, str):
        raise ValueError("forcing restriction must be text")
    normalized = restriction.lower().replace("-", " ")
    if "no residual" not in normalized or "only a,c may be fitted" not in normalized:
        raise ValueError("forcing restriction no longer forbids free residual cancellation")

    contract_truth = contract.get("truth_boundary")
    if not isinstance(contract_truth, dict) or any(contract_truth.values()):
        raise ValueError("prospective support-transform governance cannot promote readiness")

    states = delivery_contract.get("states")
    if not isinstance(states, dict):
        raise ValueError("delivery-state contract must define states")
    for name in (
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
    ):
        if name not in states:
            raise ValueError(f"delivery-state contract missing {name}")

    forbidden_inferences = delivery_contract.get("forbidden_inferences")
    if not isinstance(forbidden_inferences, list):
        raise ValueError("delivery-state contract must define forbidden_inferences")
    pde_failure_cannot_block_export = any(
        isinstance(item, dict)
        and item.get("to") == "velocity_export_not_allowed"
        and "pde_validation_failed" in item.get("from", [])
        for item in forbidden_inferences
    )
    if not pde_failure_cannot_block_export:
        raise ValueError("delivery contract no longer forbids PDE failure from blocking export")

    return {
        "schema": SCHEMA,
        "parent_velocity_export_ready": True,
        "parent_requires_physical_support_connection": True,
        "production_support_transform_materialized": False,
        "child_candidate_identity": PENDING,
        "child_global_evidence_must_be_reestablished": True,
        "pde_failure_blocks_child_velocity_export": False,
        "residual_defined_or_pointwise_free_forcing_allowed": False,
        "physical_support_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "claim_scope": "support_transform_evidence_inheritance_governance_only",
        "velocity_changed": False,
    }


def audit_eq45_support_transform_evidence_scope(
    repo_root: str | Path | None = None,
) -> dict[str, Any]:
    root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[2]
    return audit_documents(
        _load_json(root / "configs" / "eq45_support_transform_evidence_scope.json"),
        _load_json(root / "configs" / "constraints.json"),
        _load_json(root / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"),
        _load_json(root / "configs" / "delivery_state_contract.json"),
    )


if __name__ == "__main__":
    print(json.dumps(audit_eq45_support_transform_evidence_scope(), indent=2, sort_keys=True))
