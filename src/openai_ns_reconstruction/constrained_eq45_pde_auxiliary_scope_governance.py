from __future__ import annotations

import json
from pathlib import Path
from typing import Any


SCHEMA = "eq45_pde_auxiliary_scope_contract_v1"
PENDING = "pending_unknown"


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def audit_documents(
    contract: dict[str, Any],
    constraints: dict[str, Any],
    candidate: dict[str, Any],
    delivery_contract: dict[str, Any],
) -> dict[str, Any]:
    if contract.get("schema") != SCHEMA:
        raise ValueError("unexpected auxiliary-scope contract schema")

    auxiliary = contract.get("current_auxiliary_state")
    if not isinstance(auxiliary, dict):
        raise ValueError("current_auxiliary_state must be an object")
    if auxiliary.get("pressure_binding") != PENDING:
        raise ValueError("current Eq45 pressure binding must remain pending_unknown")
    if auxiliary.get("forcing_binding") != PENDING:
        raise ValueError("current Eq45 forcing binding must remain pending_unknown")

    classification = contract.get("classification")
    if not isinstance(classification, dict):
        raise ValueError("classification must be an object")
    vocabulary = delivery_contract.get("classification_vocabulary")
    if not isinstance(vocabulary, dict):
        raise ValueError("delivery-state contract must define classification_vocabulary")
    expected_classification = {
        "velocity_delivery_requirement": "user_requirement",
        "cr001_forcing_family": "autonomous_design",
        "eq45_pressure_binding": "pending_unknown",
        "eq45_forcing_binding": "pending_unknown",
    }
    if classification != expected_classification:
        raise ValueError("auxiliary-scope source classification drifted")
    if not set(classification.values()).issubset(vocabulary):
        raise ValueError("auxiliary-scope classification uses noncanonical vocabulary")

    policy = contract.get("policy")
    if not isinstance(policy, dict):
        raise ValueError("policy must be an object")
    if policy.get("pending_auxiliary_blocks_velocity_export") is not False:
        raise ValueError("pending pressure/forcing must not block velocity export")
    if policy.get("pde_validation_requires_explicit_auxiliary_binding") is not True:
        raise ValueError("PDE validation must require explicit pressure/forcing binding")
    if policy.get("residual_defined_forcing_allowed") is not False:
        raise ValueError("residual-defined forcing must remain forbidden")
    if policy.get("pressure_or_forcing_solver_success_implies_pde_validation") is not False:
        raise ValueError("pressure/forcing solver success must not imply PDE validation")

    forcing = constraints.get("forcing")
    if not isinstance(forcing, dict):
        raise ValueError("constraints.forcing must be an object")
    if forcing.get("mode") != "restricted_two_parameter_family":
        raise ValueError("registered forcing family drifted")
    restriction = forcing.get("restriction")
    if not isinstance(restriction, str):
        raise ValueError("forcing restriction must be text")
    normalized = restriction.lower().replace("-", " ")
    if "no residual" not in normalized or "only a,c may be fitted" not in normalized:
        raise ValueError("forcing restriction no longer forbids residual-dependent freedom")

    if "pressure" in candidate or "forcing" in candidate:
        raise ValueError("frozen Eq45 velocity artifact unexpectedly embeds PDE auxiliaries")
    truth = candidate.get("truth_boundary")
    if not isinstance(truth, dict):
        raise ValueError("candidate truth_boundary must be an object")
    if truth.get("velocity_export_ready") is not True:
        raise ValueError("current candidate-local velocity export readiness drifted")
    if truth.get("pde_validated") is not False:
        raise ValueError("PDE validation cannot be true while Eq45 auxiliaries are pending")

    boundary = contract.get("truth_boundary")
    if not isinstance(boundary, dict) or any(boundary.values()):
        raise ValueError("auxiliary-scope contract must not promote scientific readiness")

    states = delivery_contract.get("states")
    if not isinstance(states, dict):
        raise ValueError("delivery-state contract must define states")
    for name in ("velocity_export_ready", "pde_validated"):
        if name not in states:
            raise ValueError(f"delivery-state contract missing {name}")

    forbidden = delivery_contract.get("forbidden_inferences")
    if not isinstance(forbidden, list):
        raise ValueError("delivery-state contract must define forbidden_inferences")
    pde_failure_blocks_export = any(
        isinstance(item, dict)
        and item.get("to") == "velocity_export_not_allowed"
        and "pde_validation_failed" in item.get("from", [])
        for item in forbidden
    )
    if not pde_failure_blocks_export:
        raise ValueError("delivery contract no longer forbids PDE failure from blocking export")

    return {
        "schema": SCHEMA,
        "velocity_export_ready": True,
        "pressure_binding": PENDING,
        "forcing_binding": PENDING,
        "registered_forcing_mode": forcing["mode"],
        "residual_defined_forcing_allowed": False,
        "pde_validated": False,
        "claim_scope": "delivery_vs_pde_auxiliary_governance_only",
        "velocity_changed": False,
        "visualization_ready": bool(truth.get("visualization_ready", False)),
        "visual_correspondence_verified": bool(
            truth.get("visual_correspondence_verified", False)
        ),
        "paper_exact": bool(truth.get("paper_exact", False)),
        "openai_field_identified": bool(truth.get("openai_field_identified", False)),
    }


def audit_eq45_pde_auxiliary_scope(repo_root: str | Path | None = None) -> dict[str, Any]:
    root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[2]
    return audit_documents(
        _load_json(root / "configs" / "eq45_pde_auxiliary_scope_contract.json"),
        _load_json(root / "configs" / "constraints.json"),
        _load_json(root / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"),
        _load_json(root / "configs" / "delivery_state_contract.json"),
    )


if __name__ == "__main__":
    print(json.dumps(audit_eq45_pde_auxiliary_scope(), indent=2, sort_keys=True))
