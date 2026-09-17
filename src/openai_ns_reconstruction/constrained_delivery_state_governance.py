"""Fail-closed governance for independent constrained delivery claims.

This module governs claim semantics only. It does not set the current project
state and does not replace numerical constraints or source-provenance contracts.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "constrained_delivery_state_contract_v1"
ORIGIN_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
REQUIRED_STATES = {
    "velocity_export_ready",
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
}
FORBIDDEN_STATE_INFERENCES = {
    "velocity_export_ready": {
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
    },
    "visualization_ready": {
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
    },
    "visual_correspondence_verified": {"pde_validated", "paper_exact"},
    "pde_validated": {"visual_correspondence_verified", "paper_exact"},
}
PDE_EVIDENCE = {
    "independent_validator",
    "preregistered_thresholds",
    "held_out_validation",
}
PAPER_EXACT_EVIDENCE = {"independent_public_identity_evidence"}


def default_contract_path() -> Path:
    return Path(__file__).resolve().parents[2] / "configs" / "delivery_state_contract.json"


def load_delivery_state_contract(path: str | Path | None = None) -> dict[str, Any]:
    target = Path(path) if path is not None else default_contract_path()
    with target.open("r", encoding="utf-8") as handle:
        contract = json.load(handle)
    return validate_delivery_state_contract(contract)


def _string_set(value: Any, field: str) -> set[str]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError(f"{field} must be a list of strings")
    return set(value)


def validate_delivery_state_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(contract, Mapping):
        raise ValueError("delivery-state contract must be an object")
    if contract.get("schema") != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA!r}")

    vocabulary = contract.get("classification_vocabulary")
    if not isinstance(vocabulary, Mapping) or set(vocabulary) != ORIGIN_CLASSES:
        raise ValueError(
            "classification vocabulary must contain exactly the four canonical origin classes"
        )

    states = contract.get("states")
    if not isinstance(states, Mapping) or set(states) != REQUIRED_STATES:
        raise ValueError("states must contain exactly the canonical independent delivery claims")

    for state_name, state in states.items():
        if not isinstance(state, Mapping):
            raise ValueError(f"{state_name} must be an object")
        if state.get("definition_origin") not in ORIGIN_CLASSES:
            raise ValueError(f"{state_name}.definition_origin is not a canonical origin class")
        evidence_origins = _string_set(
            state.get("allowed_evidence_origins"),
            f"{state_name}.allowed_evidence_origins",
        )
        if not evidence_origins or not evidence_origins <= ORIGIN_CLASSES:
            raise ValueError(f"{state_name}.allowed_evidence_origins contains an invalid class")
        positive_evidence = _string_set(
            state.get("positive_evidence"),
            f"{state_name}.positive_evidence",
        )
        if not positive_evidence:
            raise ValueError(f"{state_name} must name positive evidence")
        does_not_imply = _string_set(
            state.get("does_not_imply"),
            f"{state_name}.does_not_imply",
        )
        required_negative_edges = FORBIDDEN_STATE_INFERENCES.get(state_name, set())
        if not required_negative_edges <= does_not_imply:
            missing = sorted(required_negative_edges - does_not_imply)
            raise ValueError(f"{state_name} is missing forbidden implication(s): {missing}")
        if "implies" in state:
            implied = _string_set(state["implies"], f"{state_name}.implies")
            forbidden = implied & required_negative_edges
            if forbidden:
                raise ValueError(f"{state_name} illegally implies {sorted(forbidden)}")

    if not PDE_EVIDENCE <= set(states["pde_validated"]["positive_evidence"]):
        raise ValueError("pde_validated must retain independent held-out preregistered evidence")
    if not PAPER_EXACT_EVIDENCE <= set(states["paper_exact"]["positive_evidence"]):
        raise ValueError("paper_exact requires independent public identity evidence")
    if set(states["paper_exact"]["allowed_evidence_origins"]) != {"public_source_fact"}:
        raise ValueError("paper_exact positive evidence must be public-source evidence")
    if states["paper_exact"].get("unverified_evidence_class") != "pending_unknown":
        raise ValueError("unverified paper identity must remain pending_unknown")

    rules = contract.get("forbidden_inferences")
    if not isinstance(rules, list) or not rules:
        raise ValueError("forbidden_inferences must be a nonempty list")
    for index, rule in enumerate(rules):
        if not isinstance(rule, Mapping) or rule.get("origin") not in ORIGIN_CLASSES:
            raise ValueError(f"forbidden_inferences[{index}] has invalid origin")
        _string_set(rule.get("from"), f"forbidden_inferences[{index}].from")
        if not isinstance(rule.get("to"), str):
            raise ValueError(f"forbidden_inferences[{index}].to must be a string")

    examples = contract.get("orthogonality_examples")
    if not isinstance(examples, list) or not examples:
        raise ValueError("orthogonality_examples must be present")
    if not any(
        example.get("allowed") is True
        and example.get("velocity_export_ready") is True
        and example.get("pde_validated") is False
        for example in examples
        if isinstance(example, Mapping)
    ):
        raise ValueError(
            "contract must explicitly allow export-ready while PDE validation is false"
        )

    return dict(contract)
