"""Fail-closed audit for Eq. (4.1)/(4.5) source classification.

This module does not validate a velocity field. It only keeps source facts,
user requirements, autonomous design choices, and unresolved data in the same
canonical four-class vocabulary used by the delivery-state contract.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping


EXPECTED_CLASSES = {
    "public_formula": "public_source_fact",
    "callable_delivery": "user_requirement",
    "numerical_choices": "autonomous_design",
    "unidentified_profiles": "pending_unknown",
}

UNVERIFIED_CLAIMS = (
    "final_profiles_numerically_identified",
    "full_constructed_field_recovered",
    "openai_correspondence_verified",
    "pde_validated",
    "paper_exact",
)


def _load(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _canonical_classes(delivery_contract: Mapping[str, Any]) -> set[str]:
    vocabulary = delivery_contract.get("classification_vocabulary")
    if not isinstance(vocabulary, Mapping):
        raise ValueError("delivery contract is missing classification_vocabulary")
    classes = set(vocabulary)
    expected = {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    if classes != expected:
        raise ValueError("delivery contract must expose exactly the canonical four classes")
    return classes


def _find_entry(entries: list[Mapping[str, Any]], fragment: str) -> Mapping[str, Any]:
    matches = [entry for entry in entries if fragment in str(entry.get("item", ""))]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one source classification entry matching {fragment!r}")
    return matches[0]


def validate_eq45_source_classification(
    source_contract: Mapping[str, Any],
    delivery_contract: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate source/design/pending labels without promoting mathematical claims."""
    canonical = _canonical_classes(delivery_contract)

    entries_raw = source_contract.get("source_classification")
    if not isinstance(entries_raw, list) or not entries_raw:
        raise ValueError("eq45 source contract must contain source_classification entries")

    entries: list[Mapping[str, Any]] = []
    for index, entry in enumerate(entries_raw):
        if not isinstance(entry, Mapping):
            raise ValueError(f"source_classification[{index}] must be an object")
        label = entry.get("classification")
        if label not in canonical:
            raise ValueError(
                f"source_classification[{index}] uses noncanonical classification {label!r}"
            )
        if not isinstance(entry.get("item"), str) or not entry["item"].strip():
            raise ValueError(f"source_classification[{index}].item must be nonempty")
        if not isinstance(entry.get("source"), str) or not entry["source"].strip():
            raise ValueError(f"source_classification[{index}].source must be nonempty")
        entries.append(entry)

    public_formula = _find_entry(entries, "Eq. (4.1) similarity coordinates")
    callable_delivery = _find_entry(entries, "callable velocity(x,y,z,t)")
    numerical_choices = _find_entry(entries, "finite numerical choices")
    unidentified_profiles = _find_entry(entries, "final numerical profile data")

    expected_pairs = (
        (public_formula, EXPECTED_CLASSES["public_formula"]),
        (callable_delivery, EXPECTED_CLASSES["callable_delivery"]),
        (numerical_choices, EXPECTED_CLASSES["numerical_choices"]),
        (unidentified_profiles, EXPECTED_CLASSES["unidentified_profiles"]),
    )
    for entry, expected in expected_pairs:
        if entry["classification"] != expected:
            raise ValueError(f"{entry['item']!r} must remain classified as {expected!r}")

    coordinate = source_contract.get("coordinate_contract")
    if not isinstance(coordinate, Mapping):
        raise ValueError("eq45 source contract is missing coordinate_contract")
    if coordinate.get("implicit_q_relation") != "q-z^2*q^(2*h)=1-t":
        raise ValueError("Eq. (4.1) implicit q relation drifted from the sourced +2h form")
    if coordinate.get("q_exponent_sign") != "positive":
        raise ValueError("Eq. (4.1) q exponent sign must remain positive")

    velocity = source_contract.get("velocity_contract")
    if not isinstance(velocity, Mapping) or velocity.get("leading_field_only") is not True:
        raise ValueError("Eq. (4.5) contract must remain labeled leading_field_only")

    claim_status = source_contract.get("claim_status")
    if not isinstance(claim_status, Mapping):
        raise ValueError("eq45 source contract is missing claim_status")
    for claim in UNVERIFIED_CLAIMS:
        if claim_status.get(claim) is not False:
            raise ValueError(f"{claim} must remain false until independent evidence exists")

    if claim_status.get("coordinate_formula_publicly_sourced") is not True:
        raise ValueError("coordinate_formula_publicly_sourced must remain true")
    if claim_status.get("leading_velocity_formula_publicly_sourced") is not True:
        raise ValueError("leading_velocity_formula_publicly_sourced must remain true")

    return {
        "canonical_classification_pass": True,
        "q_source_relation_pass": True,
        "claim_boundary_pass": True,
        "pde_validated": False,
        "paper_exact": False,
    }


def audit_repository_contracts(
    source_path: str | Path,
    delivery_path: str | Path,
) -> dict[str, Any]:
    return validate_eq45_source_classification(_load(source_path), _load(delivery_path))
