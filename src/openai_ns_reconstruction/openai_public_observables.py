"""Fail-closed contract for public OpenAI velocity-field observables.

The contract records only qualitative statements visible in OpenAI's public
Navier--Stokes announcement. It deliberately contains no pixel target, fitted
coefficient, hidden source datum, or PDE acceptance rule.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

_SCHEMA_VERSION = 1
_TASK_ID = "CR-A9-078"
_SOURCE_URL = "https://openai.com/index/navier-stokes-solution/"
_SOURCE_DATE = "2026-09-08"
_EXTERNAL_METHOD_REPO = "scikit-image/scikit-image"
_EXTERNAL_METHOD_COMMIT = "cc0a4b16ebf655873cd6f770e897460abdade288"
_EXTERNAL_METHOD_LICENSE = "BSD-3-Clause (project default; individual vendored files may differ)"
_EXPECTED_IDS = (
    "vortex_swirl_presence",
    "inward_spiraling_trajectories",
    "axial_stretching_trajectories",
    "shrinking_accelerating_core",
    "angular_rotation_variation",
    "radius_dependent_circulation",
)
_FALSE_TRUTH_KEYS = (
    "visualization_ready",
    "visual_correspondence_verified",
    "source_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)
_DEFAULT_PATH = (
    Path(__file__).resolve().parents[2] / "configs" / "openai_public_observables_v1.json"
)


def _canonical_bytes(contract: Mapping[str, Any]) -> bytes:
    return json.dumps(
        contract,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


def contract_sha256(contract: Mapping[str, Any]) -> str:
    """Return a deterministic digest of the validated contract payload."""
    validate_openai_public_observable_contract(contract)
    return hashlib.sha256(_canonical_bytes(contract)).hexdigest()


def validate_openai_public_observable_contract(contract: Mapping[str, Any]) -> None:
    """Reject source drift, invented numeric targets, or truth-boundary promotion."""
    if not isinstance(contract, Mapping):
        raise ValueError("observable contract must be a mapping")
    if contract.get("schema_version") != _SCHEMA_VERSION:
        raise ValueError("unexpected public-observable schema version")
    if contract.get("task_id") != _TASK_ID:
        raise ValueError("unexpected public-observable task id")
    if contract.get("status") != "preregistered_public_observable_contract":
        raise ValueError("public-observable contract status drift")

    source = contract.get("source")
    if not isinstance(source, Mapping):
        raise ValueError("missing public source record")
    if source.get("kind") != "official_public_webpage":
        raise ValueError("public source must remain the official webpage")
    if source.get("publisher") != "OpenAI":
        raise ValueError("public source publisher drift")
    if source.get("url") != _SOURCE_URL or source.get("published_date") != _SOURCE_DATE:
        raise ValueError("public source identity drift")

    observables = contract.get("observables")
    if not isinstance(observables, list):
        raise ValueError("observables must be a list")
    ids = tuple(item.get("id") for item in observables if isinstance(item, Mapping))
    if ids != _EXPECTED_IDS or len(observables) != len(_EXPECTED_IDS):
        raise ValueError("public observable set/order drift")
    for item in observables:
        if not isinstance(item, Mapping):
            raise ValueError("malformed public observable")
        if item.get("numerical_target") is not None:
            raise ValueError("public observables may not invent numerical targets")
        if not isinstance(item.get("public_observation"), str) or not item["public_observation"].strip():
            raise ValueError("public observable must carry a nonempty paraphrase")
        if item.get("kind") not in {
            "qualitative_geometry",
            "qualitative_trajectory",
            "qualitative_time_evolution",
            "qualitative_velocity_structure",
        }:
            raise ValueError("unexpected observable kind")

    screening = contract.get("external_method_screening")
    if not isinstance(screening, list) or len(screening) != 1:
        raise ValueError("external screening record drift")
    method = screening[0]
    if not isinstance(method, Mapping):
        raise ValueError("malformed external screening record")
    if method.get("source_repo") != _EXTERNAL_METHOD_REPO:
        raise ValueError("external screening source drift")
    if method.get("source_commit") != _EXTERNAL_METHOD_COMMIT:
        raise ValueError("external screening commit drift")
    if method.get("license") != _EXTERNAL_METHOD_LICENSE:
        raise ValueError("external screening license drift")
    if not isinstance(method.get("candidate_scope"), str) or not method["candidate_scope"].strip():
        raise ValueError("external screening candidate scope missing")
    if method.get("classification") != "screened_not_adopted" or method.get("migration_scope") != "none":
        raise ValueError("pixel-space tooling may not be silently promoted or migrated")
    if not isinstance(method.get("reason"), str) or not method["reason"].strip():
        raise ValueError("external screening rationale missing")

    truth = contract.get("truth_boundary")
    if not isinstance(truth, Mapping) or set(truth) != set(_FALSE_TRUTH_KEYS):
        raise ValueError("truth-boundary key drift")
    if any(truth[key] is not False for key in _FALSE_TRUTH_KEYS):
        raise ValueError("public observables cannot promote scientific or identity truth flags")

    forbidden = contract.get("forbidden_inferences")
    if not isinstance(forbidden, list) or len(forbidden) < 4:
        raise ValueError("forbidden-inference boundary missing")
    allowed = contract.get("allowed_downstream_use")
    if not isinstance(allowed, list) or len(allowed) < 3:
        raise ValueError("allowed-use boundary missing")


def load_openai_public_observable_contract(path: str | Path | None = None) -> dict[str, Any]:
    """Load and validate the repository's immutable qualitative source contract."""
    contract_path = Path(path) if path is not None else _DEFAULT_PATH
    with contract_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    validate_openai_public_observable_contract(payload)
    return payload


def observable_contract_report(path: str | Path | None = None) -> dict[str, Any]:
    """Return a compact machine-readable receipt without upgrading any truth claim."""
    payload = load_openai_public_observable_contract(path)
    method = payload["external_method_screening"][0]
    return {
        "task_id": payload["task_id"],
        "source_url": payload["source"]["url"],
        "source_published_date": payload["source"]["published_date"],
        "observable_ids": [item["id"] for item in payload["observables"]],
        "external_method_source_repo": method["source_repo"],
        "external_method_source_commit": method["source_commit"],
        "external_method_license": method["license"],
        "external_method_candidate_scope": method["candidate_scope"],
        "external_method_classification": method["classification"],
        "external_method_migration_scope": method["migration_scope"],
        "external_method_reason": method["reason"],
        "contract_sha256": contract_sha256(payload),
        "numeric_visual_target_defined": False,
        **payload["truth_boundary"],
    }


def main() -> None:
    print(json.dumps(observable_contract_report(), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
