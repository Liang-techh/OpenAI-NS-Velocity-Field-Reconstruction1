"""Fail-closed inheritance of the canonical Kokuno external-migration contract.

CR-A9-111 records the source/license/migration truth boundary.  Downstream
scientific receipts must bind that boundary into their own identity rather than
merely repeating a repository/commit/blob triple in prose.  This module provides
one small reusable adapter for that purpose.

It changes no candidate field, source mathematics, pressure, forcing, residual,
optimizer, validation sample, or project threshold.
"""
from __future__ import annotations

import copy
from typing import Any, Mapping

from .kokuno_external_migration_contract import (
    canonical_sha256,
    migration_contract,
    validate_contract,
)

SCHEMA = "kokuno-external-migration-inheritance/v1"
TASK_ID = "CR-A9-117"

_ALLOWED_EVIDENCE_ROLES = {
    "runtime_authentication",
    "source_structure_reference",
    "repository_candidate_scientific_evidence",
    "correction_lineage",
    "velocity_lineage",
}

_RESERVED_DOWNSTREAM_KEYS = {
    "external_migration",
    "inheritance_sha256",
    "schema",
    "task_id",
    "truth_boundary",
}


def _json_clone(value: Mapping[str, Any]) -> dict[str, Any]:
    """Return an independent JSON-shaped mapping while rejecting NaN/non-JSON data."""

    # canonical_sha256 already performs strict JSON serialization with allow_nan=False.
    canonical_sha256(value)
    return copy.deepcopy(dict(value))


def inherit_kokuno_migration_contract(
    downstream_payload: Mapping[str, Any],
    *,
    evidence_role: str,
) -> dict[str, Any]:
    """Bind the admitted external-migration contract into a downstream receipt.

    The downstream payload remains application-specific.  This adapter freezes its
    byte-independent semantic hash together with the canonical source/license/
    migration boundary so later rehashing cannot silently promote direct migration,
    an invented license, source-exact identity, OpenAI-field identity, or PDE truth.
    """

    if evidence_role not in _ALLOWED_EVIDENCE_ROLES:
        raise ValueError("unsupported Kokuno downstream evidence role")
    if not isinstance(downstream_payload, Mapping) or not downstream_payload:
        raise ValueError("downstream payload must be a nonempty mapping")
    overlap = _RESERVED_DOWNSTREAM_KEYS.intersection(downstream_payload)
    if overlap:
        raise ValueError(f"downstream payload uses reserved inheritance keys: {sorted(overlap)}")

    contract = migration_contract()
    validate_contract(contract)
    payload = _json_clone(downstream_payload)

    external = {
        "contract_schema": contract["schema"],
        "contract_task_id": contract["task_id"],
        "contract_sha256": contract["contract_sha256"],
        "source_repository": contract["source"]["repository"],
        "source_commit": contract["source"]["commit"],
        "source_path": contract["source"]["path"],
        "source_blob_sha1": contract["source"]["blob_sha1"],
        "license_evidence_kind": contract["license_screen"]["evidence_kind"],
        "github_license_field": contract["license_screen"]["github_license_field"],
        "license_status": contract["license_screen"]["status"],
        "migration_classification": contract["migration"]["classification"],
        "direct_code_or_text_migration_approved": contract["migration"][
            "direct_code_or_text_migration_approved"
        ],
    }

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "evidence_role": evidence_role,
        "downstream_payload": payload,
        "downstream_payload_sha256": canonical_sha256(payload),
        "external_migration": external,
        "truth_boundary": {
            "direct_third_party_code_or_text_migration": False,
            "source_exact_velocity_claimed": False,
            "source_exact_pressure_claimed": False,
            "source_exact_forcing_claimed": False,
            "openai_hidden_data_claimed": False,
            "openai_field_identified": False,
            "paper_exact": False,
            "visualization_ready_promoted_by_provenance": False,
            "pde_validated_promoted_by_provenance": False,
        },
    }
    receipt["inheritance_sha256"] = canonical_sha256(receipt)
    validate_inheritance(receipt)
    return receipt


def validate_inheritance(receipt: Mapping[str, Any]) -> None:
    """Validate one inheritance receipt against the currently admitted contract."""

    if receipt.get("schema") != SCHEMA or receipt.get("task_id") != TASK_ID:
        raise ValueError("unexpected Kokuno inheritance schema/task")
    role = receipt.get("evidence_role")
    if role not in _ALLOWED_EVIDENCE_ROLES:
        raise ValueError("unsupported Kokuno downstream evidence role")

    payload = receipt.get("downstream_payload")
    if not isinstance(payload, Mapping) or not payload:
        raise ValueError("missing downstream payload")
    overlap = _RESERVED_DOWNSTREAM_KEYS.intersection(payload)
    if overlap:
        raise ValueError("downstream payload contains reserved inheritance keys")
    if receipt.get("downstream_payload_sha256") != canonical_sha256(payload):
        raise ValueError("downstream payload checksum mismatch")

    contract = migration_contract()
    validate_contract(contract)
    expected_external = {
        "contract_schema": contract["schema"],
        "contract_task_id": contract["task_id"],
        "contract_sha256": contract["contract_sha256"],
        "source_repository": contract["source"]["repository"],
        "source_commit": contract["source"]["commit"],
        "source_path": contract["source"]["path"],
        "source_blob_sha1": contract["source"]["blob_sha1"],
        "license_evidence_kind": contract["license_screen"]["evidence_kind"],
        "github_license_field": contract["license_screen"]["github_license_field"],
        "license_status": contract["license_screen"]["status"],
        "migration_classification": contract["migration"]["classification"],
        "direct_code_or_text_migration_approved": contract["migration"][
            "direct_code_or_text_migration_approved"
        ],
    }
    if receipt.get("external_migration") != expected_external:
        raise ValueError("Kokuno external-migration inheritance drift")

    truth = receipt.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("missing Kokuno inheritance truth boundary")
    expected_truth_keys = {
        "direct_third_party_code_or_text_migration",
        "source_exact_velocity_claimed",
        "source_exact_pressure_claimed",
        "source_exact_forcing_claimed",
        "openai_hidden_data_claimed",
        "openai_field_identified",
        "paper_exact",
        "visualization_ready_promoted_by_provenance",
        "pde_validated_promoted_by_provenance",
    }
    if set(truth) != expected_truth_keys:
        raise ValueError("unexpected Kokuno inheritance truth-boundary fields")
    for key in expected_truth_keys:
        if truth.get(key) is not False:
            raise ValueError(f"Kokuno inheritance truth promotion: {key}")

    recorded = receipt.get("inheritance_sha256")
    if not isinstance(recorded, str) or len(recorded) != 64:
        raise ValueError("missing Kokuno inheritance checksum")
    unsigned = dict(receipt)
    unsigned.pop("inheritance_sha256", None)
    if recorded != canonical_sha256(unsigned):
        raise ValueError("Kokuno inheritance checksum mismatch")
