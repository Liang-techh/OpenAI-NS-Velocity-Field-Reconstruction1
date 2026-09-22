"""Canonical Agent-9 external-migration contract for the corrected Kokuno reader.

This module records what may be reused from the public Kokuno workbench and what
may not be promoted.  It is deliberately a governance/provenance object: no
third-party source code or prose is copied, and no candidate velocity, pressure,
forcing, residual, or threshold is changed.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

SCHEMA = "kokuno-external-migration-contract/v1"
TASK_ID = "CR-A9-110"

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_BLOB_SHA1 = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_DATE = "2026-09-09"
LICENSE_SCREEN_DATE = "2026-09-22"


def canonical_sha256(value: Any) -> str:
    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def migration_contract() -> dict[str, Any]:
    contract: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "date": SOURCE_DATE,
            "path": SOURCE_PATH,
            "blob_sha1": SOURCE_BLOB_SHA1,
        },
        "license_screen": {
            "screened_at_utc_date": LICENSE_SCREEN_DATE,
            "evidence_kind": "github_rest_repository_metadata",
            "github_license_field": None,
            "status": "no_explicit_repository_license_detected_by_github_metadata",
            "public_repository_does_not_imply_copy_permission": True,
        },
        "migration": {
            "classification": "reimplement_math_only",
            "direct_code_or_text_migration_approved": False,
            "allowed_scope": [
                "independent reimplementation of publicly stated mathematical structure",
                "normalized Haar probability measure convention on the auxiliary two-torus",
                "RF30-style quadratic covariance and slow-derivative identities",
                "qualitative coordinate/support/correction relationships needed to construct a candidate",
            ],
            "excluded_scope": [
                "copying upstream implementation code or substantial prose",
                "claiming recovery of hidden parameters or source-exact numerical fields",
                "claiming OpenAI velocity-field identity from the public reader",
                "transferring source pressure, forcing, or PDE validity without same-candidate verification",
            ],
            "implementation_difference": (
                "Repository modules independently implement typed NumPy/Python operators and candidate-specific providers; "
                "they are not source-code ports. The current RF30 auxiliary-T2 provider remains a separately governed "
                "same-identity dependency and is not supplied by this contract."
            ),
        },
        "truth_boundary": {
            "third_party_code_copied": False,
            "third_party_text_copied": False,
            "explicit_upstream_license_confirmed": False,
            "source_exact_velocity_recovered": False,
            "source_exact_pressure_recovered": False,
            "source_exact_forcing_recovered": False,
            "openai_hidden_data_recovered": False,
            "openai_field_identified": False,
            "paper_exact": False,
            "pde_validated": False,
            "visualization_ready_promoted": False,
        },
        "direct_contribution": (
            "gives RF30/Haar/correction adapters one reusable, checksum-bindable external-source provenance payload so "
            "independent mathematical reimplementation can advance toward a callable correction velocity without "
            "silently becoming unlicensed direct migration or source-exact evidence"
        ),
    }
    contract["contract_sha256"] = canonical_sha256(contract)
    validate_contract(contract)
    return contract


def validate_contract(contract: dict[str, Any]) -> None:
    if contract.get("schema") != SCHEMA or contract.get("task_id") != TASK_ID:
        raise ValueError("unexpected Kokuno migration contract schema/task")
    source = contract.get("source")
    license_screen = contract.get("license_screen")
    migration = contract.get("migration")
    truth = contract.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (source, license_screen, migration, truth)):
        raise ValueError("malformed Kokuno migration contract")
    expected_source = {
        "repository": SOURCE_REPOSITORY,
        "commit": SOURCE_COMMIT,
        "date": SOURCE_DATE,
        "path": SOURCE_PATH,
        "blob_sha1": SOURCE_BLOB_SHA1,
    }
    if source != expected_source:
        raise ValueError("Kokuno source identity drift")
    if license_screen.get("screened_at_utc_date") != LICENSE_SCREEN_DATE:
        raise ValueError("Kokuno license screen date drift")
    if license_screen.get("evidence_kind") != "github_rest_repository_metadata":
        raise ValueError("Kokuno license evidence-kind drift")
    if license_screen.get("github_license_field") is not None:
        raise ValueError("this frozen screen must record GitHub license metadata as null")
    if license_screen.get("status") != "no_explicit_repository_license_detected_by_github_metadata":
        raise ValueError("Kokuno license-status drift")
    if license_screen.get("public_repository_does_not_imply_copy_permission") is not True:
        raise ValueError("public-repository license boundary was relaxed")
    if migration.get("classification") != "reimplement_math_only":
        raise ValueError("Kokuno migration classification drift")
    if migration.get("direct_code_or_text_migration_approved") is not False:
        raise ValueError("direct upstream code/text migration must remain unapproved")
    allowed = migration.get("allowed_scope")
    excluded = migration.get("excluded_scope")
    if not isinstance(allowed, list) or not allowed or not isinstance(excluded, list) or not excluded:
        raise ValueError("Kokuno migration scope must be explicit")
    if not isinstance(migration.get("implementation_difference"), str) or not migration["implementation_difference"]:
        raise ValueError("Kokuno implementation difference must be explicit")
    for key in (
        "third_party_code_copied",
        "third_party_text_copied",
        "explicit_upstream_license_confirmed",
        "source_exact_velocity_recovered",
        "source_exact_pressure_recovered",
        "source_exact_forcing_recovered",
        "openai_hidden_data_recovered",
        "openai_field_identified",
        "paper_exact",
        "pde_validated",
        "visualization_ready_promoted",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"truth boundary promoted: {key}")
    recorded = contract.get("contract_sha256")
    if not isinstance(recorded, str) or len(recorded) != 64:
        raise ValueError("missing Kokuno migration contract checksum")
    unsigned = dict(contract)
    unsigned.pop("contract_sha256")
    if recorded != canonical_sha256(unsigned):
        raise ValueError("Kokuno migration contract checksum mismatch")
