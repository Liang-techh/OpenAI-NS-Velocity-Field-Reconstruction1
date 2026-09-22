from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_external_migration_contract import (
    SOURCE_BLOB_SHA1,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
    canonical_sha256,
    migration_contract,
    validate_contract,
)


def _rehash(contract: dict) -> dict:
    out = copy.deepcopy(contract)
    out.pop("contract_sha256", None)
    out["contract_sha256"] = canonical_sha256(out)
    return out


def test_canonical_contract_is_checksum_bound_and_fail_closed() -> None:
    contract = migration_contract()
    validate_contract(contract)

    assert contract["source"] == {
        "repository": SOURCE_REPOSITORY,
        "commit": SOURCE_COMMIT,
        "date": "2026-09-09",
        "path": SOURCE_PATH,
        "blob_sha1": SOURCE_BLOB_SHA1,
    }
    license_screen = contract["license_screen"]
    assert license_screen["github_license_field"] is None
    assert license_screen["status"] == "no_explicit_repository_license_detected_by_github_metadata"
    assert license_screen["public_repository_does_not_imply_copy_permission"] is True

    migration = contract["migration"]
    assert migration["classification"] == "reimplement_math_only"
    assert migration["direct_code_or_text_migration_approved"] is False
    assert migration["allowed_scope"]
    assert migration["excluded_scope"]
    assert migration["implementation_difference"]

    truth = contract["truth_boundary"]
    assert truth["third_party_code_copied"] is False
    assert truth["third_party_text_copied"] is False
    assert truth["explicit_upstream_license_confirmed"] is False
    assert truth["source_exact_velocity_recovered"] is False
    assert truth["source_exact_pressure_recovered"] is False
    assert truth["source_exact_forcing_recovered"] is False
    assert truth["openai_hidden_data_recovered"] is False
    assert truth["openai_field_identified"] is False
    assert truth["pde_validated"] is False
    assert truth["visualization_ready_promoted"] is False


def test_rehash_consistent_license_promotion_is_rejected() -> None:
    contract = migration_contract()
    mutated = copy.deepcopy(contract)
    mutated["license_screen"]["github_license_field"] = "MIT"
    mutated["license_screen"]["status"] = "explicit_license_detected"
    mutated = _rehash(mutated)
    assert mutated["contract_sha256"] != contract["contract_sha256"]
    with pytest.raises(ValueError, match="license"):
        validate_contract(mutated)


def test_rehash_consistent_direct_migration_promotion_is_rejected() -> None:
    contract = migration_contract()
    mutated = copy.deepcopy(contract)
    mutated["migration"]["classification"] = "direct_migration"
    mutated["migration"]["direct_code_or_text_migration_approved"] = True
    mutated = _rehash(mutated)
    assert mutated["contract_sha256"] != contract["contract_sha256"]
    with pytest.raises(ValueError, match="migration"):
        validate_contract(mutated)


def test_rehash_consistent_truth_promotion_is_rejected() -> None:
    contract = migration_contract()
    for key in ("third_party_code_copied", "pde_validated", "openai_field_identified"):
        mutated = copy.deepcopy(contract)
        mutated["truth_boundary"][key] = True
        mutated = _rehash(mutated)
        with pytest.raises(ValueError, match="truth boundary promoted"):
            validate_contract(mutated)


def test_rehash_consistent_source_identity_drift_is_rejected() -> None:
    contract = migration_contract()
    mutated = copy.deepcopy(contract)
    mutated["source"]["commit"] = "0" * 40
    mutated["source"]["blob_sha1"] = "1" * 40
    mutated = _rehash(mutated)
    with pytest.raises(ValueError, match="source identity drift"):
        validate_contract(mutated)


def test_checksum_detects_unrehashed_semantic_mutation() -> None:
    contract = migration_contract()
    mutated = copy.deepcopy(contract)
    mutated["migration"]["implementation_difference"] += " drift"
    with pytest.raises(ValueError, match="checksum mismatch"):
        validate_contract(mutated)
