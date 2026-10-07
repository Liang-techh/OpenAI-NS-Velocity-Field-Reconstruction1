from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_external_migration_contract import canonical_sha256
from openai_ns_reconstruction.kokuno_external_migration_inheritance import (
    inherit_kokuno_migration_contract,
    validate_inheritance,
)


def _payload() -> dict[str, object]:
    return {
        "candidate_semantic_identity_sha256": "2" * 64,
        "runtime_semantic_identity_sha256": "7" * 64,
        "stage": "current_q_s_to_correction_lineage",
        "repository_candidate_scientific_evidence": True,
    }


def _rehash(receipt: dict[str, object]) -> None:
    unsigned = dict(receipt)
    unsigned.pop("inheritance_sha256", None)
    receipt["inheritance_sha256"] = canonical_sha256(unsigned)


def test_canonical_inheritance_is_valid_and_reproducible() -> None:
    left = inherit_kokuno_migration_contract(_payload(), evidence_role="correction_lineage")
    right = inherit_kokuno_migration_contract(_payload(), evidence_role="correction_lineage")
    validate_inheritance(left)
    assert left == right
    assert left["external_migration"]["migration_classification"] == "reimplement_math_only"
    assert left["external_migration"]["direct_code_or_text_migration_approved"] is False
    assert left["truth_boundary"]["pde_validated_promoted_by_provenance"] is False


def test_downstream_payload_mutation_without_rehash_fails() -> None:
    receipt = inherit_kokuno_migration_contract(_payload(), evidence_role="correction_lineage")
    receipt["downstream_payload"]["stage"] = "mutated"
    with pytest.raises(ValueError, match="downstream payload checksum mismatch"):
        validate_inheritance(receipt)


def test_rehash_consistent_direct_migration_promotion_fails() -> None:
    receipt = inherit_kokuno_migration_contract(_payload(), evidence_role="correction_lineage")
    mutated = copy.deepcopy(receipt)
    mutated["external_migration"]["migration_classification"] = "direct_migration"
    mutated["external_migration"]["direct_code_or_text_migration_approved"] = True
    _rehash(mutated)
    with pytest.raises(ValueError, match="external-migration inheritance drift"):
        validate_inheritance(mutated)


def test_rehash_consistent_license_laundering_fails() -> None:
    receipt = inherit_kokuno_migration_contract(_payload(), evidence_role="runtime_authentication")
    mutated = copy.deepcopy(receipt)
    mutated["external_migration"]["github_license_field"] = "MIT"
    mutated["external_migration"]["license_status"] = "explicit_license_confirmed"
    _rehash(mutated)
    with pytest.raises(ValueError, match="external-migration inheritance drift"):
        validate_inheritance(mutated)


def test_rehash_consistent_truth_promotion_fails() -> None:
    receipt = inherit_kokuno_migration_contract(
        _payload(), evidence_role="repository_candidate_scientific_evidence"
    )
    mutated = copy.deepcopy(receipt)
    mutated["truth_boundary"]["openai_field_identified"] = True
    mutated["truth_boundary"]["pde_validated_promoted_by_provenance"] = True
    _rehash(mutated)
    with pytest.raises(ValueError, match="truth promotion"):
        validate_inheritance(mutated)


def test_reserved_downstream_keys_are_rejected() -> None:
    payload = _payload()
    payload["truth_boundary"] = {"pde_validated": True}
    with pytest.raises(ValueError, match="reserved inheritance keys"):
        inherit_kokuno_migration_contract(payload, evidence_role="velocity_lineage")


def test_unknown_role_is_rejected() -> None:
    with pytest.raises(ValueError, match="unsupported Kokuno downstream evidence role"):
        inherit_kokuno_migration_contract(_payload(), evidence_role="direct_migration")
