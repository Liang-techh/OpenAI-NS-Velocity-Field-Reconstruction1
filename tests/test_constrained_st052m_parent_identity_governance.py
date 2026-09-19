from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_st052m_parent_identity_governance import (
    CONTRACT_PATH,
    GovernanceError,
    audit,
)


ROOT = Path(__file__).resolve().parents[1]


def _contract() -> dict:
    return json.loads((ROOT / "configs" / CONTRACT_PATH.name).read_text(encoding="utf-8"))


def _write_mutation(tmp_path: Path, mutate) -> Path:
    obj = copy.deepcopy(_contract())
    mutate(obj)
    path = tmp_path / "mutated.json"
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def test_current_parent_identity_contract_passes() -> None:
    result = audit(ROOT)
    assert result["status"] == "PASS"
    assert result["transform_kernel_materialized"] is True
    assert result["metadata_only_parent_construction_observed"] is True
    assert result["full_parent_bound"] is False
    assert result["complete_child_velocity_export_ready"] is False
    assert result["canonical_velocity_export_ready"] is True
    assert result["pde_validated"] is False


def test_transform_spec_checksum_cannot_be_promoted_to_child_sha(tmp_path: Path) -> None:
    path = _write_mutation(
        tmp_path,
        lambda obj: obj["parent_binding_current_state"].__setitem__(
            "transform_spec_sha256_is_complete_child_candidate_sha256", True
        ),
    )
    with pytest.raises(GovernanceError, match="unsupported parent-binding claim promoted"):
        audit(ROOT, contract_path=path)


def test_parent_metadata_strings_cannot_be_promoted_to_crypto_proof(tmp_path: Path) -> None:
    path = _write_mutation(
        tmp_path,
        lambda obj: obj["parent_binding_current_state"].__setitem__(
            "correct_metadata_strings_are_cryptographic_parent_proof", True
        ),
    )
    with pytest.raises(GovernanceError, match="unsupported parent-binding claim promoted"):
        audit(ROOT, contract_path=path)


def test_exact_source_parity_cannot_cover_future_injected_parents(tmp_path: Path) -> None:
    path = _write_mutation(
        tmp_path,
        lambda obj: obj["parent_binding_current_state"].__setitem__(
            "exact_source_parity_of_migration_implies_all_future_injected_parents_are_equivalent",
            True,
        ),
    )
    with pytest.raises(GovernanceError, match="unsupported parent-binding claim promoted"):
        audit(ROOT, contract_path=path)


def test_transform_kernel_cannot_promote_complete_child_export(tmp_path: Path) -> None:
    def mutate(obj: dict) -> None:
        obj["migrated_transform_state"]["experimental_child_velocity_export_ready"] = True

    path = _write_mutation(tmp_path, mutate)
    with pytest.raises(GovernanceError, match="incomplete child promoted"):
        audit(ROOT, contract_path=path)


def test_parent_binding_prerequisite_cannot_be_dropped(tmp_path: Path) -> None:
    def mutate(obj: dict) -> None:
        obj["required_before_complete_child_promotion"][
            "bind_runtime_parent_to_immutable_identity_beyond_caller_metadata_strings"
        ] = False

    path = _write_mutation(tmp_path, mutate)
    with pytest.raises(GovernanceError, match="promotion prerequisite disabled"):
        audit(ROOT, contract_path=path)


def test_source_vocabulary_fails_closed(tmp_path: Path) -> None:
    def mutate(obj: dict) -> None:
        obj["source_classification"]["parent_identity_binding_policy"] = "repository_fact"

    path = _write_mutation(tmp_path, mutate)
    with pytest.raises(GovernanceError, match="invalid source classification"):
        audit(ROOT, contract_path=path)


def test_cr001_momentum_threshold_cannot_drift(tmp_path: Path) -> None:
    def mutate(obj: dict) -> None:
        obj["canonical_cr001"]["pde_residual_max"] = 0.01

    path = _write_mutation(tmp_path, mutate)
    with pytest.raises(GovernanceError, match="momentum max threshold drift"):
        audit(ROOT, contract_path=path)


def test_canonical_delivery_identity_cannot_be_relabelled(tmp_path: Path) -> None:
    def mutate(obj: dict) -> None:
        obj["canonical_delivery"]["candidate_family"] = "st052m_local_swirl_energy_static_v1"

    path = _write_mutation(tmp_path, mutate)
    with pytest.raises(GovernanceError, match="canonical candidate family drift"):
        audit(ROOT, contract_path=path)


def test_pde_truth_cannot_be_promoted_by_identity_work(tmp_path: Path) -> None:
    def mutate(obj: dict) -> None:
        obj["truth_boundary"]["pde_validated"] = True

    path = _write_mutation(tmp_path, mutate)
    with pytest.raises(GovernanceError, match="unsupported truth promotion"):
        audit(ROOT, contract_path=path)
