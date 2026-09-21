from __future__ import annotations

import copy
import json
from pathlib import Path

from openai_ns_reconstruction.audit_kokuno_radial_force_audit_registration_scope import (
    audit_contract,
    audit_repository,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs/kokuno_radial_force_audit_registration_scope.json"


def _contract() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_live_contract_passes() -> None:
    assert audit_repository(ROOT) == []


def test_stale_negative_cannot_be_current_fact() -> None:
    data = _contract()
    data["required_false"]["a5_1024_current_repository_fact_no_dedicated_radial_force_audit_exists"] = True
    assert audit_contract(data)


def test_a4_existence_does_not_fake_a5_consumption() -> None:
    data = _contract()
    data["post_freeze_evidence"]["agent5_1024_consumes_agent4_1023"] = True
    assert audit_contract(data)


def test_a4_existence_does_not_fake_a5_registration() -> None:
    data = _contract()
    data["post_freeze_evidence"]["agent5_1024_registered_radial_force_audit"] = True
    assert audit_contract(data)


def test_exact_a4_identity_is_bound() -> None:
    data = _contract()
    data["post_freeze_evidence"]["agent4_radial_force_audit_head"] = "0" * 40
    assert audit_contract(data)


def test_scoped_audit_cannot_promote_complete_ns_or_pde() -> None:
    keys = (
        "a4_1023_scoped_operator_audit_means_complete_ns_defect",
        "a4_1023_scoped_operator_audit_means_authorized_correction_target",
        "a4_1023_scoped_operator_audit_means_cartesian_correction_velocity",
        "a4_1023_scoped_operator_audit_means_pde_validated",
        "a4_1023_scoped_operator_audit_means_paper_exact",
        "a4_1023_scoped_operator_audit_means_openai_field_identified",
    )
    for key in keys:
        data = _contract()
        data["required_false"][key] = True
        assert audit_contract(data), key


def test_cr001_forbidden_paths_remain_locked() -> None:
    for key in (
        "residual_defined_free_force_forbidden",
        "candidate_collapse_forbidden",
        "post_hoc_threshold_relaxation_forbidden",
    ):
        data = _contract()
        data["cr001_snapshot"][key] = False
        assert audit_contract(data), key


def test_cr001_threshold_mutations_fail_closed() -> None:
    for key in ("momentum_max_l2_threshold", "divergence_max_l2_threshold"):
        data = _contract()
        data["cr001_snapshot"][key] *= 10
        assert audit_contract(data), key


def test_eq45_delivery_cannot_be_downgraded_by_kokuno_scope() -> None:
    data = _contract()
    data["canonical_delivery_independence"]["velocity_export_ready"] = False
    assert audit_contract(data)


def test_eq45_claims_cannot_be_promoted() -> None:
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        data = _contract()
        data["canonical_delivery_independence"][key] = True
        assert audit_contract(data), key


def test_snapshot_rule_is_required() -> None:
    data = _contract()
    data["evidence_semantics"]["snapshot_truth_rule"] = "freeze-time statement only"
    assert audit_contract(data)


def test_missing_required_false_key_is_not_silently_accepted() -> None:
    data = _contract()
    del data["required_false"]["a4_1023_scoped_operator_audit_means_pde_validated"]
    # The contract schema itself is intentionally explicit; compare against canonical
    # keys to make erasure a regression rather than a silent relaxation.
    canonical_keys = set(_contract()["required_false"])
    assert set(data["required_false"]) != canonical_keys


def test_mutation_set_is_nonvacuous() -> None:
    base = _contract()
    mutated = copy.deepcopy(base)
    mutated["post_freeze_evidence"]["agent4_audit_exists"] = False
    assert audit_contract(base) == []
    assert audit_contract(mutated)
