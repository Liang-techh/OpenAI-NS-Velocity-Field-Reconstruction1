from __future__ import annotations

import copy
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_finite_correction_physical_contract_scope import (
    audit,
    load_contract,
    mechanics_protocol_identity_witness,
    validate_contract,
)

ROOT = Path(__file__).resolve().parents[1]


def test_live_contract_audit_passes_without_scientific_promotion() -> None:
    receipt = audit(ROOT)
    assert receipt["status"] == "PASS"
    assert receipt["same_residual_protocol_sha_enforced"] is True
    assert receipt["same_physical_contract_identity_enforced"] is False
    assert receipt["mechanics_same_protocol_with_contract_drift_replayed"] is True
    assert receipt["fixed_physical_contract_residual_contraction_verified"] is False
    assert receipt["pde_validated"] is False
    assert receipt["callable_velocity_delivery_blocked"] is False


def test_same_protocol_sha_does_not_identify_one_physical_contract() -> None:
    witness = mechanics_protocol_identity_witness()
    assert witness["residual_protocol_sha_same"] is True
    assert witness["physical_contract_same"] is False
    assert witness["free_residual_defined_forcing_used"] is False
    assert witness["before_residual_abs"] == pytest.approx(0.8)
    assert witness["after_residual_abs"] == pytest.approx(0.6)
    assert witness["apparent_contraction_factor"] == pytest.approx(0.75)
    assert witness["candidate_scientific_evidence"] is False


@pytest.mark.parametrize(
    "mutate",
    [
        lambda p: p["existing_stage_semantics"].__setitem__(
            "same_physical_contract_sha_before_after_enforced", True
        ),
        lambda p: p["promotion_guards"].__setitem__(
            "same_residual_protocol_sha_implies_same_physical_contract", True
        ),
        lambda p: p["promotion_guards"].__setitem__(
            "silent_pressure_or_forcing_drift_allowed_in_velocity_only_correction_claim", True
        ),
        lambda p: p["promotion_guards"].__setitem__(
            "pde_pending_blocks_callable_velocity_delivery", True
        ),
        lambda p: p["canonical_cr001_snapshot"].__setitem__(
            "pde_residual_max", 2.0e-3
        ),
        lambda p: p["canonical_cr001_snapshot"].__setitem__(
            "divergence_L2", 2.0e-5
        ),
        lambda p: p["provenance"]["public_source_fact"].append(
            "repository mechanics witness"
        ),
        lambda p: p["truth_boundary"].__setitem__("pde_validated", True),
        lambda p: p["truth_boundary"].__setitem__(
            "fixed_physical_contract_residual_contraction_verified", True
        ),
        lambda p: p["audited_upstream"].__setitem__("agent3_head", "0" * 40),
    ],
)
def test_semantic_laundering_and_cr001_drift_fail_closed(mutate) -> None:
    payload = copy.deepcopy(load_contract(ROOT))
    mutate(payload)
    with pytest.raises(ValueError):
        validate_contract(payload, ROOT)


def test_four_provenance_classes_remain_separate() -> None:
    payload = load_contract(ROOT)
    provenance = payload["provenance"]
    assert set(provenance) == {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    assert provenance["public_source_fact"] == []
    assert any("mechanics witness" in item for item in provenance["autonomous_design"])
    assert any("physical-contract SHA" in item for item in provenance["pending_unknown"])


def test_future_fixed_contract_claim_requires_more_than_protocol_sha() -> None:
    payload = load_contract(ROOT)
    required = payload[
        "required_future_binding_for_fixed_contract_contraction_claim"
    ]
    assert "same_residual_protocol_sha256" in required
    assert "physical_contract_sha256_before_equals_after" in required
    assert (
        "restricted_forcing_family_and_parameter_identity_before_equals_after_or_explicitly_typed_joint_update"
        in required
    )
    assert len(required) >= 6
