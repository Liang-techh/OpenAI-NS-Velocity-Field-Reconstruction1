from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_relative_swirl_split_delivery_scope import (
    CONTRACT_PATH,
    UPSTREAM_PATH,
    ScopeAuditError,
    audit_repository,
    audit_upstream_source,
    autonomous_binary64_split_witness,
    validate_contract_payload,
)


ROOT = Path(__file__).resolve().parents[1]


def _contract() -> dict:
    return json.loads((ROOT / CONTRACT_PATH).read_text(encoding="utf-8"))


def test_repository_scope_replays_exact_parent_and_canonical_gates() -> None:
    receipt = audit_repository(ROOT)
    assert receipt["scope_id"] == "cr002-kokuno-relative-swirl-split-delivery-scope-143"
    assert receipt["exact_dependency_head"] == "657818dd83e119e6091bd2490d3804c04c3ef723"
    assert receipt["upstream_scope"] == {
        "profile_delta_channel_materialized": True,
        "cartesian_delta_channel_materialized": True,
        "velocity_split_method_present": True,
        "single_array_velocity_method_defined_on_split_class": False,
        "binary64_total_sum_resolved": False,
        "current_cartesian_relative_swirl_composed": False,
    }
    assert receipt["scientific_state_promoted"] is False
    assert receipt["candidate_bytes_changed_by_this_audit"] is False
    assert receipt["cr001_threshold_changed_by_this_audit"] is False


def test_four_way_provenance_and_delivery_states_stay_separate() -> None:
    contract = _contract()
    validate_contract_payload(contract)
    assert set(contract["provenance_classes"]) == {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    state = contract["materialization_state"]
    assert state["relative_swirl_cartesian_delta_channel_materialized"] is True
    assert state["split_configuration_save_load_materialized"] is True
    assert state["binary64_total_relative_swirl_sum_is_resolved"] is False
    assert state["single_array_total_velocity_materialized"] is False
    assert state["current_cartesian_relative_swirl_composed"] is False
    assert state["kokuno_unified_global_cartesian_velocity_export_ready"] is False
    assert state["visual_correspondence_verified"] is False
    assert state["pde_validated"] is False
    assert state["paper_exact"] is False
    assert state["openai_field_identified"] is False


def test_autonomous_witness_is_nonzero_but_binary64_total_collapses() -> None:
    witness = autonomous_binary64_split_witness()
    assert witness["role"] == "autonomous_mechanics_only"
    assert witness["delta_nonzero"] is True
    assert witness["delta"] > 0.0
    assert witness["ordinary_binary64_total"] == witness["base"]
    assert witness["ordinary_total_equals_base"] is True
    assert witness["separate_channel_roundtrip_preserves_delta"] is True
    assert witness["source_evidence"] is False
    assert witness["openai_numerical_evidence"] is False
    assert witness["pde_evidence"] is False


@pytest.mark.parametrize(
    "key",
    [
        "binary64_total_relative_swirl_sum_is_resolved",
        "single_array_total_velocity_materialized",
        "current_cartesian_relative_swirl_composed",
        "kokuno_unified_global_cartesian_velocity_export_ready",
        "visual_correspondence_verified",
        "heldout_complete_ns_residual_assessed",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ],
)
def test_false_delivery_scientific_states_cannot_be_promoted(key: str) -> None:
    contract = _contract()
    mutated = copy.deepcopy(contract)
    mutated["materialization_state"][key] = True
    with pytest.raises(ScopeAuditError):
        validate_contract_payload(mutated)


@pytest.mark.parametrize(
    "key",
    [
        "split_channels_do_not_imply_single_array_total",
        "split_configuration_save_load_does_not_imply_composed_velocity_roundtrip",
        "nonzero_split_delta_does_not_imply_nonzero_binary64_total_minus_base",
        "binary64_total_equals_base_does_not_imply_mathematical_delta_is_zero",
        "single_array_api_name_cannot_discard_delta_and_inherit_split_identity",
        "future_composition_requires_new_delivery_representation_identity",
        "future_composition_requires_exact_evaluate_save_load_export_survival_evidence",
        "future_composition_evidence_does_not_promote_visual_pde_paper_or_openai_truth",
    ],
)
def test_evidence_transfer_firewall_cannot_be_disabled(key: str) -> None:
    contract = _contract()
    mutated = copy.deepcopy(contract)
    mutated["evidence_transfer_rules"][key] = False
    with pytest.raises(ScopeAuditError):
        validate_contract_payload(mutated)


def test_upstream_truth_promotion_requires_new_governance_identity() -> None:
    source = (ROOT / UPSTREAM_PATH).read_text(encoding="utf-8")
    assert '"current_cartesian_relative_swirl_composed": False' in source
    mutated = source.replace(
        '"current_cartesian_relative_swirl_composed": False',
        '"current_cartesian_relative_swirl_composed": True',
        1,
    )
    with pytest.raises(ScopeAuditError):
        audit_upstream_source(mutated)


def test_upstream_single_array_velocity_surface_cannot_appear_silently() -> None:
    source = (ROOT / UPSTREAM_PATH).read_text(encoding="utf-8")
    needle = "    def velocity_split(\n"
    assert needle in source
    inserted = (
        "    def velocity(self, x, y, z, t):\n"
        "        return self.velocity_split(x, y, z, t)[0]\n\n"
        + needle
    )
    mutated = source.replace(needle, inserted, 1)
    with pytest.raises(ScopeAuditError):
        audit_upstream_source(mutated)


def test_canonical_constraints_cannot_be_relabelled_in_contract() -> None:
    contract = _contract()
    mutated = copy.deepcopy(contract)
    mutated["canonical_cr001"]["momentum_max"] = 0.01
    with pytest.raises(ScopeAuditError):
        validate_contract_payload(mutated)


def test_mechanics_witness_cannot_be_promoted_to_source_evidence() -> None:
    contract = _contract()
    mutated = copy.deepcopy(contract)
    mutated["autonomous_mechanics_witness"]["not_source_evidence"] = False
    with pytest.raises(ScopeAuditError):
        validate_contract_payload(mutated)
