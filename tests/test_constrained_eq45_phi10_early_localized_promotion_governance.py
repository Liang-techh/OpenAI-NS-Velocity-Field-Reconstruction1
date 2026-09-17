import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_phi10_early_localized_promotion_governance import (
    audit_phi10_early_localized_promotion_gate,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs" / "eq45_phi10_early_localized_promotion_gate.json"


def _contract():
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def test_early_localized_gate_keeps_export_visual_energy_and_pde_states_separate():
    result = audit_phi10_early_localized_promotion_gate(_contract(), repo_root=ROOT)

    assert result["contract_pass"] is True
    assert len(result["candidate_sha256"]) == 64
    assert len(result["base_supported_sha256"]) == 64
    assert result["velocity_export_ready"] is True
    assert result["visualization_candidate_only"] is True
    assert result["visualization_ready"] is False
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False
    assert result["paper_exact"] is False
    assert result["openai_field_identified"] is False
    assert result["candidate_energy_status"] == "pending_candidate_specific_revalidation"
    assert result["intermediate_morphology_status"] == "open_unconsumed_sibling_evidence"
    assert result["pde_evidence_scope"] == "pressure_free_diagnostics_not_formal_gate"


def test_target_free_temporal_shape_cannot_be_reclassified_as_public_correspondence():
    contract = copy.deepcopy(_contract())
    contract["classification"]["phi10_early_localized_temporal_shape"] = "public_source_fact"
    with pytest.raises(ValueError, match="classification"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["target_free_selection"]["may_support_public_visual_correspondence"] = True
    with pytest.raises(ValueError, match="public correspondence"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)


def test_open_intermediate_morphology_sibling_cannot_be_silently_consumed_or_promoted():
    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["intermediate_time_morphology"]["status"] = "accepted"
    with pytest.raises(ValueError, match="unconsumed"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["intermediate_time_morphology"]["may_promote_current_ancestry_state"] = True
    with pytest.raises(ValueError, match="cannot promote"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["policy"]["intermediate_resolution_stability_implies_public_visual_correspondence"] = True
    with pytest.raises(ValueError, match="forbidden promotion inference"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)


def test_pressure_free_pde_diagnostics_cannot_promote_pde_or_block_export():
    contract = copy.deepcopy(_contract())
    contract["current_state"]["pde_validated"] = True
    with pytest.raises(ValueError, match="current state"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["fixed_holdout_pde_tradeoff"]["formal_pde_gate_assessed"] = True
    with pytest.raises(ValueError, match="not the formal PDE gate"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["policy"]["pde_failure_blocks_velocity_export"] = True
    with pytest.raises(ValueError, match="forbidden promotion inference"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)


def test_pending_candidate_energy_cannot_be_inherited_or_turned_into_export_blocker():
    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["candidate_energy"]["may_inherit_parent_energy_acceptance"] = True
    with pytest.raises(ValueError, match="cannot inherit"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["candidate_energy"]["failure_blocks_velocity_export"] = True
    with pytest.raises(ValueError, match="cannot block"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["policy"]["candidate_energy_failure_blocks_velocity_export"] = True
    with pytest.raises(ValueError, match="forbidden promotion inference"):
        audit_phi10_early_localized_promotion_gate(contract, repo_root=ROOT)
