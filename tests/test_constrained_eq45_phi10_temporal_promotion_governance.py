import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_phi10_temporal_promotion_governance import (
    audit_phi10_temporal_promotion_gate,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs" / "eq45_phi10_temporal_promotion_gate.json"


def _contract():
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def test_phi10_temporal_promotion_gate_keeps_delivery_visual_and_pde_states_separate():
    result = audit_phi10_temporal_promotion_gate(_contract(), repo_root=ROOT)

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
    assert result["whole_domain_morphology_status"] == "open_unconsumed_sibling_evidence"
    assert result["pde_tradeoff"] == "mixed_seed_sensitive_not_formal_gate"


def test_target_free_morphology_selection_cannot_be_reclassified_as_public_correspondence():
    contract = copy.deepcopy(_contract())
    contract["classification"]["phi10_temporal_mode_and_slope_selection"] = "public_source_fact"

    with pytest.raises(ValueError, match="classification"):
        audit_phi10_temporal_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["target_free_selection"]["may_support_public_visual_correspondence"] = True
    with pytest.raises(ValueError, match="public correspondence"):
        audit_phi10_temporal_promotion_gate(contract, repo_root=ROOT)


def test_mixed_pde_diagnostics_cannot_promote_pde_validity_or_block_export():
    contract = copy.deepcopy(_contract())
    contract["current_state"]["pde_validated"] = True
    with pytest.raises(ValueError, match="current state"):
        audit_phi10_temporal_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["policy"]["pde_failure_blocks_velocity_export"] = True
    with pytest.raises(ValueError, match="forbidden promotion inference"):
        audit_phi10_temporal_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["fresh_seed_pde_generalization"]["fresh_seed_fractional_changes"] = [-0.01, -0.02, -0.03]
    with pytest.raises(ValueError, match="mixed"):
        audit_phi10_temporal_promotion_gate(contract, repo_root=ROOT)


def test_open_whole_domain_sibling_cannot_be_silently_consumed_or_promoted():
    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["whole_domain_morphology"]["status"] = "accepted"
    with pytest.raises(ValueError, match="unconsumed"):
        audit_phi10_temporal_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["whole_domain_morphology"]["may_promote_current_ancestry_state"] = True
    with pytest.raises(ValueError, match="cannot promote"):
        audit_phi10_temporal_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["policy"]["production_slope_is_promoted"] = True
    with pytest.raises(ValueError, match="forbidden promotion inference"):
        audit_phi10_temporal_promotion_gate(contract, repo_root=ROOT)
