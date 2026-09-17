import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_phi10_cubic_promotion_governance import (
    audit_phi10_cubic_promotion_gate,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs" / "eq45_phi10_cubic_promotion_gate.json"


def _contract():
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def test_cubic_gate_keeps_export_snapshot_derivative_visual_and_pde_states_separate():
    result = audit_phi10_cubic_promotion_gate(_contract(), repo_root=ROOT)

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
    assert result["static_return_times"] == [0.5, 0.625, 0.75]
    assert result["snapshot_identity_not_derivative_identity"] is True
    assert result["candidate_energy_status"] == "pending_candidate_specific_revalidation"
    assert result["force_sibling_status"] == "open_unconsumed_sibling_evidence"


def test_target_free_cubic_shape_cannot_be_reclassified_as_public_source_or_correspondence():
    contract = copy.deepcopy(_contract())
    contract["classification"]["phi10_cubic_temporal_shape"] = "public_source_fact"
    with pytest.raises(ValueError, match="classification"):
        audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["target_free_cubic_selection"]["may_support_public_visual_correspondence"] = True
    with pytest.raises(ValueError, match="public correspondence"):
        audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)


def test_static_return_snapshots_cannot_inherit_parent_temporal_or_scientific_acceptance():
    for key in (
        "may_inherit_parent_time_derivative_evidence",
        "may_inherit_parent_pde_acceptance",
        "may_inherit_parent_energy_acceptance",
        "may_inherit_parent_visual_acceptance",
        "may_inherit_parent_canonical_status",
    ):
        contract = copy.deepcopy(_contract())
        contract["evidence_scope"]["static_return_keyframes"][key] = True
        with pytest.raises(ValueError, match="snapshot identity cannot enable"):
            audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["policy"]["snapshot_identity_implies_time_derivative_identity"] = True
    with pytest.raises(ValueError, match="forbidden promotion inference"):
        audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)


def test_temporal_derivative_convergence_cannot_promote_formal_pde_validation():
    contract = copy.deepcopy(_contract())
    contract["current_state"]["pde_validated"] = True
    with pytest.raises(ValueError, match="current state"):
        audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["public_output_temporal_derivative"]["formal_pde_gate_assessed"] = True
    with pytest.raises(ValueError, match="not formal PDE acceptance"):
        audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["policy"]["temporal_derivative_convergence_implies_pde_validation"] = True
    with pytest.raises(ValueError, match="forbidden promotion inference"):
        audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)


def test_pending_energy_and_unconsumed_force_sibling_cannot_block_export_or_promote_pde():
    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["candidate_energy"]["may_inherit_parent_energy_acceptance"] = True
    with pytest.raises(ValueError, match="cannot inherit"):
        audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["candidate_energy"]["failure_blocks_velocity_export"] = True
    with pytest.raises(ValueError, match="cannot block"):
        audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["evidence_scope"]["restricted_force_sibling"]["may_promote_pde_validated"] = True
    with pytest.raises(ValueError, match="cannot promote PDE validity"):
        audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)

    contract = copy.deepcopy(_contract())
    contract["policy"]["pde_failure_blocks_velocity_export"] = True
    with pytest.raises(ValueError, match="forbidden promotion inference"):
        audit_phi10_cubic_promotion_gate(contract, repo_root=ROOT)
