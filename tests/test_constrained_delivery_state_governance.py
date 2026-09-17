import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_delivery_state_governance import (
    load_delivery_state_contract,
    validate_delivery_state_contract,
)


CONTRACT = Path(__file__).resolve().parents[1] / "configs" / "delivery_state_contract.json"


def raw_contract():
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_repository_contract_is_valid():
    contract = load_delivery_state_contract(CONTRACT)
    assert contract["schema"] == "constrained_delivery_state_contract_v1"
    assert set(contract["classification_vocabulary"]) == {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }


def test_export_ready_does_not_promote_pde_or_paper_identity():
    contract = raw_contract()
    export = contract["states"]["velocity_export_ready"]
    assert "pde_validated" in export["does_not_imply"]
    assert "paper_exact" in export["does_not_imply"]

    mutated = copy.deepcopy(contract)
    mutated["states"]["velocity_export_ready"]["does_not_imply"].remove("pde_validated")
    with pytest.raises(ValueError, match="missing forbidden implication"):
        validate_delivery_state_contract(mutated)


def test_visual_correspondence_cannot_promote_pde():
    mutated = raw_contract()
    mutated["states"]["visual_correspondence_verified"]["implies"] = ["pde_validated"]
    with pytest.raises(ValueError, match="illegally implies"):
        validate_delivery_state_contract(mutated)


def test_pde_claim_cannot_drop_preregistered_held_out_evidence():
    mutated = raw_contract()
    mutated["states"]["pde_validated"]["positive_evidence"].remove("preregistered_thresholds")
    with pytest.raises(ValueError, match="pde_validated"):
        validate_delivery_state_contract(mutated)


def test_paper_exact_requires_public_identity_evidence_and_pending_default():
    mutated = raw_contract()
    mutated["states"]["paper_exact"]["positive_evidence"] = ["visual_similarity"]
    with pytest.raises(ValueError, match="paper_exact"):
        validate_delivery_state_contract(mutated)

    mutated = raw_contract()
    mutated["states"]["paper_exact"]["allowed_evidence_origins"] = ["autonomous_design"]
    with pytest.raises(ValueError, match="public-source"):
        validate_delivery_state_contract(mutated)


def test_unknown_origin_class_fails_closed():
    mutated = raw_contract()
    mutated["states"]["velocity_export_ready"]["definition_origin"] = "inferred_fact"
    with pytest.raises(ValueError, match="origin class"):
        validate_delivery_state_contract(mutated)
