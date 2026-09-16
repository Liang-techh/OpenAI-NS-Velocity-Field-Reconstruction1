import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_source_classification import (
    audit_repository_contracts,
    validate_eq45_source_classification,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "configs" / "eq45_source_contract.json"
DELIVERY = ROOT / "configs" / "delivery_state_contract.json"


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_repository_contracts_use_one_canonical_vocabulary():
    result = audit_repository_contracts(SOURCE, DELIVERY)
    assert result == {
        "canonical_classification_pass": True,
        "q_source_relation_pass": True,
        "claim_boundary_pass": True,
        "pde_validated": False,
        "paper_exact": False,
    }


def test_legacy_aliases_fail_closed():
    source = load(SOURCE)
    delivery = load(DELIVERY)
    source["source_classification"][0]["classification"] = "public_source"
    with pytest.raises(ValueError, match="noncanonical classification"):
        validate_eq45_source_classification(source, delivery)


def test_pending_numerical_profiles_cannot_be_promoted_to_public_source():
    source = load(SOURCE)
    delivery = load(DELIVERY)
    source["source_classification"][3]["classification"] = "public_source_fact"
    with pytest.raises(ValueError, match="pending_unknown"):
        validate_eq45_source_classification(source, delivery)


def test_negative_q_exponent_transcription_fails_closed():
    source = load(SOURCE)
    delivery = load(DELIVERY)
    source["coordinate_contract"]["implicit_q_relation"] = "q-z^2*q^(-2*h)=1-t"
    source["coordinate_contract"]["q_exponent_sign"] = "negative"
    with pytest.raises(ValueError, match=r"\+2h"):
        validate_eq45_source_classification(source, delivery)


def test_visual_or_paper_claim_cannot_be_promoted_by_source_contract():
    source = load(SOURCE)
    delivery = load(DELIVERY)
    source["claim_status"]["openai_correspondence_verified"] = True
    with pytest.raises(ValueError, match="openai_correspondence_verified"):
        validate_eq45_source_classification(source, delivery)

    source = load(SOURCE)
    source["claim_status"]["paper_exact"] = True
    with pytest.raises(ValueError, match="paper_exact"):
        validate_eq45_source_classification(source, delivery)
