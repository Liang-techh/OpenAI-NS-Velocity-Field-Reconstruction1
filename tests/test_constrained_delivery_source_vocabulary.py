import copy
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
DELIVERY_CONTRACT = ROOT / "configs" / "velocity_delivery_contract.json"
STATE_CONTRACT = ROOT / "configs" / "delivery_state_contract.json"
CANONICAL_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _assert_canonical(delivery: dict, state: dict) -> None:
    declared = set(state["classification_vocabulary"])
    assert declared == CANONICAL_CLASSES

    invalid = {
        entry["item"]: entry["classification"]
        for entry in delivery["source_classification"]
        if entry["classification"] not in declared
    }
    assert invalid == {}, f"noncanonical source classifications: {invalid}"


def test_velocity_delivery_contract_uses_canonical_source_vocabulary() -> None:
    _assert_canonical(_load(DELIVERY_CONTRACT), _load(STATE_CONTRACT))


@pytest.mark.parametrize(
    "bad_classification",
    ["user_required", "public_source", "autonomous", "pending", "source_guess"],
)
def test_velocity_delivery_contract_rejects_legacy_or_unknown_source_classes(
    bad_classification: str,
) -> None:
    delivery = _load(DELIVERY_CONTRACT)
    state = _load(STATE_CONTRACT)
    mutated = copy.deepcopy(delivery)
    mutated["source_classification"][0]["classification"] = bad_classification

    with pytest.raises(AssertionError, match="noncanonical source classifications"):
        _assert_canonical(mutated, state)
