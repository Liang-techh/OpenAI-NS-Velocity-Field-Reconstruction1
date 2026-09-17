from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_h_regime_governance import (
    audit_checked_eq45_h_regime,
    audit_eq45_h_regime,
)


ROOT = Path(__file__).resolve().parents[1]


def _payloads():
    contract = json.loads((ROOT / "configs" / "eq45_h_regime_contract.json").read_text())
    source = json.loads((ROOT / "configs" / "eq45_source_contract.json").read_text())
    candidate = json.loads(
        (ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json").read_text()
    )
    return contract, source, candidate


def test_checked_seed_is_in_coordinate_and_theorem_h_regimes():
    report = audit_checked_eq45_h_regime(ROOT)
    assert report["contract_pass"] is True
    assert report["candidate_h"] == pytest.approx(0.005)
    assert report["coordinate_regime_valid"] is True
    assert report["theorem_regime_matched"] is True
    assert report["velocity_export_ready"] is True
    assert report["pde_validated"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["paper_exact"] is False


def test_coordinate_valid_non_theorem_h_remains_exportable_without_claim_promotion():
    contract, source, candidate = _payloads()
    candidate = deepcopy(candidate)
    candidate["h"] = 0.1

    report = audit_eq45_h_regime(contract, source_contract=source, candidate=candidate)
    assert report["coordinate_regime_valid"] is True
    assert report["theorem_regime_matched"] is False
    assert report["velocity_export_ready"] is True
    assert report["theorem_regime_required_for_velocity_export"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False


def test_h_origin_and_scientific_claim_promotions_fail_closed():
    contract, source, candidate = _payloads()

    wrong_origin = deepcopy(candidate)
    wrong_origin["classification"]["h"] = "public_source_fact"
    with pytest.raises(ValueError, match="autonomous_design"):
        audit_eq45_h_regime(contract, source_contract=source, candidate=wrong_origin)

    promoted = deepcopy(candidate)
    promoted["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="paper_exact"):
        audit_eq45_h_regime(contract, source_contract=source, candidate=promoted)

    bad_semantics = deepcopy(contract)
    bad_semantics["delivery_semantics"]["theorem_regime_required_for_velocity_export"] = True
    with pytest.raises(ValueError, match="export gate"):
        audit_eq45_h_regime(bad_semantics, source_contract=source, candidate=candidate)


def test_coordinate_h_bounds_and_source_drift_fail_closed():
    contract, source, candidate = _payloads()

    for invalid_h in (0.0, -0.01, 0.5, 0.7):
        mutated = deepcopy(candidate)
        mutated["h"] = invalid_h
        with pytest.raises(ValueError, match="coordinate regime"):
            audit_eq45_h_regime(contract, source_contract=source, candidate=mutated)

    drifted_source = deepcopy(source)
    drifted_source["coordinate_contract"]["theorem_h_upper_bound_strict"] = 0.02
    with pytest.raises(ValueError, match="theorem h bound drift"):
        audit_eq45_h_regime(contract, source_contract=drifted_source, candidate=candidate)
