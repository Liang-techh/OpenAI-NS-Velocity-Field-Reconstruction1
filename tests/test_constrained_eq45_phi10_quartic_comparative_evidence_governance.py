from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_phi10_quartic_comparative_evidence_governance import (
    audit_quartic_comparative_evidence_gate,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs/eq45_phi10_quartic_comparative_evidence_gate.json"
CONSTRAINTS = ROOT / "configs/constraints.json"
DELIVERY = ROOT / "configs/delivery_state_contract.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write(tmp_path: Path, name: str, value: dict) -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    return path


def _audit_mutated(tmp_path: Path, contract: dict, *, constraints: dict | None = None):
    return audit_quartic_comparative_evidence_gate(
        contract_path=_write(tmp_path, "contract.json", contract),
        constraints_path=_write(
            tmp_path,
            "constraints.json",
            constraints if constraints is not None else _load(CONSTRAINTS),
        ),
        delivery_state_path=DELIVERY,
    )


def test_current_quartic_comparative_evidence_gate_passes():
    report = audit_quartic_comparative_evidence_gate()
    assert report["quartic_vs_cubic_pde_ordering"] == "unresolved_seed_sensitive"
    assert report["fresh_seed_quartic_better_than_cubic_count"] == 1
    assert report["fresh_seed_count"] == 3
    assert report["fixed_holdout_quartic_vs_cubic_fractional_change"] < 0.0
    assert report["fresh_seed_quartic_vs_cubic_mean_fractional_change"] > 0.0
    assert report["fresh_seed_quartic_vs_static_mean_fractional_change"] < 0.0
    assert report["velocity_export_ready"] is True
    assert report["pde_validated"] is False
    assert report["visual_correspondence_verified"] is False


def test_rejects_formal_pde_promotion_of_pressure_free_diagnostic(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["evidence"]["fresh_seed_generalization"]["formal_pde_gate_assessed"] = True
    with pytest.raises(ValueError, match="formal PDE gate"):
        _audit_mutated(tmp_path, contract)


def test_rejects_single_holdout_as_pde_superiority(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["governed_conclusions"][
        "single_fixed_holdout_may_establish_pde_superiority"
    ] = True
    with pytest.raises(ValueError, match="must remain false"):
        _audit_mutated(tmp_path, contract)


def test_rejects_erasing_seed_sensitive_ordering(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    rows = contract["evidence"]["fresh_seed_generalization"]["rows"]
    rows[0]["quartic_zero_rms"] = rows[0]["cubic_zero_rms"] * 0.99
    rows[1]["quartic_zero_rms"] = rows[1]["cubic_zero_rms"] * 0.99
    with pytest.raises(ValueError, match="seed-sensitive"):
        _audit_mutated(tmp_path, contract)


def test_rejects_posthoc_pde_threshold_relaxation(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    constraints = copy.deepcopy(_load(CONSTRAINTS))
    constraints["validation"]["thresholds"]["pde_residual_L2"] = 0.01
    with pytest.raises(ValueError, match="pde_residual_L2 threshold drifted"):
        _audit_mutated(tmp_path, contract, constraints=constraints)


def test_rejects_reclassifying_internal_diagnostic_as_public_fact(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["source_classification"][
        "fresh_seed_pressure_free_vorticity_generalization"
    ] = "public_source_fact"
    with pytest.raises(ValueError, match="source classification drifted"):
        _audit_mutated(tmp_path, contract)


def test_rejects_consuming_open_sibling_without_replay(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["evidence"]["late_morphology"]["ancestry_status"] = "current_ancestry_evidence"
    with pytest.raises(ValueError, match="unconsumed sibling evidence"):
        _audit_mutated(tmp_path, contract)


def test_rejects_scientific_state_promotion_or_export_blocker(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["truth_states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-state boundary"):
        _audit_mutated(tmp_path, contract)

    contract = copy.deepcopy(_load(CONTRACT))
    contract["truth_states"]["velocity_export_ready"] = False
    with pytest.raises(ValueError, match="truth-state boundary"):
        _audit_mutated(tmp_path, contract)
