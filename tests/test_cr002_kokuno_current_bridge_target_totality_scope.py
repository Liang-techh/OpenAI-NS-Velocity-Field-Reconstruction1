from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_current_bridge_target_totality_scope import (
    _mechanics_full_eta_totality_witness,
    audit,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs/kokuno_current_bridge_target_totality_scope.json"


def _contract() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_cr002_current_bridge_target_totality_audit_passes() -> None:
    report = audit(ROOT)
    assert report["dependency_pr"] == 1233
    assert report["dependency_exact_head"] == "6ae44fdac17483b1fd042dca6f5bdf06701ae068"
    assert report["pointwise_diagnostic_mechanics"] == {
        "root_expansions": 12,
        "root_bisections": 96,
    }
    truth = report["truth_state"]
    assert truth["current_l_minus_h_q_s_transport_materialized"] is True
    assert truth["eta_resolved_pointwise_target_time_diagnostic_materialized"] is True
    assert truth["full_eta_target_hit_totality_established"] is False
    assert truth["smooth_eta_target_time_map_materialized"] is False
    assert truth["eta_dependent_cartesian_matching_boundary_materialized"] is False
    assert truth["kokuno_unified_global_velocity_export_ready"] is False
    assert report["canonical_states"]["velocity_export_ready"] is True
    assert report["canonical_states"]["pde_validated"] is False


def test_finite_hit_witness_does_not_certify_full_eta_totality() -> None:
    report = _mechanics_full_eta_totality_witness(_contract())
    assert report["finite_entry_margin_min"] > 0.0
    assert all(value >= 0.0 for value in report["finite_target_hit_times"])
    assert report["hidden_entry_margin"] < 0.0
    assert report["hidden_formal_hit_time"] < 0.0
    assert report["full_eta_target_totality_established"] is False


def test_contract_is_four_class_and_all_transfer_shortcuts_fail_closed() -> None:
    contract = _contract()
    assert set(contract["source_classification"]) == {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    assert contract["evidence_transfer"]
    assert all(value is False for value in contract["evidence_transfer"].values())
    assert contract["truth_state"]["pde_validated"] is False
    assert contract["truth_state"]["paper_exact"] is False
    assert contract["truth_state"]["openai_field_identified"] is False


def test_mechanics_witness_mutation_without_hidden_gap_fails_closed() -> None:
    contract = copy.deepcopy(_contract())
    contract["mechanics_witness"]["depth"] = 0.05
    with pytest.raises(AssertionError):
        _mechanics_full_eta_totality_witness(contract)
