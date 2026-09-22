from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction import audit_kokuno_postpulse_tf_representation_scope as m

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs" / "kokuno_postpulse_tf_representation_scope.json"
CONSTRAINTS = ROOT / "configs" / "constraints.json"
CANDIDATE = ROOT / "src" / "openai_ns_reconstruction" / "kokuno_pa16_current_cartesian_postpulse_eta_flattening.py"
PROJECT_STATUS = ROOT / "project_status.json"


def _contract() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_autonomous_tf_witness_is_nonvacuous_and_not_candidate_evidence():
    witness = _contract()["mechanics_witness"]
    assert witness["classification"] == "autonomous_mechanics_only"
    receipt = m.tf_sensitivity_witness(witness["T_f_a"], witness["T_f_b"], witness["y_probe"], witness["eta_probe"])
    assert receipt["different_flattening_weight"] is True
    assert receipt["different_logE_correction"] is True
    assert receipt["absolute_correction_delta"] > 0.0


def test_witness_rejects_degenerate_inputs():
    with pytest.raises(ValueError):
        m.tf_sensitivity_witness(100.0, 100.0, 50.0, 0.5)
    with pytest.raises(ValueError):
        m.tf_sensitivity_witness(100.0, 120.0, 0.0, 0.5)
    with pytest.raises(ValueError):
        m.tf_sensitivity_witness(100.0, 120.0, 50.0, 1.0)


def test_live_contract_and_exact_parent_are_fail_closed():
    assert m.audit_contract(CONTRACT, CONSTRAINTS, CANDIDATE, PROJECT_STATUS) == []


@pytest.mark.parametrize(
    ("path", "value"),
    (
        (("tf_representation_scope", "source_supplies_unique_numeric_T_f"), True),
        (("tf_representation_scope", "source_exact_T_f_recovered"), True),
        (("tf_representation_scope", "current_autonomous_T_f"), 120.0),
        (("numeric_tf_change_promotion_rule", "reuse_old_ell_or_divergence_receipts_for_new_T_f"), True),
        (("numeric_tf_change_promotion_rule", "numeric_T_f_change_requires_new_semantic_candidate_identity"), False),
        (("independent_delivery_states", "kokuno_global_velocity_export_ready"), True),
        (("independent_delivery_states", "kokuno_pde_validated"), True),
        (("cr001_invariants", "threshold_relaxation_this_increment"), True),
        (("cr001_invariants", "forcing_family_change_this_increment"), True),
    ),
)
def test_truth_or_constraint_mutations_fail_closed(tmp_path, path, value):
    payload = copy.deepcopy(_contract())
    node = payload
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value
    mutated = tmp_path / "contract.json"
    mutated.write_text(json.dumps(payload), encoding="utf-8")
    assert m.audit_contract(mutated, CONSTRAINTS, CANDIDATE, PROJECT_STATUS)


def test_four_way_provenance_classes_are_exact(tmp_path):
    payload = copy.deepcopy(_contract())
    payload["four_way_provenance"]["source_exact"] = ["not a canonical class"]
    mutated = tmp_path / "contract.json"
    mutated.write_text(json.dumps(payload), encoding="utf-8")
    errors = m.audit_contract(mutated, CONSTRAINTS, CANDIDATE, PROJECT_STATUS)
    assert any("four_way provenance classes" in error for error in errors)


def test_required_evidence_transfer_firewall_is_machine_locked(tmp_path):
    payload = copy.deepcopy(_contract())
    payload["evidence_transfer_firewall"]["forbidden_implications"] = ["stage callable does not imply global"]
    mutated = tmp_path / "contract.json"
    mutated.write_text(json.dumps(payload), encoding="utf-8")
    errors = m.audit_contract(mutated, CONSTRAINTS, CANDIDATE, PROJECT_STATUS)
    assert any("forbidden_implications" in error for error in errors)
