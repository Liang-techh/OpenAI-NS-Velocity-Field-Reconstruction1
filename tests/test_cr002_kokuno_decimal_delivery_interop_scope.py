from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_decimal_delivery_interop_scope import (
    GovernanceError,
    PARENT_SOURCE_PATH,
    audit_canonical_constraints,
    audit_delivery_representation_claim,
    audit_parent_source_text,
    autonomous_transport_witness,
    load_contract,
    run_repo_audit,
)

ROOT = Path(__file__).resolve().parents[1]


def _contract():
    return load_contract(ROOT / "configs/kokuno_decimal_delivery_interop_scope.json")


def _constraints():
    return json.loads((ROOT / "configs/constraints.json").read_text())


def _parent_text():
    return (ROOT / PARENT_SOURCE_PATH).read_text()


def _parent_claim():
    contract = _contract()
    ident = contract["delivery_representation_identity"]
    return {
        "parent_semantic_sha256": "synthetic-parent-semantic-id",
        "scalar_type": ident["parent_scalar_type"],
        "precision_digits": ident["parent_precision_digits"],
        "rounding": ident["parent_rounding"],
        "binary64_embedding": ident["parent_binary64_embedding"],
        "export_encoding": ident["parent_export_encoding"],
        "correction_survival_replayed": True,
        "correction_survives_transport": True,
        "save_load_replayed": True,
        "velocity_export_ready": False,
        "binary64_total_velocity_export_ready": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def test_exact_parent_and_canonical_constraints_audit():
    contract = _contract()
    parent = audit_parent_source_text(_parent_text(), contract)
    constraints = audit_canonical_constraints(_constraints(), contract)

    assert parent["decimal_digits"] == 96
    assert parent["binary64_total_velocity_export_ready"] is False
    assert parent["unified_global_cartesian_velocity_export_ready"] is False
    assert parent["pde_validated"] is False

    assert constraints["nu"] == 0.01
    assert constraints["evaluation_box"] == [[-2, 2], [-2, 2], [-2, 2]]
    assert constraints["momentum_max"] == 0.001
    assert constraints["momentum_L2"] == 0.001
    assert constraints["divergence_max"] == 1e-5
    assert constraints["divergence_L2"] == 1e-5


def test_four_provenance_classes_and_pending_interop_are_explicit():
    contract = _contract()
    assert set(contract["provenance"]) == {
        "user_requirements",
        "public_source_facts",
        "autonomous_design",
        "pending_unknown",
    }
    pending = "\n".join(contract["provenance"]["pending_unknown"]).lower()
    assert "binary64" in pending
    assert "matlab" in pending
    assert "global" in pending
    assert "visual correspondence" in pending
    assert "paper-exact" in pending

    rules = contract["promotion_rules"]
    assert rules["decimal_stage_callable_is_delivery_surface"] is True
    assert rules["decimal_stage_callable_is_binary64_export"] is False
    assert rules["decimal_string_stage_export_is_global_export"] is False
    assert rules["downcast_that_erases_delta_may_inherit_parent_survival_evidence"] is False


def test_autonomous_witness_distinguishes_decimal_survival_from_binary64_downcast():
    witness = autonomous_transport_witness()
    assert witness["classification"] == "autonomous_mechanics_only"
    assert witness["delta_binary64"] != 0.0
    assert witness["binary64_base_plus_delta_equals_base"] is True
    assert witness["decimal_total_differs_from_base"] is True
    assert witness["decimal_total_minus_base_replays_delta"] is True
    assert witness["decimal_string_roundtrip_exact"] is True
    assert witness["binary64_downcast_equals_base"] is True
    assert witness["binary64_downcast_erases_nonzero_delta"] is True


def test_parent_truth_promotion_is_rejected():
    contract = _contract()
    mutated = _parent_text().replace(
        '"binary64_total_velocity_export_ready": False',
        '"binary64_total_velocity_export_ready": True',
        1,
    )
    with pytest.raises(GovernanceError, match="binary64_total_velocity_export_ready"):
        audit_parent_source_text(mutated, contract)


def test_parent_export_string_encoding_drift_is_rejected():
    contract = _contract()
    mutated = _parent_text().replace(
        '"values_row_major": [str(v) for v in velocity.reshape(-1)]',
        '"values_row_major": [float(v) for v in velocity.reshape(-1)]',
        1,
    )
    with pytest.raises(GovernanceError, match="strings"):
        audit_parent_source_text(mutated, contract)


def test_representation_change_requires_survival_and_save_load_replay():
    contract = _contract()
    claim = _parent_claim()
    claim.update(
        {
            "scalar_type": "binary64",
            "precision_digits": 53,
            "export_encoding": "MATLAB native double",
            "correction_survival_replayed": False,
            "correction_survives_transport": False,
            "save_load_replayed": False,
        }
    )
    with pytest.raises(GovernanceError, match="cannot inherit"):
        audit_delivery_representation_claim(claim, contract)

    claim["correction_survival_replayed"] = True
    with pytest.raises(GovernanceError, match="save/load"):
        audit_delivery_representation_claim(claim, contract)


def test_erasing_transport_cannot_promote_velocity_export_ready():
    contract = _contract()
    claim = _parent_claim()
    claim.update(
        {
            "scalar_type": "binary64",
            "precision_digits": 53,
            "export_encoding": "MATLAB native double",
            "correction_survival_replayed": True,
            "correction_survives_transport": False,
            "save_load_replayed": True,
            "velocity_export_ready": True,
            "binary64_total_velocity_export_ready": True,
        }
    )
    with pytest.raises(GovernanceError, match="velocity_export_ready"):
        audit_delivery_representation_claim(claim, contract)


def test_representation_receipt_cannot_promote_scientific_truth():
    contract = _contract()
    for key in ("pde_validated", "paper_exact", "openai_field_identified", "blowup_proved"):
        claim = _parent_claim()
        claim[key] = True
        with pytest.raises(GovernanceError, match=key):
            audit_delivery_representation_claim(claim, contract)


def test_constraint_threshold_relaxation_is_rejected():
    contract = _contract()
    constraints = copy.deepcopy(_constraints())
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(GovernanceError, match="momentum_max"):
        audit_canonical_constraints(constraints, contract)


def test_repo_audit_is_governance_only():
    receipt = run_repo_audit(ROOT)
    assert receipt["contract_schema"] == "cr002-kokuno-decimal-delivery-interop-v1"
    assert receipt["scientific_state_promoted"] is False
    assert receipt["mechanics_witness"]["binary64_downcast_erases_nonzero_delta"] is True
