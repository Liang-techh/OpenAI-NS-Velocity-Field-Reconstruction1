import copy
import json

import pytest

from openai_ns_reconstruction.constrained_st052m_tip_basis_selection_data_audit import (
    CONTRACT_PATH,
    audit,
)


def _contract():
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def _write(tmp_path, value):
    path = tmp_path / "contract.json"
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def test_live_tip_basis_selection_data_audit_passes():
    result = audit()
    assert result["passed"] is True
    assert result["upstream"]["parent_positions_consumed"] is True
    assert result["upstream"]["same_protocol_replayed"] is True
    assert result["selection_data_scope"].startswith("development/model-selection evidence")


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        (
            "selection_data_accounting",
            "calibration_consumes_frozen_parent_tip_path_positions",
            False,
        ),
        (
            "selection_data_accounting",
            "path_screen_is_independent_held_out_validation",
            True,
        ),
        (
            "selection_data_accounting",
            "post_calibration_path_screen_reuses_same_frozen_seed_protocol",
            False,
        ),
        ("promotion_boundary", "new_spatial_basis_promoted", True),
        ("promotion_boundary", "production_candidate_selected", True),
        ("promotion_boundary", "velocity_export_ready_promoted_by_this_screen", True),
        ("promotion_boundary", "independent_held_out_visual_validation_performed", True),
        ("promotion_boundary", "visual_correspondence_verified", True),
        ("promotion_boundary", "held_out_temporal_child_pde_residual_evaluated", True),
        ("promotion_boundary", "pde_validated", True),
        ("promotion_boundary", "paper_exact", True),
        ("promotion_boundary", "openai_field_identified", True),
        (
            "fresh_evidence_requirement",
            "fresh_path_seeds_or_other_explicitly_held_out_path_protocol_required_for_independent_path_claim",
            False,
        ),
    ],
)
def test_selection_or_promotion_mutations_fail_closed(tmp_path, section, key, value):
    contract = _contract()
    contract[section][key] = value
    with pytest.raises(ValueError):
        audit(contract_path=_write(tmp_path, contract))


def test_cr001_threshold_relaxation_fails_closed(tmp_path):
    contract = _contract()
    contract["cr001_freeze"]["momentum_residual_max"] = 0.01
    with pytest.raises(ValueError):
        audit(contract_path=_write(tmp_path, contract))


def test_cr001_free_force_route_cannot_be_relabelled_allowed(tmp_path):
    contract = _contract()
    contract["cr001_freeze"]["residual_defined_or_pointwise_free_force_forbidden"] = False
    with pytest.raises(ValueError):
        audit(contract_path=_write(tmp_path, contract))


def test_four_way_source_classification_is_required(tmp_path):
    contract = _contract()
    contract["source_classification"]["public_source_fact"] = []
    with pytest.raises(ValueError):
        audit(contract_path=_write(tmp_path, contract))
