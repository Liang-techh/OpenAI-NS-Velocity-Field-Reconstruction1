from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.constrained_st052m_reservoir_axial_capacity_scope_audit import (
    audit,
    load_default_inputs,
)


def _inputs():
    contract, source_text, constraints = load_default_inputs()
    return copy.deepcopy(contract), source_text, copy.deepcopy(constraints)


def test_live_axial_capacity_scope_audit_passes() -> None:
    receipt = audit()
    assert receipt["scope"] == "correction_channel_unit_response_only"
    assert receipt["total_child_axial_velocity_evaluated"] is False
    assert receipt["public_numeric_core_target"] is None
    assert receipt["visual_correspondence_verified"] is False
    assert receipt["pde_validated"] is False
    assert receipt["cr001_unchanged"] is True


def test_rejects_total_child_core_radius_laundering() -> None:
    contract, source_text, constraints = _inputs()
    contract["channel_vs_total_field"]["channel_core_radius_may_be_called_total_child_core_radius"] = True
    with pytest.raises(AssertionError, match="claim laundering"):
        audit(contract, source_text, constraints)


def test_rejects_total_child_flux_laundering() -> None:
    contract, source_text, constraints = _inputs()
    contract["channel_vs_total_field"]["channel_zero_axial_flux_may_be_called_total_child_zero_axial_flux"] = True
    with pytest.raises(AssertionError, match="claim laundering"):
        audit(contract, source_text, constraints)


def test_rejects_invented_public_numeric_core_target() -> None:
    contract, source_text, constraints = _inputs()
    contract["frozen_channel_geometry"]["public_numeric_core_target"] = 1.0585450412096833
    with pytest.raises(AssertionError, match="invented public numerical core target"):
        audit(contract, source_text, constraints)


def test_rejects_source_provenance_laundering() -> None:
    contract, source_text, constraints = _inputs()
    contract["frozen_channel_geometry"]["geometry_provenance"] = "public_source_fact"
    with pytest.raises(AssertionError, match="laundered into source provenance"):
        audit(contract, source_text, constraints)


def test_rejects_parent_or_frozen_alpha_being_claimed_in_unit_shape_factor() -> None:
    contract, source_text, constraints = _inputs()
    contract["channel_vs_total_field"]["parent_velocity_contribution_included_in_axial_shape_factor"] = True
    with pytest.raises(AssertionError, match="claim laundering"):
        audit(contract, source_text, constraints)

    contract, source_text, constraints = _inputs()
    contract["channel_vs_total_field"]["frozen_child_coefficient_included_in_axial_shape_factor"] = True
    with pytest.raises(AssertionError, match="claim laundering"):
        audit(contract, source_text, constraints)


def test_rejects_live_pr785_truth_promotion() -> None:
    contract, source_text, constraints = _inputs()
    mutated = source_text.replace('"candidate_velocity_changed": False', '"candidate_velocity_changed": True', 1)
    assert mutated != source_text
    with pytest.raises(AssertionError, match="truth bit promoted"):
        audit(contract, mutated, constraints)


def test_rejects_visual_or_pde_promotion() -> None:
    contract, source_text, constraints = _inputs()
    contract["truth_boundary"]["visual_correspondence_verified"] = True
    with pytest.raises(AssertionError, match="truth-boundary promotion"):
        audit(contract, source_text, constraints)

    contract, source_text, constraints = _inputs()
    contract["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError, match="truth-boundary promotion"):
        audit(contract, source_text, constraints)


def test_rejects_free_force_shortcut() -> None:
    contract, source_text, constraints = _inputs()
    constraints["forcing"]["restriction"] = "pointwise free force allowed"
    with pytest.raises(AssertionError, match="free/residual-defined forcing shortcut"):
        audit(contract, source_text, constraints)


def test_rejects_cr001_threshold_relaxation() -> None:
    contract, source_text, constraints = _inputs()
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.002
    with pytest.raises(AssertionError, match="CR001 threshold drift"):
        audit(contract, source_text, constraints)


def test_rejects_cr001_amplitude_collapse_guard_drift() -> None:
    contract, source_text, constraints = _inputs()
    constraints["nontriviality"]["enforcement"] = "normalize if convenient"
    with pytest.raises(AssertionError, match="amplitude-collapse guard drift"):
        audit(contract, source_text, constraints)
