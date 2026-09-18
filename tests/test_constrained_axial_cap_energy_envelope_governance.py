from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_axial_cap_energy_envelope_governance import (
    audit_axial_cap_energy_envelope_scope,
    load_constraints,
    load_project_status,
    load_scope,
)


def test_axial_cap_energy_envelope_scope_passes() -> None:
    report = audit_axial_cap_energy_envelope_scope()
    assert report["governed_upstream_pr"] == 327
    assert report["mode_name"] == "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL"
    assert 0.06 < report["fixed_parent_diagnostic_half_width"] < 0.07
    assert report["gross_q90_q99_reach_moved_inside_envelope"] is False
    assert report["live_next_task_preserved"] is True
    assert report["canonical_cr001_unchanged"] is True
    assert report["truth_boundary_preserved"] is True


def test_source_class_cannot_be_laundered_into_public_fact() -> None:
    scope = load_scope()
    scope["source_classification"]["axial_cap_band_and_envelope"] = "public_source_fact"
    with pytest.raises(ValueError, match="source class axial_cap_band_and_envelope"):
        audit_axial_cap_energy_envelope_scope(scope=scope)


def test_diagnostic_envelope_cannot_select_bound_or_coefficient() -> None:
    for key in (
        "diagnostic_envelope_is_selected_materialization_bound",
        "diagnostic_envelope_selects_coefficient",
        "diagnostic_exact_normalization_root_selects_coefficient",
    ):
        scope = load_scope()
        scope["energy_envelope_evidence_only"][key] = True
        with pytest.raises(ValueError, match=key):
            audit_axial_cap_energy_envelope_scope(scope=scope)


def test_no_q90_q99_motion_cannot_be_promoted_to_pde_or_visual_claim() -> None:
    for key in (
        "zero_q90_q99_movement_proves_no_pde_leverage",
        "zero_q90_q99_movement_proves_mode_useless",
        "relative_tail_gain_establishes_visual_correspondence",
        "morphology_screen_is_momentum_residual_jacobian",
    ):
        scope = load_scope()
        scope["morphology_semantics"][key] = True
        with pytest.raises(ValueError, match=key):
            audit_axial_cap_energy_envelope_scope(scope=scope)


def test_experimental_axial_cap_result_cannot_replace_live_materialization_route() -> None:
    scope = load_scope()
    scope["integration_routing_scope"][
        "axial_cap_energy_envelope_may_replace_live_next_task_without_explicit_integration_decision"
    ] = True
    with pytest.raises(ValueError, match="routing promotion guard"):
        audit_axial_cap_energy_envelope_scope(scope=scope)

    status = load_project_status()
    status["next_integration_task"] = "promote axial cap immediately"
    with pytest.raises(ValueError, match="project next task"):
        audit_axial_cap_energy_envelope_scope(project_status=status)


def test_cr001_threshold_or_force_drift_is_rejected() -> None:
    constraints = load_constraints()
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="pde_residual_max"):
        audit_axial_cap_energy_envelope_scope(constraints=constraints)

    constraints = load_constraints()
    constraints["forcing"]["restriction"] = "Allow residual-dependent basis"
    with pytest.raises(ValueError, match="free residual-dependent forcing"):
        audit_axial_cap_energy_envelope_scope(constraints=constraints)


def test_truth_states_remain_independent_and_false() -> None:
    for key in (
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        scope = load_scope()
        scope["truth_boundary"][key] = True
        with pytest.raises(ValueError, match=f"truth state {key}"):
            audit_axial_cap_energy_envelope_scope(scope=scope)


def test_joint_renormalization_requires_new_candidate_validation() -> None:
    scope = load_scope()
    scope["future_nonzero_materialization"]["joint_renormalization_requires_fresh_nonlinear_pde_validation"] = False
    with pytest.raises(ValueError, match="joint_renormalization_requires_fresh_nonlinear_pde_validation"):
        audit_axial_cap_energy_envelope_scope(scope=scope)

    scope = deepcopy(load_scope())
    scope["future_nonzero_materialization"]["new_candidate_sha_required"] = False
    with pytest.raises(ValueError, match="new_candidate_sha_required"):
        audit_axial_cap_energy_envelope_scope(scope=scope)
