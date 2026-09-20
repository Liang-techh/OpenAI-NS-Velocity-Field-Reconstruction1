from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st052m_tip_return_flow_scope import (
    GovernanceError,
    audit_contract,
    audit_live,
    load_contract,
    mutated_contract,
)


def test_live_tip_return_flow_scope_audit_passes() -> None:
    report = audit_live()
    assert report["status"] == "PASS"
    assert report["audited_pr"] == 723
    assert report["independent_physical_abs_z_sign_change"] == pytest.approx(
        1.46844135228393, abs=1e-12
    )
    assert report["fresh_714_tip_seed_in_inward_subband"] is True
    assert report["single_channel_return_flow_obstruction"] is True
    assert report["obstruction_applies_to_total_701_child"] is False
    assert report["obstruction_applies_to_arbitrary_divergence_free_fields"] is False
    assert report["whole_tip_support_inwardness_verified"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("scope_limits", "obstruction_proves_total_701_child_radial_velocity_has_return_flow"), True),
        (("scope_limits", "obstruction_proves_any_divergence_free_compact_field_requires_same_return_flow_pattern"), True),
        (("scope_limits", "fresh_714_pass_may_support_whole_tip_support_inwardness"), True),
        (("scope_limits", "fresh_714_pass_may_support_global_visual_correspondence"), True),
        (("scope_limits", "return_flow_obstruction_is_a_pde_acceptance_result"), True),
        (("scope_limits", "return_flow_obstruction_is_a_public_openai_field_identity_result"), True),
        (("claim_states", "full_candidate_tip_support_inwardness_verified"), True),
        (("claim_states", "full_candidate_return_flow_required"), True),
        (("claim_states", "arbitrary_divergence_free_return_flow_theorem_proved"), True),
        (("claim_states", "production_candidate_selected"), True),
        (("claim_states", "velocity_export_ready_promoted_by_this_audit"), True),
        (("claim_states", "visual_correspondence_verified"), True),
        (("claim_states", "pde_validated"), True),
        (("claim_states", "paper_exact"), True),
        (("claim_states", "openai_field_identified"), True),
        (("claim_states", "blowup_proved"), True),
    ],
)
def test_scope_and_claim_laundering_fail_closed(path: tuple[str, ...], value: object) -> None:
    contract = mutated_contract(load_contract(), path, value)
    with pytest.raises(GovernanceError):
        audit_contract(contract)


def test_local_fresh_path_evidence_cannot_be_erased_or_promoted_to_global() -> None:
    contract = load_contract()
    local_erased = mutated_contract(
        contract,
        ("scope_limits", "fresh_714_pass_may_support_local_post_freeze_directional_generalization"),
        False,
    )
    with pytest.raises(GovernanceError):
        audit_contract(local_erased)

    global_promoted = mutated_contract(
        contract,
        ("scope_limits", "fresh_714_pass_may_support_whole_tip_support_inwardness"),
        True,
    )
    with pytest.raises(GovernanceError):
        audit_contract(global_promoted)


def test_autonomous_tip_basis_cannot_be_relabelled_as_public_source_fact() -> None:
    contract = deepcopy(load_contract())
    items = contract["source_classification"]["items"]
    target = next(
        item
        for item in items
        if item["item"]
        == "701_tip_local_compact_odd_poloidal_channel_formula_support_orientation_and_alpha_rule"
    )
    target["class"] = "public_source_fact"
    with pytest.raises(GovernanceError):
        audit_contract(contract)


def test_whole_lobe_inwardness_cannot_be_relabelled_as_known_public_fact() -> None:
    contract = deepcopy(load_contract())
    items = contract["source_classification"]["items"]
    target = next(
        item
        for item in items
        if item["item"]
        == "claim_that_every_point_of_the_compact_tip_lobe_has_inward_radial_velocity"
    )
    target["class"] = "public_source_fact"
    with pytest.raises(GovernanceError):
        audit_contract(contract)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("cr001_lock", "nu"), 0.02),
        (("cr001_lock", "time_interval"), [0.2, 0.8]),
        (("cr001_lock", "forcing_bounds"), {"a": [0.0, 20.0], "c": [0.0, 10.0]}),
        (("cr001_lock", "reference_energy_abs_tolerance"), 0.01),
        (("cr001_lock", "validation_seed"), 914028),
        (("cr001_lock", "held_out_points"), 2048),
        (("cr001_lock", "derivative_steps"), [0.04, 0.02, 0.01]),
        (("cr001_lock", "quadrature_orders_per_axis"), [16, 32, 64]),
        (("cr001_lock", "divergence_max"), 1e-4),
        (("cr001_lock", "pde_residual_max"), 1e-2),
        (("cr001_lock", "residual_defined_or_pointwise_free_force_allowed"), True),
        (("cr001_lock", "amplitude_collapse_success_allowed"), True),
        (("cr001_lock", "threshold_relaxation_allowed"), True),
    ],
)
def test_cr001_drift_and_shortcuts_fail_closed(path: tuple[str, ...], value: object) -> None:
    contract = mutated_contract(load_contract(), path, value)
    with pytest.raises(GovernanceError):
        audit_contract(contract)
