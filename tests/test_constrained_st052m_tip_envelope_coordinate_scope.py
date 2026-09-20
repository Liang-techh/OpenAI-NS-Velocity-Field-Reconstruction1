from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st052m_tip_envelope_coordinate_scope import (
    GovernanceError,
    audit,
    independent_unit_response_witness,
    load_contract,
)


def test_exact_732_envelope_coordinate_scope_audits_cleanly():
    result = audit()
    assert result["status"] == "pass"
    assert result["selected_inner_exponent_p"] == 9
    assert result["physical_abs_z_sign_change"] >= 1.60
    assert result["old_701_alpha_transfer_validated"] is False
    assert result["next_child_requires_explicit_frozen_alpha_semantics"] is True


def test_same_dimension_peak_normalization_does_not_preserve_unit_response():
    witness = independent_unit_response_witness()
    assert witness["baseline_axial_z_factor"] > 0.0
    assert witness["reshaped_axial_z_factor"] > 0.0
    assert witness["reshaped_to_baseline_magnitude_ratio"] < 0.1
    assert witness["reshaped_to_baseline_magnitude_ratio"] != pytest.approx(1.0)


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("representation_scope", "spatial_basis_function_identical", True),
        ("representation_scope", "unit_velocity_response_identical", True),
        ("representation_scope", "coefficient_coordinate_identity_preserved", True),
        ("representation_scope", "old_701_alpha_transfer_validated", True),
        ("representation_scope", "new_alpha_calibration_performed_by_732", True),
        ("representation_scope", "nonlinear_child_evaluated_by_732", True),
        ("promotion_rules", "732_pass_may_claim_old_alpha_is_invariant_under_reshape", True),
        ("promotion_rules", "next_child_must_freeze_alpha_semantics_before_result_evaluation", False),
        (
            "promotion_rules",
            "rederiving_alpha_from_701_historical_tip_paths_remains_development_selection_data",
            False,
        ),
        (
            "promotion_rules",
            "fresh_714_paths_must_not_be_used_to_retune_alpha_if_they_are_to_remain_post_freeze_evidence",
            False,
        ),
        ("promotion_rules", "same_dimension_must_not_be_relabelled_as_same_basis_vector", False),
        ("promotion_rules", "same_support_must_not_be_relabelled_as_same_unit_response", False),
        ("promotion_rules", "visual_correspondence_requires_independent_visual_evidence", False),
        ("promotion_rules", "pde_validation_requires_cr001_held_out_acceptance", False),
        ("claim_states", "production_candidate_selected", True),
        ("claim_states", "velocity_changed_by_732", True),
        ("claim_states", "velocity_export_ready_promoted_by_this_preflight", True),
        ("claim_states", "visual_correspondence_verified", True),
        ("claim_states", "pde_validated", True),
        ("claim_states", "paper_exact", True),
        ("claim_states", "openai_field_identified", True),
    ],
)
def test_coordinate_scope_mutations_fail_closed(section, key, value):
    contract = deepcopy(load_contract())
    contract[section][key] = value
    with pytest.raises(GovernanceError):
        audit(contract=contract)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("pde_residual_max", 0.01),
        ("pde_residual_L2", 0.01),
        ("divergence_max", 0.001),
        ("divergence_L2", 0.001),
        ("validation_seed", 9177001),
        ("held_out_points", 1024),
        ("residual_defined_or_pointwise_free_force_allowed", True),
        ("amplitude_collapse_success_allowed", True),
        ("threshold_relaxation_allowed", True),
    ],
)
def test_cr001_mutations_fail_closed(key, value):
    contract = deepcopy(load_contract())
    contract["cr001_lock"][key] = value
    with pytest.raises(GovernanceError):
        audit(contract=contract)


def test_source_classification_cannot_drop_pending_unknown():
    contract = deepcopy(load_contract())
    contract["source_classification"]["items"] = [
        item for item in contract["source_classification"]["items"] if item["class"] != "pending_unknown"
    ]
    with pytest.raises(GovernanceError):
        audit(contract=contract)


def test_732_replay_permission_is_preserved_but_not_alpha_identity():
    contract = load_contract()
    rules = contract["promotion_rules"]
    scope = contract["representation_scope"]
    assert rules["732_pass_may_justify_one_same_dimension_nonlinear_replay"] is True
    assert rules["732_pass_may_claim_old_alpha_is_invariant_under_reshape"] is False
    assert scope["same_basis_dimension"] is True
    assert scope["spatial_basis_function_identical"] is False
    assert scope["old_701_alpha_transfer_validated"] is False
