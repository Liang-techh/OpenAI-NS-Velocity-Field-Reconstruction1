from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st052m_tip_holdout_protocol_scope import (
    GovernanceError,
    audit_contract,
    audit_live,
    load_contract,
)


def _reject(mutator) -> None:
    contract = deepcopy(load_contract())
    mutator(contract)
    with pytest.raises(GovernanceError):
        audit_contract(contract)


def test_live_tip_holdout_protocol_scope_audit_passes() -> None:
    result = audit_live()
    assert result["status"] == "PASS"
    assert result["audited_pr"] == 714
    assert result["fresh_path_count"] == 64
    assert result["exact_seed_overlap_required"] == 0
    assert result["post_freeze_directional_generalization_eligible_if_upstream_gate_passes"] is True
    assert result["blind_global_visual_validation"] is False
    assert result["cr001_pde_acceptance_sample"] is False
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False


def test_rejects_fresh_data_feedback_or_seed_overlap() -> None:
    _reject(lambda c: c["selection_evaluation_separation"].__setitem__("fresh_results_feed_back_into_alpha", True))
    _reject(lambda c: c["selection_evaluation_separation"].__setitem__("exact_seed_overlap_with_development_required", 1))
    _reject(lambda c: c["selection_evaluation_separation"].__setitem__("candidate_parameter_frozen_before_fresh_path_evaluation", False))


def test_rejects_disjoint_holdout_laundered_into_blind_global_validation() -> None:
    _reject(lambda c: c["protocol_independence_scope"].__setitem__("protocol_is_blind_to_prior_band_diagnosis", True))
    _reject(lambda c: c["protocol_independence_scope"].__setitem__("protocol_is_random_statistical_holdout", True))
    _reject(lambda c: c["protocol_independence_scope"].__setitem__("protocol_is_global_visual_validation", True))
    _reject(lambda c: c["protocol_independence_scope"].__setitem__("fresh_seed_disjointness_alone_establishes_visual_correspondence", True))


def test_rejects_path_holdout_laundered_into_cr001_pde_acceptance() -> None:
    _reject(lambda c: c["protocol_independence_scope"].__setitem__("protocol_is_cr001_pde_acceptance_sample", True))
    _reject(lambda c: c["protocol_independence_scope"].__setitem__("fresh_seed_disjointness_alone_establishes_pde_validation", True))
    _reject(lambda c: c["claim_states"].__setitem__("pde_validated", True))


def test_rejects_unreviewed_upstream_result_or_delivery_promotion() -> None:
    _reject(lambda c: c["claim_states"].__setitem__("fresh_post_freeze_path_directional_generalization_verified", True))
    _reject(lambda c: c["claim_states"].__setitem__("visual_correspondence_verified", True))
    _reject(lambda c: c["claim_states"].__setitem__("openai_field_identified", True))
    _reject(lambda c: c["candidate_identity_scope"].__setitem__("velocity_export_ready_promoted_by_714", True))
    _reject(lambda c: c["candidate_identity_scope"].__setitem__("standalone_serialized_714_candidate_identity_created", True))


def test_rejects_source_class_laundering() -> None:
    def promote_holdout_protocol(c):
        for item in c["source_classification"]["items"]:
            if item["item"] == "fresh_radii_angles_z_bands_solver_tolerances_and_directional_gates":
                item["class"] = "public_source_fact"

    def promote_hidden_data(c):
        for item in c["source_classification"]["items"]:
            if item["item"] == "hidden_openai_path_seeds_camera_frame_time_mapping_coefficients_and_numeric_field":
                item["class"] = "public_source_fact"

    _reject(promote_holdout_protocol)
    _reject(promote_hidden_data)


def test_rejects_cr001_threshold_or_force_shortcut_drift() -> None:
    _reject(lambda c: c["cr001_lock"].__setitem__("pde_residual_max", 0.01))
    _reject(lambda c: c["cr001_lock"].__setitem__("divergence_max", 1e-4))
    _reject(lambda c: c["cr001_lock"].__setitem__("residual_defined_or_pointwise_free_force_allowed", True))
    _reject(lambda c: c["cr001_lock"].__setitem__("amplitude_collapse_success_allowed", True))
    _reject(lambda c: c["cr001_lock"].__setitem__("threshold_relaxation_allowed", True))
