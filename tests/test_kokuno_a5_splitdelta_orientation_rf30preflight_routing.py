from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_splitdelta_orientation_rf30preflight_routing import (
    AGENT1_RELATIVE_SWIRL_SPLIT,
    AGENT2_AZIMUTHAL_FRAME_ORIENTATION,
    AGENT3_AUXILIARY_T2_HAAR_BRIDGE,
    FROZEN_GATES,
    TASK,
    build_registration,
    load_registration,
    registration_sha256,
    validate_registration,
    write_registration,
)


def _rehash(reg):
    reg["registration_sha256"] = registration_sha256(reg)
    return reg


def test_registration_harvests_split_delta_orientation_and_haar_bridge_without_promotion():
    reg = build_registration()
    frontiers = reg["frontiers"]
    truth = reg["truth_boundary"]

    assert reg["task"] == TASK
    assert reg["parent_a5"]["pr"] == 1169
    assert reg["parent_a5"]["head"] == "a035784fcba6d373f5dffc062872b632c6a26309"

    a1 = frontiers["relative_swirl_split_delta"]
    assert a1 == AGENT1_RELATIVE_SWIRL_SPLIT
    assert a1["pr"] == 1171
    assert a1["relative_swirl_profile_delta_channel_materialized"] is True
    assert a1["relative_swirl_cartesian_delta_channel_materialized"] is True
    assert a1["binary64_total_relative_swirl_sum_is_resolved"] is False
    assert a1["current_cartesian_relative_swirl_composed"] is False
    assert a1["scientifically_admitted"] is False

    a2 = frontiers["azimuthal_frame_oriented_multiharmonic"]
    assert a2 == AGENT2_AZIMUTHAL_FRAME_ORIENTATION
    assert a2["pr"] == 1170
    assert a2["bounded_azimuthal_frame_orientation_materialized"] is True
    assert a2["same_rigid_rotation_applied_to_vector_potential_and_complete_curl"] is True
    assert a2["orientation_classification"] == "repository_autonomous_candidate_design"
    assert a2["source_exact_orientation_recovered"] is False
    assert a2["self_contained_velocity_xyzt_provider"] is False
    assert a2["complete_ns_residual_assessed"] is False

    a3 = frontiers["auxiliary_t2_haar_rf30_bridge"]
    assert a3 == AGENT3_AUXILIARY_T2_HAAR_BRIDGE
    assert a3["pr"] == 1172
    assert a3["normalized_source_auxiliary_t2_haar_operator_executable"] is True
    assert a3["haar_measure_total_mass_one"] is True
    assert a3["covariance_formed_in_bridge_from_raw_wave_samples"] is True
    assert a3["preaveraged_covariance_input_exposed"] is False
    assert a3["repository_provider_blob_pinned"] is False
    assert a3["rf30_repository_candidate_state_authorized"] is False
    assert a3["rf30_defect_materialized"] is False

    assert frontiers["latest_self_contained_project_composite"]["pr"] == 1117
    assert truth["relative_swirl_cartesian_delta_channel_materialized"] is True
    assert truth["current_cartesian_relative_swirl_composed"] is False
    assert truth["current_i4_source_auxiliary_t2_haar_operator_executable"] is True
    assert truth["current_i4_source_auxiliary_t2_provider_blob_pinned"] is False
    assert truth["current_i4_rf30_repository_candidate_state_authorized"] is False
    assert truth["agent4_matching_agent1_1171_split_audit_present"] is False
    assert truth["agent4_matching_agent2_1170_orientation_audit_present"] is False
    assert truth["agent4_matching_agent3_1172_bridge_audit_present"] is False
    assert truth["pde_validated"] is False


def test_a4_remains_scoped_and_all_new_identities_are_firewalled():
    reg = build_registration()
    frontiers = reg["frontiers"]
    firewall = reg["identity_firewall"]

    a4 = frontiers["relative_swirl_compensator_validator"]
    assert a4["pr"] == 1160
    assert a4["audited_agent1_pr"] == 1154
    assert a4["scoped_gate_passed"] is None
    assert a4["scientifically_admitted"] is False

    assert firewall["agent1_1171_split_delta_not_unified_float64_cartesian_composition"] is True
    assert firewall["agent4_1160_not_evidence_for_agent1_1171_split_representation"] is True
    assert firewall["agent2_1170_orientation_is_repository_autonomous_not_source_exact"] is True
    assert firewall["agent2_1170_provider_family_not_self_contained_project_candidate"] is True
    assert firewall["agent3_1172_haar_operator_not_repository_candidate_without_pinned_provider"] is True
    assert firewall["agent3_1172_raw_wave_bridge_not_rf30_defect_or_correction"] is True
    assert firewall["agent4_1160_not_evidence_for_agent3_1172_bridge"] is True
    assert firewall["cross_identity_evidence_transfer_allowed"] is False


def test_core_state_and_frozen_gates_unchanged():
    reg = build_registration()
    assert reg["core_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert reg["frozen_gates"] == FROZEN_GATES
    assert FROZEN_GATES["st006_momentum_sampled_max"] == 0.1082289305112118
    assert FROZEN_GATES["st006_momentum_volume_l2"] == 0.10758432876230622
    assert FROZEN_GATES["normalized_momentum_sampled_max"] == 1e-3
    assert FROZEN_GATES["normalized_momentum_volume_l2"] == 1e-3
    assert FROZEN_GATES["normalized_divergence_sampled_max"] == 1e-5
    assert FROZEN_GATES["normalized_divergence_volume_l2"] == 1e-5
    assert FROZEN_GATES["canonical_quadrature"] == [24, 48, 96]
    assert FROZEN_GATES["residual_defined_free_forcing_allowed"] is False


def test_roundtrip_is_deterministic(tmp_path):
    path = tmp_path / "routing.json"
    first = write_registration(path)
    second = load_registration(path)
    assert first == second
    assert json.loads(path.read_text()) == first


def test_rejects_split_delta_laundered_into_unified_float64_field():
    reg = build_registration()
    a1 = reg["frontiers"]["relative_swirl_split_delta"]
    a1["binary64_total_relative_swirl_sum_is_resolved"] = True
    a1["current_cartesian_relative_swirl_composed"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_autonomous_orientation_laundered_into_source_exact_or_self_contained():
    reg = build_registration()
    a2 = reg["frontiers"]["azimuthal_frame_oriented_multiharmonic"]
    a2["source_exact_orientation_recovered"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))

    reg = build_registration()
    a2 = reg["frontiers"]["azimuthal_frame_oriented_multiharmonic"]
    a2["source_provider_self_contained"] = True
    a2["self_contained_velocity_xyzt_provider"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_unpinned_haar_bridge_promoted_to_rf30_candidate_or_defect():
    reg = build_registration()
    a3 = reg["frontiers"]["auxiliary_t2_haar_rf30_bridge"]
    a3["repository_provider_blob_pinned"] = True
    a3["rf30_repository_candidate_state_authorized"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))

    reg = build_registration()
    a3 = reg["frontiers"]["auxiliary_t2_haar_rf30_bridge"]
    a3["rf30_defect_materialized"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_a4_cross_identity_transfer():
    reg = build_registration()
    reg["frontiers"]["relative_swirl_compensator_validator"]["audited_agent1_pr"] = 1171
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_final_truth_or_gate_promotion():
    reg = build_registration()
    reg["truth_boundary"]["pde_validated"] = True
    reg["core_state"]["pde_validated"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))

    relaxed = copy.deepcopy(build_registration())
    relaxed["frozen_gates"]["normalized_momentum_volume_l2"] = 2e-3
    with pytest.raises(ValueError):
        validate_registration(_rehash(relaxed))
