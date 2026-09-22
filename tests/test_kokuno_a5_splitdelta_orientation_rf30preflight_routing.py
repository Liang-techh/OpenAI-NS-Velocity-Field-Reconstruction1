from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_splitdelta_orientation_rf30preflight_routing import (
    AGENT1_RELATIVE_SWIRL_SPLIT,
    AGENT2_AZIMUTHAL_FRAME_ORIENTATION,
    AGENT3_AUXILIARY_T2_HAAR_BRIDGE,
    AGENT4_RELATIVE_SWIRL_SPLIT_AUDIT,
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


def test_registration_harvests_fresh_a1_a2_a3_a4_without_promotion():
    reg = build_registration()
    f = reg["frontiers"]
    t = reg["truth_boundary"]
    assert reg["task"] == TASK
    assert reg["parent_a5"]["pr"] == 1169

    a1 = f["relative_swirl_split_delta"]
    assert a1 == AGENT1_RELATIVE_SWIRL_SPLIT
    assert a1["pr"] == 1171
    assert a1["relative_swirl_cartesian_delta_channel_materialized"] is True
    assert a1["binary64_total_relative_swirl_sum_is_resolved"] is False
    assert a1["current_cartesian_relative_swirl_composed"] is False

    a2 = f["azimuthal_frame_oriented_multiharmonic"]
    assert a2 == AGENT2_AZIMUTHAL_FRAME_ORIENTATION
    assert a2["pr"] == 1170
    assert a2["bounded_azimuthal_frame_orientation_materialized"] is True
    assert a2["source_exact_orientation_recovered"] is False
    assert a2["self_contained_velocity_xyzt_provider"] is False

    a3 = f["auxiliary_t2_haar_rf30_bridge"]
    assert a3 == AGENT3_AUXILIARY_T2_HAAR_BRIDGE
    assert a3["pr"] == 1172
    assert a3["normalized_source_auxiliary_t2_haar_operator_executable"] is True
    assert a3["repository_provider_blob_pinned"] is False
    assert a3["rf30_repository_candidate_state_authorized"] is False
    assert a3["rf30_defect_materialized"] is False

    a4 = f["relative_swirl_split_validator"]
    assert a4 == AGENT4_RELATIVE_SWIRL_SPLIT_AUDIT
    assert a4["pr"] == 1174
    assert a4["audited_agent1_pr"] == 1171
    assert a4["scientific_path_consumes_public_velocity_split_only_after_save_load"] is True
    assert a4["split_relative_swirl_cartesian_channel_independently_audited"] is True
    assert a4["scoped_gate_passed"] is None
    assert a4["complete_ns_residual_assessed"] is False
    assert a4["pde_validated"] is False

    assert f["latest_self_contained_project_composite"]["pr"] == 1117
    assert t["agent4_matching_agent1_1171_split_audit_present"] is True
    assert t["agent4_matching_agent1_1171_split_scoped_gate_passed"] is None
    assert t["current_i4_source_auxiliary_t2_haar_operator_executable"] is True
    assert t["current_i4_source_auxiliary_t2_provider_blob_pinned"] is False
    assert t["pde_validated"] is False


def test_identity_firewalls_keep_scoped_evidence_scoped():
    reg = build_registration()
    fw = reg["identity_firewall"]
    assert fw["agent1_1171_split_delta_not_unified_float64_cartesian_composition"] is True
    assert fw["agent4_1174_evidence_bound_to_exact_agent1_1171_split_identity"] is True
    assert fw["agent4_1174_relative_fd_ladder_not_replacement_for_canonical_project_fd"] is True
    assert fw["agent4_1174_split_divergence_evidence_not_complete_ns_validation"] is True
    assert fw["agent2_1170_orientation_is_repository_autonomous_not_source_exact"] is True
    assert fw["agent3_1172_haar_operator_not_repository_candidate_without_pinned_provider"] is True
    assert fw["agent4_1174_not_evidence_for_agent2_1170_or_agent3_1172"] is True
    assert fw["cross_identity_evidence_transfer_allowed"] is False


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
    assert load_registration(path) == first
    assert json.loads(path.read_text()) == first


def test_rejects_split_delta_laundered_into_unified_field():
    reg = build_registration()
    reg["frontiers"]["relative_swirl_split_delta"]["current_cartesian_relative_swirl_composed"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_orientation_promotion():
    reg = build_registration()
    reg["frontiers"]["azimuthal_frame_oriented_multiharmonic"]["source_exact_orientation_recovered"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_unpinned_haar_bridge_promoted_to_rf30():
    reg = build_registration()
    a3 = reg["frontiers"]["auxiliary_t2_haar_rf30_bridge"]
    a3["repository_provider_blob_pinned"] = True
    a3["rf30_repository_candidate_state_authorized"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_unresolved_a4_gate_or_full_ns_promotion():
    reg = build_registration()
    a4 = reg["frontiers"]["relative_swirl_split_validator"]
    a4["scoped_gate_passed"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))

    reg = build_registration()
    a4 = reg["frontiers"]["relative_swirl_split_validator"]
    a4["complete_ns_residual_assessed"] = True
    a4["pde_validated"] = True
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
