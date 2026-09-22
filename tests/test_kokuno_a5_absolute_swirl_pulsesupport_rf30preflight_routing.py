from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_absolute_swirl_pulsesupport_rf30preflight_routing import (
    AGENT1_ABSOLUTE_RELATIVE_SWIRL_TARGET,
    AGENT2_PULSE_COORDINATE_SUPPORT,
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


def test_registration_harvests_absolute_target_and_pulse_support_without_promotion():
    reg = build_registration()
    frontiers = reg["frontiers"]
    truth = reg["truth_boundary"]

    assert reg["task"] == TASK
    assert reg["parent_a5"]["pr"] == 1162
    assert reg["parent_a5"]["head"] == "ab8501a3f2e89e2f2e6f47c2ff67b53ead8eaa8f"

    a1 = frontiers["absolute_relative_swirl_target"]
    assert a1 == AGENT1_ABSOLUTE_RELATIVE_SWIRL_TARGET
    assert a1["pr"] == 1168
    assert a1["exact_compare_ahead"] == 3
    assert a1["exact_compare_behind"] == 0
    assert a1["open_pr_mergeable_at_registration"] is False
    assert a1["absolute_current_r_I_at_flattening_endpoint_materialized"] is True
    assert a1["current_relative_swirl_total_target_materialized"] is True
    assert a1["imported_base_absolute_angular_target_materialized"] is True
    assert a1["current_cartesian_relative_swirl_composed"] is False
    assert a1["scientifically_admitted"] is False

    a2 = frontiers["pulse_coordinate_support_multiharmonic"]
    assert a2 == AGENT2_PULSE_COORDINATE_SUPPORT
    assert a2["pr"] == 1164
    assert a2["provider_certified_pulse_coordinate_support_materialized"] is True
    assert a2["pulse_coordinate_mapping_identity_bound"] is True
    assert a2["source_exact_pulse_coordinate_mapping_recovered"] is False
    assert a2["source_pulse_coordinate_identified_with_chart_T"] is False
    assert a2["source_pulse_coordinate_identified_with_physical_t"] is False
    assert a2["physical_time_support_certificate_materialized"] is False
    assert a2["source_provider_self_contained"] is False
    assert a2["self_contained_velocity_xyzt_provider"] is False

    assert frontiers["latest_cartesian_leading"]["pr"] == 1148
    assert frontiers["latest_self_contained_project_composite"]["pr"] == 1117
    assert truth["absolute_current_r_I_at_flattening_endpoint_materialized"] is True
    assert truth["current_relative_swirl_total_target_materialized"] is True
    assert truth["current_cartesian_relative_swirl_composed"] is False
    assert truth["agent4_1160_audits_agent1_1168_absolute_target"] is False
    assert truth["pde_validated"] is False


def test_a3_a4_are_retained_with_identity_firewalls():
    reg = build_registration()
    frontiers = reg["frontiers"]
    firewall = reg["identity_firewall"]

    a3 = frontiers["current_i4_rf30_fixedq_preflight"]
    assert a3["pr"] == 1161
    assert a3["normalized_physical_azimuthal_mean_used"] is True
    assert a3["source_auxiliary_t2_provider_available"] is False
    assert a3["normalized_source_auxiliary_t2_haar_mean_used"] is False
    assert a3["rf30_repository_candidate_state_authorized"] is False
    assert a3["current_i4_rf30_defect_materialized"] is False

    a4 = frontiers["relative_swirl_compensator_validator"]
    assert a4["pr"] == 1160
    assert a4["audited_agent1_pr"] == 1154
    assert a4["audits_agent1_1159_current_correction"] is False
    assert a4["scoped_gate_passed"] is None
    assert a4["scientifically_admitted"] is False

    assert firewall["agent1_1168_absolute_target_not_cartesian_relative_swirl_composition"] is True
    assert firewall["agent1_1168_nonmergeable_open_pr_not_scientific_admission"] is True
    assert firewall["agent4_1160_not_evidence_for_agent1_1168_absolute_target"] is True
    assert firewall["agent2_1164_pulse_coordinate_not_chart_T_or_physical_t"] is True
    assert firewall["agent3_1161_physical_theta_mean_not_source_auxiliary_t2_haar_mean"] is True
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


def test_rejects_mergeability_or_scientific_admission_laundering():
    reg = build_registration()
    reg["frontiers"]["absolute_relative_swirl_target"]["open_pr_mergeable_at_registration"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))

    reg = build_registration()
    reg["frontiers"]["absolute_relative_swirl_target"]["scientifically_admitted"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_absolute_target_laundered_into_cartesian_composition():
    reg = build_registration()
    reg["frontiers"]["absolute_relative_swirl_target"]["current_cartesian_relative_swirl_composed"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_pulse_coordinate_laundered_into_physical_time_or_self_contained_candidate():
    reg = build_registration()
    a2 = reg["frontiers"]["pulse_coordinate_support_multiharmonic"]
    a2["source_pulse_coordinate_identified_with_physical_t"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))

    reg = build_registration()
    a2 = reg["frontiers"]["pulse_coordinate_support_multiharmonic"]
    a2["source_provider_self_contained"] = True
    a2["self_contained_velocity_xyzt_provider"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_a4_cross_identity_transfer_or_rf30_theta_substitution():
    reg = build_registration()
    reg["frontiers"]["relative_swirl_compensator_validator"]["audited_agent1_pr"] = 1168
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))

    reg = build_registration()
    pre = reg["frontiers"]["current_i4_rf30_fixedq_preflight"]
    pre["source_auxiliary_t2_provider_available"] = True
    pre["normalized_source_auxiliary_t2_haar_mean_used"] = True
    pre["rf30_repository_candidate_state_authorized"] = True
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
