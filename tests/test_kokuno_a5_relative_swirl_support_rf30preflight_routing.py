from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_relative_swirl_support_rf30preflight_routing import (
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


def test_registration_binds_fresh_frontiers_without_promotion():
    reg = build_registration()
    f = reg["frontiers"]
    truth = reg["truth_boundary"]

    assert reg["task"] == TASK
    assert reg["parent_a5"]["pr"] == 1153
    assert reg["parent_a5"]["head"] == "eb15089a312b9ea1a3c4632145d501d244ec1324"

    assert f["latest_cartesian_leading"]["pr"] == 1148
    assert f["current_angular_target_correction"]["pr"] == 1159
    assert f["current_angular_target_correction"]["current_lineage_angular_target_correction_materialized"] is True
    assert f["current_angular_target_correction"]["imported_base_absolute_angular_target_materialized"] is False
    assert f["current_angular_target_correction"]["current_cartesian_relative_swirl_composed"] is False

    assert f["support_certified_multiharmonic"]["pr"] == 1158
    assert f["support_certified_multiharmonic"]["provider_certified_outer_spatial_support_materialized"] is True
    assert f["support_certified_multiharmonic"]["source_provider_self_contained"] is False
    assert f["latest_self_contained_project_composite"]["pr"] == 1117

    pre = f["current_i4_rf30_fixedq_preflight"]
    assert pre["pr"] == 1161
    assert pre["fixed_q_coordinate_map_materialized"] is True
    assert pre["normalized_physical_azimuthal_mean_used"] is True
    assert pre["source_auxiliary_t2_provider_available"] is False
    assert pre["normalized_source_auxiliary_t2_haar_mean_used"] is False
    assert pre["rf30_repository_candidate_state_authorized"] is False
    assert pre["current_i4_rf30_defect_materialized"] is False

    assert f["holdprefix_leading_validator"]["pr"] == 1152
    assert f["holdprefix_leading_validator"]["audited_agent1_pr"] == 1148
    assert f["holdprefix_leading_validator"]["scoped_gate_passed"] is None
    assert truth["agent4_1152_audits_agent1_1159"] is False
    assert truth["pde_validated"] is False


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


def test_rejects_laundered_absolute_angular_target():
    reg = build_registration()
    reg["frontiers"]["current_angular_target_correction"]["imported_base_absolute_angular_target_materialized"] = True
    reg["truth_boundary"]["imported_base_absolute_angular_target_materialized"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_provider_family_laundered_as_self_contained():
    reg = build_registration()
    reg["frontiers"]["support_certified_multiharmonic"]["source_provider_self_contained"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_physical_theta_mean_laundered_as_source_haar():
    reg = build_registration()
    pre = reg["frontiers"]["current_i4_rf30_fixedq_preflight"]
    pre["source_auxiliary_t2_provider_available"] = True
    pre["normalized_source_auxiliary_t2_haar_mean_used"] = True
    pre["rf30_repository_candidate_state_authorized"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_cross_identity_a4_transfer():
    reg = build_registration()
    reg["frontiers"]["holdprefix_leading_validator"]["audited_agent1_pr"] = 1159
    with pytest.raises(ValueError):
        validate_registration(_rehash(reg))


def test_rejects_pde_or_gate_promotion():
    reg = build_registration()
    promoted = copy.deepcopy(reg)
    promoted["truth_boundary"]["pde_validated"] = True
    promoted["core_state"]["pde_validated"] = True
    with pytest.raises(ValueError):
        validate_registration(_rehash(promoted))

    relaxed = copy.deepcopy(reg)
    relaxed["frozen_gates"]["normalized_momentum_sampled_max"] = 2e-3
    with pytest.raises(ValueError):
        validate_registration(_rehash(relaxed))
