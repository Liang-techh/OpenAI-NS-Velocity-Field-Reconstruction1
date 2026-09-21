from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_current_i4_routing_i2_stress_ingest as a5


def test_registration_materializes_exact_split_frontiers_and_roundtrips(tmp_path):
    registration = a5.build_registration()
    assert registration["schema"] == a5.SCHEMA
    assert registration["agent2_current_i4_composite"]["head"] == "c40d8ddecd2971544a6e07dab093436b423cf326"
    assert registration["agent4_current_i4_composite_audit"]["head"] == "0ff663fcabaa52af984d273ba87356a72481c85c"
    assert registration["agent3_current_i2_stress"]["head"] == "cb9ee25d0c466d42ec7ca30ee89cc10517d4b5c4"
    assert registration["truth_boundary"]["velocity_frontier_stage"] == "current-I4"
    assert registration["truth_boundary"]["correction_frontier_stage"] == "current-I2"
    assert registration["truth_boundary"]["current_i2_stress_consumed_by_current_i4_candidate"] is False
    assert registration["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }

    path = tmp_path / "registration.json"
    written = a5.write_registration(path)
    loaded = a5.load_registration(path)
    assert loaded == written == registration
    assert json.loads(path.read_text(encoding="utf-8"))["digest"] == registration["digest"]


def test_i2_stress_cannot_be_transplanted_to_i4_candidate():
    registration = a5.build_registration()
    mutated = copy.deepcopy(registration)
    mutated["agent3_current_i2_stress"]["consumed_by_current_i4_candidate"] = True
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)


def test_i4_divergence_audit_cannot_be_laundered_into_stress_or_full_ns_evidence():
    registration = a5.build_registration()
    mutated = copy.deepcopy(registration)
    mutated["agent4_current_i4_composite_audit"]["audits_agent3_i2_stress"] = True
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(registration)
    mutated["agent4_current_i4_composite_audit"]["complete_ns_residual_evidence"] = True
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)


def test_radial_force_and_correction_authorization_remain_fail_closed():
    registration = a5.build_registration()
    assert registration["agent3_current_i2_stress"]["radial_force_partial_z_sigma1_materialized"] is False
    assert registration["agent3_current_i2_stress"]["authorized_as_complete_ns_correction_target"] is False
    assert registration["truth_boundary"]["complete_identity_bound_ns_defect_materialized"] is False
    assert registration["truth_boundary"]["current_cartesian_correction_velocity_materialized"] is False
    assert registration["truth_boundary"]["finite_correction_cycle_admitted"] is False
    assert registration["truth_boundary"]["heldout_normalized_full_ns_residual_assessed"] is False
    assert registration["truth_boundary"]["pde_validated"] is False


def test_source_role_and_global_boundaries_remain_false():
    registration = a5.build_registration()
    truth = registration["truth_boundary"]
    assert truth["source_positive_order_i3_correction_materialized"] is False
    assert truth["source_i4_mean_correction_materialized"] is False
    assert truth["terminal_global_leading_completion_materialized"] is False
    assert truth["global_leading_plus_oscillatory_velocity_materialized"] is False
    assert truth["matched_cartesian_pressure_materialized"] is False
    assert truth["preregistered_restricted_forcing_materialized"] is False


def test_project_gates_and_free_forcing_cannot_drift():
    registration = a5.build_registration()
    assert registration["final_gates"]["normalized_momentum_sampled_max"] == 1.0e-3
    assert registration["final_gates"]["normalized_momentum_volume_l2"] == 1.0e-3
    assert registration["final_gates"]["normalized_divergence_sampled_max"] == 1.0e-5
    assert registration["final_gates"]["normalized_divergence_volume_l2"] == 1.0e-5
    assert registration["final_gates"]["canonical_quadrature"] == [24, 48, 96]
    assert registration["final_gates"]["residual_defined_free_forcing_allowed"] is False

    mutated = copy.deepcopy(registration)
    mutated["final_gates"]["normalized_momentum_sampled_max"] = 2.0e-3
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(registration)
    mutated["final_gates"]["residual_defined_free_forcing_allowed"] = True
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)


def test_digest_and_exact_identity_mutations_fail_closed():
    registration = a5.build_registration()
    mutated = copy.deepcopy(registration)
    mutated["digest"] = "0" * 64
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(registration)
    mutated["agent2_current_i4_composite"]["head"] = "0" * 40
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)
