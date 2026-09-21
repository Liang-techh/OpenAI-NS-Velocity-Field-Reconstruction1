from __future__ import annotations

import copy
import inspect

import pytest

from openai_ns_reconstruction import kokuno_a5_current_i2_composite_mean_ingest as a5


def _redigest(registration: dict[str, object]) -> None:
    payload = {
        key: copy.deepcopy(value)
        for key, value in registration.items()
        if key not in {"schema", "task", "digest"}
    }
    registration["digest"] = a5._sha256(payload)


def test_registration_binds_exact_current_i2_composite_and_mean_without_a4_laundering() -> None:
    registration = a5.build_registration()
    a5.validate_registration(registration)

    assert registration["parent_a5"]["head"] == "2d2cd726107a221b88b3bdfeaa0adcce83c3d57f"
    assert registration["agent1_current_i2"]["head"] == "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3"
    assert registration["agent2_current_i2_composite"]["head"] == "48d69e37e78e7f7f0e4e9936f28ff7974719288d"
    assert registration["agent2_current_i2_composite"]["consumes_agent1_head"] == registration["agent1_current_i2"]["head"]
    assert registration["agent2_current_i2_composite"]["identity_preserving_save_load"] is True
    assert registration["agent3_current_i2_mean"]["head"] == "a798056d7ddabe0bf4020082861807bc5ac1efbd"
    assert registration["agent3_current_i2_mean"]["consumes_agent2_head"] == registration["agent2_current_i2_composite"]["head"]
    assert registration["truth_boundary"]["current_i2_matching_leading_plus_oscillatory_composite_materialized"] is True
    assert registration["truth_boundary"]["current_i2_nonlinear_m0_mean_attribution_materialized"] is True

    assert registration["agent1_current_i2"]["public_float64_i2_delta_verified"] is False
    assert registration["agent2_current_i2_composite"]["public_float64_i2_delta_verified"] is False
    assert registration["agent3_current_i2_mean"]["public_float64_i2_delta_verified"] is False
    assert registration["truth_boundary"]["current_i2_public_float64_delta_verified"] is False

    assert registration["agent4_current_i1_only_evidence"]["head"] == "33626e55f301bbb5a28f75684de49d020bc6039d"
    assert registration["agent4_current_i1_only_evidence"]["audits_agent2_1062_current_i1_composite"] is True
    assert registration["agent4_current_i1_only_evidence"]["audits_agent2_1071_current_i2_composite"] is False
    assert registration["agent4_current_i1_only_evidence"]["audits_agent3_1073_current_i2_mean"] is False
    assert registration["agent4_current_i2_audit_status"] == {
        "dedicated_composite_audit_present": False,
        "independent_composite_audit_registered": False,
        "independent_composite_audit_admitted": False,
        "independent_nonlinear_mean_audit_present": False,
        "complete_ns_residual_evidence": False,
    }

    assert registration["agent1_current_i3_sibling"]["head"] == "5fb7b583062b4a991db86e39fdcdb231022b70c9"
    assert registration["agent1_current_i3_sibling"]["leading_velocity_through_i3_materialized"] is True
    assert registration["agent1_current_i3_sibling"]["source_positive_order_i3_profiles_materialized"] is False
    assert registration["agent1_current_i3_sibling"]["matching_agent2_i3_composite_materialized"] is False
    assert registration["agent1_current_i3_sibling"]["consumed_by_this_i2_checkpoint"] is False

    assert registration["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert registration["truth_boundary"]["complete_identity_bound_ns_defect_materialized"] is False
    assert registration["truth_boundary"]["heldout_normalized_full_ns_residual_assessed"] is False
    assert registration["truth_boundary"]["pde_validated"] is False


def test_i1_a4_i3_sibling_and_subulp_evidence_cannot_be_laundered() -> None:
    base = a5.build_registration()

    mutated = copy.deepcopy(base)
    mutated["agent4_current_i1_only_evidence"]["audits_agent2_1071_current_i2_composite"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="agent4_current_i1_only_evidence identity/truth drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["agent4_current_i2_audit_status"]["dedicated_composite_audit_present"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="current-I2 A4 audit status drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["agent1_current_i3_sibling"]["consumed_by_this_i2_checkpoint"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="agent1_current_i3_sibling identity/truth drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["agent1_current_i3_sibling"]["source_positive_order_i3_profiles_materialized"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="agent1_current_i3_sibling identity/truth drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["agent2_current_i2_composite"]["public_float64_i2_delta_verified"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="agent2_current_i2_composite identity/truth drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["truth_boundary"]["current_i2_public_float64_delta_verified"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="truth boundary drifted"):
        a5.validate_registration(mutated)


def test_project_gates_and_pde_state_fail_closed() -> None:
    base = a5.build_registration()

    mutated = copy.deepcopy(base)
    mutated["final_gates"]["normalized_momentum_sampled_max"] = 2.0e-3
    _redigest(mutated)
    with pytest.raises(ValueError, match="fixed project gates drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["final_gates"]["residual_defined_free_forcing_allowed"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="fixed project gates drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["truth_boundary"]["pde_validated"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="truth boundary drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["truth_boundary"]["current_i2_radial_force_materialized"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="truth boundary drifted"):
        a5.validate_registration(mutated)


def test_registration_roundtrip_is_deterministic_and_public_builder_has_no_tuning_knobs(tmp_path) -> None:
    path = tmp_path / "registration.json"
    first = a5.write_registration(path)
    loaded = a5.load_registration(path)
    second = a5.build_registration()
    assert first == loaded == second

    params = inspect.signature(a5.build_registration).parameters
    assert not params
    forbidden = {
        "threshold",
        "momentum_gate",
        "divergence_gate",
        "forcing",
        "residual",
        "correction",
        "pressure",
        "gain",
        "public_float64_delta",
    }
    assert forbidden.isdisjoint(params)
