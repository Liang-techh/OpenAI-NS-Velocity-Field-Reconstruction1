from __future__ import annotations

import copy
import inspect

import pytest

from openai_ns_reconstruction import kokuno_a5_current_i1_composite_mean_ingest as a5


def _redigest(registration: dict[str, object]) -> None:
    payload = {
        key: copy.deepcopy(value)
        for key, value in registration.items()
        if key not in {"schema", "task", "digest"}
    }
    registration["digest"] = a5._sha256(payload)


def test_registration_binds_exact_current_i1_composite_mean_and_matching_a4_audit() -> None:
    registration = a5.build_registration()
    a5.validate_registration(registration)

    assert registration["parent_a5"]["head"] == "ee6e1909e06edabbdf44e204fab1bc072c37d4b6"
    assert registration["agent2_current_i1_composite"]["head"] == "109527f520abb29bbe10372b0517eda44bcad0b6"
    assert registration["agent2_current_i1_composite"]["identity_preserving_save_load"] is True
    assert registration["agent3_current_i1_mean"]["head"] == "e6ac8bed0ef775232e18e81e87f6a7572dcba8fb"
    assert registration["agent3_current_i1_mean"]["consumes_agent2_head"] == registration["agent2_current_i1_composite"]["head"]
    assert registration["truth_boundary"]["current_i1_leading_plus_frozen_complete_curl_composite_materialized"] is True
    assert registration["truth_boundary"]["current_i1_nonlinear_m0_mean_attribution_materialized"] is True

    assert registration["agent4_current_i1_leading_audit"]["audits_agent1_1051_leading"] is True
    assert registration["agent4_current_i1_leading_audit"]["audits_agent2_1062_composite"] is False
    assert registration["agent4_current_i1_leading_audit"]["audits_agent3_1063_nonlinear_mean"] is False
    assert registration["agent4_current_i1_composite_audit"]["head"] == "33626e55f301bbb5a28f75684de49d020bc6039d"
    assert registration["agent4_current_i1_composite_audit"]["audited_agent2_head"] == registration["agent2_current_i1_composite"]["head"]
    assert registration["agent4_current_i1_composite_audit"]["scientific_admission"] is False
    assert registration["truth_boundary"]["agent4_matching_current_i1_composite_audit_present"] is True
    assert registration["truth_boundary"]["agent4_matching_current_i1_composite_audit_registered"] is True
    assert registration["truth_boundary"]["agent4_matching_current_i1_composite_audit_admitted"] is False
    assert registration["truth_boundary"]["agent4_matching_current_i1_nonlinear_mean_audit_registered"] is False

    assert registration["agent1_current_i2_sibling"]["leading_velocity_through_i2_materialized"] is True
    assert registration["agent1_current_i2_sibling"]["matching_agent2_i2_composite_materialized"] is False
    assert registration["agent1_current_i2_sibling"]["consumed_by_this_i1_checkpoint"] is False
    assert registration["truth_boundary"]["current_i2_matching_leading_plus_oscillatory_composite_materialized"] is False

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


def test_a4_evidence_and_i2_sibling_cannot_be_laundered() -> None:
    base = a5.build_registration()

    mutated = copy.deepcopy(base)
    mutated["agent4_current_i1_leading_audit"]["audits_agent2_1062_composite"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="agent4_current_i1_leading_audit identity/truth drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["agent4_current_i1_composite_audit"]["scientific_admission"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="agent4_current_i1_composite_audit identity/truth drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["truth_boundary"]["agent4_matching_current_i1_composite_audit_admitted"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="truth boundary drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["truth_boundary"]["agent4_matching_current_i1_nonlinear_mean_audit_registered"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="truth boundary drifted"):
        a5.validate_registration(mutated)

    mutated = copy.deepcopy(base)
    mutated["agent1_current_i2_sibling"]["consumed_by_this_i1_checkpoint"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="agent1_current_i2_sibling identity/truth drifted"):
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
    }
    assert forbidden.isdisjoint(params)
