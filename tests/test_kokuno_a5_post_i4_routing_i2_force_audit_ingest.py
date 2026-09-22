from __future__ import annotations

import copy
import inspect
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_post_i4_routing_i2_force_audit_ingest as a5


def test_registration_binds_new_frontiers_without_transplanting_evidence(tmp_path):
    receipt = a5.build_registration()
    assert receipt["schema"] == a5.SCHEMA
    assert receipt["task"] == a5.TASK

    truth = receipt["truth_boundary"]
    assert truth["leading_frontier_stage"] == "post-I4/pre-pulse-leading-only"
    assert truth["latest_self_contained_composite_frontier_stage"] == "current-I4"
    assert truth["correction_frontier_stage"] == "current-I2-radial-force"
    assert truth["post_i4_pre_pulse_leading_materialized"] is True
    assert truth["matching_post_i4_agent2_project_cartesian_composite_materialized"] is False
    assert truth["agent2_source_localized_harmonic_adapter_present"] is True
    assert truth["agent2_source_localized_harmonic_adapter_self_contained_project_velocity"] is False
    assert truth["current_i2_radial_force_materialized"] is True
    assert truth["agent4_matching_current_i2_radial_force_audit_registered"] is True
    assert truth["agent4_matching_current_i2_radial_force_audit_admitted"] is False
    assert truth["current_i2_radial_force_consumed_by_current_i4_candidate"] is False
    assert truth["current_i2_radial_force_consumed_by_post_i4_candidate"] is False
    assert truth["pde_validated"] is False

    path = a5.save_registration(tmp_path / "routing.json")
    assert a5.load_registration(path) == receipt


def test_exact_lineage_and_independent_audit_are_frozen():
    receipt = a5.build_registration()
    assert receipt["parent_a5"]["head"] == "2886d888cafa80a30e9b21ffee5011e079d944a0"
    assert receipt["agent1_post_i4_leading"]["head"] == "bd85ed48323feb3b33b9c404a078630245dfe783"
    assert receipt["agent2_source_localized_harmonic"]["head"] == "4bc9408d10b2f3c636abd4f6ce33cefab391093d"
    assert receipt["agent3_current_i2_radial_force"]["head"] == "5e8512392df683375fabd132f35577e1b94773ff"
    audit = receipt["agent4_current_i2_radial_force_audit"]
    assert audit["head"] == "a01bc1f088b6b0d6556d81cf417c2b73bd7c99e8"
    assert audit["audited_agent3_head"] == receipt["agent3_current_i2_radial_force"]["head"]
    assert audit["implementation_distinct"] is True
    assert "five-point FD4" in audit["operator"]
    assert audit["scientific_admission"] is False


def test_fixed_project_gates_and_readiness_remain_unchanged():
    receipt = a5.build_registration()
    gates = receipt["final_gates"]
    assert gates["normalized_momentum_sampled_max"] == 1e-3
    assert gates["normalized_momentum_volume_l2"] == 1e-3
    assert gates["normalized_divergence_sampled_max"] == 1e-5
    assert gates["normalized_divergence_volume_l2"] == 1e-5
    assert gates["canonical_quadrature"] == [24, 48, 96]
    assert gates["residual_defined_free_forcing_allowed"] is False
    assert receipt["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_mutations_fail_closed():
    base = a5.build_registration()

    cases = []
    x = copy.deepcopy(base)
    x["truth_boundary"]["matching_post_i4_agent2_project_cartesian_composite_materialized"] = True
    cases.append(x)
    x = copy.deepcopy(base)
    x["agent2_source_localized_harmonic"]["self_contained_velocity_xyzt_provider"] = True
    cases.append(x)
    x = copy.deepcopy(base)
    x["agent3_current_i2_radial_force"]["consumed_by_post_i4_candidate"] = True
    cases.append(x)
    x = copy.deepcopy(base)
    x["agent3_current_i2_radial_force"]["authorized_as_complete_ns_correction_target"] = True
    cases.append(x)
    x = copy.deepcopy(base)
    x["agent4_current_i2_radial_force_audit"]["scientific_admission"] = True
    cases.append(x)
    x = copy.deepcopy(base)
    x["truth_boundary"]["pde_validated"] = True
    cases.append(x)
    x = copy.deepcopy(base)
    x["final_gates"]["normalized_momentum_sampled_max"] = 2e-3
    cases.append(x)

    for mutated in cases:
        # Recomputing the digest must not make a truth/provenance mutation legal.
        payload = {k: copy.deepcopy(v) for k, v in mutated.items() if k not in {"schema", "task", "digest"}}
        mutated["digest"] = a5._sha256(payload)
        with pytest.raises(ValueError):
            a5.validate_registration(mutated)


def test_registration_surface_exposes_no_scientific_tuning_knobs():
    assert set(inspect.signature(a5.build_registration).parameters) == set()
    assert set(inspect.signature(a5.validate_registration).parameters) == {"registration"}
    source = inspect.getsource(a5.build_registration)
    for forbidden in ("residual=", "forcing=", "pressure=", "viscosity=", "threshold=", "gain=", "damping="):
        assert forbidden not in source


def test_json_digest_detects_manual_corruption(tmp_path):
    path = a5.save_registration(tmp_path / "routing.json")
    raw = json.loads(path.read_text())
    raw["pipeline_position"]["stage_firewall"] = "launder evidence"
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError):
        a5.load_registration(path)
