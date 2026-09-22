from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_xi11_composite_i4_force_routing as a5


def test_registration_is_deterministic_and_roundtrips(tmp_path):
    first = a5.build_registration()
    second = a5.build_registration()
    assert first == second
    assert len(first["digest"]) == 64

    path = tmp_path / "routing.json"
    saved = a5.save_registration(path)
    loaded = a5.load_registration(path)
    assert saved == loaded == first
    assert json.loads(path.read_text()) == first


def test_latest_frontiers_are_identity_separated():
    reg = a5.build_registration()
    truth = reg["truth_boundary"]

    assert truth["matching_xi11_agent2_project_composite_materialized"] is True
    assert truth["matching_xi11_agent2_project_composite_consumes_exact_agent1_1107"] is True
    assert truth["matching_xi11_agent2_project_composite_consumes_agent1_1116"] is False
    assert truth["current_i4_radial_force_materialized"] is True

    assert reg["agent2_xi11_logx_composite"]["head"] == "27741d9c0a27262f7fabf61eebaa2fbd507e9f03"
    assert reg["agent3_current_i4_radial_force"]["head"] == "922f7aa10460ded44af212d313eceff33a2ac647"
    assert reg["agent2_xi11_logx_composite"]["consumed_agent1_head"] == "45da043dd2b4cd067f005a72c7e21fd0f2bcf309"
    assert reg["agent3_current_i4_radial_force"]["agent2_composite_head"] == "c40d8ddecd2971544a6e07dab093436b423cf326"


def test_agent1_1116_is_registered_as_non_consumed_algebra_sibling():
    reg = a5.build_registration()
    a1 = reg["agent1_pulse_end_algebra"]
    assert a1["corrected_two_row_c1_c2_system_materialized"] is True
    assert a1["source_exact_amplitude_root_materialized"] is False
    assert a1["current_lineage_J_materialized"] is False
    assert a1["cartesian_c1_c2_composition_materialized"] is False
    assert a1["consumed_by_latest_agent2_composite"] is False


def test_a4_1110_is_not_laundered_onto_new_frontiers():
    reg = a5.build_registration()
    a4 = reg["agent4_finite_prefix_divergence_audit"]
    truth = reg["truth_boundary"]

    assert a4["audited_agent2_pr"] == 1108
    assert a4["audited_agent2_head"] == "9f221dd57ef4b5e2d2e80e24c8cf531e991a5387"
    assert a4["audited_agent2_1117"] is False
    assert a4["audited_agent3_1118_radial_force"] is False
    assert truth["agent4_matching_xi11_composite_audit_present"] is False
    assert truth["agent4_matching_current_i4_radial_force_audit_present"] is False
    assert truth["agent4_1110_finite_prefix_audit_transferred_to_xi11_identity"] is False


def test_project_scientific_gates_and_readiness_are_inherited_unchanged():
    reg = a5.build_registration()
    assert reg["final_gates"] == a5.parent_a5.FINAL_GATES
    assert reg["st006_baseline"] == a5.parent_a5.ST006_BASELINE
    assert reg["readiness"] == a5.parent_a5.READINESS
    assert reg["st006_baseline"]["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert reg["st006_baseline"]["momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert reg["truth_boundary"]["pde_validated"] is False
    assert reg["truth_boundary"]["heldout_normalized_full_ns_residual_assessed"] is False


@pytest.mark.parametrize(
    "key",
    [
        "agent4_matching_xi11_composite_audit_present",
        "agent4_matching_current_i4_radial_force_audit_present",
        "current_i4_radial_force_authorized_as_complete_ns_correction_target",
        "global_leading_plus_oscillatory_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "preregistered_restricted_forcing_materialized",
        "complete_identity_bound_ns_defect_materialized",
        "finite_correction_cycle_admitted",
        "heldout_normalized_full_ns_residual_assessed",
        "scientific_admission",
        "pde_validated",
    ],
)
def test_unsupported_truth_promotions_fail_closed(key):
    reg = copy.deepcopy(a5.build_registration())
    reg["truth_boundary"][key] = True
    with pytest.raises(ValueError):
        a5.validate_registration(reg)


@pytest.mark.parametrize(
    "key",
    [
        "agent4_1110_finite_prefix_audit_transferred_to_xi11_identity",
        "current_i4_force_evidence_transferred_to_xi11_composite",
        "xi11_composite_evidence_transferred_to_current_i4_force",
    ],
)
def test_cross_identity_evidence_transfer_fails_closed(key):
    reg = copy.deepcopy(a5.build_registration())
    reg["truth_boundary"][key] = True
    with pytest.raises(ValueError):
        a5.validate_registration(reg)


def test_exact_provenance_mutation_fails_closed():
    reg = copy.deepcopy(a5.build_registration())
    reg["agent2_xi11_logx_composite"]["head"] = "0" * 40
    with pytest.raises(ValueError):
        a5.validate_registration(reg)


def test_radial_force_cannot_be_authorized_by_scoped_operator_receipt():
    reg = copy.deepcopy(a5.build_registration())
    reg["agent3_current_i4_radial_force"]["authorized_as_complete_ns_correction_target"] = True
    with pytest.raises(ValueError):
        a5.validate_registration(reg)
