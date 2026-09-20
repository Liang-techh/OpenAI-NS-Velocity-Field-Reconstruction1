from __future__ import annotations

import math

from openai_ns_reconstruction.kokuno_pa10_source_pair_mk import (
    KokunoPA10SourcePairMK,
)


def test_pair_consumes_complete_agent1_post_j_self_envelopes() -> None:
    calc = KokunoPA10SourcePairMK()
    values = calc.full_post_j_inputs()
    r1 = values["post_J2_R1"]
    r2 = values["post_J1_R2"]

    assert r1 == calc.r1_source.full_r1_after_j2()
    assert r2 == calc.r2_source.full_r2_after_j1()
    assert r1.norm > 0.0 and r1.lipschitz > 0.0
    assert r2.norm > 0.0 and r2.lipschitz > 0.0
    assert math.isfinite(r1.norm) and math.isfinite(r1.lipschitz)
    assert math.isfinite(r2.norm) and math.isfinite(r2.lipschitz)


def test_pair_uses_existing_source_inverse_and_never_reapplies_j() -> None:
    calc = KokunoPA10SourcePairMK()
    source_inverse = calc.phi_source.phi_ball_certificate()[
        "inverse_one_plus_T_absolute_series_upper"
    ]
    pair = calc.pair_mk_self_certificate()
    report = calc.report()
    cert = report["operator_certificate"]

    assert pair["inverse_one_plus_T_absolute_series_upper"] == source_inverse
    assert cert["R1_input_is_already_post_J2"] is True
    assert cert["R2_input_is_already_post_J1"] is True
    assert cert["J2_reapplied_in_this_increment"] is False
    assert cert["J1_reapplied_in_this_increment"] is False
    assert cert["pair_factor_one_half_applied_exactly_once"] is True
    assert cert["diagnostic_unit_ordinary_fixture_used"] is False


def test_pair_mk_matches_existing_post_j_arithmetic_exactly() -> None:
    calc = KokunoPA10SourcePairMK()
    pair = calc.pair_mk_self_certificate()
    replay = calc.bridge_post_j_replay()

    assert calc.replay_matches_self_certificate() is True
    assert pair["inverse_one_plus_T_absolute_series_upper"] == replay[
        "inverse_one_plus_T_absolute_series_upper"
    ]
    assert pair["M_phi_component"] == replay["diagnostic_M_phi_component"]
    assert pair["M_u_component"] == replay["diagnostic_M_u_component"]
    assert pair["M_pair_max"] == replay["diagnostic_M_pair_max"]
    assert pair["K_phi_component"] == replay["diagnostic_K_phi_component"]
    assert pair["K_u_component"] == replay["diagnostic_K_u_component"]
    assert pair["K_pair_max"] == replay["diagnostic_K_pair_max"]


def test_pair_mk_is_finite_nontrivial_but_not_contraction_promotion() -> None:
    calc = KokunoPA10SourcePairMK()
    pair = calc.pair_mk_self_certificate()
    for name in (
        "inverse_one_plus_T_absolute_series_upper",
        "M_phi_component",
        "M_u_component",
        "M_pair_max",
        "K_phi_component",
        "K_u_component",
        "K_pair_max",
    ):
        assert math.isfinite(pair[name])
        assert pair[name] > 0.0

    report = calc.report()
    progress = report["progress"]
    assert progress["pair_M_agent1_self_certificate_materialized"] is True
    assert progress["pair_K_agent1_self_certificate_materialized"] is True
    assert progress["M_K_independent_admission_still_open"] is True
    assert progress["contraction_not_promoted"] is True
    assert progress["fixed_point_not_promoted"] is True


def test_report_is_deterministic_and_keeps_science_gates_closed() -> None:
    calc = KokunoPA10SourcePairMK()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert first["receipt_sha256"] == calc.sha256

    post_j = first["full_post_J_source_compatible_agent1_self_envelopes"]
    assert set(post_j) == {"post_J2_R1", "post_J1_R2"}
    assert "diagnostic_unit_ordinary_fixture" not in first

    truth = first["truth_boundary"]
    assert truth[
        "source_compatible_pair_M_agent1_self_certificate_machine_bound"
    ] is True
    assert truth[
        "source_compatible_pair_K_agent1_self_certificate_machine_bound"
    ] is True
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["source_operator_constant_M_independent_agent4_admission"] is False
    assert truth["source_operator_constant_K_independent_agent4_admission"] is False
    assert truth["source_contraction_invariant_ball_machine_verified"] is False
    assert truth["source_fixed_point_distance_machine_bound"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
