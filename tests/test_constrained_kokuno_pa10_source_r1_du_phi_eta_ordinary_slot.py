from __future__ import annotations

import math

from openai_ns_reconstruction.kokuno_pa10_post_j_remainder_bridge import (
    R1_ORDINARY_TERMS,
)
from openai_ns_reconstruction.kokuno_pa10_source_r1_du_phi_eta_ordinary_slot import (
    R1_ALL_ORDINARY_POST_J2_SOURCE_BOUND_TERMS,
    R1_DU_PHI_ETA_ORDINARY_TERMS,
    R1_REMAINING_AFTER_DU_PHI_ETA_TERMS,
    KokunoPA10SourceR1DuPhiEtaOrdinarySlot,
)


def test_du_phi_eta_increment_is_exactly_final_registered_r1_slot() -> None:
    calc = KokunoPA10SourceR1DuPhiEtaOrdinarySlot()
    current = calc.r1_du_phi_eta_after_j2()
    assert tuple(current) == R1_DU_PHI_ETA_ORDINARY_TERMS
    assert tuple(current) == ("lambda_inv_d_u_times_Phi_eta",)
    assert set(current).issubset(set(R1_ORDINARY_TERMS))

    all_ordinary = calc.all_r1_ordinary_after_j2()
    assert tuple(all_ordinary) == R1_ALL_ORDINARY_POST_J2_SOURCE_BOUND_TERMS
    assert tuple(all_ordinary) == R1_ORDINARY_TERMS
    assert len(all_ordinary) == 11
    assert R1_REMAINING_AFTER_DU_PHI_ETA_TERMS == ()


def test_source_single_eta_product_post_j2_factor_is_executable_and_frozen() -> None:
    calc = KokunoPA10SourceR1DuPhiEtaOrdinarySlot()
    factor = calc.source_single_eta_after_j2_factor_upper()
    expected = 80.0 / calc.rho
    assert math.isfinite(factor)
    assert factor >= expected
    assert factor <= math.nextafter(expected, math.inf)


def test_du_phi_eta_uses_source_u_phi_balls_and_fixed_eta_only_multipliers() -> None:
    calc = KokunoPA10SourceR1DuPhiEtaOrdinarySlot()
    inputs = calc.operator_inputs()
    assert tuple(inputs) == ("d", "u", "Phi", "L_inverse")
    assert inputs["d"].norm > 0.0 and inputs["d"].lipschitz == 0.0
    assert inputs["u"].norm > 1.0 and inputs["u"].lipschitz > 0.0
    assert inputs["Phi"].norm > 1.0 and inputs["Phi"].lipschitz > 0.0
    assert inputs["L_inverse"].norm > 0.0
    assert inputs["L_inverse"].lipschitz == 0.0


def test_du_phi_eta_produces_finite_nontrivial_post_j2_bound() -> None:
    calc = KokunoPA10SourceR1DuPhiEtaOrdinarySlot()
    core = calc.j2_u_times_phi_eta()
    slot = calc.r1_du_phi_eta_after_j2()["lambda_inv_d_u_times_Phi_eta"]

    assert math.isfinite(core.norm) and core.norm > 0.0
    assert math.isfinite(core.lipschitz) and core.lipschitz > 0.0
    assert math.isfinite(slot.norm) and slot.norm > 0.0
    assert math.isfinite(slot.lipschitz) and slot.lipschitz > 0.0

    report = calc.report()
    certificate = report["operator_certificate"]
    assert certificate["post_J2_single_eta_factor_upper"] >= 80.0 / calc.rho
    assert certificate["single_undifferentiated_radial_factor_inside_J2_rule"] == "u"
    assert certificate["eta_only_multipliers_outside_J2_rule"] == [
        "L_inverse",
        "d",
    ]
    assert certificate["standalone_Phi_eta_norm_used"] is False
    assert certificate["J2_already_in_source_single_eta_product_factor"] is True
    assert certificate["ordinary_bridge_J2_reapplied"] is False
    assert certificate["Lambda_inverse_uniform_upper"] == 1.0


def test_all_ordinary_preserves_prior_ten_slots_exactly() -> None:
    calc = KokunoPA10SourceR1DuPhiEtaOrdinarySlot()
    prior = calc.parent.combined_r1_after_j2_subset()
    all_ordinary = calc.all_r1_ordinary_after_j2()
    assert len(prior) == 10
    for name, bound in prior.items():
        assert all_ordinary[name] == bound
    assert "lambda_inv_d_u_times_Phi_eta" not in prior
    assert "lambda_inv_d_u_times_Phi_eta" in all_ordinary


def test_full_r1_adds_existing_mixed_seam_once_and_is_nontrivial() -> None:
    calc = KokunoPA10SourceR1DuPhiEtaOrdinarySlot()
    all_ordinary = calc.all_r1_ordinary_after_j2()
    mixed = calc.bridge.full_mixed_terms()[
        "lambda_inv_Linv_d_detaAu_Y_Phi_Y_after_J2"
    ]
    full = calc.full_r1_after_j2()

    ordinary_norm_sum = sum(bound.norm for bound in all_ordinary.values())
    ordinary_lip_sum = sum(bound.lipschitz for bound in all_ordinary.values())
    assert full.norm >= ordinary_norm_sum
    assert full.norm >= mixed.norm
    assert full.lipschitz >= ordinary_lip_sum
    assert full.lipschitz >= mixed.lipschitz
    assert math.isfinite(full.norm) and full.norm > 0.0
    assert math.isfinite(full.lipschitz) and full.lipschitz > 0.0


def test_report_is_deterministic_closes_ordinary_r1_but_keeps_science_gates_closed() -> None:
    calc = KokunoPA10SourceR1DuPhiEtaOrdinarySlot()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert first["receipt_sha256"] == calc.sha256

    progress = first["ledger_progress"]
    assert progress["R1_post_J2_source_bound_ordinary_slots_count"] == 11
    assert progress["R1_total_ordinary_slots"] == 11
    assert progress["R1_ordinary_slots_still_required"] == 0
    assert progress["R1_remaining_ordinary_slots"] == []
    assert progress["R1_numerator_space_source_bound_count_remains_six"] is True
    assert progress["R1_mixed_slot_handled_upstream_by_A1_668"] is True
    assert progress["full_R1_agent1_self_certificate_materialized"] is True
    assert progress["full_R1_independent_admission_still_open"] is True
    assert progress["M_K_not_promoted"] is True

    truth = first["truth_boundary"]
    assert truth["source_single_eta_product_post_J2_rule_executable"] is True
    assert truth[
        "source_compatible_R1_du_Phi_eta_ordinary_post_J2_slot_machine_bound"
    ] is True
    assert truth[
        "source_compatible_all_R1_ordinary_post_J2_slots_machine_bound"
    ] is True
    assert truth[
        "source_compatible_full_post_J2_R1_radius_one_ball_bound_agent1_self_certificate"
    ] is True
    assert truth["standalone_Phi_eta_coefficient_norm_invented"] is False
    assert truth["ordinary_bridge_J2_reapplied"] is False
    assert truth["full_R1_independent_agent4_admission"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
