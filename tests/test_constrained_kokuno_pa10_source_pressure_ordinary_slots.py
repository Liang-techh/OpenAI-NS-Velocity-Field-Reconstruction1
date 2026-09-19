from __future__ import annotations

from dataclasses import asdict
import math

import pytest

from openai_ns_reconstruction.kokuno_pa10_post_j_remainder_bridge import (
    R2_ORDINARY_TERMS,
)
from openai_ns_reconstruction.kokuno_pa10_source_pressure_ordinary_slots import (
    R2_PRESSURE_ORDINARY_TERMS,
    KokunoPA10SourcePressureOrdinarySlots,
)


def test_eta_coefficient_ball_matches_source_weight() -> None:
    calc = KokunoPA10SourcePressureOrdinarySlots()
    eta = calc.eta_coefficient_ball()
    expected = max(
        1.0 + calc.domain.enlarged_real_margin,
        4.0 * calc.domain.coefficient_rho,
    )
    assert eta.norm == pytest.approx(expected, rel=2e-16)
    assert eta.lipschitz == 0.0


def test_pressure_slots_fill_only_registered_r2_pressure_subset() -> None:
    calc = KokunoPA10SourcePressureOrdinarySlots()
    slots = calc.r2_pressure_numerator_slots()
    assert tuple(slots) == R2_PRESSURE_ORDINARY_TERMS
    assert set(slots).issubset(set(R2_ORDINARY_TERMS))
    assert len(slots) == 3
    assert len(R2_ORDINARY_TERMS) == 11
    for bound in slots.values():
        assert math.isfinite(bound.norm) and bound.norm > 0.0
        assert math.isfinite(bound.lipschitz) and bound.lipschitz > 0.0


def test_pressure_slots_reuse_source_pressure_ball_without_selected_local_values() -> None:
    calc = KokunoPA10SourcePressureOrdinarySlots()
    pressure = calc.source_pressure_ball()
    assert set(pressure) == {"integrand_g2_Phi2", "p", "p_eta", "Y_p_Y"}
    assert pressure["p_eta"].norm > pressure["p"].norm
    assert pressure["p_eta"].lipschitz > pressure["p"].lipschitz


def test_post_j1_pressure_slots_apply_bridge_exactly_once() -> None:
    calc = KokunoPA10SourcePressureOrdinarySlots()
    numerator = calc.r2_pressure_numerator_slots()
    after = calc.r2_pressure_after_j1()
    assert tuple(after) == R2_PRESSURE_ORDINARY_TERMS
    for name in R2_PRESSURE_ORDINARY_TERMS:
        direct = calc.bridge.ordinary_term_after_jnu(numerator[name], nu=1)
        assert asdict(after[name]) == asdict(direct)
        assert after[name].norm > numerator[name].norm
        assert after[name].lipschitz > numerator[name].lipschitz


def test_report_is_fail_closed_for_remaining_ordinary_slots_and_pde() -> None:
    calc = KokunoPA10SourcePressureOrdinarySlots()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert first["receipt_sha256"] == calc.sha256
    assert first["ledger_progress"]["R2_pressure_slots_filled_count"] == 3
    assert first["ledger_progress"]["R2_other_ordinary_slots_still_required"] == 8
    assert first["ledger_progress"]["R1_ordinary_slots_still_required"] == 11
    truth = first["truth_boundary"]
    assert truth["source_compatible_R2_pressure_ordinary_post_J1_slots_machine_bound"] is True
    assert truth["all_R2_ordinary_slots_source_bound"] is False
    assert truth["source_full_post_J1_R2_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert first["independent_audit_boundary"][
        "A4_source_axis_admission_required_before_downstream_promotion"
    ] is True
