from __future__ import annotations

import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_post_j_remainder_bridge import (
    R1_ORDINARY_TERMS,
    R2_ORDINARY_TERMS,
    KokunoPA10PostJRemainderBridge,
)
from openai_ns_reconstruction.kokuno_pa10_remainder_ball_bounds import BallFactorBound


def test_l_inverse_coefficient_ball_is_finite_and_covers_real_beta_zero() -> None:
    calc = KokunoPA10PostJRemainderBridge()
    bound = calc.l_inverse_coefficient_ball()
    assert math.isfinite(bound.norm) and bound.norm > 1.0
    assert bound.lipschitz == 0.0

    eta = np.linspace(-1.0, 1.0, 2001)
    L = 1.0 - 2.0 * calc.domain.h * eta * eta
    sampled = float(np.max(np.abs(1.0 / L)))
    assert bound.norm >= sampled


def test_term_ledgers_are_exact_and_fail_closed() -> None:
    calc = KokunoPA10PostJRemainderBridge()
    r1, r2 = calc.diagnostic_unit_ordinary_terms()
    assert tuple(r1) == R1_ORDINARY_TERMS
    assert tuple(r2) == R2_ORDINARY_TERMS

    broken = dict(r1)
    broken.pop(R1_ORDINARY_TERMS[-1])
    with pytest.raises(ValueError, match="term ledger mismatch"):
        calc.assemble(
            r1_ordinary_numerator_terms=broken,
            r2_ordinary_numerator_terms=r2,
        )

    extra = dict(r2)
    extra["not_a_source_term"] = BallFactorBound(1.0, 1.0)
    with pytest.raises(ValueError, match="term ledger mismatch"):
        calc.assemble(
            r1_ordinary_numerator_terms=r1,
            r2_ordinary_numerator_terms=extra,
        )


def test_ordinary_term_bridge_applies_L_inverse_and_Jnu_once() -> None:
    calc = KokunoPA10PostJRemainderBridge()
    unit = BallFactorBound(1.0, 2.0)
    j2 = calc.ordinary_term_after_jnu(unit, nu=2)
    j1 = calc.ordinary_term_after_jnu(unit, nu=1)
    assert math.isfinite(j2.norm) and j2.norm > 0.0
    assert math.isfinite(j1.norm) and j1.norm > j2.norm
    assert j2.lipschitz / j2.norm == pytest.approx(2.0, rel=3e-15)
    assert j1.lipschitz / j1.norm == pytest.approx(2.0, rel=3e-15)

    # At vanishing order zero the source factors are J2<=40 and J1<=80.
    assert j1.norm / j2.norm == pytest.approx(2.0, rel=3e-15)


def test_full_mixed_terms_include_L_inverse_without_selecting_Lambda() -> None:
    calc = KokunoPA10PostJRemainderBridge()
    full = calc.full_mixed_terms()
    old = calc.u_source.mixed_handoff()

    r1 = full["lambda_inv_Linv_d_detaAu_Y_Phi_Y_after_J2"]
    r2 = full["lambda_inv_Linv_d_detaAu_Y_u_Y_after_J1"]
    assert math.isfinite(r1.norm) and r1.norm > old["post_J2_d_detaAu_YPhiY"]["norm"]
    assert math.isfinite(r2.norm) and r2.norm > old["post_J1_d_detaAu_YuY"]["norm"]
    assert r1.lipschitz > 0.0
    assert r2.lipschitz > 0.0

    report = calc.report()
    assert report["source_compatible_fixed_inputs"][
        "Lambda_inverse_uniform_upper_for_Lambda_ge_1"
    ] == 1.0
    assert report["integration_boundary"]["full_mixed_terms_now_include_L_inverse"] is True


def test_full_post_j_assembly_runs_every_slot_but_stays_conditional() -> None:
    calc = KokunoPA10PostJRemainderBridge()
    r1, r2 = calc.diagnostic_unit_ordinary_terms()
    assembled = calc.assemble(
        r1_ordinary_numerator_terms=r1,
        r2_ordinary_numerator_terms=r2,
    )
    assert set(assembled["R1_ordinary_after_J2"]) == set(R1_ORDINARY_TERMS)
    assert set(assembled["R2_ordinary_after_J1"]) == set(R2_ORDINARY_TERMS)
    assert assembled["post_J2_R1"]["norm"] > 0.0
    assert assembled["post_J1_R2"]["norm"] > 0.0
    assert assembled["ordinary_inputs_source_certified"] is False

    pair = calc.diagnostic_pair_envelope_from_post_j(assembled)
    assert math.isfinite(pair["diagnostic_M_pair_max"])
    assert math.isfinite(pair["diagnostic_K_pair_max"])
    assert pair["diagnostic_M_pair_max"] > 0.0
    assert pair["diagnostic_K_pair_max"] > 0.0

    truth = calc.truth_boundary
    assert truth["termwise_post_J_remainder_bridge_executable"] is True
    assert truth["source_compatible_full_R1_mixed_post_J2_term_executable"] is True
    assert truth["source_compatible_full_R2_mixed_post_J1_term_executable"] is True
    assert truth["ordinary_term_inputs_are_source_certified"] is False
    assert truth["source_full_post_J2_R1_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_full_post_J1_R2_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["pde_validated"] is False


def test_report_is_deterministic_saveable_and_truth_safe(tmp_path) -> None:
    calc = KokunoPA10PostJRemainderBridge()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert len(first["receipt_sha256"]) == 64
    assert first["diagnostic_unit_ordinary_fixture"]["status"].startswith("arithmetic-only")
    assert first["heldout_ns_residual_assessed"] is False
    assert first["pde_validated"] is False

    target = tmp_path / "post_j_bridge.json"
    saved = calc.save_report(target)
    assert json.loads(target.read_text()) == saved
    assert saved["receipt_sha256"] == calc.sha256
