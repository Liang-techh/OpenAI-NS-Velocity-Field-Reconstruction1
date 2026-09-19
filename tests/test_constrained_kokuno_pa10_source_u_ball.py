from __future__ import annotations

import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_source_u_ball import (
    KokunoPA10SourceUBallBounds,
)


def test_source_u0_center_is_vectorized_and_radially_exact() -> None:
    calc = KokunoPA10SourceUBallBounds()
    Y = np.array([[0.0], [0.7], [4.1]])
    eta = np.array([[-0.8, -0.13, 0.0, 0.31, 0.91]])
    values = calc.center_values(Y, eta)
    assert values["u_0"].shape == (3, 5)
    np.testing.assert_allclose(
        values["u_0"], Y * values["u_0_Y"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        values["Y_u_0_Y"], values["u_0"], rtol=0.0, atol=0.0
    )
    assert np.max(np.abs(values["u_0"][0])) == 0.0
    assert np.max(np.abs(values["u_0"][1:])) > 0.0


def test_source_u0_coefficient_certificate_covers_beta_zero_slice() -> None:
    calc = KokunoPA10SourceUBallBounds()
    cert = calc.u0_certificate()
    assert cert["cauchy_L_abs_lower"] > 0.0
    assert cert["cauchy_one_plus_z2_abs_lower"] > 0.0
    assert cert["rho_over_cauchy_radius"] == pytest.approx(0.1)
    assert cert["alpha1_cauchy_weight_sum_upper"] >= 1.0 / 0.9**2
    assert math.isfinite(cert["u0_coefficient_norm_upper"])
    assert cert["u0_coefficient_norm_upper"] > 0.0

    # For alpha=1,beta=0 the exact source weight is a_10=1/80.
    # A dense real slice is only a regression lower bound, not the certificate.
    eta = np.linspace(-1.0, 1.0, 1001)
    slope = calc.center_values(np.ones_like(eta), eta)["u_0_Y"]
    sampled_beta0_ratio = 80.0 * float(np.max(np.abs(slope)))
    assert cert["u0_coefficient_norm_upper"] >= sampled_beta0_ratio


def test_source_radius_one_u_ball_closes_660_unit_u_placeholder() -> None:
    calc = KokunoPA10SourceUBallBounds()
    cert = calc.u0_certificate()
    u = calc.u_ball()
    assert math.isfinite(u.norm) and u.norm >= cert["u0_coefficient_norm_upper"] + 1.0
    assert u.lipschitz == 1.0

    mixed = calc.mixed_handoff()
    for key in ("post_J2_d_detaAu_YPhiY", "post_J1_d_detaAu_YuY"):
        assert math.isfinite(mixed[key]["norm"]) and mixed[key]["norm"] > 0.0
        assert math.isfinite(mixed[key]["lipschitz"]) and mixed[key]["lipschitz"] > 0.0

    report = calc.report()
    assert report["integration_boundary"][
        "diagnostic_unit_u_fixture_replaced_for_these_two_mixed_terms"
    ] is True
    assert report["integration_boundary"]["termwise_post_J_R1_R2_bridge_still_required"] is True
    truth = report["truth_boundary"]
    assert truth["source_compatible_u_radius_one_ball_norm_machine_bound"] is True
    assert truth["source_compatible_post_J2_R1_mixed_numeric_ball_bound_executable"] is True
    assert truth["source_compatible_post_J1_R2_mixed_numeric_ball_bound_executable"] is True
    assert truth["source_R1_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_R2_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert report["heldout_ns_residual_assessed"] is False
    assert report["pde_validated"] is False


def test_source_u_ball_report_is_deterministic_and_saveable(tmp_path) -> None:
    calc = KokunoPA10SourceUBallBounds()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert len(first["receipt_sha256"]) == 64
    target = tmp_path / "u_ball.json"
    saved = calc.save_report(target)
    assert json.loads(target.read_text()) == saved
    assert saved["receipt_sha256"] == calc.sha256


def test_source_u0_rejects_outside_physical_inner_rectangle() -> None:
    calc = KokunoPA10SourceUBallBounds()
    with pytest.raises(ValueError):
        calc.center_values(-1e-6, 0.0)
    with pytest.raises(ValueError):
        calc.center_values(4.100001, 0.0)
    with pytest.raises(ValueError):
        calc.center_values(1.0, 1.0 + calc.domain.enlarged_real_margin + 1e-6)
