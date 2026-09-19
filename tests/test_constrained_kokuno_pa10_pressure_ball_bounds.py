import hashlib
import json
import math

import pytest

from openai_ns_reconstruction.kokuno_pa10_pressure_ball_bounds import (
    KokunoPA10PressureBallBounds,
)
from openai_ns_reconstruction.kokuno_pa10_remainder_ball_bounds import (
    BallFactorBound,
    KokunoPA10RemainderBallBounds,
    KokunoPA10RemainderPrimitiveBounds,
)


def test_source_operator_factors_follow_displayed_weight_ratios():
    bounds = KokunoPA10PressureBallBounds()
    assert bounds.radial_integral_factor() == 80.0
    assert bounds.multiply_Y_factor() == 80.0
    assert bounds.eta_radial_integral_factor(0.5) >= 160.0
    assert math.isclose(bounds.eta_radial_integral_factor(0.5), 160.0, rel_tol=2e-16)
    assert bounds.product_constant_upper >= bounds.operators.product_constant


def test_pressure_map_propagates_norm_and_one_factor_lipschitz():
    bounds = KokunoPA10PressureBallBounds()
    g_norm = 1.2
    phi = BallFactorBound(norm=2.0, lipschitz=0.5)
    rho = 0.25
    result = bounds.pressure_bounds(rho=rho, g_norm=g_norm, Phi=phi)

    product_constant = bounds.product_constant_upper
    raw_integrand_norm = product_constant**3 * g_norm**2 * phi.norm**2
    raw_integrand_lip = (
        product_constant**3 * g_norm**2 * 2.0 * phi.norm * phi.lipschitz
    )
    integrand = result["integrand_g2_Phi2"]
    assert integrand.norm >= raw_integrand_norm
    assert math.isclose(integrand.norm, raw_integrand_norm, rel_tol=2e-15)
    assert integrand.lipschitz >= raw_integrand_lip
    assert math.isclose(integrand.lipschitz, raw_integrand_lip, rel_tol=2e-15)

    assert result["p"].norm >= 80.0 * integrand.norm
    assert result["Y_p_Y"].norm >= 80.0 * integrand.norm
    assert result["p_eta"].norm >= (80.0 / rho) * integrand.norm
    assert result["p"].lipschitz > 0.0
    assert result["p_eta"].lipschitz > result["p"].lipschitz


def test_pressure_injection_changes_only_pressure_primitives():
    bounds = KokunoPA10PressureBallBounds()
    base = KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture()
    pressure = bounds.pressure_bounds(
        rho=0.1,
        g_norm=1.1,
        Phi=BallFactorBound(norm=2.0, lipschitz=0.75),
    )
    patched = bounds.inject_pressure_bounds(base, pressure)

    assert patched.p == pressure["p"]
    assert patched.p_eta == pressure["p_eta"]
    assert patched.Y_p_Y == pressure["Y_p_Y"]
    for name in base.__dataclass_fields__:
        if name not in {"p", "p_eta", "Y_p_Y"}:
            assert getattr(patched, name) == getattr(base, name)


def test_conditional_bridge_executes_existing_R1_R2_algebra():
    bounds = KokunoPA10PressureBallBounds()
    base = KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture()
    remainder = KokunoPA10RemainderBallBounds(
        h=0.01,
        A=0.51,
        D=0.49,
        rescaling_lambda=2.0,
    )
    bridge = bounds.conditional_remainder_bridge(
        rho=5e-4,
        g_norm=1.25,
        Phi=BallFactorBound(norm=2.5, lipschitz=1.0),
        primitives=base,
        remainder=remainder,
    )
    for key in ("R1", "R2"):
        total = bridge["remainder"][key]
        assert total["norm"] > 0.0
        assert total["lipschitz"] > 0.0
        assert math.isfinite(total["norm"])
        assert math.isfinite(total["lipschitz"])
    assert bridge["pressure"]["p_eta"].norm > bridge["pressure"]["p"].norm


def test_truth_boundary_refuses_source_and_pde_promotion():
    truth = KokunoPA10PressureBallBounds().truth_boundary
    assert truth["source_pressure_operator_algebra_executable"] is True
    assert truth["conditional_pressure_to_R1_R2_bridge_executable"] is True
    assert truth["diagnostic_inputs_are_source_bounds"] is False
    assert truth["source_rho_machine_bound"] is False
    assert truth["source_g_coefficient_norm_machine_bound"] is False
    assert truth["source_pressure_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_pressure_radius_one_ball_lipschitz_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_report_is_deterministic_and_hash_bound(tmp_path):
    bounds = KokunoPA10PressureBallBounds()
    first = bounds.report()
    second = bounds.report()
    assert first == second
    assert first["diagnostic_inputs_not_source_bounds"]["rho"] == 5e-4
    assert first["conditional_pressure_bounds"]["p"]["norm"] > 0.0
    assert first["conditional_pressure_bounds"]["p_eta"]["norm"] > 0.0
    assert first["conditional_remainder_totals_after_pressure_injection"]["R2"]["norm"] > 0.0

    identity = dict(first)
    expected = identity.pop("receipt_sha256")
    raw = json.dumps(identity, sort_keys=True, separators=(",", ":"), allow_nan=False)
    assert hashlib.sha256(raw.encode("utf-8")).hexdigest() == expected

    path = tmp_path / "pressure.json"
    saved = bounds.save_report(path)
    assert json.loads(path.read_text()) == saved


def test_invalid_inputs_fail_closed():
    bounds = KokunoPA10PressureBallBounds()
    phi = BallFactorBound(norm=1.0, lipschitz=1.0)
    for bad_rho in (0.0, -1.0, math.inf, math.nan):
        with pytest.raises(ValueError):
            bounds.pressure_bounds(rho=bad_rho, g_norm=1.0, Phi=phi)
    with pytest.raises(ValueError):
        bounds.pressure_bounds(rho=0.1, g_norm=-1.0, Phi=phi)
    with pytest.raises(TypeError):
        bounds.pressure_bounds(rho=0.1, g_norm=1.0, Phi=object())
    with pytest.raises(ValueError):
        bounds.inject_pressure_bounds(
            KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture(),
            {"p": phi},
        )
