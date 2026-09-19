import hashlib
import json
import math

import pytest

from openai_ns_reconstruction.kokuno_pa10_mixed_jnu_seams import (
    KokunoPA10MixedJnuSeamBounds,
)
from openai_ns_reconstruction.kokuno_pa10_remainder_ball_bounds import BallFactorBound


def test_source_compatible_fixed_inputs_are_bound_without_standalone_radial_derivative():
    seam = KokunoPA10MixedJnuSeamBounds()
    assert seam.rho == 1.0e-6
    assert seam.product_constant_upper > 170.0
    assert seam.mixed_operator_factor() > 1.0e10
    assert seam.mixed_operator_factor(1) > seam.mixed_operator_factor()
    d = seam.d_coefficient_ball()
    assert d.norm == 1.0
    assert d.lipschitz == 0.0
    phi = seam.source_phi_ball()
    assert math.isfinite(phi.norm) and phi.norm > 1.0
    assert phi.lipschitz == 1.0


def test_generic_post_jnu_bound_uses_one_factor_at_a_time_lipschitz():
    seam = KokunoPA10MixedJnuSeamBounds()
    u = BallFactorBound(2.0, 0.5)
    g = BallFactorBound(3.0, 0.25)
    result = seam.averaged_eta_logradial_after_jnu(u=u, G=g, nu=2)
    c = seam.mixed_operator_factor()
    assert result.norm >= c * 6.0
    assert result.norm <= math.nextafter(c * 6.0, math.inf) * (1.0 + 1e-15)
    assert result.lipschitz >= c * 2.0
    assert result.lipschitz <= math.nextafter(c * 2.0, math.inf) * (1.0 + 1e-15)


def test_actual_two_mixed_operator_shapes_are_executable_but_u_conditional():
    seam = KokunoPA10MixedJnuSeamBounds()
    unit_u = BallFactorBound(1.0, 1.0)
    r1 = seam.r1_mixed_after_j2(u=unit_u)
    r2 = seam.r2_mixed_after_j1(u=unit_u)
    for bound in (r1, r2):
        assert math.isfinite(bound.norm) and bound.norm > 0.0
        assert math.isfinite(bound.lipschitz) and bound.lipschitz > 0.0
    assert r1.norm > r2.norm
    assert r2.lipschitz == pytest.approx(2.0 * r2.norm, rel=1e-15)


def test_extra_fixed_factor_adds_product_convolution_and_not_lipschitz():
    seam = KokunoPA10MixedJnuSeamBounds()
    u = BallFactorBound(1.0, 1.0)
    g = BallFactorBound(1.0, 0.0)
    plain = seam.averaged_eta_logradial_after_jnu(u=u, G=g, nu=1)
    with_d = seam.averaged_eta_logradial_after_jnu(
        u=u, G=g, nu=1, extra_factors=(seam.d_coefficient_ball(),)
    )
    assert with_d.norm > plain.norm
    assert with_d.lipschitz > plain.lipschitz
    assert with_d.norm / plain.norm == pytest.approx(seam.product_constant_upper, rel=2e-15)


def test_truth_boundary_closes_only_operator_algebra_not_source_u_or_R1_R2():
    seam = KokunoPA10MixedJnuSeamBounds()
    truth = seam.truth_boundary
    assert truth["source_post_Jnu_mixed_derivative_bound_executable"] is True
    assert truth["source_averaged_factor_cancellation_executable"] is True
    assert truth["source_compatible_d_coefficient_norm_machine_bound"] is True
    assert truth["standalone_Y_Phi_Y_coefficient_norm_invented"] is False
    assert truth["standalone_Y_u_Y_coefficient_norm_invented"] is False
    assert truth["source_u_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_post_J2_R1_mixed_numeric_ball_bound_complete"] is False
    assert truth["source_post_J1_R2_mixed_numeric_ball_bound_complete"] is False
    assert truth["source_R1_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_R2_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_report_is_deterministic_hash_bound_and_saveable(tmp_path):
    seam = KokunoPA10MixedJnuSeamBounds()
    first = seam.report()
    second = seam.report()
    assert first == second
    identity = dict(first)
    expected = identity.pop("receipt_sha256")
    raw = json.dumps(identity, sort_keys=True, separators=(",", ":"), allow_nan=False)
    assert hashlib.sha256(raw.encode()).hexdigest() == expected
    assert first["integration_boundary"]["standalone_Y_Phi_Y_bound_required"] is False
    assert first["integration_boundary"]["source_u_ball_still_required"] is True
    path = tmp_path / "mixed_jnu.json"
    saved = seam.save_report(path)
    assert json.loads(path.read_text()) == saved


def test_rejects_invalid_nu_and_non_ball_inputs():
    seam = KokunoPA10MixedJnuSeamBounds()
    u = BallFactorBound(1.0, 1.0)
    with pytest.raises(ValueError):
        seam.averaged_eta_logradial_after_jnu(u=u, G=u, nu=3)
    with pytest.raises(TypeError):
        seam.averaged_eta_logradial_after_jnu(u=u, G=object(), nu=1)
