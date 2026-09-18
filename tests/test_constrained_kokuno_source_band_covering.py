import math

import pytest

from openai_ns_reconstruction.kokuno_physical_evaluation import (
    SOURCE_KAPPA_S,
    SOURCE_LAMBDA_G,
    SOURCE_T_G,
)
from openai_ns_reconstruction.kokuno_source_band_covering import KokunoSourceBandCovering


def test_displayed_band_schedule_matches_independent_formulas():
    h = 0.004
    for ell in (5, 12, 20, 40):
        band = KokunoSourceBandCovering(ell=ell, h=h)
        q = 2.0 ** (-ell)
        eps = q**h
        s_star = float(ell**2)
        rho = math.log(SOURCE_LAMBDA_G) / math.log(SOURCE_T_G)
        d_r = 2.0 * ((1.0 + h) * rho - h * SOURCE_KAPPA_S)
        i = math.floor(math.log(q ** (-1.0 - h) / s_star) / math.log(SOURCE_T_G))
        c_i = SOURCE_T_G**i * q ** (1.0 + h)
        m_i = SOURCE_LAMBDA_G**i * q ** (d_r / 2.0)

        assert band.Q == pytest.approx(q, rel=2e-15)
        assert band.epsilon == pytest.approx(eps, rel=2e-15)
        assert band.S_star == s_star
        assert band.product_grid_mesh == pytest.approx(ell ** -6)
        assert band.covering_level == i
        assert band.c_i == pytest.approx(c_i, rel=3e-15)
        assert band.M_i == pytest.approx(m_i, rel=3e-15)


def test_floor_bounds_hold_across_multiple_source_bands():
    h = 0.003
    for ell in range(5, 61):
        band = KokunoSourceBandCovering(ell=ell, h=h)
        d = band.bound_diagnostics()
        assert d["c_i_bound_passed"] is True
        assert d["M_i_bound_passed"] is True
        assert 1.0 / SOURCE_T_G < d["c_i_scaled"] <= 1.0
        assert SOURCE_T_G ** (-band.rho_g) < d["M_i_scaled"] <= 1.0


def test_interacting_band_covering_level_bound_is_executable():
    h = 0.004
    bands = [KokunoSourceBandCovering(ell=e, h=h) for e in range(12, 21)]
    checked = 0
    for left in bands:
        for right in bands:
            if abs(left.ell - right.ell) <= 4:
                receipt = left.interaction_receipt(right, ell0=12)
                assert receipt["bound_passed"] is True
                assert receipt["level_difference"] <= receipt["source_upper_bound"]
                checked += 1
    assert checked > 20


def test_pulse_length_keeps_nonunique_r0_truth_boundary():
    band = KokunoSourceBandCovering(ell=20, h=0.004)
    pulse = band.pulse_length_from_r0(0.03125)
    assert pulse["r0_binding"] == "caller_supplied_not_source_recovered"
    assert pulse["L_s"] == pytest.approx(2.0 * 0.03125 / band.c_i)

    receipt = band.receipt()
    truth = receipt["truth_boundary"]
    assert truth["source_band_covering_schedule_executable"] is True
    assert truth["source_rectangle_centers_recovered"] is False
    assert truth["source_rectangle_radius_r0_recovered"] is False
    assert truth["source_partition_grid_origin_recovered"] is False
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["actual_auxiliary_torus_mode_family_bound"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["paper_exact"] is False
    assert len(receipt["sha256"]) == 64


def test_fail_closed_band_and_interaction_guards():
    with pytest.raises(ValueError, match="integer"):
        KokunoSourceBandCovering(ell=12.0, h=0.004)
    with pytest.raises(ValueError, match="exceed 4"):
        KokunoSourceBandCovering(ell=4, h=0.004)
    with pytest.raises(ValueError, match="0<h<1/100"):
        KokunoSourceBandCovering(ell=12, h=0.0)

    a = KokunoSourceBandCovering(ell=12, h=0.004)
    b = KokunoSourceBandCovering(ell=17, h=0.004)
    with pytest.raises(ValueError, match="restricted"):
        a.interaction_receipt(b)
    with pytest.raises(ValueError, match="same h"):
        a.interaction_receipt(KokunoSourceBandCovering(ell=13, h=0.003))
    with pytest.raises(ValueError, match="ell0 cannot exceed"):
        a.interaction_receipt(KokunoSourceBandCovering(ell=13, h=0.004), ell0=13)
    with pytest.raises(ValueError, match="positive"):
        a.pulse_length_from_r0(0.0)
