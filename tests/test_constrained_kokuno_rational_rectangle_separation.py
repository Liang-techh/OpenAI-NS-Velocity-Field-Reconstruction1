from fractions import Fraction
import math

import pytest

from openai_ns_reconstruction.kokuno_source_band_covering import KokunoSourceBandCovering
from openai_ns_reconstruction.kokuno_rational_rectangle_separation import (
    KokunoRationalRectangleSeparation,
)


def test_autonomous_witness_exactly_separates_all_source_constraints():
    cert = KokunoRationalRectangleSeparation.autonomous_for_band_window(
        h=0.004,
        ell0=12,
        color_count=4,
        denominator=17,
        safety=0.5,
    )
    assert cert.delta_max == 4
    rows = cert.constraint_rows()
    assert len(rows) == 4 * 4 * 5 - 4
    assert all(row["separated"] is True for row in rows)
    exact_min = min(
        Fraction(row["distance_squared_numerator"], row["distance_squared_denominator"])
        for row in rows
    )
    assert cert.d_c_squared_exact == exact_min
    assert cert.d_c > 0.0
    assert cert.injectivity_guard_passed is True
    assert cert.separation_guard_passed is True
    assert 4.0 * cert.r0 * cert.C_v < 1.0
    assert 2.0 * cert.r0 * cert.C_v * (cert.C_J + 1.0) < cert.d_c


def test_autonomous_witness_is_deterministic_and_truthful():
    kwargs = dict(color_count=5, denominator=17, delta_max=3, safety=0.4)
    left = KokunoRationalRectangleSeparation.autonomous_rational_witness(**kwargs)
    right = KokunoRationalRectangleSeparation.autonomous_rational_witness(**kwargs)
    assert left.centers == right.centers
    assert left.r0 == right.r0
    assert left.receipt()["sha256"] == right.receipt()["sha256"]

    truth = left.receipt()["truth_boundary"]
    assert truth["source_finite_rectangle_separation_structure_executable"] is True
    assert truth["source_small_r0_inequalities_executable"] is True
    assert truth["autonomous_rational_center_witness_available"] is True
    assert truth["source_rectangle_centers_recovered"] is False
    assert truth["source_rectangle_radius_r0_recovered"] is False
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["actual_auxiliary_torus_mode_family_bound"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["paper_exact"] is False


def test_source_norm_bounds_and_band_pulse_lengths():
    cert = KokunoRationalRectangleSeparation.autonomous_for_band_window(
        h=0.004,
        ell0=12,
        color_count=3,
        denominator=17,
    )
    b_g = math.sqrt(2.0) - 1.0
    expected_cv = 2.0 * math.sqrt(1.0 + b_g * b_g)
    expected_cj = (4.0 + math.sqrt(2.0)) ** cert.delta_max
    assert cert.C_v == pytest.approx(expected_cv, rel=2e-15)
    assert cert.C_J == pytest.approx(expected_cj, rel=2e-15)

    bands = [
        KokunoSourceBandCovering(ell=12, h=0.004),
        KokunoSourceBandCovering(ell=13, h=0.004),
    ]
    rows = cert.pulse_lengths(bands)
    assert [row["ell"] for row in rows] == [12, 13]
    for row, band in zip(rows, bands):
        assert row["L_s"] == pytest.approx(2.0 * cert.r0 / band.c_i)
        assert row["r0_binding"] == "repository_autonomous_source_compatible_safety_fraction"


def test_collision_and_radius_guards_fail_closed():
    with pytest.raises(ValueError, match="separation"):
        KokunoRationalRectangleSeparation(
            centers=(
                (Fraction(0), Fraction(0)),
                (Fraction(0), Fraction(1, 3)),
            ),
            delta_max=2,
            r0=1e-8,
        )

    cert = KokunoRationalRectangleSeparation.autonomous_rational_witness(
        color_count=3,
        denominator=17,
        delta_max=2,
        safety=0.5,
    )
    with pytest.raises(ValueError, match="separation guard|injectivity guard"):
        KokunoRationalRectangleSeparation(
            centers=cert.centers,
            delta_max=cert.delta_max,
            r0=1.01 * cert.strict_radius_upper,
        )


def test_invalid_autonomous_choices_fail_closed():
    with pytest.raises(ValueError, match="color_count"):
        KokunoRationalRectangleSeparation.autonomous_rational_witness(
            color_count=1, denominator=17, delta_max=2
        )
    with pytest.raises(ValueError, match="denominator"):
        KokunoRationalRectangleSeparation.autonomous_rational_witness(
            color_count=3, denominator=1, delta_max=2
        )
    with pytest.raises(ValueError, match="delta_max"):
        KokunoRationalRectangleSeparation.autonomous_rational_witness(
            color_count=3, denominator=17, delta_max=-1
        )
    with pytest.raises(ValueError, match="safety"):
        KokunoRationalRectangleSeparation.autonomous_rational_witness(
            color_count=3, denominator=17, delta_max=2, safety=1.0
        )
