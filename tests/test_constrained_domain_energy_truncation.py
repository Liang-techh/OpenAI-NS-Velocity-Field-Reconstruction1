import math

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_domain_energy_truncation import (
    diagnose_finite_domain_energy_truncation,
    integrate_cube_kinetic_energy,
)


def gaussian_velocity(points, time):
    points = np.asarray(points, dtype=float)
    amplitude = np.exp(-0.5 * np.sum(points * points, axis=-1))
    out = np.zeros_like(points)
    out[:, 0] = amplitude
    return out


def constant_velocity(points, time):
    points = np.asarray(points, dtype=float)
    out = np.zeros_like(points)
    out[:, 0] = 1.0
    return out


def analytic_gaussian_energy(half_width):
    one_dimensional = math.sqrt(math.pi) * math.erf(half_width)
    return 0.5 * one_dimensional**3


def test_gaussian_energy_matches_analytic_nested_cube_and_reports_tail_unknown():
    widths = (1.0, 2.0, 3.0)
    report = diagnose_finite_domain_energy_truncation(
        gaussian_velocity,
        0.5,
        half_widths=widths,
        quadrature_orders=(12, 24, 48),
    )

    finest = report["finest_energies"]
    expected = [analytic_gaussian_energy(width) for width in widths]
    assert np.allclose(finest, expected, rtol=2e-12, atol=2e-12)
    assert report["finest_captured_fractions_of_largest_cube"][-1] == pytest.approx(1.0)
    assert 0.0 < report["finest_outermost_added_shell_fraction"] < 0.02
    assert report["tail_beyond_largest_cube_resolved"] is False
    assert report["visualization_ready"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False


def test_constant_field_keeps_large_outer_shell_instead_of_claiming_domain_convergence():
    report = diagnose_finite_domain_energy_truncation(
        constant_velocity,
        0.5,
        half_widths=(1.0, 2.0, 3.0),
        quadrature_orders=(8, 16, 32),
    )

    assert report["finest_energies"] == pytest.approx([4.0, 32.0, 108.0], rel=1e-13)
    assert report["finest_captured_fractions_of_largest_cube"] == pytest.approx(
        [1.0 / 27.0, 8.0 / 27.0, 1.0], rel=1e-13
    )
    assert report["finest_outermost_added_shell_fraction"] == pytest.approx(
        19.0 / 27.0, rel=1e-13
    )
    assert all(row["energy_monotone_with_domain"] for row in report["energy_rows"])


def test_single_box_integrator_and_fail_closed_inputs():
    energy = integrate_cube_kinetic_energy(
        constant_velocity,
        0.25,
        half_width=2.0,
        quadrature_order=8,
    )
    assert energy == pytest.approx(32.0, rel=1e-13)

    with pytest.raises(ValueError, match="at least three"):
        diagnose_finite_domain_energy_truncation(
            constant_velocity,
            0.5,
            half_widths=(1.0, 2.0),
            quadrature_orders=(8, 16, 32),
        )
    with pytest.raises(ValueError, match="strictly increasing"):
        diagnose_finite_domain_energy_truncation(
            constant_velocity,
            0.5,
            half_widths=(1.0, 3.0, 2.0),
            quadrature_orders=(8, 16, 32),
        )
    with pytest.raises(ValueError, match="finite integers"):
        diagnose_finite_domain_energy_truncation(
            constant_velocity,
            0.5,
            half_widths=(1.0, 2.0, 3.0),
            quadrature_orders=(8, 16.5, 32),
        )

    def bad_shape(points, time):
        return np.zeros((len(points), 2))

    with pytest.raises(ValueError, match="shape"):
        integrate_cube_kinetic_energy(bad_shape, 0.5, 1.0, 8)

    def bad_values(points, time):
        out = np.zeros_like(points)
        out[:, 0] = np.nan
        return out

    with pytest.raises(ValueError, match="non-finite"):
        integrate_cube_kinetic_energy(bad_values, 0.5, 1.0, 8)
