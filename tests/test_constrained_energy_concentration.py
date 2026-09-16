import math
import numpy as np
import pytest

from openai_ns_reconstruction.constrained_energy_concentration import (
    diagnose_axisymmetric_energy_concentration,
)


def compact_swirl(points, time):
    points = np.asarray(points, dtype=float)
    x, y, z = np.moveaxis(points, -1, 0)
    r2 = x*x + y*y
    radial = np.maximum(1.0-r2, 0.0)**3
    axial = np.maximum(1.0-z*z, 0.0)**3
    amp = radial*axial*(1.0+0.1*time)
    return np.stack((-y*amp, x*amp, np.zeros_like(x)), axis=-1)


def test_compact_swirl_energy_and_concentration_are_resolution_stable():
    report = diagnose_axisymmetric_energy_concentration(
        compact_swirl, 0.5,
        support_radius=1.0,
        support_half_height=1.0,
        levels=(0.25,0.5,0.75,1.0),
        orders=(16,32,64),
    )
    exact_base = 128.0*math.pi/21021.0
    exact = exact_base*(1.0+0.05)**2
    assert report.finest_total_energy == pytest.approx(exact, rel=2e-13, abs=2e-15)
    assert report.max_adjacent_order_fraction_delta < 2e-13
    assert report.concentration_scale_spread["q50"] < 2e-13
    assert report.concentration_scale_spread["q90"] < 2e-13
    assert all(a <= b for a,b in zip(report.finest_fraction_profile, report.finest_fraction_profile[1:]))
    assert report.finest_fraction_profile[-1] == 1.0
    assert "not evidence of blow-up" in report.scope


def test_time_amplitude_changes_total_energy_but_not_normalised_profile():
    a = diagnose_axisymmetric_energy_concentration(
        compact_swirl, 0.0, support_radius=1.0, support_half_height=1.0,
        orders=(12,24,48),
    )
    b = diagnose_axisymmetric_energy_concentration(
        compact_swirl, 1.0, support_radius=1.0, support_half_height=1.0,
        orders=(12,24,48),
    )
    assert b.finest_total_energy/a.finest_total_energy == pytest.approx(1.21, rel=2e-13)
    assert b.finest_fraction_profile == pytest.approx(a.finest_fraction_profile, abs=2e-13)


def test_fail_closed_inputs_and_velocity_contract():
    with pytest.raises(ValueError):
        diagnose_axisymmetric_energy_concentration(compact_swirl, 0.0, support_radius=0.0, support_half_height=1.0)
    with pytest.raises(ValueError):
        diagnose_axisymmetric_energy_concentration(compact_swirl, 0.0, support_radius=1.0, support_half_height=1.0, levels=(0.5,0.4,1.0))
    with pytest.raises(ValueError):
        diagnose_axisymmetric_energy_concentration(compact_swirl, 0.0, support_radius=1.0, support_half_height=1.0, orders=(16,))
    def bad_shape(points, time):
        return np.zeros(points.shape[:-1])
    with pytest.raises(ValueError):
        diagnose_axisymmetric_energy_concentration(bad_shape, 0.0, support_radius=1.0, support_half_height=1.0)
    def nonfinite(points, time):
        out=np.zeros_like(points); out[...,0]=np.nan; return out
    with pytest.raises(ValueError):
        diagnose_axisymmetric_energy_concentration(nonfinite, 0.0, support_radius=1.0, support_half_height=1.0)
