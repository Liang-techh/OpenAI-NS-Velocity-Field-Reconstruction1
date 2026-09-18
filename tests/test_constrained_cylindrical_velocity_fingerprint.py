import math

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_cylindrical_velocity_fingerprint import (
    audit_cylindrical_velocity_resolution,
    diagnose_cylindrical_velocity_fingerprint,
)


def mixed_linear(points, time):
    p = np.asarray(points, dtype=float)
    x, y, z = np.moveaxis(p, -1, 0)
    del z
    a = 2.0 + 0.0 * time
    b = 1.0
    return np.stack((b * x - a * y, b * y + a * x, np.zeros_like(x)), axis=-1)


def counter_rotating(points, time):
    p = np.asarray(points, dtype=float)
    x, y, z = np.moveaxis(p, -1, 0)
    del time
    s = np.where(z >= 0.0, 1.0, -1.0)
    return np.stack((-s * y, s * x, np.zeros_like(x)), axis=-1)


def test_mixed_linear_exact_swirl_poloidal_split_and_pitch():
    fp = diagnose_cylindrical_velocity_fingerprint(
        mixed_linear, 0.5, box=((-1, 1), (-1, 1), (-1, 1)), grid_size=17
    )
    assert fp.swirl_energy_fraction == pytest.approx(4.0 / 5.0, abs=2e-15)
    assert fp.poloidal_energy_fraction == pytest.approx(1.0 / 5.0, abs=2e-15)
    assert fp.swirl_to_poloidal_ratio == pytest.approx(4.0, abs=2e-14)
    assert fp.radial_fraction_of_poloidal == pytest.approx(1.0, abs=2e-15)
    expected_pitch = math.degrees(math.atan2(1.0, 2.0))
    assert fp.pitch_deg_q10 == pytest.approx(expected_pitch, abs=2e-13)
    assert fp.pitch_deg_q50 == pytest.approx(expected_pitch, abs=2e-13)
    assert fp.pitch_deg_q90 == pytest.approx(expected_pitch, abs=2e-13)
    assert fp.swirl_sign_coherence == pytest.approx(1.0, abs=2e-15)
    assert fp.axis_transverse_energy_fraction == pytest.approx(0.0, abs=1e-15)


def test_three_resolution_levels_stay_stable_for_linear_field():
    audit = audit_cylindrical_velocity_resolution(
        mixed_linear,
        0.5,
        box=((-1, 1), (-1, 1), (-1, 1)),
        grid_sizes=(17, 25, 33),
    )
    levels = audit["levels"]
    assert [row["grid_size"] for row in levels] == [17, 25, 33]
    for row in levels:
        assert row["swirl_energy_fraction"] == pytest.approx(0.8, abs=2e-15)
        assert row["pitch_deg_q50"] == pytest.approx(
            math.degrees(math.atan2(1.0, 2.0)), abs=2e-13
        )
    exact_speed2_integral = 80.0 / 3.0
    errors = [abs(row["speed2_integral"] - exact_speed2_integral) for row in levels]
    assert errors[2] < errors[1] < errors[0]
    assert audit["absolute_change_to_finest"]["swirl_energy_fraction"][-1] < 2e-15
    assert audit["absolute_change_to_finest"]["pitch_deg_q50"][-1] < 2e-13


def test_counter_rotation_reduces_sign_coherence_without_changing_swirl_fraction():
    fp = diagnose_cylindrical_velocity_fingerprint(
        counter_rotating, 0.0, box=((-1, 1), (-1, 1), (-1, 1)), grid_size=17
    )
    assert fp.swirl_energy_fraction == pytest.approx(1.0, abs=2e-15)
    assert fp.poloidal_energy_fraction == pytest.approx(0.0, abs=2e-15)
    assert fp.swirl_to_poloidal_ratio is None
    assert fp.swirl_sign_coherence < 0.08
    assert fp.pitch_deg_q50 == pytest.approx(0.0, abs=1e-15)


def test_fail_closed_inputs():
    with pytest.raises(ValueError):
        diagnose_cylindrical_velocity_fingerprint(mixed_linear, 0.0, grid_size=6)
    with pytest.raises(ValueError):
        diagnose_cylindrical_velocity_fingerprint(
            mixed_linear, 0.0, box=((1, 2), (-1, 1), (-1, 1)), grid_size=9
        )
    with pytest.raises(ValueError):
        audit_cylindrical_velocity_resolution(
            mixed_linear, 0.0, grid_sizes=(9, 9, 17)
        )

    def zero(points, time):
        del time
        return np.zeros_like(points, dtype=float)

    with pytest.raises(ValueError):
        diagnose_cylindrical_velocity_fingerprint(zero, 0.0, grid_size=9)
