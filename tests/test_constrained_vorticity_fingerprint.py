import numpy as np
import pytest

from openai_ns_reconstruction.constrained_vorticity_fingerprint import (
    audit_vorticity_resolution,
    diagnose_vorticity_fingerprint,
)


def solid_rotation(points, time):
    p = np.asarray(points, dtype=float)
    x, y, _ = np.moveaxis(p, -1, 0)
    amp = 1.0 + 0.25 * time
    return np.stack((-amp * y, amp * x, np.zeros_like(x)), axis=-1)


def smooth_axisymmetric_swirl(points, time):
    p = np.asarray(points, dtype=float)
    x, y, z = np.moveaxis(p, -1, 0)
    r2 = x*x + y*y
    a = np.exp(-2.0*r2 - 1.5*z*z) * (1.0 + 0.2*time)
    return np.stack((-y*a, x*a, np.zeros_like(x)), axis=-1)


def test_solid_rotation_has_exact_constant_curl():
    fp = diagnose_vorticity_fingerprint(
        solid_rotation, 0.4, box=((-1,1),(-1,1),(-1,1)), grid_size=17
    )
    expected = 2.0 * (1.0 + 0.25*0.4)
    assert fp.omega_peak == pytest.approx(expected, rel=2e-13, abs=2e-13)
    assert fp.omega_rms == pytest.approx(expected, rel=2e-13, abs=2e-13)
    assert fp.half_peak_volume_fraction == 1.0


def test_three_level_smooth_fingerprint_refines():
    audit = audit_vorticity_resolution(
        smooth_axisymmetric_swirl,
        0.5,
        box=((-1.5,1.5),(-1.5,1.5),(-1.5,1.5)),
        grid_sizes=(17,25,33),
    )
    levels = audit["levels"]
    assert [row["grid_size"] for row in levels] == [17,25,33]
    assert levels[-1]["omega_peak"] > 0
    assert 0 < levels[-1]["enstrophy_r50"] < levels[-1]["enstrophy_r90"] < 1.5*np.sqrt(2)
    assert 0 <= levels[-1]["enstrophy_absz50"] <= levels[-1]["enstrophy_absz90"] < 1.5
    changes = audit["relative_change_to_finest"]
    assert changes["omega_peak"][-1] < 0.03
    assert changes["omega_rms"][-1] < 0.03
    assert changes["enstrophy_r90"][-1] < 0.08


def test_fail_closed_inputs():
    with pytest.raises(ValueError):
        diagnose_vorticity_fingerprint(solid_rotation, 0.2, grid_size=4)
    with pytest.raises(ValueError):
        diagnose_vorticity_fingerprint(
            solid_rotation, 0.2, box=((1,0),(-1,1),(-1,1)), grid_size=9
        )
    with pytest.raises(ValueError):
        audit_vorticity_resolution(solid_rotation, 0.2, grid_sizes=(9,9,17))

    def zero(points, time):
        return np.zeros_like(points, dtype=float)
    with pytest.raises(ValueError):
        diagnose_vorticity_fingerprint(zero, 0.2, grid_size=9)
