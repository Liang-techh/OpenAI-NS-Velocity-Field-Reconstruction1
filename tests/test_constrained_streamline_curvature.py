import numpy as np
import pytest

from openai_ns_reconstruction.constrained_streamline_curvature import measure_streamline_curvature_profile


def test_straight_line_has_zero_curvature():
    z = np.linspace(-1.0, 1.0, 9)
    pts = np.column_stack((0.2 * z, -0.1 * z, z))
    r = measure_streamline_curvature_profile(pts, provenance="manufactured-straight")
    assert r.max_curvature < 1e-10
    assert r.total_turning < 1e-10
    assert r.tortuosity == pytest.approx(1.0, abs=1e-12)
    assert not r.metadata["pde_validated"]
    assert not r.points.flags.writeable


def test_helix_curvature_and_tortuosity_are_nonzero():
    th = np.linspace(0.0, 4 * np.pi, 121)
    a, b = 1.25, 0.35
    pts = np.column_stack((a * np.cos(th), a * np.sin(th), b * th))
    r = measure_streamline_curvature_profile(
        pts, sample_count=257, provenance="manufactured-helix"
    )
    expected = a / (a * a + b * b)
    core = r.curvature[20:-20]
    assert np.median(core) == pytest.approx(expected, rel=0.01)
    assert r.tortuosity > 2.0
    assert r.total_turning > 5.0
    assert r.metadata["smoothing_applied"] is False


def test_scale_changes_curvature_inverse_but_not_tortuosity():
    th = np.linspace(0.0, 2 * np.pi, 81)
    pts = np.column_stack((np.cos(th), np.sin(th), 0.2 * th))
    a = measure_streamline_curvature_profile(pts, provenance="scale-a")
    b = measure_streamline_curvature_profile(3.0 * pts, provenance="scale-b")
    assert b.mean_curvature == pytest.approx(a.mean_curvature / 3.0, rel=2e-3)
    assert b.tortuosity == pytest.approx(a.tortuosity, rel=2e-3)


def test_fail_closed_inputs():
    good = np.column_stack((np.linspace(0, 1, 6), np.zeros(6), np.linspace(0, 2, 6)))
    with pytest.raises(ValueError):
        measure_streamline_curvature_profile(good[:4], provenance="x")
    bad = good.copy()
    bad[2] = bad[1]
    with pytest.raises(ValueError):
        measure_streamline_curvature_profile(bad, provenance="x")
    bad = good.copy()
    bad[1, 0] = np.nan
    with pytest.raises(ValueError):
        measure_streamline_curvature_profile(bad, provenance="x")
    with pytest.raises(ValueError):
        measure_streamline_curvature_profile(good, sample_count=32, provenance="x")
    with pytest.raises(ValueError):
        measure_streamline_curvature_profile(good, provenance="")
