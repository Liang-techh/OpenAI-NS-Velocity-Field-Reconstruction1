import numpy as np
import pytest

from openai_ns_reconstruction.constrained_streamline_torsion import (
    measure_streamline_torsion_profile,
)


def helix(n=81, radius=1.4, rise=0.65, turns=2.5):
    theta = np.linspace(-np.pi * turns, np.pi * turns, n)
    return np.column_stack((radius * np.cos(theta), radius * np.sin(theta), rise * theta))


def test_analytic_helix_torsion_and_scale():
    radius = 1.4
    rise = 0.65
    pts = helix(radius=radius, rise=rise)
    report = measure_streamline_torsion_profile(pts, provenance="analytic helix")
    expected = rise / (radius * radius + rise * rise)
    assert report.valid_arclength_fraction > 0.99
    assert report.mean_signed_torsion == pytest.approx(expected, rel=3e-3)
    assert report.rms_torsion == pytest.approx(expected, rel=4e-3)
    assert report.handedness_coherence > 0.999

    scale = 3.25
    scaled = measure_streamline_torsion_profile(scale * pts, provenance="scaled helix")
    assert scaled.mean_signed_torsion == pytest.approx(report.mean_signed_torsion / scale, rel=5e-4)
    assert scaled.integrated_signed_torsion == pytest.approx(report.integrated_signed_torsion, rel=5e-4)


def test_rigid_motion_preserves_torsion_and_reflection_flips_sign():
    pts = helix()
    angle = 0.71
    c, s = np.cos(angle), np.sin(angle)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    moved = pts @ rot.T + np.array([2.0, -3.0, 0.4])
    base = measure_streamline_torsion_profile(pts, provenance="base")
    rigid = measure_streamline_torsion_profile(moved, provenance="rigid")
    assert rigid.mean_signed_torsion == pytest.approx(base.mean_signed_torsion, rel=2e-10, abs=2e-10)
    assert rigid.integrated_signed_torsion == pytest.approx(base.integrated_signed_torsion, rel=2e-10, abs=2e-10)

    reversed_report = measure_streamline_torsion_profile(pts[::-1], provenance="reversed parameterization")
    assert reversed_report.mean_signed_torsion == pytest.approx(base.mean_signed_torsion, rel=2e-10, abs=2e-10)
    assert reversed_report.integrated_signed_torsion == pytest.approx(base.integrated_signed_torsion, rel=2e-10, abs=2e-10)

    reflected = pts.copy()
    reflected[:, 0] *= -1.0
    mirror = measure_streamline_torsion_profile(reflected, provenance="mirror")
    assert mirror.mean_signed_torsion == pytest.approx(-base.mean_signed_torsion, rel=2e-10, abs=2e-10)
    assert mirror.mean_abs_torsion == pytest.approx(base.mean_abs_torsion, rel=2e-10, abs=2e-10)


def test_planar_curve_has_zero_torsion_without_false_handness():
    theta = np.linspace(-1.2, 1.4, 65)
    pts = np.column_stack((2.0 * np.cos(theta), 2.0 * np.sin(theta), np.zeros_like(theta)))
    report = measure_streamline_torsion_profile(pts, provenance="planar arc")
    assert report.valid_arclength_fraction > 0.99
    assert report.max_abs_torsion < 1e-10
    assert report.integrated_abs_torsion < 1e-10
    assert report.handedness_coherence == 0.0


def test_outputs_are_read_only_and_truth_boundary_is_false():
    report = measure_streamline_torsion_profile(helix(), provenance="frozen candidate streamline")
    for arr in (report.s_fraction, report.points, report.torsion, report.valid_mask):
        assert not arr.flags.writeable
    assert report.metadata["smoothing_applied"] is False
    assert report.metadata["camera_or_registration_fit"] is False
    assert report.metadata["velocity_changed"] is False
    assert report.metadata["visualization_ready"] is False
    assert report.metadata["visual_correspondence_verified"] is False
    assert report.metadata["pde_validated"] is False
    assert report.metadata["openai_field_identified"] is False


def test_fail_closed_on_degenerate_or_malformed_inputs():
    good = helix()
    bad_cases = [
        good[:5],
        np.zeros((8, 3)),
        np.column_stack((np.linspace(0.0, 1.0, 8), np.zeros(8), np.zeros(8))),
        good[:, :2],
    ]
    duplicate = good.copy()
    duplicate[12] = duplicate[11]
    bad_cases.append(duplicate)
    nonfinite = good.copy()
    nonfinite[4, 1] = np.nan
    bad_cases.append(nonfinite)
    for bad in bad_cases:
        with pytest.raises(ValueError):
            measure_streamline_torsion_profile(bad, provenance="bad")

    with pytest.raises(ValueError):
        measure_streamline_torsion_profile(good, sample_count=64, provenance="bad count")
    with pytest.raises(ValueError):
        measure_streamline_torsion_profile(good, provenance="")
