import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_blend_spectral_fingerprint import (
    SCHEMA,
    audit_supported_phi10_blend_spectral_fingerprint,
    diagnose_blend_spectrum_level,
)


class PeriodicMode:
    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        x = points[..., 0]
        return np.stack((np.sin(0.5 * np.pi * x), np.zeros_like(x), np.zeros_like(x)), axis=-1)


def test_periodic_mode_parseval_shell_and_seam():
    row = diagnose_blend_spectrum_level(
        PeriodicMode(), time=0.3125, grid_size=24, box_half_width=2.0, fixed_shell_cutoff=6
    )
    assert row["parseval_relative_error"] < 1e-13
    assert row["dominant_shell"] == 1
    assert row["shell_q90"] == 1
    assert row["fixed_shell_tail_fraction"] < 1e-28
    assert row["opposite_face_seam_rms"] < 1e-14


def test_materialized_supported_blend_three_resolution_regression():
    report = audit_supported_phi10_blend_spectral_fingerprint()
    assert report["schema"] == SCHEMA
    assert report["blend_weights"] == [0.0, 0.5, 1.0]
    assert report["grid_sizes"] == [24, 32, 40]
    assert len(report["rows"]) == 9
    for row in report["rows"]:
        assert row["parseval_relative_error"] < 1e-12
        assert 0.0 <= row["fixed_shell_tail_fraction"] <= 1.0
        assert 0.0 <= row["near_nyquist_energy_fraction"] <= 1.0
        assert row["seam_rms_over_velocity_rms"] < 1e-10

    summaries = report["summaries"]
    energy = np.array([row["finest_direct_energy"] for row in summaries])
    shell_rms = np.array([row["finest_shell_rms"] for row in summaries])
    fixed_tail = np.array([row["finest_fixed_shell_tail_fraction"] for row in summaries])
    near_nyquist = np.array([row["finest_near_nyquist_energy_fraction"] for row in summaries])

    np.testing.assert_allclose(
        energy,
        [0.7585189378488971, 0.7322297758443781, 0.7076837868151513],
        rtol=2e-6,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        shell_rms,
        [2.4256510007454577, 2.419252310899209, 2.4145813548907067],
        rtol=2e-6,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        fixed_tail,
        [0.01033102991654133, 0.009924795360989112, 0.009570897430566603],
        rtol=2e-6,
        atol=1e-12,
    )
    assert np.all(np.diff(energy) < 0.0)
    assert np.all(np.diff(shell_rms) < 0.0)
    assert np.all(np.diff(fixed_tail) < 0.0)
    assert [row["finest_dominant_shell"] for row in summaries] == [2, 2, 2]
    assert [row["finest_shell_q50"] for row in summaries] == [2, 2, 2]
    assert [row["finest_shell_q90"] for row in summaries] == [4, 4, 4]
    assert [row["finest_shell_q99"] for row in summaries] == [6, 5, 5]
    assert np.max(near_nyquist) < 4e-5

    for row in summaries:
        assert row["medium_to_fine_energy_relative_change"] < 2e-5
        assert row["medium_to_fine_shell_rms_relative_change"] < 5e-4
        assert row["medium_to_fine_fixed_tail_relative_change"] < 5e-3

    comparison = report["finest_endpoint_comparison"]
    assert 0.92 < comparison["compact_over_quartic_energy"] < 0.95
    assert 0.99 < comparison["compact_over_quartic_shell_rms"] < 1.0
    assert 0.91 < comparison["compact_over_quartic_fixed_tail"] < 0.95
    assert comparison["midpoint_shell_rms_between_endpoints"] is True
    assert comparison["midpoint_fixed_tail_between_endpoints"] is True


def test_fail_closed_grid_and_cutoff_contract():
    with pytest.raises(ValueError):
        audit_supported_phi10_blend_spectral_fingerprint(grid_sizes=(24, 32))
    with pytest.raises(ValueError):
        audit_supported_phi10_blend_spectral_fingerprint(grid_sizes=(24, 32, 40), fixed_shell_cutoff=12)
