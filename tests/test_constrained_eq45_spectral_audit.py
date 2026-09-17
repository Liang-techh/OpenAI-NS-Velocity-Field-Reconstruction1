from pathlib import Path

import numpy as np

from openai_ns_reconstruction.constrained_eq45_spectral_audit import (
    audit_checked_eq45_seed,
    audit_spectral_resolution,
    diagnose_spectrum_level,
)


class _PeriodicMode:
    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        out = np.zeros_like(points)
        out[..., 0] = np.sin(0.5 * np.pi * points[..., 0])
        return out


class _NonperiodicRamp:
    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        out = np.zeros_like(points)
        out[..., 0] = points[..., 0]
        return out


def test_periodic_single_mode_has_shell_one_and_zero_seam():
    row = diagnose_spectrum_level(
        _PeriodicMode(), time=0.5, grid_size=16, fixed_shell_cutoff=6
    )
    assert row["dominant_shell"] == 1
    assert row["shell_q50"] == 1
    assert row["shell_q90"] == 1
    assert row["parseval_relative_error"] < 1e-14
    assert row["opposite_face_seam_max"] < 1e-12
    assert row["near_nyquist_energy_fraction"] < 1e-20


def test_nonperiodic_ramp_exposes_fft_seam_obstruction():
    row = diagnose_spectrum_level(
        _NonperiodicRamp(), time=0.5, grid_size=16, fixed_shell_cutoff=6
    )
    assert row["opposite_face_seam_rms"] > 1.0
    assert row["seam_rms_over_velocity_rms"] > 0.5
    assert row["near_nyquist_energy_fraction"] > 0.0


def test_audit_requires_three_strictly_increasing_levels():
    try:
        audit_spectral_resolution(_PeriodicMode(), grid_sizes=(16, 16, 24))
    except ValueError as exc:
        assert "strictly increasing" in str(exc)
    else:
        raise AssertionError("duplicate grid levels must be rejected")


def test_checked_eq45_seed_reports_resolution_and_time_scale_shift():
    candidate = (
        Path(__file__).resolve().parents[1]
        / "artifacts"
        / "constrained"
        / "eq45_velocity_candidate_seed.json"
    )
    report = audit_checked_eq45_seed(candidate)
    assert report["candidate_sha256"] == (
        "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
    )
    assert report["grid_sizes"] == [16, 24, 32]
    assert report["times"] == [0.25, 0.5, 0.75]

    finest = {
        item["time"]: item for item in report["time_summaries"]
    }
    assert [finest[t]["finest_dominant_shell"] for t in report["times"]] == [1, 2, 3]
    assert [finest[t]["finest_shell_q90"] for t in report["times"]] == [3, 4, 6]
    assert finest[0.25]["finest_seam_rms_over_velocity_rms"] > 0.03
    assert finest[0.75]["finest_seam_rms_over_velocity_rms"] < 0.001
    assert max(
        finest[t]["finest_near_nyquist_energy_fraction"] for t in report["times"]
    ) < 5e-4

    rows = {(row["grid_size"], row["time"]): row for row in report["rows"]}
    for time in report["times"]:
        medium = rows[(24, time)]["direct_energy"]
        fine = rows[(32, time)]["direct_energy"]
        assert abs(medium - fine) / fine < 0.01

    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["blowup_proved"] is False
