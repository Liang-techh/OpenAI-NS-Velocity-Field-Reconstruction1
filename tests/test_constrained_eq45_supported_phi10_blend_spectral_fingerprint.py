import json

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


def test_materialized_supported_blend_three_resolution_calibration():
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
    # One sentinel calibration run records the exact governed numbers before
    # replacing this with tolerance-based regression expectations.
    raise AssertionError(json.dumps({
        "summaries": report["summaries"],
        "comparison": report["finest_endpoint_comparison"],
    }, sort_keys=True))


def test_fail_closed_grid_and_cutoff_contract():
    with pytest.raises(ValueError):
        audit_supported_phi10_blend_spectral_fingerprint(grid_sizes=(24, 32))
    with pytest.raises(ValueError):
        audit_supported_phi10_blend_spectral_fingerprint(grid_sizes=(24, 32, 40), fixed_shell_cutoff=12)
