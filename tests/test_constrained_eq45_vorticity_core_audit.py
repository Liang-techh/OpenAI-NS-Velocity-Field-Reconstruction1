import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_vorticity_core_audit import (
    _cell_centres,
    audit_eq45_vorticity_core,
    centered_vorticity,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE_PATH = ROOT / "artifacts/constrained/eq45_velocity_candidate_seed.json"
REPORT_PATH = ROOT / "artifacts/constrained/eq45_vorticity_core_audit_seed.json"
EXPECTED_SHA256 = "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"


def _grid(grid_size=17, box_half_width=2.0):
    axis, dx = _cell_centres(grid_size, box_half_width)
    return (*np.meshgrid(axis, axis, axis, indexing="ij"), dx)


def test_centered_vorticity_catches_curl_sign_and_zero_calibrations():
    x, y, z, dx = _grid()

    rigid_rotation = np.stack((-y, x, np.zeros_like(z)), axis=-1)
    omega = centered_vorticity(rigid_rotation, dx)
    np.testing.assert_allclose(omega[..., 0], 0.0, atol=1e-13)
    np.testing.assert_allclose(omega[..., 1], 0.0, atol=1e-13)
    np.testing.assert_allclose(omega[..., 2], 2.0, atol=1e-13)

    irrotational_linear = np.stack((x, y, z), axis=-1)
    np.testing.assert_allclose(centered_vorticity(irrotational_linear, dx), 0.0, atol=1e-13)


def test_vorticity_audit_requires_three_increasing_resolution_levels():
    candidate = Eq45VelocityCandidate.load_json(CANDIDATE_PATH)
    with pytest.raises(ValueError, match="at least three"):
        audit_eq45_vorticity_core(candidate, grid_sizes=(16, 24), times=(0.5,))
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_eq45_vorticity_core(candidate, grid_sizes=(16, 24, 24), times=(0.5,))


def test_frozen_eq45_report_replays_public_velocity_and_bulk_core_is_stable():
    candidate = Eq45VelocityCandidate.load_json(CANDIDATE_PATH)
    frozen = json.loads(REPORT_PATH.read_text(encoding="utf-8"))
    replay = audit_eq45_vorticity_core(candidate)

    assert candidate.sha256 == EXPECTED_SHA256
    assert replay["candidate_sha256"] == EXPECTED_SHA256
    assert replay["grid_sizes"] == [48, 64, 80]
    assert replay["times"] == [0.25, 0.5, 0.75]
    assert replay["truth_boundary"] == frozen["truth_boundary"]

    keys = (
        "omega_max",
        "enstrophy_integral",
        "enstrophy_radial_rms",
        "enstrophy_axial_rms",
        "enstrophy_rms_aspect_ratio",
        "enstrophy_radial_q90",
        "enstrophy_abs_z_q90",
        "halfmax_radial_extent",
        "halfmax_abs_z_extent",
    )
    for actual, expected in zip(replay["rows"], frozen["rows"]):
        assert actual["grid_size"] == expected["grid_size"]
        assert actual["time"] == expected["time"]
        for key in keys:
            np.testing.assert_allclose(actual[key], expected[key], rtol=2e-10, atol=2e-12)

    summaries = replay["time_summaries"]
    assert [item["time"] for item in summaries] == [0.25, 0.5, 0.75]
    # Bulk vorticity-core second moments are substantially less grid-sensitive
    # than a peak or a plane-quantized axial q90.  These are visualization
    # stability checks, not PDE acceptance thresholds.
    assert max(item["mid_to_finest_radial_rms_relative_delta"] for item in summaries) < 0.001
    assert max(item["mid_to_finest_axial_rms_relative_delta"] for item in summaries) < 0.011
    assert max(item["mid_to_finest_aspect_relative_delta"] for item in summaries) < 0.010
    assert max(item["mid_to_finest_omega_max_relative_delta"] for item in summaries) < 0.042
    assert max(item["mid_to_finest_enstrophy_relative_delta"] for item in summaries) < 0.031

    finest = [row for row in replay["rows"] if row["grid_size"] == 80]
    radial = [row["enstrophy_radial_rms"] for row in finest]
    axial = [row["enstrophy_axial_rms"] for row in finest]
    peaks = [row["omega_max"] for row in finest]
    enstrophy = [row["enstrophy_integral"] for row in finest]
    assert radial[0] > radial[1] > radial[2]
    assert axial[0] > axial[1] > axial[2]
    assert peaks[0] < peaks[1] < peaks[2]
    assert enstrophy[0] < enstrophy[1] < enstrophy[2]
    assert all(0.24 < row["enstrophy_rms_aspect_ratio"] < 0.26 for row in finest)

    assert replay["truth_boundary"]["velocity_changed"] is False
    assert replay["truth_boundary"]["visualization_ready"] is False
    assert replay["truth_boundary"]["pde_validated"] is False
    assert replay["truth_boundary"]["paper_exact"] is False
    assert replay["truth_boundary"]["openai_field_identified"] is False
