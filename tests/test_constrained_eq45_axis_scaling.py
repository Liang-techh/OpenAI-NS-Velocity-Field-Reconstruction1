import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.constrained_eq45_axis_scaling import (
    audit_eq45_near_axis_scaling,
)
from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE_PATH = ROOT / "artifacts/constrained/eq45_velocity_candidate_seed.json"
AUDIT_PATH = ROOT / "artifacts/constrained/eq45_near_axis_scaling_audit_seed.json"


def _compare_rows(actual, expected):
    assert len(actual) == len(expected)
    for got, want in zip(actual, expected):
        assert got["time"] == want["time"]
        assert got["z"] == want["z"]
        assert np.allclose(got["axis_velocity"], want["axis_velocity"], rtol=5e-10, atol=1e-13)
        assert np.allclose(got["transverse_orders"], want["transverse_orders"], rtol=5e-10, atol=1e-13)
        assert np.allclose(
            got["axial_departure_orders"], want["axial_departure_orders"], rtol=5e-10, atol=1e-13
        )
        assert np.isclose(
            got["max_transverse_over_r_azimuthal_spread"],
            want["max_transverse_over_r_azimuthal_spread"],
            rtol=0.0,
            atol=5e-14,
        )
        assert np.isclose(
            got["max_axial_azimuthal_spread"],
            want["max_axial_azimuthal_spread"],
            rtol=0.0,
            atol=5e-14,
        )


def test_frozen_eq45_near_axis_scaling_replays_checked_artifact():
    field = Eq45VelocityCandidate.load_json(CANDIDATE_PATH)
    expected = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))

    assert expected["schema"] == "eq45_near_axis_scaling_audit_v1"
    assert expected["candidate_sha256"] == field.sha256
    sampling = expected["sampling"]
    actual = audit_eq45_near_axis_scaling(
        field.at_points,
        times=sampling["times"],
        z_levels=sampling["z_levels"],
        radii=sampling["radii"],
        azimuth_count=sampling["azimuth_count"],
    )
    _compare_rows(actual["rows"], expected["rows"])
    for key, value in expected["summary"].items():
        assert np.isclose(actual["summary"][key], value, rtol=5e-10, atol=5e-14)

    summary = actual["summary"]
    # Current-seed regression observations, not universal physical thresholds.
    assert summary["min_finest_transverse_order"] > 0.998
    assert summary["min_finest_axial_departure_order"] > 1.999
    assert summary["max_axis_transverse_speed"] < 1e-14
    assert summary["max_transverse_over_r_azimuthal_spread"] < 1e-12
    assert summary["max_axial_azimuthal_spread"] < 1e-12

    truth = expected["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["physical_support_validated"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert truth["blowup_proved"] is False


def test_near_axis_scaling_detects_constant_transverse_contamination():
    field = Eq45VelocityCandidate.load_json(CANDIDATE_PATH)

    def contaminated(points, time):
        values = field.at_points(points, time).copy()
        values[..., 0] += 0.02
        return values

    report = audit_eq45_near_axis_scaling(contaminated)
    assert report["summary"]["max_axis_transverse_speed"] > 0.019
    assert report["summary"]["min_finest_transverse_order"] < 0.5


def test_near_axis_scaling_fail_closed_inputs():
    field = Eq45VelocityCandidate.load_json(CANDIDATE_PATH)
    with np.testing.assert_raises(ValueError):
        audit_eq45_near_axis_scaling(field.at_points, radii=(0.01, 0.02, 0.04))
    with np.testing.assert_raises(ValueError):
        audit_eq45_near_axis_scaling(field.at_points, azimuth_count=3)
