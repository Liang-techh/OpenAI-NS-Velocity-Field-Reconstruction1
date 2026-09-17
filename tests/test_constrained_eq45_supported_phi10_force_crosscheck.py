from __future__ import annotations

import math

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_force_crosscheck import (
    FROZEN_SLOPE,
    compare_supported_phi10_temporal_force,
)


def test_frozen_phi10_temporal_force_crosscheck_contract() -> None:
    report = compare_supported_phi10_temporal_force()

    assert report["schema"] == "eq45_supported_phi10_temporal_force_crosscheck_v1"
    assert report["slope_frozen_before_pde_fit"] is True
    assert report["slope"] == FROZEN_SLOPE == 1.4
    assert report["mode"] == {"family": "phi", "index": [1, 0]}
    assert report["midpoint_exact"] is True
    assert report["force_family"] == "preregistered_restricted_two_parameter_family"
    assert report["force_bounds"] == [0.0, 10.0]
    assert report["pressure_fitted"] is False
    assert report["new_force_direction_added"] is False
    assert report["residual_defined_force_allowed"] is False
    assert report["formal_pde_gate_assessed"] is False

    static_fit = report["static_projection"]["fit"]
    temporal_fit = report["temporal_projection"]["fit"]
    for fit in (static_fit, temporal_fit):
        assert 0.0 <= fit["a"] <= 10.0
        assert 0.0 <= fit["c"] <= 10.0
        assert fit["probe_count"] == 16
        assert fit["design_rank"] == 2
        assert math.isfinite(fit["design_condition"])
        assert math.isfinite(fit["before"]["rms"])
        assert math.isfinite(fit["after"]["rms"])

    static_rows = report["static_projection"]["holdout_rows"]
    temporal_rows = report["temporal_projection"]["holdout_rows"]
    assert [row["spatial_step"] for row in static_rows] == [0.02, 0.01, 0.005]
    assert [row["spatial_step"] for row in temporal_rows] == [0.02, 0.01, 0.005]
    for rows in (static_rows, temporal_rows):
        for row in rows:
            assert math.isfinite(row["before"]["rms"])
            assert math.isfinite(row["after"]["rms"])
            assert set(row["by_region"]) == {
                "plateau",
                "radial_collar",
                "axial_collar",
                "corner_collar",
            }

    comparison = report["finest_holdout_comparison"]
    for key, value in comparison.items():
        assert math.isfinite(float(value)), key

    assert report["static_supported_sha256"] != report["temporal_trial_sha256"]
    assert report["visual_correspondence_verified"] is False
    assert report["visualization_ready"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False


def test_crosscheck_rejects_slope_drift() -> None:
    with pytest.raises(ValueError, match="frozen"):
        compare_supported_phi10_temporal_force(slope=1.39)


def test_crosscheck_report_has_nontrivial_temporal_difference() -> None:
    report = compare_supported_phi10_temporal_force()
    comparison = report["finest_holdout_comparison"]
    assert abs(float(comparison["zero_force_fractional_change"])) > 1.0e-6
    assert abs(float(comparison["projected_force_fractional_change"])) > 1.0e-6
    assert not np.isclose(
        report["temporal_projection"]["fit"]["c"],
        report["static_projection"]["fit"]["c"],
        rtol=0.0,
        atol=1.0e-10,
    )
