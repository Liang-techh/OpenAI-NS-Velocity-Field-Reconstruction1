from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_early_localized_force_crosscheck import (
    compare_early_localized_phi10_force,
)


@pytest.fixture(scope="module")
def report():
    return compare_early_localized_phi10_force()


def test_crosscheck_reuses_frozen_force_contract_and_known_affine_baseline(report):
    assert report["schema"] == "eq45_supported_phi10_early_localized_force_crosscheck_v1"
    assert report["schedule_frozen_before_pde_fit"] is True
    assert report["early_delta"] == pytest.approx(-1.4)
    assert report["affine_reference_slope"] == pytest.approx(1.4)
    assert report["mode"] == {"family": "phi", "index": [1, 0]}
    assert report["derivative_steps"] == pytest.approx([0.02, 0.01, 0.005])
    assert report["force_family"] == "preregistered_restricted_two_parameter_family"
    assert report["force_bounds"] == pytest.approx([0.0, 10.0])
    assert report["pressure_fitted"] is False
    assert report["new_force_direction_added"] is False
    assert report["residual_defined_force_allowed"] is False
    assert report["temporal_schedule_refit_on_pde"] is False
    assert all(report["endpoint_checks"].values())

    finest = report["finest_holdout_comparison"]
    # Lock the inherited #144 static/affine references so this new comparison
    # cannot silently change the old PDE/force operator while assessing the new shape.
    assert finest["static_zero_force_rms"] == pytest.approx(3.5872881930415295, rel=2e-10)
    assert finest["static_projected_force_rms"] == pytest.approx(3.5860556463766144, rel=2e-10)
    assert finest["affine_zero_force_rms"] == pytest.approx(3.722131026637941, rel=2e-10)
    assert finest["affine_projected_force_rms"] == pytest.approx(3.7208188786928376, rel=2e-10)


def test_quadratic_projection_is_bounded_separate_and_finite(report):
    projection = report["early_localized_projection"]
    fit = projection["fit"]
    assert fit["probe_count"] == 16
    assert projection["holdout_probe_count"] == 16
    assert fit["design_rank"] == 2
    assert math.isfinite(fit["design_condition"])
    assert 0.0 <= fit["a"] <= 10.0
    assert 0.0 <= fit["c"] <= 10.0
    assert len(projection["holdout_rows"]) == 3
    assert [row["spatial_step"] for row in projection["holdout_rows"]] == pytest.approx(
        [0.02, 0.01, 0.005]
    )
    for row in projection["holdout_rows"]:
        for state in ("before", "after"):
            assert math.isfinite(row[state]["rms"])
            assert math.isfinite(row[state]["max"])
            assert row[state]["rms"] > 0.0
            assert row[state]["max"] > 0.0

    finest = report["finest_holdout_comparison"]
    for key, value in finest.items():
        if key != "spatial_step":
            assert math.isfinite(value)
    assert finest["early_localized_zero_force_rms"] > 0.0
    assert finest["early_localized_projected_force_rms"] > 0.0


def test_crosscheck_truth_boundary_stays_fail_closed(report):
    assert report["formal_pde_gate_assessed"] is False
    assert report["visualization_candidate_only"] is True
    for key in (
        "visual_correspondence_verified",
        "visualization_ready",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert report[key] is False


def test_crosscheck_rejects_forcing_convention_drift(tmp_path):
    root = Path(__file__).resolve().parents[1]
    data = json.loads((root / "configs/constraints.json").read_text(encoding="utf-8"))
    data["forcing"]["mode"] = "free_residual_canceller"
    path = tmp_path / "constraints.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(AssertionError, match="forcing convention"):
        compare_early_localized_phi10_force(constraints_path=path)
