import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_compact_quartic_blend_force_screen import (
    FROZEN_BLEND_WEIGHTS,
    TASK_ID,
    screen_compact_quartic_blend_force,
)


def test_blend_force_screen_keeps_selection_and_force_contracts_separate():
    report = screen_compact_quartic_blend_force()

    assert report["task_id"] == TASK_ID
    assert report["blend_weights"] == list(FROZEN_BLEND_WEIGHTS)
    assert report["blend_weight_selected"] is False
    assert report["blend_weight_selection_rule"] is None
    assert report["velocity_schedule_frozen_before_pde_fit"] is True
    assert report["fit_and_holdout_separate"] is True
    assert report["holdout_force_refit"] is False
    assert report["force_bounds"] == [0.0, 10.0]
    assert report["pressure_fitted"] is False
    assert report["new_force_direction_added"] is False
    assert report["residual_defined_force_allowed"] is False
    assert report["temporal_schedule_refit_on_pde"] is False
    assert report["public_image_fitted"] is False

    assert len(report["rows"]) == 5
    assert [row["blend_weight"] for row in report["rows"]] == list(FROZEN_BLEND_WEIGHTS)
    for row in report["rows"]:
        assert 0.0 <= row["fit"]["a"] <= 10.0
        assert 0.0 <= row["fit"]["c"] <= 10.0
        assert row["fit"]["probe_count"] == 16
        assert len(row["holdout_rows"]) == 3
        assert [entry["spatial_step"] for entry in row["holdout_rows"]] == pytest.approx(
            report["derivative_steps"], rel=0.0, abs=0.0
        )
        assert np.isfinite(row["finest_zero_force_rms"])
        assert np.isfinite(row["finest_projected_force_rms"])

    # These endpoint numbers are the already-independent #169/#178 diagnostics;
    # they guard against accidentally changing the PDE/force operator while adding
    # the interpolation screen.
    assert report["rows"][0]["finest_zero_force_rms"] == pytest.approx(
        3.7016893, rel=2.0e-6
    )
    assert report["rows"][0]["finest_projected_force_rms"] == pytest.approx(
        3.7004148, rel=2.0e-6
    )
    assert report["rows"][-1]["finest_zero_force_rms"] == pytest.approx(
        4.0299396, rel=2.0e-6
    )
    assert report["rows"][-1]["finest_projected_force_rms"] == pytest.approx(
        4.0287319, rel=2.0e-6
    )

    export = report["direct_velocity_export_check"]
    assert export["checked_blend_weight"] == 0.5
    assert export["serialization_roundtrip_exact"] is True
    assert export["grid_shape"] == [3, 3, 3, 3, 3]
    assert export["grid_finite"] is True
    assert export["grid_nonzero"] is True

    assert report["formal_pde_gate_assessed"] is False
    assert report["visualization_candidate_only"] is True
    assert report["canonical_velocity_changed"] is False
    assert report["production_temporal_shape_promoted"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["visualization_ready"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
