import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_blend_axial_taper_force_screen import (
    FROZEN_AXIAL_PLATEAU_Q,
    FROZEN_BLEND_WEIGHT,
    FROZEN_IDENTITY_HALF_HEIGHTS,
    TASK_ID,
    screen_blend_axial_taper_force,
)


def test_axial_taper_force_screen_keeps_selection_and_force_contracts_separate():
    report = screen_blend_axial_taper_force()

    assert report["task_id"] == TASK_ID
    assert report["blend_weight"] == FROZEN_BLEND_WEIGHT
    assert report["blend_weight_selected_by_this_module"] is False
    assert report["axial_identity_half_heights"] == list(FROZEN_IDENTITY_HALF_HEIGHTS)
    assert report["axial_plateau_q_values"] == pytest.approx(
        list(FROZEN_AXIAL_PLATEAU_Q), rel=0.0, abs=1.0e-15
    )
    assert report["axial_taper_value_selected"] is False
    assert report["axial_taper_selection_rule"] is None
    assert report["velocity_family_frozen_before_pde_fit"] is True
    assert report["fit_and_holdout_separate"] is True
    assert report["holdout_force_refit"] is False
    assert report["force_bounds"] == [0.0, 10.0]
    assert report["pressure_fitted"] is False
    assert report["new_force_direction_added"] is False
    assert report["residual_defined_force_allowed"] is False
    assert report["taper_refit_on_pde"] is False
    assert report["public_image_fitted"] is False

    rows = report["rows"]
    assert len(rows) == 5
    assert [row["axial_identity_half_height"] for row in rows] == list(
        FROZEN_IDENTITY_HALF_HEIGHTS
    )
    assert [row["axial_plateau_q"] for row in rows] == pytest.approx(
        list(FROZEN_AXIAL_PLATEAU_Q), rel=0.0, abs=1.0e-15
    )
    for row in rows:
        assert 0.0 <= row["fit"]["a"] <= 10.0
        assert 0.0 <= row["fit"]["c"] <= 10.0
        assert row["fit"]["probe_count"] == 16
        assert len(row["holdout_rows"]) == 3
        assert [entry["spatial_step"] for entry in row["holdout_rows"]] == pytest.approx(
            report["derivative_steps"], rel=0.0, abs=0.0
        )
        assert np.isfinite(row["finest_zero_force_rms"])
        assert np.isfinite(row["finest_projected_force_rms"])
        assert set(row["finest_by_region"]) == {
            "plateau",
            "radial_collar",
            "axial_collar",
            "corner_collar",
        }

    # q=.64 must exactly replay the already-screened lambda=.5 blend PDE lane.
    assert rows[0]["finest_zero_force_rms"] == pytest.approx(3.851129051, rel=3.0e-6)
    assert rows[0]["finest_projected_force_rms"] == pytest.approx(3.849884649, rel=3.0e-6)
    assert rows[0]["public_axial_probe_velocity_delta_rms_vs_q064"] == pytest.approx(0.0)
    assert rows[0]["public_axial_probe_velocity_relative_delta_vs_q064"] == pytest.approx(0.0)
    assert rows[-1]["public_axial_probe_velocity_delta_rms_vs_q064"] > 0.0

    export = report["direct_velocity_export_check"]
    assert export["checked_axial_identity_half_height"] == pytest.approx(1.70)
    assert export["checked_axial_plateau_q"] == pytest.approx(0.7225)
    assert export["serialization_roundtrip_exact"] is True
    assert export["grid_shape"] == [3, 3, 3, 3, 3]
    assert export["grid_finite"] is True
    assert export["grid_nonzero"] is True

    assert report["formal_pde_gate_assessed"] is False
    assert report["visualization_candidate_only"] is True
    assert report["canonical_velocity_changed"] is False
    assert report["production_axial_taper_promoted"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["visualization_ready"] is False
    assert report["physical_support_validated"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
