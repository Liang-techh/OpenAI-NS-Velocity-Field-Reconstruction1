import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_c2_compact_force_crosscheck import (
    FROZEN_EARLY_DELTA,
    TASK_ID,
    compare_c2_compact_phi10_force,
)


def _report():
    return compare_c2_compact_phi10_force()


def test_c2_compact_force_crosscheck_freezes_target_free_schedule_before_pde_fit():
    report = _report()

    assert report["task_id"] == TASK_ID
    assert report["velocity_schedule_frozen_before_pde_fit"] is True
    assert report["early_delta"] == pytest.approx(FROZEN_EARLY_DELTA, rel=0.0, abs=0.0)
    assert report["window_tau_length"] == pytest.approx(0.6760634706279541, rel=1.0e-12)
    assert report["return_tau"] == pytest.approx(-0.32393652937204587, rel=1.0e-12)
    assert report["return_time"] == pytest.approx(0.4190158676569885, rel=1.0e-12)
    assert report["pde_objective_used_to_choose_temporal_schedule"] is False
    assert report["temporal_schedule_refit_on_pde"] is False
    assert all(report["public_identity_checks"].values())


def test_c2_compact_force_crosscheck_keeps_preregistered_force_and_disjoint_contract():
    report = _report()

    assert report["nu"] == pytest.approx(0.01)
    assert report["derivative_steps"] == pytest.approx([0.02, 0.01, 0.005])
    assert report["fit_and_holdout_separate"] is True
    assert report["holdout_force_refit"] is False
    assert report["force_family"] == "preregistered_restricted_two_parameter_family"
    assert report["force_bounds"] == pytest.approx([0.0, 10.0])
    assert report["pressure_fitted"] is False
    assert report["new_force_direction_added"] is False
    assert report["residual_defined_force_allowed"] is False

    projection = report["compact_projection"]
    fit = projection["fit"]
    assert fit["probe_count"] == 16
    assert fit["design_rank"] == 2
    assert 0.0 <= fit["a"] <= 10.0
    assert 0.0 <= fit["c"] <= 10.0
    assert fit["after"]["rms"] <= fit["before"]["rms"] + 1.0e-12
    assert len(projection["holdout_rows"]) == 3

    for row in projection["holdout_rows"]:
        for state in ("before", "after"):
            assert np.isfinite(row[state]["rms"])
            assert np.isfinite(row[state]["max"])
            assert row[state]["rms"] > 0.0
            assert row[state]["max"] > 0.0


def test_c2_compact_force_crosscheck_reports_comparison_without_scientific_promotion():
    report = _report()
    comparison = report["finest_holdout_comparison"]

    assert comparison["spatial_step"] == pytest.approx(0.005)
    for key, value in comparison.items():
        if key != "spatial_step":
            assert np.isfinite(value)

    assert report["static_supported_sha256"] != report["compact_screen_sha256"]
    assert report["quartic_trial_sha256"] != report["compact_screen_sha256"]
    assert report["formal_pde_gate_assessed"] is False
    assert report["visualization_candidate_only"] is True
    assert report["canonical_velocity_changed"] is False
    assert report["production_temporal_shape_promoted"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["visualization_ready"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False
