import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_force_crosscheck import (
    FROZEN_EARLY_DELTA,
    FROZEN_NULLSPACE_COEFFICIENT,
    compare_quartic_balanced_phi10_force,
)


def _report():
    return compare_quartic_balanced_phi10_force()


def test_quartic_force_crosscheck_freezes_velocity_schedule_before_pde_fit():
    report = _report()

    assert report["task_id"] == "CR005-EQ45-SUPPORTED-PHI10-QUARTIC-FORCE-CROSSCHECK-027"
    assert report["velocity_schedule_frozen_before_pde_fit"] is True
    assert report["early_delta"] == pytest.approx(FROZEN_EARLY_DELTA, rel=0.0, abs=0.0)
    assert report["nullspace_coefficient"] == pytest.approx(
        FROZEN_NULLSPACE_COEFFICIENT, rel=0.0, abs=0.0
    )
    assert report["nullspace_coefficient"] == pytest.approx(
        -0.2831460674157303, rel=0.0, abs=1.0e-14
    )
    assert report["pde_objective_used_to_choose_temporal_schedule"] is False
    assert report["temporal_schedule_refit_on_pde"] is False
    assert report["mode"] == {"family": "phi", "index": [1, 0]}
    assert all(report["public_identity_checks"].values())
    assert report["static_supported_sha256"] != report["quartic_trial_sha256"]
    assert report["cubic_trial_sha256"] != report["quartic_trial_sha256"]


def test_quartic_force_crosscheck_keeps_preregistered_fit_holdout_contract():
    report = _report()

    assert report["nu"] == pytest.approx(0.01)
    assert report["derivative_steps"] == pytest.approx([0.02, 0.01, 0.005])
    assert report["fit_and_holdout_separate"] is True
    assert report["training_randomness_used"] is False
    assert report["force_family"] == "preregistered_restricted_two_parameter_family"
    assert report["force_bounds"] == pytest.approx([0.0, 10.0])
    assert report["pressure_fitted"] is False
    assert report["new_force_direction_added"] is False
    assert report["residual_defined_force_allowed"] is False

    for key in ("cubic_projection", "quartic_projection"):
        projection = report[key]
        fit = projection["fit"]
        assert 0.0 <= fit["a"] <= 10.0
        assert 0.0 <= fit["c"] <= 10.0
        assert fit["probe_count"] == 16
        assert projection["holdout_probe_count"] == 16
        assert len(projection["holdout_rows"]) == 3
        for row, step in zip(projection["holdout_rows"], [0.02, 0.01, 0.005]):
            assert row["spatial_step"] == pytest.approx(step)
            for state in ("before", "after"):
                assert np.isfinite(row[state]["rms"])
                assert np.isfinite(row[state]["max"])
                assert row[state]["rms"] > 0.0
                assert row[state]["max"] > 0.0


def test_quartic_force_crosscheck_reports_small_force_capacity_without_pde_promotion():
    report = _report()
    comparison = report["finest_holdout_comparison"]
    quartic = report["quartic_projection"]

    # This family has repeatedly shown only tiny global capacity on the supported
    # candidate.  Keep a generous guard that catches accidental force-family or
    # operator changes without turning the diagnostic into a PDE acceptance gate.
    assert abs(comparison["quartic_force_reduction_fraction"]) < 0.01
    assert quartic["fit"]["rms_reduction_fraction"] < 0.01

    assert report["formal_pde_gate_assessed"] is False
    assert report["visualization_candidate_only"] is True
    assert report["canonical_velocity_changed"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["visualization_ready"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False
