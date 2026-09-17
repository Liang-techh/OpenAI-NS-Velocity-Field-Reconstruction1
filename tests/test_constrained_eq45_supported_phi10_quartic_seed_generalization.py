import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_seed_generalization import (
    DEFAULT_SEEDS,
    REGIONS,
    audit_quartic_phi10_seed_generalization,
    held_out_probe_set,
)


STATIC_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
CUBIC_SHA = "310fc2ad1bf3d721860369b9a58c7599202182964158708867691f4445f21576"
QUARTIC_SHA = "03fdae73170469ae1160297b489b1b83a27adddfc941119e116a21b98bd15049"


def test_fresh_probe_set_is_deterministic_stratified_and_off_grid():
    points_a, times_a, labels_a = held_out_probe_set(DEFAULT_SEEDS[0], 2)
    points_b, times_b, labels_b = held_out_probe_set(DEFAULT_SEEDS[0], 2)
    np.testing.assert_array_equal(points_a, points_b)
    np.testing.assert_array_equal(times_a, times_b)
    np.testing.assert_array_equal(labels_a, labels_b)
    assert points_a.shape == (8, 3)
    assert times_a.shape == (8,)
    assert set(labels_a.tolist()) == set(REGIONS)
    assert np.all((times_a >= 0.34) & (times_a <= 0.66))
    # Random polar draws should not accidentally collapse to the old fixed probe grid.
    assert np.all(np.abs(points_a / 0.005 - np.round(points_a / 0.005)) > 1.0e-9)


def test_seed_contract_rejects_malformed_requests():
    with pytest.raises(ValueError, match="at least three"):
        audit_quartic_phi10_seed_generalization(seeds=(1, 2), points_per_region=1)
    with pytest.raises(ValueError, match="unique"):
        audit_quartic_phi10_seed_generalization(seeds=(1, 1, 2), points_per_region=1)
    with pytest.raises(ValueError, match="positive integer"):
        held_out_probe_set(1, 0)


def test_quartic_generalization_freezes_upstream_identity_force_and_truth_boundary():
    report = audit_quartic_phi10_seed_generalization(points_per_region=1)
    assert report["task_id"] == "CR009-EQ45-SUPPORTED-PHI10-QUARTIC-SEED-GENERALIZATION-028"
    assert report["validation_seeds"] == list(DEFAULT_SEEDS)
    assert report["seed_count"] == 3
    assert report["probe_count_per_seed"] == 4
    assert report["spatial_steps"] == [0.02, 0.01, 0.005]
    assert report["fixed_time_step"] == 0.0025
    assert report["static_supported_sha256"] == STATIC_SHA
    assert report["cubic_trial_sha256"] == CUBIC_SHA
    assert report["quartic_trial_sha256"] == QUARTIC_SHA
    assert report["candidates_serialized_and_reloaded_before_validation"] is True
    assert report["schedule_frozen_before_validation"] is True
    assert report["force_coefficients_frozen_before_validation"] is True
    assert report["fresh_seed_force_refit"] is False
    assert report["pressure_fitted"] is False
    assert report["new_force_direction_added"] is False
    assert report["residual_defined_force_allowed"] is False
    assert report["formal_pde_gate_assessed"] is False
    assert report["visualization_candidate_only"] is True
    assert report["canonical_velocity_changed"] is False
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False

    np.testing.assert_allclose(
        [report["static_frozen_force"]["a"], report["static_frozen_force"]["c"]],
        [4.841316873036678e-22, 0.0766923261606361],
        rtol=1.0e-10,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        [report["cubic_frozen_force"]["a"], report["cubic_frozen_force"]["c"]],
        [5.070373823229265e-16, 0.07802699746086716],
        rtol=1.0e-10,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        [report["quartic_frozen_force"]["a"], report["quartic_frozen_force"]["c"]],
        [2.1615475834844833e-15, 0.07770084654400304],
        rtol=1.0e-10,
        atol=1.0e-12,
    )

    assert len(report["seed_reports"]) == 3
    assert len(report["finest_rows"]) == 3
    for seed_report in report["seed_reports"]:
        assert len(seed_report["levels"]) == 3
        for level in seed_report["levels"]:
            assert level["spatial_step"] in (0.02, 0.01, 0.005)
            for field in ("static", "cubic", "quartic"):
                for force_state in ("zero", "frozen_force"):
                    metrics = level[field][force_state]
                    assert metrics["max"] >= 0.0
                    assert metrics["rms"] >= 0.0
                    assert metrics["term_normalized_rms"] >= 0.0
                    assert np.all(np.isfinite(list(metrics.values())))

    summary = report["finest_cross_seed_summary"]
    assert 0 <= summary["quartic_beats_cubic_zero_seed_count"] <= 3
    assert 0 <= summary["quartic_beats_cubic_forced_seed_count"] <= 3
    for key in (
        "quartic_vs_static_zero_fractional_change",
        "quartic_vs_cubic_zero_fractional_change",
        "quartic_vs_static_forced_fractional_change",
        "quartic_vs_cubic_forced_fractional_change",
    ):
        assert np.all(np.isfinite(list(summary[key].values())))
