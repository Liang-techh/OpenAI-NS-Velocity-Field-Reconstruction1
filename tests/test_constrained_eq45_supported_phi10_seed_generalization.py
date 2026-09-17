import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_seed_generalization import (
    DEFAULT_SEEDS,
    REGIONS,
    audit_supported_phi10_temporal_seed_generalization,
    held_out_probe_set,
)


def test_stratified_probe_draw_is_deterministic_and_inside_contract():
    points, times, labels = held_out_probe_set(914401, 2)
    again_points, again_times, again_labels = held_out_probe_set(914401, 2)
    assert np.array_equal(points, again_points)
    assert np.array_equal(times, again_times)
    assert np.array_equal(labels, again_labels)
    assert points.shape == (8, 3)
    assert set(labels.tolist()) == set(REGIONS)
    assert np.all((times >= 0.34) & (times <= 0.66))
    radius = np.hypot(points[:, 0], points[:, 1])
    assert np.all(radius < 1.79)
    assert np.all(np.abs(points[:, 2]) < 1.79)


def test_phi10_temporal_tradeoff_generalizes_without_validation_refit():
    report = audit_supported_phi10_temporal_seed_generalization(
        seeds=DEFAULT_SEEDS,
        points_per_region=1,
    )
    assert report["schema"] == "eq45_supported_phi10_temporal_seed_generalization_v1"
    assert report["validation_seeds"] == list(DEFAULT_SEEDS)
    assert report["seed_count"] == 3
    assert report["probe_count_per_seed"] == 4
    assert report["spatial_steps"] == [0.02, 0.01, 0.005]
    assert report["slope"] == 1.4
    assert report["slope_frozen_before_validation"] is True
    assert report["force_coefficients_frozen_before_validation"] is True
    assert report["force_fitted_on_generalization_seeds"] is False
    assert report["pressure_fitted"] is False
    assert report["new_force_direction_added"] is False
    assert report["velocity_access"] == "public_at_points_only"
    assert report["static_supported_sha256"] != report["temporal_trial_sha256"]

    assert len(report["seed_reports"]) == 3
    assert len(report["finest_rows"]) == 3
    for seed_report in report["seed_reports"]:
        assert len(seed_report["levels"]) == 3
        assert seed_report["velocity_relative_rms_change"] > 0.0
        for level in seed_report["levels"]:
            for key in (
                "static_zero",
                "temporal_zero",
                "static_frozen_force",
                "temporal_frozen_force",
            ):
                metrics = level[key]
                assert np.isfinite(metrics["max"]) and metrics["max"] >= 0.0
                assert np.isfinite(metrics["rms"]) and metrics["rms"] >= 0.0
                assert np.isfinite(metrics["term_normalized_rms"])
                assert metrics["term_normalized_rms"] >= 0.0
            assert np.isfinite(level["zero_fractional_change"])
            assert np.isfinite(level["forced_fractional_change"])
            assert set(level["by_region"]["static_zero"]) == set(REGIONS)

    summary = report["cross_seed_summary"]
    assert 0.0 <= summary["zero_improved_seed_fraction"] <= 1.0
    assert 0.0 <= summary["forced_improved_seed_fraction"] <= 1.0
    assert summary["velocity_relative_rms_change"]["min"] > 0.0
    assert report["formal_pde_gate_assessed"] is False
    assert report["visualization_candidate_only"] is True
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False


def test_seed_contract_fails_closed():
    with pytest.raises(ValueError, match="at least three"):
        audit_supported_phi10_temporal_seed_generalization(seeds=(1, 2), points_per_region=1)
    with pytest.raises(ValueError, match="unique"):
        audit_supported_phi10_temporal_seed_generalization(seeds=(1, 1, 2), points_per_region=1)
    with pytest.raises(ValueError, match="positive integer"):
        held_out_probe_set(1, 0)
