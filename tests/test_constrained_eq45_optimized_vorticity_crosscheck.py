import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_optimized_vorticity_crosscheck import (
    DEFAULT_VALIDATION_SEEDS,
    audit_optimized_eq45_vorticity,
)


def test_optimized_profile_crosscheck_uses_frozen_artifacts_and_public_field():
    report = audit_optimized_eq45_vorticity(sample_count=4)

    assert report["schema"] == "eq45_optimized_profile_vorticity_crosscheck_v1"
    assert report["task_id"] == "CR009-EQ45-OPTIMIZED-PROFILE-VORTICITY-CROSSCHECK-015"
    assert report["source_optimization_task_id"] == "CR005-EQ45-PROFILE-FORCE-OPT-015"
    assert report["baseline_candidate_sha256"] == (
        "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
    )
    assert report["optimized_candidate_sha256"] == (
        "b34902ee28b3965f245eceb5867df7409498063cb0ea80f1758931723b303c98"
    )
    assert report["optimized_candidate_roundtrip_equal"] is True
    assert report["validation_seeds"] == list(DEFAULT_VALIDATION_SEEDS)
    assert report["spatial_steps"] == [0.02, 0.01, 0.005]
    assert report["fixed_time_step"] == 0.0025
    assert report["profile_parameters"] == {
        "Phi_02": pytest.approx(3.8740984117202455),
        "F_02": pytest.approx(3.9916795387937536),
    }
    assert report["frozen_force_choices"]["fitted_in_this_audit"] is False

    assert len(report["seed_reports"]) == 3
    for seed_report in report["seed_reports"]:
        assert len(seed_report["levels"]) == 3
        steps = [level["spatial_step"] for level in seed_report["levels"]]
        assert steps == [0.02, 0.01, 0.005]
        for level in seed_report["levels"]:
            for family in ("baseline", "optimized"):
                metrics = level[family]
                values = np.array(
                    [
                        metrics["zero_force_max"],
                        metrics["zero_force_rms"],
                        metrics["zero_force_relative_rms"],
                        metrics["restricted_force_max"],
                        metrics["restricted_force_rms"],
                        metrics["restricted_force_relative_rms"],
                    ]
                )
                assert np.all(np.isfinite(values))
                assert np.all(values >= 0.0)
            assert np.isfinite(level["profile_only_zero_force_rms_reduction_fraction"])
            assert np.isfinite(
                level["profile_plus_frozen_force_rms_reduction_fraction"]
            )

    finest = [seed_report["levels"][-1] for seed_report in report["seed_reports"]]
    observed_zero = np.array(
        [row["profile_only_zero_force_rms_reduction_fraction"] for row in finest]
    )
    observed_restricted = np.array(
        [row["profile_plus_frozen_force_rms_reduction_fraction"] for row in finest]
    )
    summary = report["finest_level_summary"]
    assert summary["profile_only_zero_force_rms_reduction_fraction"]["mean"] == pytest.approx(
        float(np.mean(observed_zero))
    )
    assert summary["profile_plus_frozen_force_rms_reduction_fraction"]["mean"] == pytest.approx(
        float(np.mean(observed_restricted))
    )

    assert report["velocity_changed_by_source_optimization"] is True
    assert report["velocity_changed_by_this_validator"] is False
    assert report["training_or_force_fit_performed_here"] is False
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False


@pytest.mark.parametrize(
    "seeds",
    [
        (1, 2),
        (1, 1, 2),
        (1, 2, -3),
        (1, 2, 3.5),
        (914117, 2, 3),
        (914211, 2, 3),
    ],
)
def test_optimized_profile_crosscheck_seed_contract_fails_closed(seeds):
    with pytest.raises(ValueError):
        audit_optimized_eq45_vorticity(seeds=seeds, sample_count=4)
