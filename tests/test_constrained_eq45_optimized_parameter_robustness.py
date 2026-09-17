import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_optimized_parameter_robustness import (
    DEFAULT_VALIDATION_SEEDS,
    DIRECTIONS,
    audit_optimized_eq45_parameter_robustness,
)


def test_optimized_parameter_robustness_is_held_out_bounded_and_public():
    report = audit_optimized_eq45_parameter_robustness(sample_count=4)

    assert report["schema"] == "eq45_optimized_parameter_robustness_v1"
    assert report["task_id"] == "CR009-EQ45-OPTIMIZED-PARAM-PERTURB-VORTICITY-016"
    assert report["source_optimization_task_id"] == "CR005-EQ45-PROFILE-FORCE-OPT-015"
    assert report["optimized_candidate_sha256"] == (
        "b34902ee28b3965f245eceb5867df7409498063cb0ea80f1758931723b303c98"
    )
    assert report["validation_seeds"] == list(DEFAULT_VALIDATION_SEEDS)
    assert report["perturbation_amplitudes"] == [0.02, 0.04]
    assert set(report["perturbation_directions"]) == set(DIRECTIONS)
    assert report["coefficient_bound"] == [-4.0, 4.0]
    assert report["spatial_steps"] == [0.02, 0.01, 0.005]
    assert report["fixed_time_step"] == 0.0025
    assert report["frozen_force"]["refit_on_validation"] is False

    assert len(report["seed_reports"]) == 3
    for seed_report in report["seed_reports"]:
        assert len(seed_report["baseline_levels"]) == 3
        assert [row["spatial_step"] for row in seed_report["baseline_levels"]] == [
            0.02,
            0.01,
            0.005,
        ]
        assert len(seed_report["perturbations"]) == 6
        identities = set()
        for perturbation in seed_report["perturbations"]:
            identities.add((perturbation["direction"], perturbation["amplitude"]))
            assert -4.0 <= perturbation["Phi_02"] <= 4.0
            assert -4.0 <= perturbation["F_02"] <= 4.0
            assert len(perturbation["candidate_sha256"]) == 64
            response = perturbation["velocity_response"]
            assert response["baseline_velocity_rms"] > 0.0
            assert response["delta_velocity_rms"] > 0.0
            assert response["delta_velocity_relative_rms"] > 0.0
            assert np.isfinite(response["delta_velocity_max"])
            assert len(perturbation["levels"]) == 3
            for level in perturbation["levels"]:
                values = np.array(
                    [
                        level["perturbed"]["zero_force_max"],
                        level["perturbed"]["zero_force_rms"],
                        level["perturbed"]["zero_force_relative_rms"],
                        level["perturbed"]["restricted_force_max"],
                        level["perturbed"]["restricted_force_rms"],
                        level["perturbed"]["restricted_force_relative_rms"],
                        level["zero_force_rms_relative_change"],
                        level["frozen_force_rms_relative_change"],
                    ],
                    dtype=float,
                )
                assert np.all(np.isfinite(values))
                assert np.all(values[:6] >= 0.0)
        assert identities == {
            ("phi_inward", 0.02),
            ("phi_inward", 0.04),
            ("F_inward", 0.02),
            ("F_inward", 0.04),
            ("joint_inward", 0.02),
            ("joint_inward", 0.04),
        }

    for summary in report["finest_level_summary"].values():
        assert set(summary) == {"min", "median", "mean", "max", "population_std"}
        assert np.all(np.isfinite(np.array(list(summary.values()), dtype=float)))

    assert report["training_performed_here"] is False
    assert report["pressure_fitted_here"] is False
    assert report["force_refit_here"] is False
    assert report["velocity_changed_by_audit"] is False
    assert report["visualization_candidate_only"] is True
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False


@pytest.mark.parametrize(
    "kwargs",
    [
        {"seeds": (1, 2)},
        {"seeds": (1, 1, 2)},
        {"seeds": (1, 2, -3)},
        {"seeds": (914307, 2, 3)},
        {"seeds": (914117, 2, 3)},
        {"amplitudes": (0.02,)},
        {"amplitudes": (0.02, 0.02)},
        {"amplitudes": (0.04, 0.02)},
        {"amplitudes": (0.02, np.nan)},
        {"amplitudes": (0.02, 20.0)},
    ],
)
def test_optimized_parameter_robustness_fails_closed(kwargs):
    with pytest.raises((ValueError, RuntimeError)):
        audit_optimized_eq45_parameter_robustness(sample_count=4, **kwargs)
