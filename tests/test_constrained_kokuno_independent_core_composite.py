import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_independent_core_composite import (
    _cylindrical_component_rms,
    generate_core_composite_independent_report,
)


def test_cylindrical_component_projection_is_exact():
    points = np.asarray([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    residual = np.asarray([[3.0, 4.0, 5.0], [-4.0, 3.0, 5.0]])
    metrics = _cylindrical_component_rms(residual, points)
    assert metrics["radial_rms"] == pytest.approx(3.0)
    assert metrics["theta_rms"] == pytest.approx(4.0)
    assert metrics["axial_rms"] == pytest.approx(5.0)


def test_core_composite_report_keeps_gate_fail_closed(tmp_path):
    report = generate_core_composite_independent_report(
        output=tmp_path / "report.json",
        annulus_count=8,
        axis_near_count=4,
        times=(0.5,),
        steps=(0.02, 0.01, 0.005),
    )
    assert report["held_out"]["seed"] == 9172941
    assert report["held_out"]["points_reused_from_agent3_training_or_holdout"] is False
    assert report["upstream_stage_replay"]["selection_rerun_here"] is False
    assert report["upstream_stage_replay"]["mean_correction_upstream_status"].startswith(
        "rejected_by_agent3"
    )
    assert report["support_and_nontriviality"]["oscillatory_rms_on_annulus"] > 0.0
    assert report["support_and_nontriviality"]["mean_correction_rms_on_annulus"] > 0.0
    assert report["support_and_nontriviality"][
        "mean_correction_outside_declared_support_max"
    ] == 0.0
    assert report["mutation_calibration"]["mean_divergence_shift"] == pytest.approx(
        0.02, abs=1.0e-10
    )
    assert report["gate"]["formal_full_domain_gate_assessed"] is False
    assert report["gate"]["formal_full_domain_gate_passed"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["gate"]["local_reference_threshold_met"] is False

    summary = report["finest_step_summary"]
    assert set(summary) == {
        "leading_only",
        "leading_plus_oscillatory",
        "after_rejected_mean_correction",
        "stage_ratios",
    }
    for stage in (
        "leading_only",
        "leading_plus_oscillatory",
        "after_rejected_mean_correction",
    ):
        assert np.isfinite(summary[stage]["mean_residual_rms_over_times"])
        assert summary[stage]["mean_residual_rms_over_times"] > 0.0


def test_report_requires_three_resolution_levels(tmp_path):
    with pytest.raises(ValueError, match="at least three resolutions"):
        generate_core_composite_independent_report(
            output=tmp_path / "report.json",
            annulus_count=8,
            axis_near_count=4,
            times=(0.5,),
            steps=(0.01, 0.005),
        )
