import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_missing_covariance_column_target import (
    build_missing_covariance_column_target,
)
from openai_ns_reconstruction.kokuno_multislice_covariance_preflight import (
    HELD_OUT_CYCLE_TIMES,
    evaluate_multislice_velocity_column_preflight,
    generate_actual_core_report,
)


def _basis(points):
    points = np.asarray(points, dtype=float)
    x, y, _ = points.T
    r = np.hypot(x, y)
    er = np.column_stack((x / r, y / r, np.zeros_like(r)))
    et = np.column_stack((-y / r, x / r, np.zeros_like(r)))
    ez = np.column_stack((np.zeros_like(r), np.zeros_like(r), np.ones_like(r)))
    return r, er, et, ez


def _theta_unit_column(points, time, phase):
    r, er, et, _ = _basis(points)
    c = np.cos(phase)
    return (r * c)[:, None] * er + (2.0 * r * c)[:, None] * et


def _axial_unit_column(points, time, phase):
    r, er, _, ez = _basis(points)
    c = np.cos(phase)
    return (r * c)[:, None] * er + (2.0 * r * c)[:, None] * ez


def _synthetic_slices(amplitude=0.25):
    radii = np.linspace(0.1, 0.4, 9)
    current = np.column_stack((2.0 * amplitude * radii**2, np.zeros_like(radii)))
    second = np.column_stack((np.zeros_like(radii), 2.0 * amplitude * radii**2))
    target = 0.01 * current + 0.02 * second
    receipt1 = build_missing_covariance_column_target(radii, target, current)
    receipt2 = build_missing_covariance_column_target(radii, 1.1 * target, current)
    return [(0.4, receipt1), (0.6, receipt2)]


def test_multislice_independent_velocity_column_passes_every_slice():
    audit = evaluate_multislice_velocity_column_preflight(
        _synthetic_slices(),
        _axial_unit_column,
        amplitude=0.25,
        z=0.08,
        angular_count=8,
        phase_count=16,
        coefficient_budget=0.05,
    )
    assert audit["slice_count"] == 2
    assert audit["all_slices_bounded_inverse_preflight_passed"] is True
    assert audit["finite_cycle_rerun_allowed"] is True
    assert audit["rank2_required_nodes_total"] == audit["required_nodes_total"]
    for item in audit["slices"]:
        bounded = item["bounded_preflight"]["bounded_inverse"]
        assert bounded["all_required_nodes_rank2"] is True
        assert bounded["both_signed_updates_within_budget_on_all_required"] is True
        assert bounded["two_column_relative_stress_residual_rms"] < 1.0e-12


def test_multislice_duplicate_velocity_column_is_rejected_on_every_slice():
    audit = evaluate_multislice_velocity_column_preflight(
        _synthetic_slices(),
        _theta_unit_column,
        amplitude=0.25,
        z=0.08,
        angular_count=8,
        phase_count=16,
        coefficient_budget=0.05,
    )
    assert audit["finite_cycle_rerun_allowed"] is False
    assert audit["rank2_required_nodes_total"] == 0
    assert audit["required_nodes_total"] > 0
    assert all(
        item["finite_cycle_rerun_allowed_for_slice"] is False
        for item in audit["slices"]
    )


def test_multislice_preflight_requires_common_grid_and_frozen_budget():
    slices = _synthetic_slices()
    bad_radii = np.linspace(0.11, 0.41, 9)
    current = np.column_stack((bad_radii**2, np.zeros_like(bad_radii)))
    target = current + np.column_stack((np.zeros_like(bad_radii), bad_radii**2))
    bad = build_missing_covariance_column_target(bad_radii, target, current)

    with pytest.raises(ValueError, match="same radial grid"):
        evaluate_multislice_velocity_column_preflight(
            [slices[0], (0.6, bad)],
            _axial_unit_column,
            amplitude=0.25,
            coefficient_budget=0.05,
        )
    with pytest.raises(ValueError, match="coefficient_budget"):
        evaluate_multislice_velocity_column_preflight(
            slices,
            _axial_unit_column,
            amplitude=0.25,
            coefficient_budget=0.0,
        )
    with pytest.raises(ValueError, match="strictly increasing"):
        evaluate_multislice_velocity_column_preflight(
            [slices[1], slices[0]],
            _axial_unit_column,
            amplitude=0.25,
            coefficient_budget=0.05,
        )


def test_real_multislice_report_uses_real_defects_and_blocks_duplicate_column(tmp_path):
    report = generate_actual_core_report(output=tmp_path / "report.json")
    assert report["inputs"]["times"] == list(HELD_OUT_CYCLE_TIMES)
    assert report["inputs"]["surrogate_defect_used"] is False
    assert report["inputs"]["coefficient_budget"] > 0.0
    assert report["truth_boundary"]["coefficient_budget_changed"] is False
    assert report["truth_boundary"]["finite_correction_cycle_run"] is False
    assert report["routing"]["finite_correction_cycle_rerun_allowed"] is False

    target = report["real_multislice_target"]
    assert target["slice_count"] == 3
    assert target["missing_relative_vector_rms_min"] > 0.0
    assert target["required_transverse_over_current_max_across_slices"] > 0.0
    for item in target["slices"]:
        assert item["target_metrics"]["nodes_requiring_second_direction"] > 0
        assert item["required_second_column_envelope"][
            "coefficient_budget"
        ] == report["inputs"]["coefficient_budget"]

    duplicate = report["duplicate_existing_column_negative_control"]
    assert duplicate["slice_count"] == 3
    assert duplicate["rank2_required_nodes_total"] == 0
    assert duplicate["required_nodes_total"] > 0
    assert duplicate["finite_cycle_rerun_allowed"] is False
