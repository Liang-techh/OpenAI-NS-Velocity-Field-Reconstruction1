from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_support_localized_axis_refinement import (
    BASE_HEAD,
    CASES,
    FD4_STEPS,
    LOCAL_DIVERGENCE_FINE_ABS_GUARD,
    TASK_ID,
    _cartesian_curl_fd4,
    _public_cartesian_velocity,
    build_report,
    write_report,
)


@pytest.fixture(scope="module")
def report():
    return build_report()


def test_fd4_cartesian_oracle_agrees_with_public_mode_on_one_probe():
    case = CASES[0]
    radius = 0.19
    theta = 0.93
    xyz = np.asarray(
        (radius * np.cos(theta), radius * np.sin(theta), 0.17), dtype=float
    )
    oracle = _cartesian_curl_fd4(case, xyz, 0.0005)
    public = _public_cartesian_velocity(case, xyz)
    relative = np.linalg.norm(public - oracle) / np.linalg.norm(oracle)
    # Broad software sanity only; report guards remain the scientific local checks.
    assert relative < 1.0e-2


def test_report_uses_fresh_three_resolution_small_radius_ladder(report):
    assert report["task_id"] == TASK_ID
    assert report["dependency"]["exact_head"] == BASE_HEAD
    assert [row["step"] for row in report["resolution_ladder"]] == list(FD4_STEPS)
    assert report["sampling"]["regions"] == {"off_grid": 15, "small_radius": 21}
    assert report["sampling"]["total_points"] == 36
    assert report["sampling"]["fresh_seed_relative_to_parent"] is True
    assert len(report["parameter_cases"]) == 3
    for row in report["resolution_ladder"]:
        assert np.isfinite(row["curl_relative_rms"])
        assert np.isfinite(row["divergence_max_abs"])
        assert row["sample_count"] == 36
        assert set(row["regions"]) == {"off_grid", "small_radius"}
        assert set(row["parameter_cases"]) == {case.label for case in CASES}


def test_parent_rejection_threshold_is_not_relaxed(report):
    thresholds = report["local_guards"]["thresholds"]
    assert LOCAL_DIVERGENCE_FINE_ABS_GUARD == 2.0e-4
    assert thresholds["finest_divergence_max_abs_max"] == 2.0e-4
    assert thresholds["finest_curl_relative_rms_max"] == 2.0e-4
    assert report["local_guards"]["same_thresholds_as_parent_rejection"] is True


def test_missing_support_gradient_mutation_remains_detectable(report):
    finest = report["resolution_ladder"][-1]
    mutation = report["mutation"]
    assert mutation["curl_point_relative_rms"] > 10.0 * finest["curl_relative_rms"]
    assert mutation["mutated_divergence_max_abs"] > finest["divergence_max_abs"]
    assert report["local_guards"]["checks"]["missing_support_gradient_curl_detected"]
    assert report["local_guards"]["checks"]["missing_support_gradient_divergence_detected"]


def test_truth_boundary_keeps_formal_pde_gate_unassessed(report):
    operator = report["operator"]
    formal = report["formal_project_gates"]
    truth = report["truth_boundary"]
    assert operator["construction_cylindrical_fd4_oracle_reused"] is False
    assert operator["previous_agent4_cartesian_fd2_validation_operator_reused"] is False
    assert operator["training_tensor_or_loss_read"] is False
    assert operator["free_forcing_used"] is False
    assert formal["normalized_momentum_max_l2"] == 1.0e-3
    assert formal["divergence_max_l2"] == 1.0e-5
    assert formal["formal_full_domain_pde_gate_assessed"] is False
    assert truth["actual_source_public_oscillatory_xyz_t_velocity_available"] is False
    assert truth["global_leading_velocity_pressure_available"] is False
    assert truth["after_correction_global_velocity_available"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False


def test_report_write_refuses_overwrite(tmp_path):
    path = tmp_path / "axis_refinement.json"
    write_report(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["task_id"] == TASK_ID
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_report(path)
