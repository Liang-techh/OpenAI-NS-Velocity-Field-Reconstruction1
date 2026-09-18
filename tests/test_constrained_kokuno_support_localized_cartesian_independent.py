from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_support_localized_cartesian_independent import (
    BASE_HEAD,
    CASES,
    FD_STEPS,
    TASK_ID,
    _cartesian_curl_fd2,
    _public_cartesian_velocity,
    build_report,
    write_report,
)


@pytest.fixture(scope="module")
def report():
    return build_report()


def test_cartesian_oracle_agrees_with_public_localized_mode_on_one_probe():
    case = CASES[0]
    R = 0.73
    theta = 0.81
    xyz = np.asarray((R * np.cos(theta), R * np.sin(theta), 0.21), dtype=float)
    oracle = _cartesian_curl_fd2(case, xyz, 0.001)
    public = _public_cartesian_velocity(case, xyz)
    relative = np.linalg.norm(public - oracle) / np.linalg.norm(oracle)
    # Broad software sanity guard only.  The report contains the frozen stricter
    # scientific/local implementation guards and preserves a rejection if they fail.
    assert relative < 1.0e-2


def test_report_is_three_resolution_parameter_and_small_radius_audit(report):
    assert report["task_id"] == TASK_ID
    assert report["dependency"]["exact_head"] == BASE_HEAD
    assert [row["step"] for row in report["resolution_ladder"]] == list(FD_STEPS)
    assert report["sampling"]["regions"] == {"off_grid": 18, "small_radius": 15}
    assert report["sampling"]["total_points"] == 33
    assert len(report["parameter_cases"]) == 3
    for row in report["resolution_ladder"]:
        assert np.isfinite(row["curl_relative_rms"])
        assert np.isfinite(row["divergence_max_abs"])
        assert row["sample_count"] == 33
        assert set(row["regions"]) == {"off_grid", "small_radius"}
        assert set(row["parameter_cases"]) == {case.label for case in CASES}


def test_missing_support_gradient_mutation_is_detected(report):
    finest = report["resolution_ladder"][-1]
    mutation = report["mutation"]
    # This is a calibration of the independent oracle, not a PDE gate: removing
    # the product-rule support derivatives must be visibly worse than the correct path.
    assert mutation["curl_point_relative_rms"] > 10.0 * finest["curl_relative_rms"]
    assert mutation["mutated_divergence_max_abs"] > finest["divergence_max_abs"]
    assert report["local_guards"]["checks"]["missing_support_gradient_curl_detected"]
    assert report["local_guards"]["checks"]["missing_support_gradient_divergence_detected"]


def test_truth_boundary_keeps_formal_pde_gate_unassessed(report):
    formal = report["formal_project_gates"]
    truth = report["truth_boundary"]
    assert formal["normalized_momentum_max_l2"] == 1.0e-3
    assert formal["divergence_max_l2"] == 1.0e-5
    assert formal["formal_full_domain_pde_gate_assessed"] is False
    assert truth["actual_source_public_oscillatory_xyz_t_velocity_available"] is False
    assert truth["global_leading_velocity_pressure_available"] is False
    assert truth["after_correction_global_velocity_available"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert report["operator"]["construction_cylindrical_fd4_oracle_reused"] is False
    assert report["operator"]["training_tensor_or_loss_read"] is False
    assert report["operator"]["free_forcing_used"] is False


def test_report_write_refuses_overwrite(tmp_path):
    path = tmp_path / "audit.json"
    write_report(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["task_id"] == TASK_ID
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_report(path)
