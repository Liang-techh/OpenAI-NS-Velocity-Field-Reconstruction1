from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_wave_shell_cartesian_independent import (
    BASE_HEAD,
    CASES,
    FD4_STEPS,
    LOCAL_DIVERGENCE_FINE_ABS_GUARD,
    TASK_ID,
    _cartesian_curl_fd4,
    _public_cartesian_velocity,
    _sample_probes,
    build_report,
    write_report,
)


@pytest.fixture(scope="module")
def report():
    return build_report()


def test_source_grid_cartesian_oracle_agrees_on_one_held_out_probe():
    probe = _sample_probes()[0]
    case = probe["case"]
    oracle = _cartesian_curl_fd4(case, probe["xyz"], 0.0005)
    public = _public_cartesian_velocity(case, probe["s_Q"], probe["xyz"])
    relative = np.linalg.norm(public - oracle) / np.linalg.norm(oracle)
    # Broad software sanity only. The frozen report guards own scientific routing.
    assert relative < 1.0e-2


def test_report_uses_source_grid_shell_and_fresh_three_resolution_ladder(report):
    assert report["task_id"] == TASK_ID
    assert report["dependency"]["exact_head"] == BASE_HEAD
    assert [row["step"] for row in report["resolution_ladder"]] == list(FD4_STEPS)
    assert report["sampling"]["regions"] == {
        "source_shell_interior": 18,
        "near_inner_shell_edge": 18,
    }
    assert report["sampling"]["total_points"] == 36
    assert report["sampling"]["fresh_seed_relative_to_prior_agent4_audits"] is True
    assert report["sampling"]["stencil_source_domain"][
        "all_fd4_stencil_points_source_shell_valid"
    ] is True
    assert report["sampling"]["stencil_source_domain"][
        "minimum_X_over_X_a_across_stencils"
    ] > 1.0
    assert len(report["parameter_cases"]) == 3
    for case_row in report["parameter_cases"]:
        assert case_row["physical_azimuthal_mode_mj"] != 0
        assert case_row["p"] * case_row["k"] == pytest.approx(case_row["j"])
    for row in report["resolution_ladder"]:
        assert np.isfinite(row["curl_relative_rms"])
        assert np.isfinite(row["divergence_max_abs"])
        assert row["sample_count"] == 36
        assert set(row["regions"]) == {
            "source_shell_interior",
            "near_inner_shell_edge",
        }
        assert set(row["parameter_cases"]) == {case.label for case in CASES}


def test_frozen_local_and_formal_thresholds_are_not_relaxed(report):
    thresholds = report["local_guards"]["thresholds"]
    formal = report["formal_project_gates"]
    assert LOCAL_DIVERGENCE_FINE_ABS_GUARD == 2.0e-4
    assert thresholds["finest_divergence_max_abs_max"] == 2.0e-4
    assert thresholds["finest_curl_relative_rms_max"] == 2.0e-4
    assert thresholds["minimum_refinement_ratio"] == 2.5
    assert report["local_guards"]["same_thresholds_as_agent4_prs_353_360"] is True
    assert formal["normalized_momentum_max_l2"] == 1.0e-3
    assert formal["divergence_max_l2"] == 1.0e-5
    assert formal["formal_full_domain_pde_gate_assessed"] is False


def test_source_domain_negative_controls_are_retained(report):
    controls = report["negative_controls"]
    checks = report["local_guards"]["checks"]
    assert controls["off_source_angular_grid_phase_rejected"] is True
    assert controls["below_inner_source_wave_shell_rejected"] is True
    assert checks["off_source_angular_grid_phase_rejected"] is True
    assert checks["below_inner_source_wave_shell_rejected"] is True


def test_missing_support_gradient_mutation_is_detectable(report):
    finest = report["resolution_ladder"][-1]
    mutation = report["mutation"]
    assert mutation["curl_point_relative_rms"] > 10.0 * finest["curl_relative_rms"]
    assert mutation["mutated_divergence_max_abs"] > finest["divergence_max_abs"]
    assert report["local_guards"]["checks"]["missing_support_gradient_curl_detected"]
    assert report["local_guards"]["checks"][
        "missing_support_gradient_divergence_detected"
    ]


def test_truth_boundary_keeps_full_pde_validation_false(report):
    operator = report["operator"]
    truth = report["truth_boundary"]
    partition = report["manufactured_partition_family"]
    assert operator["agent2_cylindrical_validation_oracle_reused"] is False
    assert operator["previous_agent4_generic_p_phase_reused"] is False
    assert operator["training_tensor_or_loss_read"] is False
    assert operator["construction_residual_operator_reused"] is False
    assert operator["free_forcing_used"] is False
    assert partition["sum_eta_beta_squared_exact_algebraically"] is True
    assert partition["concrete_source_chi_bumps_claimed"] is False
    assert truth["actual_source_Lambda_recovered"] is False
    assert truth["concrete_source_partition_bumps_reconstructed"] is False
    assert truth["actual_source_public_oscillatory_xyz_t_velocity_available"] is False
    assert truth["global_leading_velocity_pressure_available"] is False
    assert truth["genuinely_independent_second_covariance_column_available"] is False
    assert truth["after_correction_global_velocity_available"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False


def test_report_write_refuses_overwrite(tmp_path):
    path = tmp_path / "source_wave_shell_cartesian_audit.json"
    write_report(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["task_id"] == TASK_ID
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_report(path)
