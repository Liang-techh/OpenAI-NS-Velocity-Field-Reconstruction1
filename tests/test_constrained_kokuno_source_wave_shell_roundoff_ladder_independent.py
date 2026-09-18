from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_wave_shell_roundoff_ladder_independent import (
    BASE_HEAD,
    FD4_STEPS,
    FORMAL_DIVERGENCE_GATE,
    FORMAL_MOMENTUM_GATE,
    LOCAL_CURL_FINE_REL_GUARD,
    LOCAL_DIVERGENCE_FINE_ABS_GUARD,
    LOCAL_REFINEMENT_RATIO_GUARD,
    PRIOR_REJECTION,
    SEED,
    TASK_ID,
    _sample_probes,
    build_report,
    write_report,
)


@pytest.fixture(scope="module")
def report():
    return build_report()


def test_roundoff_ladder_is_fresh_and_preregistered(report):
    assert report["task_id"] == TASK_ID
    assert report["dependency"]["exact_head"] == BASE_HEAD
    assert SEED == 9173111
    assert FD4_STEPS == (0.004, 0.002, 0.001)
    assert report["seed"] == SEED
    assert [row["step"] for row in report["resolution_ladder"]] == list(FD4_STEPS)
    prereg = report["preregistration"]
    assert prereg["historical_rejection_overwritten"] is False
    assert prereg["unchanged_scientific_guards"] is True
    assert prereg["roundoff_hypothesis_is_not_an_acceptance_override"] is True
    assert prereg["changed_before_this_run"] == {
        "seed": [9173101, 9173111],
        "fd4_steps": [[0.001, 0.0005, 0.00025], [0.004, 0.002, 0.001]],
    }


def test_pr371_rejection_is_preserved_verbatim(report):
    prior = report["preregistration"]["prior_rejection"]
    assert prior == PRIOR_REJECTION
    assert prior["pr"] == 371
    assert prior["exact_head"] == BASE_HEAD
    assert prior["source_wave_shell_independent_cartesian_preflight_passed"] is False
    assert prior["divergence_rms_refinement_ratios"][-1] == pytest.approx(
        1.193958216011321
    )


def test_fresh_source_domain_sampling_and_stencils(report):
    probes = _sample_probes()
    assert len(probes) == 36
    assert report["sampling"]["regions"] == {
        "source_shell_interior": 18,
        "near_inner_shell_edge": 18,
    }
    assert report["sampling"]["fresh_seed_relative_to_pr371"] is True
    stencil = report["sampling"]["stencil_source_domain"]
    assert stencil["all_fd4_stencil_points_source_shell_valid"] is True
    assert stencil["minimum_X_over_X_a_across_stencils"] > 1.0
    assert stencil["minimum_R_over_axis_lower_bound_across_stencils"] > 1.0
    for row in report["resolution_ladder"]:
        assert row["sample_count"] == 36
        assert np.isfinite(row["curl_relative_rms"])
        assert np.isfinite(row["divergence_max_abs"])
        assert np.isfinite(row["divergence_rms_abs"])


def test_scientific_thresholds_are_exactly_unchanged(report):
    thresholds = report["local_guards"]["thresholds"]
    assert LOCAL_CURL_FINE_REL_GUARD == 2.0e-4
    assert LOCAL_DIVERGENCE_FINE_ABS_GUARD == 2.0e-4
    assert LOCAL_REFINEMENT_RATIO_GUARD == 2.5
    assert thresholds["finest_curl_relative_rms_max"] == 2.0e-4
    assert thresholds["finest_divergence_max_abs_max"] == 2.0e-4
    assert thresholds["minimum_refinement_ratio"] == 2.5
    assert report["local_guards"][
        "same_thresholds_as_agent4_prs_353_360_371"
    ] is True
    assert FORMAL_MOMENTUM_GATE == 1.0e-3
    assert FORMAL_DIVERGENCE_GATE == 1.0e-5
    assert report["formal_project_gates"]["normalized_momentum_max_l2"] == 1.0e-3
    assert report["formal_project_gates"]["divergence_max_l2"] == 1.0e-5


def test_scientific_verdict_is_derived_without_forcing_a_pass(report):
    checks = report["local_guards"]["checks"]
    verdict = report["local_guards"][
        "source_wave_shell_roundoff_ladder_preflight_passed"
    ]
    assert verdict is all(checks.values())
    # CI success means the scientific verdict was reproduced, not that it passed.
    assert isinstance(verdict, bool)


def test_mutation_and_negative_controls_remain_live(report):
    mutation = report["mutation"]
    finest = report["resolution_ladder"][-1]
    assert mutation["curl_point_relative_rms"] > 10.0 * finest["curl_relative_rms"]
    assert mutation["mutated_divergence_max_abs"] > finest["divergence_max_abs"]
    controls = report["negative_controls"]
    assert controls["off_source_angular_grid_phase_rejected"] is True
    assert controls["below_inner_source_wave_shell_rejected"] is True


def test_truth_boundary_keeps_formal_pde_gate_false(report):
    operator = report["operator"]
    truth = report["truth_boundary"]
    assert operator["agent2_cylindrical_validation_oracle_reused"] is False
    assert operator["training_tensor_or_loss_read"] is False
    assert operator["construction_residual_operator_reused"] is False
    assert operator["free_forcing_used"] is False
    assert truth["actual_source_public_oscillatory_xyz_t_velocity_available"] is False
    assert truth["global_leading_velocity_pressure_available"] is False
    assert truth["genuinely_independent_second_covariance_column_available"] is False
    assert truth["after_correction_global_velocity_available"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False


def test_report_write_refuses_overwrite(tmp_path):
    path = tmp_path / "source_wave_shell_roundoff_ladder_audit.json"
    write_report(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["task_id"] == TASK_ID
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_report(path)
