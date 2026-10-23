import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_refinement_audit import audit_refinement


THRESHOLDS = {
    "pde_residual_max": 1e-3,
    "pde_residual_L2": 1e-3,
    "divergence_max": 1e-5,
    "divergence_L2": 1e-5,
}


def _rows(values):
    rows = []
    for time in (0.25, 0.75):
        for step, value in values.items():
            rows.append(
                {
                    "time": time,
                    "step": step,
                    "residual_sampled_max": value[0],
                    "residual_L2_estimate": value[1],
                    "divergence_sampled_max": value[2],
                    "divergence_L2_estimate": value[3],
                }
            )
    return rows


def test_coupled_joint_refinement_separates_momentum_plateau_from_divergence():
    config = json.loads(Path("configs/constraints.json").read_text())
    validation = json.loads(
        Path("artifacts/constrained/coupled_joint/validation.json").read_text()
    )
    result = audit_refinement(validation, config["validation"]["thresholds"])

    assert result["source_candidate"] == "artifacts/constrained/coupled_joint/candidate.json"
    assert result["validation_seed"] == 914027
    assert result["validation_points"] == 4096
    assert result["refinement_steps"] == [0.005, 0.01, 0.02]
    assert result["momentum_plateau_all_times"]
    assert result["divergence_contracting_all_times"]
    assert result["classification"] == (
        "momentum_mismatch_not_explained_by_current_fd_refinement"
    )
    assert result["finest_residual_sampled_max_range"] == pytest.approx(
        [0.6042252486915877, 1.0022183043210762]
    )
    assert result["finest_residual_L2_range"] == pytest.approx(
        [1.169647947564117, 2.2215324099175593]
    )
    assert max(
        row["momentum_relative_change_medium_to_finest"]["residual_sampled_max"]
        for row in result["rows"]
    ) < 0.075
    assert min(
        row["momentum_threshold_multiples_at_finest"]["residual_sampled_max"]
        for row in result["rows"]
    ) > 600
    assert min(
        row["momentum_threshold_multiples_at_finest"]["residual_L2_estimate"]
        for row in result["rows"]
    ) > 1100
    assert max(
        max(row["divergence_fine_over_medium"].values()) for row in result["rows"]
    ) < 0.32


def test_discretization_only_case_is_not_labeled_plateau():
    values = {}
    for step in (0.02, 0.01, 0.005):
        error = 40 * step * step
        values[step] = (error, error, error, error)
    result = audit_refinement({"rows": _rows(values)}, THRESHOLDS)
    assert not result["momentum_plateau_all_times"]
    assert result["classification"] == "refinement_interpretation_unresolved"


def test_refinement_audit_fails_closed_on_incomplete_or_invalid_inputs():
    two_levels = {
        0.01: (1.0, 1.0, 0.1, 0.1),
        0.005: (1.0, 1.0, 0.01, 0.01),
    }
    with pytest.raises(ValueError, match="at least three"):
        audit_refinement({"rows": _rows(two_levels)}, THRESHOLDS)

    three_levels = {0.02: (1.0, 1.0, 0.2, 0.2), **two_levels}
    with pytest.raises(ValueError, match="relative tolerances"):
        audit_refinement(
            {"rows": _rows(three_levels)},
            THRESHOLDS,
            plateau_relative_tolerance=1.0,
        )
