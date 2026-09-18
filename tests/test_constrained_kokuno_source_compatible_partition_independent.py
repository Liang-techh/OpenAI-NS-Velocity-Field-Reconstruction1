import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_compatible_partition_independent import (
    BASE_HEAD,
    BASE_PR,
    DROPPED_DERIVATIVE_RELATIVE_ERROR_FLOOR,
    DROPPED_LABEL_PARTITION_ERROR_FLOOR,
    FD6_STEPS,
    SEED,
    _align,
    build_report,
    write_report,
)


def test_label_alignment_zero_fills_changed_active_sets():
    labels = ((5, (0, 0, 0)), (6, (1, 0, 0)))
    union = ((5, (0, 0, 0)), (5, (0, 1, 0)), (6, (1, 0, 0)))
    values = np.array([[2.0, 3.0], [4.0, 5.0]])
    aligned = _align(values, labels, union)
    np.testing.assert_array_equal(aligned, np.array([[2.0, 0.0, 3.0], [4.0, 0.0, 5.0]]))


def test_independent_partition_report_is_black_box_and_mutation_sensitive():
    report = build_report()
    assert report["dependency"]["agent2_parent_pr"] == BASE_PR
    assert report["dependency"]["agent2_exact_head"] == BASE_HEAD
    assert report["preregistration"]["seed"] == SEED
    assert report["preregistration"]["fd6_steps"] == list(FD6_STEPS)
    assert report["preregistration"]["thresholds_changed_after_results"] is False

    operator = report["independent_operator"]
    assert operator["construction_raw_bump_helper_used"] is False
    assert operator["construction_normalization_helper_used"] is False
    assert operator["construction_analytic_derivative_helper_used"] is False
    assert operator["agent2_fd4_test_operator_reused"] is False
    assert operator["training_tensor_or_loss_read"] is False
    assert operator["free_forcing_used"] is False

    for origin in report["origin_reports"]:
        assert origin["held_out_points"] == 24
        assert origin["maximum_active_ell_span"] <= 2
        assert origin["mutation"]["dropped_label_partition_error_max_abs"] >= DROPPED_LABEL_PARTITION_ERROR_FLOOR
        assert origin["mutation"]["dropped_derivative_relative_error_vs_finest_fd6"] >= DROPPED_DERIVATIVE_RELATIVE_ERROR_FLOOR
        for direction in ("r", "z"):
            ladder = origin["fd6_ladders"][direction]
            assert [row["step"] for row in ladder] == list(FD6_STEPS)
            assert all(np.isfinite(row["relative_rms"]) for row in ladder)
            assert all(row["aligned_beta_count"] >= origin["base_beta_count"] for row in ladder)
        assert origin["local_preflight_passed"] == all(origin["checks"].values())

    expected = all(
        value
        for origin in report["origin_reports"]
        for value in origin["checks"].values()
    )
    assert report["local_guards"]["source_compatible_partition_independent_preflight_passed"] == expected
    assert report["formal_project_gates"]["formal_full_domain_pde_gate_assessed"] is False
    assert report["truth_boundary"]["pde_validated"] is False


def test_report_round_trip_and_no_overwrite(tmp_path):
    path = tmp_path / "audit.json"
    write_report(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["schema"].startswith("kokuno-agent4-source-compatible-partition")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_report(path)
