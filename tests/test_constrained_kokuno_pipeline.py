from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_pipeline import (
    UPSTREAM_HEADS,
    _receipt_sha,
    build_current_pipeline_receipt,
    load_pipeline_receipt,
    validate_pipeline_receipt,
)


def _resign(receipt: dict) -> dict:
    receipt["sha256"] = _receipt_sha(receipt)
    return receipt


def test_kokuno_agent1_to_agent4_preflight_pipeline_is_replayable_and_fail_closed(tmp_path: Path) -> None:
    receipt = build_current_pipeline_receipt(output_dir=tmp_path, mean_defect_points=8)

    assert receipt["upstream_heads"] == UPSTREAM_HEADS
    assert receipt["status"] == {
        "leading_axis_profile_ready": True,
        "leading_ready": False,
        "oscillatory_ready": True,
        "mean_defect_ready": True,
        "correction_ready": False,
        "independent_structural_preflight_ready": True,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert receipt["leading"]["coordinate_mapping_status"] == "pending_explicit_convention_reconciliation"
    assert receipt["bridge_base"]["may_be_promoted_to_kokuno_leading"] is False
    assert receipt["registered_pde_gate"] == {
        "held_out_normalized_ns_residual_threshold": 1.0e-3,
        "threshold_changed": False,
        "assessed_on_full_kokuno_candidate": False,
        "passed": False,
    }

    preflight = receipt["oscillatory"]["structural_preflight"]
    assert preflight["structural_preflight_passed"] is True
    assert preflight["finest_curl_error_rms_max_over_times"] < 1.0e-7
    assert preflight["finest_divergence_rms_max_over_times"] < 1.0e-7
    assert preflight["outside_support_max"] == 0.0

    mean_rows = receipt["mean_radial_correction"]["mean_defect_rows"]
    assert len(mean_rows) == 3
    assert min(float(row["mean_defect_increment_rms"]) for row in mean_rows) > 1.0e-3
    assert receipt["mean_radial_correction"]["radial_stress_inverse_applied"] is False
    assert receipt["mean_radial_correction"]["finite_correction_cycle_run"] is False

    loaded = load_pipeline_receipt(tmp_path / "pipeline_receipt.json")
    assert loaded["sha256"] == receipt["sha256"]
    assert validate_pipeline_receipt(loaded) is True

    promoted = deepcopy(receipt)
    promoted["status"]["pde_validated"] = True
    _resign(promoted)
    with pytest.raises(ValueError, match="pipeline status drift"):
        validate_pipeline_receipt(promoted)

    relaxed = deepcopy(receipt)
    relaxed["registered_pde_gate"]["held_out_normalized_ns_residual_threshold"] = 2.0e-3
    _resign(relaxed)
    with pytest.raises(ValueError, match="1e-3 PDE gate changed"):
        validate_pipeline_receipt(relaxed)

    laundered = deepcopy(receipt)
    laundered["bridge_base"]["may_be_promoted_to_kokuno_leading"] = True
    _resign(laundered)
    with pytest.raises(ValueError, match="relabeled as Kokuno leading"):
        validate_pipeline_receipt(laundered)
