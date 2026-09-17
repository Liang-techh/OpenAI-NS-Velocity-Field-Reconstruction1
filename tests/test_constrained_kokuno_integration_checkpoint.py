from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_integration_checkpoint import (
    BRIDGE_CANDIDATE_SHA256,
    OPEN_UNCONSUMED_SIBLINGS,
    REGISTERED_PDE_THRESHOLD,
    UPSTREAM_HEADS,
    build_integration_checkpoint,
    checkpoint_sha256,
    load_integration_checkpoint,
)


def _write_rehashed(path: Path, payload: dict) -> None:
    payload = copy.deepcopy(payload)
    payload["receipt_sha256"] = checkpoint_sha256(payload)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def test_agent5_checkpoint_runs_current_agent1_to_agent4_chain(tmp_path: Path) -> None:
    payload = build_integration_checkpoint(
        output_dir=tmp_path,
        radial_count=33,
        angular_count=16,
    )

    assert payload["upstream_heads"] == UPSTREAM_HEADS
    assert payload["open_unconsumed_siblings"] == OPEN_UNCONSUMED_SIBLINGS
    assert payload["engineering_bridge"]["candidate_sha256"] == BRIDGE_CANDIDATE_SHA256
    assert payload["engineering_bridge"]["role"] == "temporary_engineering_bridge_not_kokuno_leading"

    source = payload["source_bundle"]
    assert Path(source["leading_axis_profile_path"]).is_file()
    assert Path(source["native_coordinates_path"]).is_file()
    assert Path(source["agent3_report_path"]).is_file()
    assert Path(source["agent4_report_path"]).is_file()
    assert source["native_relation_probe_max_abs"] < 1.0e-11

    state = payload["stage_state"]
    assert state == {
        "leading_axis_profile_ready": True,
        "native_coordinates_ready": True,
        "leading_ready": False,
        "oscillatory_ready": True,
        "mean_defect_ready": True,
        "radial_inverse_ready": True,
        "radial_inverse_independent_validation_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }

    metrics = payload["stage_metrics"]
    assert metrics["full_sampled_mean_defect_increment_rms"] > 0.0
    assert metrics["ring_theta_raw_rms"] > 0.0
    assert 0.0 < metrics["theta_gate_capture_rms_ratio"] < 1.0
    assert metrics["theta_radial_stress_rms"] > 0.0
    assert metrics["theta_radial_independent_stress_disagreement_normalized_rms"] < 5.0e-3
    assert metrics["theta_radial_independent_identity_rms"] < 5.0e-6
    assert metrics["theta_radial_independent_mutation_ratio"] > 20.0
    assert metrics["full_stage_residual_comparison_available"] is False

    gate = payload["pde_gate"]
    assert gate["threshold"] == REGISTERED_PDE_THRESHOLD == 1.0e-3
    assert gate["assessed"] is False
    assert gate["passed"] is False
    assert payload["truth_boundary"]["free_forcing_used"] is False
    assert payload["truth_boundary"]["pde_validated"] is False

    receipt_path = tmp_path / "integration_checkpoint.json"
    loaded = load_integration_checkpoint(receipt_path)
    assert loaded["receipt_sha256"] == payload["receipt_sha256"]


def test_agent5_checkpoint_fail_closes_scientific_promotion(tmp_path: Path) -> None:
    payload = build_integration_checkpoint(
        output_dir=tmp_path / "source",
        radial_count=33,
        angular_count=16,
    )

    mutations = []

    relaxed_gate = copy.deepcopy(payload)
    relaxed_gate["pde_gate"]["threshold"] = 2.0e-3
    mutations.append(relaxed_gate)

    assessed_gate = copy.deepcopy(payload)
    assessed_gate["pde_gate"]["assessed"] = True
    assessed_gate["pde_gate"]["passed"] = True
    assessed_gate["stage_state"]["pde_validated"] = True
    assessed_gate["truth_boundary"]["pde_validated"] = True
    mutations.append(assessed_gate)

    laundered_bridge = copy.deepcopy(payload)
    laundered_bridge["engineering_bridge"]["role"] = "kokuno_leading_field"
    laundered_bridge["truth_boundary"]["engineering_bridge_is_kokuno_leading"] = True
    mutations.append(laundered_bridge)

    premature_export = copy.deepcopy(payload)
    premature_export["stage_state"]["leading_ready"] = True
    premature_export["stage_state"]["correction_ready"] = True
    premature_export["stage_state"]["velocity_export_ready"] = True
    mutations.append(premature_export)

    for index, mutation in enumerate(mutations):
        path = tmp_path / f"mutation_{index}.json"
        _write_rehashed(path, mutation)
        with pytest.raises(ValueError):
            load_integration_checkpoint(path)
