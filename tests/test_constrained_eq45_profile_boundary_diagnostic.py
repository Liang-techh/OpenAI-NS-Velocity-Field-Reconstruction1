from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_profile_boundary_diagnostic import (
    TASK_ID,
    audit_eq45_profile_boundary_pressure,
    default_receipt_path,
    load_and_validate_receipt,
)


def test_eq45_profile_boundary_diagnostic_replays_registered_problem(capsys):
    report = audit_eq45_profile_boundary_pressure()

    assert report.task_id == TASK_ID
    assert report.receipt_candidate_sha256 == report.reconstructed_candidate_sha256
    assert report.fit_seed != report.holdout_seed
    assert report.derivative_steps == (0.02, 0.01, 0.005)
    assert report.profile_bounds == (-4.0, 4.0)
    assert report.force_bounds == (0.0, 10.0)
    assert report.best_boundary_probe in {"phi_upper", "F_upper", "upper_corner"}
    assert np.isfinite(report.local_gradient_phi)
    assert np.isfinite(report.local_gradient_F)

    probes = {probe.name: probe for probe in report.boundary_probes}
    assert set(probes) == {
        "receipt",
        "phi_upper",
        "F_upper",
        "upper_corner",
        "phi_inward",
        "F_inward",
    }
    receipt = probes["receipt"]
    assert receipt.training_rms_after_force == pytest.approx(
        report.receipt_training_rms_after_force, rel=0.0, abs=1e-12
    )
    for probe in report.boundary_probes:
        assert report.profile_bounds[0] <= probe.phi_02 <= report.profile_bounds[1]
        assert report.profile_bounds[0] <= probe.F_02 <= report.profile_bounds[1]
        assert report.force_bounds[0] <= probe.force_a <= report.force_bounds[1]
        assert report.force_bounds[0] <= probe.force_c <= report.force_bounds[1]
        assert probe.training_rms_after_force > 0.0
        assert probe.training_max_after_force >= probe.training_rms_after_force

    restart = report.alternate_restart
    assert restart.function_evaluations > 0
    assert restart.iterations >= 0
    assert report.profile_bounds[0] <= restart.phi_02 <= report.profile_bounds[1]
    assert report.profile_bounds[0] <= restart.F_02 <= report.profile_bounds[1]
    assert report.force_bounds[0] <= restart.force_a <= report.force_bounds[1]
    assert report.force_bounds[0] <= restart.force_c <= report.force_bounds[1]
    assert restart.training_rms_after_force > 0.0

    assert tuple(row.spatial_step for row in report.holdout_levels_for_best_boundary) == (
        0.02,
        0.01,
        0.005,
    )
    for row in report.holdout_levels_for_best_boundary:
        assert row.receipt_rms_after_force > 0.0
        assert row.probe_rms_after_force > 0.0
        assert row.probe_max_after_force >= row.probe_rms_after_force
        assert np.isfinite(row.relative_change_vs_receipt)

    checked_path = default_receipt_path().with_name("eq45_profile_boundary_diagnostic.json")
    checked = json.loads(checked_path.read_text(encoding="utf-8"))
    assert checked["schema"] == "eq45_profile_boundary_diagnostic_v1"
    assert checked["task_id"] == TASK_ID
    assert checked["dependency"]["candidate_sha256"] == report.receipt_candidate_sha256
    assert checked["diagnosis"]["one_start_local_basin_evidence"] is False
    assert checked["diagnosis"]["exact_upper_bound_improves_registered_objective"] is False
    assert checked["diagnosis"]["two_parameter_subspace_appears_locally_exhausted"] is True
    assert checked["receipt_replay"]["training_rms_after_force"] == pytest.approx(
        report.receipt_training_rms_after_force, rel=5e-10, abs=5e-10
    )
    assert checked["alternate_restart"]["training_rms_after_force"] == pytest.approx(
        restart.training_rms_after_force, rel=5e-10, abs=5e-10
    )
    assert checked["best_boundary_holdout"]["probe"] == report.best_boundary_probe
    for expected, measured in zip(
        checked["best_boundary_holdout"]["levels"],
        report.holdout_levels_for_best_boundary,
    ):
        assert expected["step"] == measured.spatial_step
        assert expected["probe_rms_after_force"] == pytest.approx(
            measured.probe_rms_after_force, rel=5e-10, abs=5e-10
        )

    assert report.velocity_changed is False
    assert report.pressure_fitted is False
    assert report.forcing_family_changed is False
    assert report.holdout_force_refit is False
    assert report.pde_validated is False
    assert report.visualization_ready_promoted is False
    assert report.visual_correspondence_verified is False
    assert report.paper_exact is False
    assert report.openai_field_identified is False

    # Keep the exact deterministic production numerics visible in Actions logs without
    # turning them into a PDE acceptance threshold.
    with capsys.disabled():
        print("EQ45_PROFILE_BOUNDARY_DIAGNOSTIC=" + json.dumps(report.to_dict(), sort_keys=True))


def test_eq45_profile_boundary_receipt_fails_closed_on_contract_drift(tmp_path):
    payload = json.loads(default_receipt_path().read_text(encoding="utf-8"))
    payload["contract"]["force_family"] = "residual_defined_free_force"
    mutated = tmp_path / "mutated_receipt.json"
    mutated.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="forcing-family contract drift"):
        load_and_validate_receipt(mutated)


def test_eq45_profile_boundary_diagnostic_rejects_invalid_restart():
    with pytest.raises(ValueError, match="outside registered profile bounds"):
        audit_eq45_profile_boundary_pressure(alternate_start=(4.1, 0.0))
