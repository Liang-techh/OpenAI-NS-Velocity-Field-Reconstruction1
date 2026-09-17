import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_phi10_profile_force_optimization import (
    PROFILE_BOUNDS,
    candidate_with_phi10,
    optimize_eq45_phi10_profile_force,
)
from openai_ns_reconstruction.constrained_eq45_profile_force_optimization import (
    candidate_with_profile_pair,
)
from openai_ns_reconstruction.constrained_restricted_force_curl_capacity import FORCE_BOUNDS


ARTIFACT = (
    Path(__file__).resolve().parents[1]
    / "artifacts"
    / "constrained"
    / "eq45_phi10_profile_force_optimization.json"
)


def _payload():
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def _parent_candidate(data):
    seed = Eq45VelocityCandidate.seed()
    frozen = data["contract"]["frozen_profile_parameters"]
    return candidate_with_profile_pair(seed, frozen["Phi(0,2)"], frozen["F(0,2)"])


def test_checked_phi10_receipt_preserves_registered_contract_and_truth_boundary():
    data = _payload()
    assert data["schema"] == "eq45_phi10_profile_force_optimization_v1"
    assert data["task_id"] == "CR005-EQ45-PHI10-PROFILE-FORCE-OPT-017"
    contract = data["contract"]
    assert contract["varied_profile_parameter"] == "Phi(1,0)"
    assert contract["profile_bounds"] == list(PROFILE_BOUNDS)
    assert contract["force_bounds"] == list(FORCE_BOUNDS)
    assert contract["force_family"] == "preregistered_restricted_two_parameter_family"
    assert contract["training_holdout_separate"] is True
    assert contract["holdout_force_refit"] is False
    assert contract["max_function_evaluations"] == 48

    fit = data["fit"]
    assert fit["optimizer_success"] is True
    assert fit["optimizer_function_evaluations"] <= contract["max_function_evaluations"]
    assert PROFILE_BOUNDS[0] <= fit["Phi(1,0)"] <= PROFILE_BOUNDS[1]
    assert fit["Phi(1,0)"] > 3.99
    assert FORCE_BOUNDS[0] <= fit["force_a"] <= FORCE_BOUNDS[1]
    assert FORCE_BOUNDS[0] <= fit["force_c"] <= FORCE_BOUNDS[1]
    assert fit["curl_rms_after_force"] < data["parent"]["training_curl_rms_after_force"]
    assert fit["training_improvement_vs_parent"] > 0.15

    levels = data["holdout"]
    assert [level["step"] for level in levels] == [0.02, 0.01, 0.005]
    assert min(level["improvement_vs_parent"] for level in levels) > 0.16
    assert levels[-1]["optimized_curl_rms_after_force"] > 1.0

    status = data["status"]
    assert status["new_basis_added"] is False
    assert status["forcing_family_changed"] is False
    assert status["pressure_fitted"] is False
    assert status["holdout_force_refit"] is False
    assert status["pde_validated"] is False
    assert status["visualization_ready_promoted"] is False
    assert status["visual_correspondence_verified"] is False
    assert status["paper_exact"] is False
    assert status["openai_field_identified"] is False


def test_phi10_update_changes_only_existing_phi10_and_roundtrips_public_velocity():
    data = _payload()
    parent = _parent_candidate(data)
    assert parent.sha256 == data["parent"]["candidate_sha256"]
    optimized = candidate_with_phi10(parent, data["fit"]["Phi(1,0)"])
    mode_index = parent.profile_basis.mode_indices.index((1, 0))

    for index in range(parent.profile_basis.mode_count):
        if index == mode_index:
            assert optimized.profile_basis.phi_coefficients[index] == data["fit"]["Phi(1,0)"]
        else:
            assert optimized.profile_basis.phi_coefficients[index] == parent.profile_basis.phi_coefficients[index]
        assert optimized.profile_basis.swirl_coefficients[index] == parent.profile_basis.swirl_coefficients[index]

    with pytest.raises(ValueError):
        candidate_with_phi10(parent, PROFILE_BOUNDS[1] + 1e-6)

    roundtrip = Eq45VelocityCandidate.from_dict(optimized.to_dict())
    assert roundtrip.sha256 == optimized.sha256
    assert optimized.sha256 == data["optimized_candidate"]["sha256"]
    axis = np.linspace(-0.2, 0.2, 3)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((x, y, z), axis=-1).reshape(-1, 3)
    values = roundtrip.at_points(points, np.full(points.shape[0], 0.5))
    assert values.shape == (27, 3)
    assert np.all(np.isfinite(values))
    assert np.linalg.norm(values) > 0.0

    change = data["velocity_change_vs_parent_on_holdout_probes"]
    assert 0.05 < change["relative_delta_rms"] < 0.10
    assert data["optimized_candidate"]["canonical_candidate_promoted"] is False


def test_default_phi10_optimizer_replays_checked_result_without_holdout_refit():
    data = _payload()
    optimized, report = optimize_eq45_phi10_profile_force()

    assert report.optimizer_success is True
    assert report.optimizer_function_evaluations <= data["contract"]["max_function_evaluations"]
    assert abs(report.optimized_phi_10 - data["fit"]["Phi(1,0)"]) < 5e-3
    assert abs(report.force_c - data["fit"]["force_c"]) < 5e-3
    assert report.force_a < 1e-6
    assert abs(
        report.optimized_training_rms_after_force - data["fit"]["curl_rms_after_force"]
    ) < 5e-3
    assert report.training_improvement_vs_parent > 0.15
    assert min(level.improvement_vs_parent for level in report.holdout_levels) > 0.16
    assert abs(
        report.holdout_levels[-1].optimized_rms_after_force
        - data["holdout"][-1]["optimized_curl_rms_after_force"]
    ) < 5e-3
    assert 0.05 < report.velocity_relative_delta_rms < 0.10
    assert report.new_basis_added is False
    assert report.pressure_fitted is False
    assert report.forcing_family_changed is False
    assert report.holdout_force_refit is False
    assert report.canonical_candidate_promoted is False
    assert report.pde_validated is False
    assert report.visualization_ready_promoted is False
    assert report.visual_correspondence_verified is False
    assert report.paper_exact is False
    assert report.openai_field_identified is False
    assert optimized.sha256 == report.optimized_candidate_sha256
    assert report.optimized_candidate_roundtrip_equal is True
