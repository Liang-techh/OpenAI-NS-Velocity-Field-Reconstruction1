import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_profile_force_optimization import (
    PROFILE_BOUNDS,
    candidate_with_profile_pair,
    optimize_eq45_profile_force,
)
from openai_ns_reconstruction.constrained_restricted_force_curl_capacity import FORCE_BOUNDS


ARTIFACT = (
    Path(__file__).resolve().parents[1]
    / "artifacts"
    / "constrained"
    / "eq45_profile_force_optimization.json"
)


def _payload():
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_checked_profile_force_result_is_bounded_and_truthful():
    data = _payload()
    assert data["schema"] == "eq45_profile_force_optimization_v1"
    assert data["task_id"] == "CR005-EQ45-PROFILE-FORCE-OPT-015"
    assert data["contract"]["profile_bounds"] == list(PROFILE_BOUNDS)
    assert data["contract"]["force_bounds"] == list(FORCE_BOUNDS)
    assert data["contract"]["training_holdout_separate"] is True
    assert data["contract"]["max_function_evaluations"] == 80

    fit = data["fit"]
    assert fit["optimizer_success"] is True
    assert fit["optimizer_function_evaluations"] <= 80
    assert PROFILE_BOUNDS[0] <= fit["Phi_02"] <= PROFILE_BOUNDS[1]
    assert PROFILE_BOUNDS[0] <= fit["F_02"] <= PROFILE_BOUNDS[1]
    assert FORCE_BOUNDS[0] <= fit["force_a"] <= FORCE_BOUNDS[1]
    assert FORCE_BOUNDS[0] <= fit["force_c"] <= FORCE_BOUNDS[1]
    assert fit["curl_rms_after_force"] < data["seed"]["training_curl_rms_after_force"]

    levels = data["holdout"]
    assert [level["step"] for level in levels] == [0.02, 0.01, 0.005]
    assert min(level["improvement_vs_seed"] for level in levels) > 0.50
    assert levels[-1]["optimized_curl_rms_after_force"] > 1.0

    status = data["status"]
    assert status["forcing_family_changed"] is False
    assert status["pressure_fitted"] is False
    assert status["pde_validated"] is False
    assert status["visualization_ready_promoted"] is False
    assert status["visual_correspondence_verified"] is False
    assert status["paper_exact"] is False
    assert status["openai_field_identified"] is False


def test_optimized_pair_roundtrips_and_exports_public_velocity_grid_smoke():
    data = _payload()
    fit = data["fit"]
    seed = Eq45VelocityCandidate.seed()
    optimized = candidate_with_profile_pair(seed, fit["Phi_02"], fit["F_02"])

    mode_index = seed.profile_basis.mode_indices.index((0, 2))
    for index in range(seed.profile_basis.mode_count):
        if index == mode_index:
            assert optimized.profile_basis.phi_coefficients[index] == fit["Phi_02"]
            assert optimized.profile_basis.swirl_coefficients[index] == fit["F_02"]
        else:
            assert optimized.profile_basis.phi_coefficients[index] == seed.profile_basis.phi_coefficients[index]
            assert optimized.profile_basis.swirl_coefficients[index] == seed.profile_basis.swirl_coefficients[index]

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

    velocity_change = data["velocity_change_on_holdout_probes"]
    assert 0.15 < velocity_change["relative_delta_rms"] < 0.25
    assert data["optimized_candidate"]["canonical_seed_promoted"] is False


def test_default_optimizer_replays_checked_capacity_result():
    data = _payload()
    fit = data["fit"]
    optimized, report = optimize_eq45_profile_force()

    assert report.optimizer_success is True
    assert report.optimizer_function_evaluations <= data["contract"]["max_function_evaluations"]
    assert abs(report.phi_02 - fit["Phi_02"]) < 5e-3
    assert abs(report.F_02 - fit["F_02"]) < 5e-3
    assert abs(report.force_c - fit["force_c"]) < 5e-3
    assert report.force_a < 1e-6
    assert abs(report.optimized_training_rms_after_force - fit["curl_rms_after_force"]) < 5e-3
    assert min(level.improvement_vs_seed for level in report.holdout_levels) > 0.53
    assert report.holdout_levels[-1].optimized_rms_after_force > 9.0
    assert 0.15 < report.velocity_relative_delta_rms < 0.25
    assert report.pde_validated is False
    assert report.visualization_ready_promoted is False
    assert report.paper_exact is False
    assert report.openai_field_identified is False
    assert optimized.sha256 == report.optimized_candidate_sha256
    assert report.optimized_candidate_roundtrip_equal is True
