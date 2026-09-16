import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_force import RestrictedForce
from openai_ns_reconstruction.constrained_public_velocity_time_refinement import (
    audit_public_velocity_time_refinement,
    temporal_refinement_residual,
)
from openai_ns_reconstruction.velocity_components import VelocityField


THRESHOLDS = {
    "pde_residual_max": 1e-3,
    "pde_residual_L2": 1e-3,
    "divergence_max": 1e-5,
    "divergence_L2": 1e-5,
}


class CubicShearField:
    @property
    def candidate(self):
        raise AssertionError("validation must use public at_points, not candidate internals")

    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        result = np.zeros_like(points)
        result[:, 0] = float(time) ** 3 * np.sin(points[:, 1])
        return result


def zero_pressure(points, time):
    return np.zeros(len(points), dtype=float)


def exact_force(points, time, nu=0.01):
    points = np.asarray(points, dtype=float)
    result = np.zeros_like(points)
    t = float(time)
    result[:, 0] = (3.0 * t * t + nu * t**3) * np.sin(points[:, 1])
    return result


def test_time_step_only_refinement_is_second_order_at_interior_and_endpoints():
    rng = np.random.default_rng(914053)
    points = rng.uniform(-1.0, 1.0, size=(257, 3))
    report = audit_public_velocity_time_refinement(
        CubicShearField(),
        zero_pressure,
        exact_force,
        points,
        times=(0.25, 0.5, 0.75),
        time_steps=(0.08, 0.04, 0.02),
        nu=0.01,
        spatial_step=0.0025,
        time_bounds=(0.25, 0.75),
        volume=8.0,
        thresholds=THRESHOLDS,
    )

    assert report["velocity_source"] == "public field.at_points only"
    assert report["pde_validated"] is False
    for block in report["time_refinement"]:
        maxima = [row["residual_sampled_max"] for row in block["rows"]]
        assert maxima[0] / maxima[1] == pytest.approx(4.0, rel=2e-6)
        assert maxima[1] / maxima[2] == pytest.approx(4.0, rel=8e-6)
        assert block["finest_to_coarsest_residual_max_ratio"] == pytest.approx(1.0 / 16.0, rel=1e-5)
        assert max(row["divergence_sampled_max"] for row in block["rows"]) < 1e-11
        deltas = block["adjacent_time_step_deltas"]
        assert deltas[0]["momentum_delta_rms"] > deltas[1]["momentum_delta_rms"]


def test_wrong_force_creates_refinement_plateau_instead_of_fake_convergence():
    rng = np.random.default_rng(914054)
    points = rng.uniform(-1.0, 1.0, size=(193, 3))

    def wrong_force(points, time):
        return -exact_force(points, time)

    report = audit_public_velocity_time_refinement(
        CubicShearField(),
        zero_pressure,
        wrong_force,
        points,
        times=(0.5,),
        time_steps=(0.08, 0.04, 0.02),
        nu=0.01,
        spatial_step=0.0025,
        time_bounds=(0.25, 0.75),
        volume=8.0,
        thresholds=THRESHOLDS,
    )
    rows = report["time_refinement"][0]["rows"]
    assert rows[-1]["residual_sampled_max"] > 1.0
    assert report["time_refinement"][0]["finest_to_coarsest_residual_max_ratio"] > 0.99
    assert rows[-1]["residual_max_over_threshold"] > 1000.0


def test_public_temporal_residual_fail_closes_on_bad_inputs():
    points = np.zeros((4, 3))
    with pytest.raises(ValueError):
        audit_public_velocity_time_refinement(
            CubicShearField(),
            zero_pressure,
            exact_force,
            points,
            times=(0.5,),
            time_steps=(0.02, 0.01),
            nu=0.01,
            spatial_step=0.005,
            time_bounds=(0.25, 0.75),
            volume=64.0,
            thresholds=THRESHOLDS,
        )

    class BadField:
        def at_points(self, points, time):
            return np.zeros((len(points), 2))

    with pytest.raises(ValueError):
        temporal_refinement_residual(
            BadField(),
            zero_pressure,
            exact_force,
            points,
            0.5,
            time_step=0.01,
        )

    with pytest.raises(ValueError):
        temporal_refinement_residual(
            CubicShearField(),
            zero_pressure,
            exact_force,
            points,
            0.25,
            time_step=0.3,
        )


def test_packaged_velocity_smoke_uses_frozen_public_artifact():
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root / "configs" / "constraints.json").read_text())
    training = json.loads((root / "artifacts" / "constrained" / "coupled_joint" / "training.json").read_text())
    field = VelocityField()
    force = RestrictedForce(**training["force"])
    rng = np.random.default_rng(914059)
    box = np.asarray(config["domain"]["evaluation_box"], dtype=float)
    points = rng.uniform(box[:, 0], box[:, 1], size=(8, 3))

    report = audit_public_velocity_time_refinement(
        field,
        field.candidate.pressure,
        force,
        points,
        times=(0.5,),
        time_steps=config["validation"]["derivative_steps"],
        nu=config["nu"],
        spatial_step=min(config["validation"]["derivative_steps"]),
        time_bounds=config["domain"]["time_interval"],
        volume=float(np.prod(box[:, 1] - box[:, 0])),
        thresholds=config["validation"]["thresholds"],
    )

    assert report["candidate_sha256"] == field.sha256
    assert report["pde_validated"] is False
    assert report["visualization_candidate_only_allowed"] is True
    assert len(report["rows"]) == 3
    assert all(np.isfinite(row["residual_sampled_max"]) for row in report["rows"])
