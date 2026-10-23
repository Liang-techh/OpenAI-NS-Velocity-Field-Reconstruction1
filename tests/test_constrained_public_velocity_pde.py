import numpy as np
import pytest

from openai_ns_reconstruction.constrained_public_velocity_pde import (
    audit_packaged_candidate_smoke,
    audit_public_velocity_refinement,
)


THRESHOLDS = {
    "pde_residual_max": 1e-5,
    "pde_residual_L2": 1e-5,
    "divergence_max": 1e-5,
    "divergence_L2": 1e-5,
}


def zero_scalar(points, time):
    points = np.asarray(points, dtype=float)
    return np.zeros(points.shape[:-1], dtype=float)


def zero_vector(points, time):
    return np.zeros_like(np.asarray(points, dtype=float))


class TrackingShearField:
    def __init__(self):
        self.calls = 0
        self.sha256 = "synthetic-shear"

    @property
    def candidate(self):
        raise AssertionError("public audit must not access an internal candidate velocity")

    def at_points(self, points, time):
        self.calls += 1
        points = np.asarray(points, dtype=float)
        result = np.zeros_like(points)
        result[..., 0] = points[..., 1]
        return result


class DivergentMutationField:
    sha256 = "synthetic-divergent-mutation"

    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        result = np.zeros_like(points)
        result[..., 0] = points[..., 0]
        return result


def test_public_api_only_shear_passes_three_level_refinement():
    rng = np.random.default_rng(771)
    points = rng.uniform(-0.8, 0.8, size=(31, 3))
    field = TrackingShearField()
    report = audit_public_velocity_refinement(
        field,
        zero_scalar,
        zero_vector,
        points,
        (0.3, 0.5, 0.7),
        (0.04, 0.02, 0.01),
        nu=0.01,
        time_bounds=(0.25, 0.75),
        volume=8.0,
        thresholds=THRESHOLDS,
    )
    assert field.calls > 0
    assert report["velocity_source"] == "public field.at_points only"
    assert report["sampled_finest_thresholds_passed"] is True
    assert report["pde_validated"] is False
    assert len(report["rows"]) == 9
    assert max(row["residual_sampled_max"] for row in report["rows"]) < 1e-10
    assert max(row["divergence_sampled_max"] for row in report["rows"]) < 1e-10


def test_divergent_public_velocity_mutation_is_detected():
    points = np.array(
        [
            [-0.7, -0.2, 0.1],
            [-0.2, 0.4, -0.5],
            [0.3, -0.6, 0.2],
            [0.8, 0.1, -0.3],
        ],
        dtype=float,
    )
    report = audit_public_velocity_refinement(
        DivergentMutationField(),
        zero_scalar,
        zero_vector,
        points,
        (0.5,),
        (0.04, 0.02, 0.01),
        nu=0.01,
        time_bounds=(0.25, 0.75),
        volume=8.0,
        thresholds=THRESHOLDS,
    )
    assert report["sampled_finest_thresholds_passed"] is False
    finest = [row for row in report["rows"] if row["step"] == 0.01][0]
    assert finest["divergence_sampled_max"] == pytest.approx(1.0, abs=1e-11)
    assert finest["divergence_max_over_threshold"] > 9e4
    assert finest["momentum_component_rms"][0] > 0.1


def test_packaged_candidate_smoke_loads_artifact_and_keeps_truth_boundary():
    report = audit_packaged_candidate_smoke(
        seed=914033,
        point_count=8,
        times=(0.5,),
    )
    assert isinstance(report["candidate_sha256"], str)
    assert len(report["candidate_sha256"]) == 64
    assert report["seed"] != report["training_seed"]
    assert report["velocity_source"] == "public field.at_points only"
    assert report["pde_validated"] is False
    assert len(report["rows"]) == 3
    assert all(np.isfinite(row["residual_sampled_max"]) for row in report["rows"])
    assert all(np.isfinite(row["residual_L2_over_threshold"]) for row in report["rows"])


def test_invalid_inputs_fail_closed():
    field = TrackingShearField()
    points = np.zeros((8, 3))
    with pytest.raises(ValueError, match="three"):
        audit_public_velocity_refinement(
            field,
            zero_scalar,
            zero_vector,
            points,
            (0.5,),
            (0.02, 0.01),
            nu=0.01,
            time_bounds=(0.25, 0.75),
            volume=8.0,
            thresholds=THRESHOLDS,
        )
    broken = dict(THRESHOLDS)
    broken["pde_residual_max"] = 0.0
    with pytest.raises(ValueError, match="positive"):
        audit_public_velocity_refinement(
            field,
            zero_scalar,
            zero_vector,
            points,
            (0.5,),
            (0.04, 0.02, 0.01),
            nu=0.01,
            time_bounds=(0.25, 0.75),
            volume=8.0,
            thresholds=broken,
        )
