import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_force import RestrictedForce
from openai_ns_reconstruction.constrained_public_velocity_alt_pde import (
    audit_public_velocity_operator_crosscheck,
    second_order_public_residual,
)
from openai_ns_reconstruction.velocity_components import VelocityField


class PublicField:
    def __init__(self, fn):
        self.fn = fn
        self.sha256 = "synthetic"

    @property
    def candidate(self):
        raise AssertionError("validation must not read the field's internal candidate")

    def at_points(self, points, time):
        return self.fn(np.asarray(points, dtype=float), float(time))


def zero_pressure(points, time):
    return np.zeros(len(points), dtype=float)


def test_second_order_operator_converges_on_manufactured_sine_shear():
    nu = 0.01

    def velocity(points, time):
        out = np.zeros_like(points)
        out[:, 0] = np.sin(points[:, 1])
        return out

    def force(points, time):
        out = np.zeros_like(points)
        out[:, 0] = nu * np.sin(points[:, 1])
        return out

    field = PublicField(velocity)
    rng = np.random.default_rng(914041)
    points = rng.uniform(-0.9, 0.9, size=(257, 3))
    errors = []
    for step in (0.08, 0.04, 0.02):
        result = second_order_public_residual(
            field,
            zero_pressure,
            force,
            points,
            0.5,
            nu=nu,
            step=step,
            time_bounds=(0.25, 0.75),
        )
        errors.append(np.max(np.linalg.norm(result["momentum"], axis=1)))
        assert np.max(np.abs(result["divergence"])) < 1e-12
    orders = [np.log(errors[i] / errors[i + 1]) / np.log(2.0) for i in range(2)]
    assert min(orders) > 1.99
    assert max(orders) < 2.01


def test_crosscheck_uses_public_interface_and_preserves_thresholds():
    nu = 0.01

    def velocity(points, time):
        out = np.zeros_like(points)
        out[:, 0] = np.sin(points[:, 1])
        return out

    def force(points, time):
        out = np.zeros_like(points)
        out[:, 0] = nu * np.sin(points[:, 1])
        return out

    thresholds = {
        "pde_residual_max": 1e-3,
        "pde_residual_L2": 1e-3,
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
    }
    points = np.random.default_rng(914043).uniform(-0.8, 0.8, size=(64, 3))
    report = audit_public_velocity_operator_crosscheck(
        PublicField(velocity),
        zero_pressure,
        force,
        points,
        (0.5,),
        (0.08, 0.04, 0.02),
        nu=nu,
        time_bounds=(0.25, 0.75),
        volume=8.0,
        thresholds=thresholds,
    )
    assert report["velocity_source"] == "public field.at_points only"
    assert report["thresholds_unchanged"] == thresholds
    assert report["pde_validated"] is False
    assert len(report["rows"]) == 3
    deltas = [row["operator_delta_residual_rms"] for row in report["rows"]]
    assert deltas[2] < deltas[1] < deltas[0]


def test_force_sign_mutation_is_detected_by_both_operators():
    nu = 0.01

    def velocity(points, time):
        out = np.zeros_like(points)
        out[:, 0] = np.sin(points[:, 1])
        return out

    def wrong_force(points, time):
        out = np.zeros_like(points)
        out[:, 0] = -nu * np.sin(points[:, 1])
        return out

    points = np.column_stack((np.zeros(25), np.linspace(-1.2, 1.2, 25), np.zeros(25)))
    report = audit_public_velocity_operator_crosscheck(
        PublicField(velocity),
        zero_pressure,
        wrong_force,
        points,
        (0.5,),
        (0.08, 0.04, 0.02),
        nu=nu,
        time_bounds=(0.25, 0.75),
        volume=8.0,
        thresholds={
            "pde_residual_max": 1e-3,
            "pde_residual_L2": 1e-3,
            "divergence_max": 1e-5,
            "divergence_L2": 1e-5,
        },
    )
    fine = report["rows"][-1]
    assert fine["second_order_residual_max_over_threshold"] > 10.0
    assert fine["fourth_order_residual_max_over_threshold"] > 10.0


def test_packaged_candidate_smoke_uses_public_velocity_artifact():
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root / "configs" / "constraints.json").read_text(encoding="utf-8"))
    training = json.loads(
        (root / "artifacts" / "constrained" / "coupled_joint" / "training.json").read_text(
            encoding="utf-8"
        )
    )
    seed = 914047
    assert seed != int(config["optimization"]["seed"])
    assert seed != int(config["validation"]["seed"])
    box = np.asarray(config["domain"]["evaluation_box"], dtype=float)
    points = np.random.default_rng(seed).uniform(box[:, 0], box[:, 1], size=(16, 3))
    field = VelocityField()
    force = RestrictedForce(**training["force"])
    report = audit_public_velocity_operator_crosscheck(
        field,
        field.candidate.pressure,
        force,
        points,
        (0.5,),
        config["validation"]["derivative_steps"],
        nu=config["nu"],
        time_bounds=config["domain"]["time_interval"],
        volume=float(np.prod(box[:, 1] - box[:, 0])),
        thresholds=config["validation"]["thresholds"],
    )
    assert report["candidate_sha256"] == field.sha256
    assert report["velocity_source"] == "public field.at_points only"
    assert report["pde_validated"] is False
    assert len(report["rows"]) == 3
    for row in report["rows"]:
        assert np.isfinite(row["operator_delta_residual_rms"])
        assert np.isfinite(row["operator_delta_divergence_rms"])
        assert np.isfinite(row["second_order_residual_max_over_threshold"])


def test_fail_closed_on_too_few_levels_or_malformed_public_velocity():
    good = PublicField(lambda points, time: np.zeros_like(points))
    points = np.zeros((8, 3))
    common = dict(
        pressure=zero_pressure,
        force=lambda p, t: np.zeros_like(p),
        points=points,
        times=(0.5,),
        nu=0.01,
        time_bounds=(0.25, 0.75),
        volume=8.0,
        thresholds={
            "pde_residual_max": 1e-3,
            "pde_residual_L2": 1e-3,
            "divergence_max": 1e-5,
            "divergence_L2": 1e-5,
        },
    )
    with pytest.raises(ValueError, match="three"):
        audit_public_velocity_operator_crosscheck(good, steps=(0.02, 0.01), **common)

    bad = PublicField(lambda points, time: np.zeros((len(points), 2)))
    with pytest.raises(ValueError, match="public velocity"):
        audit_public_velocity_operator_crosscheck(bad, steps=(0.04, 0.02, 0.01), **common)
