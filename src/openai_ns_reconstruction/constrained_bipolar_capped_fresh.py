"""Fresh held-out refinement audit for the energy-capped bipolar candidate.

This is an independent validation consumer: it reloads the saved candidate,
reads velocity only through the public ``at_points`` interface, and never uses
optimizer/training tensors. The full vector diagnostic fixes p=0 because no
pressure artifact exists for this candidate; the azimuthal projection is a
necessary pressure-independent condition for axisymmetric pressure.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .constrained_force import RestrictedForce
from .constrained_validation import residual
from .eq45_supported_delivery import Eq45SupportedDeliveryField


VOLUME = 64.0
DEFAULT_SEED = 9172631
DEFAULT_STEPS = (0.02, 0.01, 0.005)
DEFAULT_TIMES = (0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75)


def _theta(points: np.ndarray) -> np.ndarray:
    radius = np.hypot(points[:, 0], points[:, 1])
    if np.any(radius <= 1e-14):
        raise ValueError("theta direction undefined on the axis")
    return np.column_stack(
        (-points[:, 1] / radius, points[:, 0] / radius, np.zeros(len(points)))
    )


def _volume_l2(values: np.ndarray) -> float:
    return float(np.sqrt(VOLUME * np.mean(np.asarray(values, dtype=float) ** 2)))


def _row(field, force, points, time: float, step: float, nu: float) -> dict:
    zero_pressure = lambda x, t: np.zeros(len(x), dtype=float)
    result = residual(
        field.at_points,
        zero_pressure,
        force,
        points,
        time,
        nu=nu,
        step=step,
        time_bounds=(0.25, 0.75),
    )
    momentum = result["momentum"]
    divergence = result["divergence"]
    momentum_norm = np.linalg.norm(momentum, axis=1)
    theta_residual = np.sum(momentum * _theta(points), axis=1)
    momentum_l2 = _volume_l2(momentum_norm)
    theta_l2 = _volume_l2(theta_residual)
    return {
        "time": float(time),
        "step": float(step),
        "divergence_sampled_max": float(np.max(np.abs(divergence))),
        "divergence_volume_L2": _volume_l2(divergence),
        "zero_pressure_momentum_sampled_max": float(np.max(momentum_norm)),
        "zero_pressure_momentum_volume_L2": momentum_l2,
        "theta_momentum_sampled_max": float(np.max(np.abs(theta_residual))),
        "theta_momentum_volume_L2": theta_l2,
        "theta_fraction_of_zero_pressure_L2": float(theta_l2 / momentum_l2)
        if momentum_l2 > 0
        else float("inf"),
    }


def _refinement(rows: list[dict], key: str) -> list[dict]:
    out = []
    times = sorted({row["time"] for row in rows})
    for time in times:
        selected = sorted(
            (row for row in rows if row["time"] == time),
            key=lambda row: row["step"],
            reverse=True,
        )
        if len(selected) != 3:
            raise ValueError("refinement audit requires exactly three derivative steps")
        values = [row[key] for row in selected]
        steps = [row["step"] for row in selected]
        coarse_to_mid = abs(values[1] - values[0]) / max(abs(values[1]), 1e-300)
        mid_to_fine = abs(values[2] - values[1]) / max(abs(values[2]), 1e-300)
        out.append(
            {
                "time": time,
                "steps": steps,
                "values": values,
                "coarse_to_mid_relative_change": float(coarse_to_mid),
                "mid_to_fine_relative_change": float(mid_to_fine),
                "fine_to_coarse_ratio": float(values[2] / values[0])
                if values[0] != 0
                else float("inf"),
            }
        )
    return out


def run(
    *,
    candidate_path: str | Path = "artifacts/bipolar_joint_capped/candidate.json",
    training_report_path: str | Path = "artifacts/bipolar_joint_capped/report.json",
    config_path: str | Path = "configs/constraints.json",
    output: str | Path = "artifacts/bipolar_joint_capped/fresh_refinement.json",
    seed: int = DEFAULT_SEED,
    points: int = 2048,
    times=DEFAULT_TIMES,
    steps=DEFAULT_STEPS,
) -> dict:
    cfg = json.loads(Path(config_path).read_text())
    training = json.loads(Path(training_report_path).read_text())
    field = Eq45SupportedDeliveryField.load_candidate(candidate_path)
    force_parameters = training["parameters"][-2:]
    force = RestrictedForce(a=float(force_parameters[0]), c=float(force_parameters[1]))

    rng = np.random.default_rng(seed)
    box = np.asarray(cfg["domain"]["evaluation_box"], dtype=float)
    sample = rng.uniform(box[:, 0], box[:, 1], (points, 3))
    rows = [
        _row(field, force, sample, float(time), float(step), float(cfg["nu"]))
        for step in steps
        for time in times
    ]

    fine_step = min(float(step) for step in steps)
    fine = [row for row in rows if row["step"] == fine_step]
    limits = cfg["validation"]["thresholds"]
    worst_div_max = max(row["divergence_sampled_max"] for row in fine)
    worst_div_l2 = max(row["divergence_volume_L2"] for row in fine)
    worst_theta_max = max(row["theta_momentum_sampled_max"] for row in fine)
    worst_theta_l2 = max(row["theta_momentum_volume_L2"] for row in fine)

    mutation_strength = 0.03

    def mutated_velocity(x, t):
        values = np.array(field.at_points(x, t), dtype=float, copy=True)
        values[:, 0] += mutation_strength * np.asarray(x)[:, 0]
        return values

    zero_pressure = lambda x, t: np.zeros(len(x), dtype=float)
    baseline_mutation_time = residual(
        field.at_points,
        zero_pressure,
        force,
        sample,
        0.5,
        nu=float(cfg["nu"]),
        step=fine_step,
        time_bounds=(0.25, 0.75),
    )
    mutation_result = residual(
        mutated_velocity,
        zero_pressure,
        force,
        sample,
        0.5,
        nu=float(cfg["nu"]),
        step=fine_step,
        time_bounds=(0.25, 0.75),
    )
    baseline_div_l2 = _volume_l2(baseline_mutation_time["divergence"])
    mutation_div = mutation_result["divergence"]
    mutation_div_l2 = _volume_l2(mutation_div)
    mutation = {
        "definition": "u_x <- u_x + 0.03*x",
        "expected_added_divergence": mutation_strength,
        "time": 0.5,
        "step": fine_step,
        "divergence_sampled_max": float(np.max(np.abs(mutation_div))),
        "divergence_volume_L2": mutation_div_l2,
        "baseline_divergence_volume_L2": baseline_div_l2,
        "caught": bool(mutation_div_l2 > 100.0 * baseline_div_l2),
    }

    report = {
        "candidate": str(candidate_path),
        "candidate_sha256": field.sha256,
        "training_report": str(training_report_path),
        "seed": int(seed),
        "points": int(points),
        "sampling": "fresh uniform Cartesian held-out points in [-2,2]^3",
        "times": [float(t) for t in times],
        "steps": [float(h) for h in steps],
        "nu": float(cfg["nu"]),
        "force": {"a": float(force.a), "c": float(force.c)},
        "rows": rows,
        "refinement": {
            "divergence_volume_L2": _refinement(rows, "divergence_volume_L2"),
            "theta_momentum_volume_L2": _refinement(rows, "theta_momentum_volume_L2"),
            "zero_pressure_momentum_volume_L2": _refinement(
                rows, "zero_pressure_momentum_volume_L2"
            ),
        },
        "finest_worst": {
            "step": fine_step,
            "divergence_sampled_max": worst_div_max,
            "divergence_volume_L2": worst_div_l2,
            "theta_momentum_sampled_max": worst_theta_max,
            "theta_momentum_volume_L2": worst_theta_l2,
            "divergence_max_over_threshold": worst_div_max / limits["divergence_max"],
            "divergence_L2_over_threshold": worst_div_l2 / limits["divergence_L2"],
            "theta_max_over_full_PDE_threshold": worst_theta_max
            / limits["pde_residual_max"],
            "theta_L2_over_full_PDE_threshold": worst_theta_l2
            / limits["pde_residual_L2"],
        },
        "preregistered_thresholds_applied_without_change": True,
        "fresh_sample_divergence_thresholds_passed": bool(
            worst_div_max <= limits["divergence_max"]
            and worst_div_l2 <= limits["divergence_L2"]
        ),
        "formal_registered_divergence_gate_assessed": False,
        "necessary_axisymmetric_pressure_PDE_condition_assessed": True,
        "necessary_axisymmetric_pressure_PDE_condition_passed": bool(
            worst_theta_max <= limits["pde_residual_max"]
            and worst_theta_l2 <= limits["pde_residual_L2"]
        ),
        "formal_full_momentum_gate_assessed": False,
        "mutation": mutation,
        "classification": {
            "visualization_candidate_only": True,
            "pde_validated": False,
            "visual_correspondence_verified": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
        "scope": (
            "Candidate is reloaded from its artifact and all velocity values are obtained "
            "through public at_points. The preregistered thresholds are reused unchanged, "
            "but this fresh seed/2048-point generalization audit is not the registered "
            "4096-point gate. Zero-pressure full-vector momentum is diagnostic only because "
            "no pressure artifact is supplied. For axisymmetric pressure, the theta projection "
            "is pressure-independent and therefore a necessary PDE condition, not a sufficient "
            "full Navier-Stokes validation. Monte Carlo L2 values are sampled volume estimates, "
            "not certified global bounds."
        ),
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["finest_worst"], indent=2))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", default="artifacts/bipolar_joint_capped/candidate.json")
    parser.add_argument("--training-report", default="artifacts/bipolar_joint_capped/report.json")
    parser.add_argument("--config", default="configs/constraints.json")
    parser.add_argument("--output", default="artifacts/bipolar_joint_capped/fresh_refinement.json")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--points", type=int, default=2048)
    args = parser.parse_args()
    run(
        candidate_path=args.candidate,
        training_report_path=args.training_report,
        config_path=args.config,
        output=args.output,
        seed=args.seed,
        points=args.points,
    )


if __name__ == "__main__":
    main()
