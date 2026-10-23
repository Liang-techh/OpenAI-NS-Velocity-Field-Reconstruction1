"""Nonlinear manufactured-solution calibration for the independent PDE verifier.

This is a verifier calibration, not a candidate-quality test.  The manufactured
field is divergence-free, has nonzero convection, pressure gradient, viscosity,
and time derivative, and is unrelated to the optimization/training loss.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .constrained_validation import residual


DEFAULT_STEPS = (0.08, 0.04, 0.02, 0.01)
DEFAULT_TIMES = (0.25, 0.5, 0.75)
NU = 0.01


def manufactured_velocity(points: np.ndarray, time: float) -> np.ndarray:
    x = np.asarray(points, dtype=float)
    xx, yy, _ = x.T
    a = 1.0 + float(time)
    return np.column_stack((
        a * np.sin(xx) * np.cos(yy),
        -a * np.cos(xx) * np.sin(yy),
        np.zeros(len(x)),
    ))


def manufactured_pressure(points: np.ndarray, time: float) -> np.ndarray:
    x = np.asarray(points, dtype=float)
    return float(time) * np.sin(np.sum(x, axis=1))


def manufactured_terms(points: np.ndarray, time: float, *, nu: float = NU) -> dict[str, np.ndarray]:
    """Return exact PDE terms for the manufactured field."""
    x = np.asarray(points, dtype=float)
    xx, yy, zz = x.T
    t = float(time)
    a = 1.0 + t
    sx, cx = np.sin(xx), np.cos(xx)
    sy, cy = np.sin(yy), np.cos(yy)
    u = manufactured_velocity(x, t)
    ut = np.column_stack((sx * cy, -cx * sy, np.zeros(len(x))))
    convection = np.column_stack((
        a * a * sx * cx,
        a * a * sy * cy,
        np.zeros(len(x)),
    ))
    gp_scalar = t * np.cos(xx + yy + zz)
    grad_pressure = np.column_stack((gp_scalar, gp_scalar, gp_scalar))
    laplacian = -2.0 * u
    viscous = -float(nu) * laplacian
    return {
        "time_derivative": ut,
        "convection": convection,
        "grad_pressure": grad_pressure,
        "laplacian": laplacian,
        "viscous": viscous,
    }


def manufactured_force(
    points: np.ndarray,
    time: float,
    *,
    nu: float = NU,
    mutation: str | None = None,
) -> np.ndarray:
    """Exact force, optionally with one deliberate sign mutation.

    Mutations emulate common verifier/reference bugs while the finite-difference
    verifier remains untouched.
    """
    terms = manufactured_terms(points, time, nu=nu)
    ut = terms["time_derivative"]
    conv = terms["convection"]
    gp = terms["grad_pressure"]
    visc = terms["viscous"]
    if mutation is None:
        return ut + conv + gp + visc
    if mutation == "convection_sign":
        return ut - conv + gp + visc
    if mutation == "pressure_sign":
        return ut + conv - gp + visc
    if mutation == "viscosity_sign":
        return ut + conv + gp - visc
    raise ValueError(f"unknown mutation: {mutation}")


def _row(points: np.ndarray, time: float, step: float, *, mutation: str | None = None) -> dict[str, float]:
    force = lambda p, t: manufactured_force(p, t, nu=NU, mutation=mutation)
    out = residual(
        manufactured_velocity,
        manufactured_pressure,
        force,
        points,
        time,
        nu=NU,
        step=step,
        time_bounds=(0.25, 0.75),
    )
    r = np.linalg.norm(out["momentum"], axis=1)
    return {
        "time": float(time),
        "step": float(step),
        "residual_max": float(np.max(r)),
        "residual_rms": float(np.sqrt(np.mean(r * r))),
        "divergence_max": float(np.max(np.abs(out["divergence"]))),
    }


def run_calibration(
    *,
    seed: int = 38117,
    point_count: int = 257,
    steps: tuple[float, ...] = DEFAULT_STEPS,
    times: tuple[float, ...] = DEFAULT_TIMES,
) -> dict:
    if point_count < 8:
        raise ValueError("point_count must be at least 8")
    if len(steps) < 3 or any(not np.isfinite(s) or s <= 0 for s in steps):
        raise ValueError("need at least three positive finite steps")
    if any(steps[i + 1] >= steps[i] for i in range(len(steps) - 1)):
        raise ValueError("steps must be strictly decreasing")

    points = np.random.default_rng(seed).uniform(-0.8, 0.8, (point_count, 3))
    rows: list[dict[str, float]] = []
    observed_orders: dict[str, list[float]] = {}
    for t in times:
        local = [_row(points, t, h) for h in steps]
        rows.extend(local)
        observed_orders[str(t)] = [
            float(np.log(local[i]["residual_max"] / local[i + 1]["residual_max"]) /
                  np.log(local[i]["step"] / local[i + 1]["step"]))
            for i in range(len(local) - 1)
        ]

    mutation_rows = {
        name: _row(points, 0.5, steps[-1], mutation=name)
        for name in ("convection_sign", "pressure_sign", "viscosity_sign")
    }
    min_order = min(min(v) for v in observed_orders.values())
    return {
        "schema_version": 1,
        "task_id": "CR006-VERIFIER-NONLINEAR-CALIBRATION-001",
        "scope": (
            "independent verifier calibration only; nonlinear manufactured "
            "divergence-free field with nonzero convection/pressure/viscosity/time terms"
        ),
        "seed": int(seed),
        "point_count": int(point_count),
        "times": [float(t) for t in times],
        "steps": [float(s) for s in steps],
        "rows": rows,
        "observed_orders": observed_orders,
        "minimum_observed_order": float(min_order),
        "mutation_rows": mutation_rows,
        "interpretation": (
            "Approximately fourth-order step convergence demonstrates the spatial "
            "finite-difference path is behaving as designed on this manufactured case. "
            "Large mutation residuals show common sign errors are not silently accepted. "
            "This does not validate any optimized candidate."
        ),
        "claims": {
            "candidate_pde_passed": False,
            "project_acceptance": False,
            "blow_up_claim": False,
        },
    }


def write_report(path: str | Path) -> dict:
    report = run_calibration()
    Path(path).write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    target = Path("artifacts/constrained/verifier_nonlinear_calibration.json")
    result = write_report(target)
    print(json.dumps(result, indent=2))
