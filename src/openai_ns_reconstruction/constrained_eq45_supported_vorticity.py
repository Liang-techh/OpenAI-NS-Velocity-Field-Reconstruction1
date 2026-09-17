"""Independent pressure-free vorticity-equation audit for the supported Eq45 child.

The support transform changes public ``[u,v,w]`` in the exterior collar, so PDE
evidence from the untapered parent cannot be inherited.  This module treats the
supported child as a black box: it serializes/reloads the child, samples only
``at_points(points, time)``, and reconstructs Cartesian derivatives with centered
finite differences outside the candidate implementation.

The pressure-free vorticity equation is the necessary condition

    omega_t + (u . grad)omega - (omega . grad)u - nu Delta omega = curl(f),
    omega = curl(u).

This increment evaluates the zero-force operator only.  The active project allows
only the preregistered two-parameter RestrictedForce family, but coefficients fitted
against the untapered parent are not silently transferred after the support transform.
A supported-child force fit therefore remains pending and ``pde_validated`` stays
false regardless of the numbers reported here.
"""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable, Protocol

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


class VelocityLike(Protocol):
    def at_points(self, points, time): ...


def _points(points) -> np.ndarray:
    value = np.asarray(points, dtype=float)
    if value.ndim != 2 or value.shape[1] != 3 or len(value) == 0:
        raise ValueError("points must have shape (n,3) with n>0")
    if not np.all(np.isfinite(value)):
        raise ValueError("points must be finite")
    return value


def _times(time, n: int) -> np.ndarray:
    value = np.asarray(time, dtype=float)
    if not np.all(np.isfinite(value)):
        raise ValueError("time must be finite")
    try:
        return np.broadcast_to(value, (n,)).astype(float, copy=False)
    except ValueError as exc:
        raise ValueError("time must be scalar or broadcast to point count") from exc


def _vector_values(function: Callable, points: np.ndarray, time: np.ndarray) -> np.ndarray:
    value = np.asarray(function(points, time), dtype=float)
    if value.shape != points.shape:
        raise ValueError("vector evaluator must return shape (n,3)")
    if not np.all(np.isfinite(value)):
        raise FloatingPointError("vector evaluator returned nonfinite values")
    return value


def _velocity(field: VelocityLike, points: np.ndarray, time: np.ndarray) -> np.ndarray:
    return _vector_values(field.at_points, points, time)


def _jacobian(function: Callable, points: np.ndarray, time: np.ndarray, step: float) -> np.ndarray:
    if not np.isfinite(step) or step <= 0:
        raise ValueError("finite-difference step must be positive and finite")
    out = np.empty((len(points), 3, 3), dtype=float)
    for axis in range(3):
        shift = np.zeros_like(points)
        shift[:, axis] = step
        plus = _vector_values(function, points + shift, time)
        minus = _vector_values(function, points - shift, time)
        out[:, :, axis] = (plus - minus) / (2.0 * step)
    return out


def _curl_from_jacobian(jacobian: np.ndarray) -> np.ndarray:
    return np.column_stack(
        (
            jacobian[:, 2, 1] - jacobian[:, 1, 2],
            jacobian[:, 0, 2] - jacobian[:, 2, 0],
            jacobian[:, 1, 0] - jacobian[:, 0, 1],
        )
    )


def _curl(function: Callable, points: np.ndarray, time: np.ndarray, step: float) -> np.ndarray:
    return _curl_from_jacobian(_jacobian(function, points, time, step))


def _omega(field: VelocityLike, points: np.ndarray, time: np.ndarray, spatial_step: float) -> np.ndarray:
    return _curl(lambda p, t: _velocity(field, p, t), points, time, spatial_step)


def vorticity_equation_terms(
    field: VelocityLike,
    points,
    time,
    *,
    spatial_step: float,
    time_step: float,
    nu: float,
) -> dict[str, np.ndarray]:
    """Return independently reconstructed zero-force vorticity-equation terms."""
    p = _points(points)
    t = _times(time, len(p))
    if not np.isfinite(time_step) or time_step <= 0:
        raise ValueError("time_step must be positive and finite")
    if not np.isfinite(nu) or nu < 0:
        raise ValueError("nu must be finite and nonnegative")

    velocity = _velocity(field, p, t)
    velocity_jacobian = _jacobian(lambda q, s: _velocity(field, q, s), p, t, spatial_step)
    omega = _curl_from_jacobian(velocity_jacobian)

    omega_t = (
        _omega(field, p, t + time_step, spatial_step)
        - _omega(field, p, t - time_step, spatial_step)
    ) / (2.0 * time_step)

    omega_gradient = np.empty((len(p), 3, 3), dtype=float)
    omega_laplacian = np.zeros((len(p), 3), dtype=float)
    for axis in range(3):
        shift = np.zeros_like(p)
        shift[:, axis] = spatial_step
        plus = _omega(field, p + shift, t, spatial_step)
        minus = _omega(field, p - shift, t, spatial_step)
        omega_gradient[:, :, axis] = (plus - minus) / (2.0 * spatial_step)
        omega_laplacian += (plus - 2.0 * omega + minus) / (spatial_step**2)

    advection = np.einsum("nj,nij->ni", velocity, omega_gradient)
    stretching = np.einsum("nj,nij->ni", omega, velocity_jacobian)
    diffusion = float(nu) * omega_laplacian
    residual = omega_t + advection - stretching - diffusion
    if not np.all(np.isfinite(residual)):
        raise FloatingPointError("vorticity-equation residual became nonfinite")

    return {
        "velocity": velocity,
        "omega": omega,
        "omega_t": omega_t,
        "advection": advection,
        "stretching": stretching,
        "diffusion": diffusion,
        "residual": residual,
    }


def _metrics(terms: dict[str, np.ndarray], mask: np.ndarray | None = None) -> dict[str, float]:
    residual = np.asarray(terms["residual"], dtype=float)
    if mask is not None:
        residual = residual[mask]
    norms = np.linalg.norm(residual, axis=1)
    rms = float(np.sqrt(np.mean(norms**2)))

    scale_squared = np.zeros(len(terms["residual"]), dtype=float)
    for key in ("omega_t", "advection", "stretching", "diffusion"):
        value = np.asarray(terms[key], dtype=float)
        scale_squared += np.sum(value * value, axis=1)
    if mask is not None:
        scale_squared = scale_squared[mask]
    scale = float(np.sqrt(np.mean(scale_squared)))
    return {
        "max": float(np.max(norms)),
        "rms": rms,
        "term_normalized_rms": rms / max(scale, np.finfo(float).tiny),
    }


def _probe_set() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fixed off-axis probes spanning plateau and support-transform collars."""
    points = np.asarray(
        [
            [0.35, 0.17, 0.25],
            [-0.42, 0.23, -0.30],
            [0.55, -0.31, 0.40],
            [-0.28, -0.48, -0.45],
            [1.70, 0.25, 0.25],
            [-1.65, 0.40, -0.35],
            [1.20, 1.25, 0.45],
            [-1.30, -1.15, -0.50],
            [0.45, 0.15, 1.72],
            [-0.40, 0.30, -1.70],
            [0.65, -0.20, 1.76],
            [-0.25, -0.55, -1.74],
            [1.55, 0.65, 1.72],
            [-1.50, 0.75, -1.70],
            [1.25, -1.10, 1.75],
            [-1.20, -1.20, -1.73],
        ],
        dtype=float,
    )
    labels = np.asarray(
        ["plateau"] * 4
        + ["radial_collar"] * 4
        + ["axial_collar"] * 4
        + ["corner_collar"] * 4,
        dtype=object,
    )
    times = np.asarray([0.35, 0.45, 0.55, 0.65] * 4, dtype=float)
    return points, times, labels


def audit_supported_eq45(
    *,
    parent_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
    time_step: float = 0.0025,
) -> dict:
    """Revalidate the support-connected child against its untapered parent."""
    root = Path(__file__).resolve().parents[2]
    parent_path = Path(parent_path or root / "artifacts/constrained/eq45_velocity_candidate_seed.json")
    constraints_path = Path(constraints_path or root / "configs/constraints.json")

    parent = Eq45VelocityCandidate.load_json(parent_path)
    constructed_child = Eq45SupportedVelocityCandidate(parent=parent)
    with TemporaryDirectory() as tmp:
        child_path = Path(tmp) / "supported.json"
        constructed_child.save_json(child_path)
        child = Eq45SupportedVelocityCandidate.load_json(child_path)

    if child.sha256 != constructed_child.sha256:
        raise AssertionError("supported-child serialization changed candidate identity")

    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    nu = float(constraints["nu"])
    steps = tuple(float(step) for step in constraints["validation"]["derivative_steps"])
    if len(steps) < 3 or any(right >= left for left, right in zip(steps, steps[1:])):
        raise ValueError("registered derivative ladder must contain >=3 decreasing levels")

    points, times, labels = _probe_set()
    regions = tuple(dict.fromkeys(labels.tolist()))
    rows = []
    for step in steps:
        parent_terms = vorticity_equation_terms(
            parent,
            points,
            times,
            spatial_step=step,
            time_step=time_step,
            nu=nu,
        )
        child_terms = vorticity_equation_terms(
            child,
            points,
            times,
            spatial_step=step,
            time_step=time_step,
            nu=nu,
        )
        parent_metrics = _metrics(parent_terms)
        child_metrics = _metrics(child_terms)
        by_region = {}
        for region in regions:
            mask = labels == region
            by_region[region] = {
                "parent": _metrics(parent_terms, mask),
                "supported_child": _metrics(child_terms, mask),
            }
        rows.append(
            {
                "spatial_step": step,
                "parent": parent_metrics,
                "supported_child": child_metrics,
                "child_vs_parent_rms_fractional_change": (
                    child_metrics["rms"] - parent_metrics["rms"]
                ) / max(parent_metrics["rms"], np.finfo(float).tiny),
                "by_region": by_region,
            }
        )

    finest = rows[-1]
    plateau_parent = finest["by_region"]["plateau"]["parent"]["rms"]
    plateau_child = finest["by_region"]["plateau"]["supported_child"]["rms"]
    plateau_rel = abs(plateau_child - plateau_parent) / max(
        plateau_parent, np.finfo(float).tiny
    )

    return {
        "schema": "eq45_supported_vorticity_revalidation_v1",
        "claim_scope": "independent pressure_free_vorticity_equation_on_supported_child",
        "parent_sha256": parent.sha256,
        "supported_child_sha256": child.sha256,
        "serialized_reloaded_child": True,
        "velocity_access": "public_at_points_only",
        "nu": nu,
        "spatial_steps": list(steps),
        "fixed_time_step": float(time_step),
        "probe_count": int(len(points)),
        "regions": list(regions),
        "rows": rows,
        "finest_supported_zero_force_rms": finest["supported_child"]["rms"],
        "finest_supported_zero_force_max": finest["supported_child"]["max"],
        "finest_supported_term_normalized_rms": finest["supported_child"]["term_normalized_rms"],
        "finest_child_vs_parent_rms_fractional_change": finest["child_vs_parent_rms_fractional_change"],
        "finest_plateau_parent_child_rms_relative_difference": plateau_rel,
        "forcing_family": constraints["forcing"]["mode"],
        "supported_force_coefficients_status": "pending_refit_and_independent_revalidation",
        "untapered_fitted_force_transferred": False,
        "pde_thresholds": {
            "max": float(constraints["validation"]["thresholds"]["pde_residual_max"]),
            "L2": float(constraints["validation"]["thresholds"]["pde_residual_L2"]),
        },
        "formal_pde_gate_assessed": False,
        "reason_formal_pde_gate_unassessed": (
            "pressure and supported-child restricted-force coefficients are not bound here; "
            "probe RMS is not the preregistered volume-weighted spatial L2"
        ),
        "velocity_changed": False,
        "physical_support_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    print(json.dumps(audit_supported_eq45(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
