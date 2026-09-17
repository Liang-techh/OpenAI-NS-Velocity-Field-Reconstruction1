"""Independent pressure-free vorticity-equation audit for callable velocity fields.

The audit deliberately differentiates only public ``at_points(points, time)`` velocity
samples.  It does not reuse an optimizer loss, candidate profile derivatives, pressure
fit, or the momentum-residual implementation.  For incompressible flow, taking curl
of the momentum equation removes pressure and gives the necessary condition

    omega_t + (u . grad) omega - (omega . grad) u - nu Delta omega = curl(f),
    omega = curl(u).

A small value is not sufficient for full Navier--Stokes validation.  A large, stable
value is useful independent evidence that pressure alone cannot close the candidate.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Callable, Protocol

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_force import RestrictedForce


class VelocityLike(Protocol):
    def at_points(self, points, time): ...


@dataclass(frozen=True)
class VorticityResidualRow:
    spatial_step: float
    zero_force_max: float
    zero_force_rms: float
    zero_force_relative_rms: float
    restricted_force_max: float
    restricted_force_rms: float
    restricted_force_relative_rms: float
    force_rms_reduction_fraction: float

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


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
        raise ValueError("time must be scalar or broadcast to the point count") from exc


def _vector_values(function: Callable, points: np.ndarray, time: np.ndarray) -> np.ndarray:
    value = np.asarray(function(points, time), dtype=float)
    if value.shape != points.shape:
        raise ValueError("vector evaluator must return shape (n,3)")
    if not np.all(np.isfinite(value)):
        raise FloatingPointError("vector evaluator returned nonfinite values")
    return value


def _velocity(field: VelocityLike, points: np.ndarray, time: np.ndarray) -> np.ndarray:
    return _vector_values(field.at_points, points, time)


def _jacobian(
    function: Callable,
    points: np.ndarray,
    time: np.ndarray,
    step: float,
) -> np.ndarray:
    """Return J[n, component, axis] by centered Cartesian differences."""
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


def _omega(
    field: VelocityLike,
    points: np.ndarray,
    time: np.ndarray,
    spatial_step: float,
) -> np.ndarray:
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
    """Compute the pressure-free vorticity equation before subtracting curl(force)."""
    p = _points(points)
    t = _times(time, len(p))
    if not np.isfinite(nu) or nu < 0:
        raise ValueError("nu must be finite and nonnegative")
    if not np.isfinite(time_step) or time_step <= 0:
        raise ValueError("time_step must be positive and finite")

    velocity = _velocity(field, p, t)
    velocity_jacobian = _jacobian(
        lambda q, s: _velocity(field, q, s), p, t, spatial_step
    )
    omega = _curl_from_jacobian(velocity_jacobian)

    omega_plus_t = _omega(field, p, t + time_step, spatial_step)
    omega_minus_t = _omega(field, p, t - time_step, spatial_step)
    omega_t = (omega_plus_t - omega_minus_t) / (2.0 * time_step)

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
    residual_without_force = omega_t + advection - stretching - diffusion
    if not np.all(np.isfinite(residual_without_force)):
        raise FloatingPointError("vorticity-equation terms became nonfinite")

    return {
        "velocity": velocity,
        "omega": omega,
        "omega_t": omega_t,
        "advection": advection,
        "stretching": stretching,
        "diffusion": diffusion,
        "residual_without_force": residual_without_force,
    }


def vorticity_equation_residual(
    field: VelocityLike,
    points,
    time,
    *,
    spatial_step: float,
    time_step: float,
    nu: float,
    force: Callable | None = None,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Return ``omega_t+u.grad(omega)-omega.grad(u)-nu*Delta(omega)-curl(f)``."""
    p = _points(points)
    t = _times(time, len(p))
    terms = vorticity_equation_terms(
        field,
        p,
        t,
        spatial_step=spatial_step,
        time_step=time_step,
        nu=nu,
    )
    if force is None:
        force_curl = np.zeros_like(p)
    else:
        force_curl = _curl(force, p, t, spatial_step)
    residual = terms["residual_without_force"] - force_curl
    if not np.all(np.isfinite(residual)):
        raise FloatingPointError("vorticity-equation residual became nonfinite")
    terms = dict(terms)
    terms["force_curl"] = force_curl
    return residual, terms


def _metrics(residual: np.ndarray, terms: dict[str, np.ndarray]) -> tuple[float, float, float]:
    norms = np.linalg.norm(residual, axis=1)
    rms = float(np.sqrt(np.mean(norms**2)))
    scale_squared = np.zeros(len(residual), dtype=float)
    for key in ("omega_t", "advection", "stretching", "diffusion", "force_curl"):
        value = np.asarray(terms[key], dtype=float)
        scale_squared += np.sum(value * value, axis=1)
    scale = float(np.sqrt(np.mean(scale_squared)))
    return float(np.max(norms)), rms, rms / max(scale, np.finfo(float).tiny)


def audit_vorticity_equation(
    field: VelocityLike,
    points,
    time,
    *,
    spatial_steps=(0.02, 0.01, 0.005),
    time_step: float = 0.0025,
    nu: float = 0.01,
    restricted_force: RestrictedForce | None = None,
) -> list[VorticityResidualRow]:
    """Compare zero forcing with one frozen restricted-force choice at >=3 spatial levels."""
    p = _points(points)
    t = _times(time, len(p))
    steps = tuple(float(step) for step in spatial_steps)
    if len(steps) < 3 or any(not np.isfinite(step) or step <= 0 for step in steps):
        raise ValueError("spatial_steps must contain at least three positive finite levels")
    if any(right >= left for left, right in zip(steps, steps[1:])):
        raise ValueError("spatial_steps must be strictly decreasing")
    if restricted_force is None:
        restricted_force = RestrictedForce(a=0.0, c=2.19114231)

    rows: list[VorticityResidualRow] = []
    for step in steps:
        zero_residual, zero_terms = vorticity_equation_residual(
            field,
            p,
            t,
            spatial_step=step,
            time_step=time_step,
            nu=nu,
            force=None,
        )
        zero_max, zero_rms, zero_relative = _metrics(zero_residual, zero_terms)

        force_curl = _curl(restricted_force, p, t, step)
        restricted_residual = zero_terms["residual_without_force"] - force_curl
        restricted_terms = dict(zero_terms)
        restricted_terms["force_curl"] = force_curl
        force_max, force_rms, force_relative = _metrics(
            restricted_residual, restricted_terms
        )
        reduction = (zero_rms - force_rms) / max(zero_rms, np.finfo(float).tiny)
        rows.append(
            VorticityResidualRow(
                spatial_step=step,
                zero_force_max=zero_max,
                zero_force_rms=zero_rms,
                zero_force_relative_rms=zero_relative,
                restricted_force_max=force_max,
                restricted_force_rms=force_rms,
                restricted_force_relative_rms=force_relative,
                force_rms_reduction_fraction=float(reduction),
            )
        )
    return rows


def audit_frozen_eq45(
    *,
    candidate_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
    sample_count: int = 12,
    seed: int = 914131,
) -> dict:
    """Run the independent audit on the checked frozen Eq45 artifact."""
    if not isinstance(sample_count, int) or sample_count < 4:
        raise ValueError("sample_count must be an integer >=4")
    root = Path(__file__).resolve().parents[2]
    candidate_path = Path(candidate_path or root / "artifacts/constrained/eq45_velocity_candidate_seed.json")
    constraints_path = Path(constraints_path or root / "configs/constraints.json")

    candidate = Eq45VelocityCandidate.load_json(candidate_path)
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    nu = float(constraints["nu"])
    registered_steps = tuple(float(x) for x in constraints["validation"]["derivative_steps"])

    rng = np.random.default_rng(seed)
    points = rng.uniform(-0.22, 0.22, size=(sample_count, 3))
    times = rng.uniform(0.40, 0.60, size=sample_count)
    frozen_force = RestrictedForce(a=0.0, c=2.19114231)
    rows = audit_vorticity_equation(
        candidate,
        points,
        times,
        spatial_steps=registered_steps,
        time_step=0.0025,
        nu=nu,
        restricted_force=frozen_force,
    )
    finest = rows[-1]
    return {
        "schema": "eq45_vorticity_equation_crosscheck_v1",
        "candidate_sha256": candidate.sha256,
        "validation_seed": seed,
        "sample_count": sample_count,
        "sample_box": [[-0.22, 0.22], [-0.22, 0.22], [-0.22, 0.22]],
        "sample_time_range": [0.40, 0.60],
        "nu": nu,
        "spatial_steps": list(registered_steps),
        "fixed_time_step": 0.0025,
        "restricted_force_comparison": {
            "a": frozen_force.a,
            "c": frozen_force.c,
            "classification": "frozen comparison copied from upstream diagnostic PR #82; not refit here",
        },
        "rows": [row.as_dict() for row in rows],
        "finest_zero_force_rms": finest.zero_force_rms,
        "finest_restricted_force_rms": finest.restricted_force_rms,
        "finest_force_rms_reduction_fraction": finest.force_rms_reduction_fraction,
        "claim_scope": "pressure-free vorticity-equation necessary-condition cross-check",
        "threshold_interpretation": "no new vorticity-residual acceptance threshold is introduced",
        "velocity_changed": False,
        "pressure_fitted": False,
        "forcing_fitted_in_this_audit": False,
        "physical_support_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    print(json.dumps(audit_frozen_eq45(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
