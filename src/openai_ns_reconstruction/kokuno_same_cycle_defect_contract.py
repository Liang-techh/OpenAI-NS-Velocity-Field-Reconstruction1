"""Typed same-cycle Navier--Stokes momentum-defect contract for Kokuno Agent 3.

This module is deliberately narrower than a correction cycle. It turns raw,
same-cycle field providers into the actual momentum defect

    R = u_t + (u . grad)u + grad(p) - nu Delta(u) - f

with fixed centered fourth-order spatial differences. Callers cannot inject a
precomputed residual, stress, correction target, pressure gradient, or gain.

The contract is an engineering interface, not independent PDE validation.
Provider provenance and restricted-forcing semantics still require upstream
admission before the resulting defect may drive the full Agent-3 correction
machinery.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Sequence

import numpy as np

VectorEvaluator = Callable[[float, float, float, float], Sequence[float]]
ScalarEvaluator = Callable[[float, float, float, float], float]


@dataclass(frozen=True)
class CycleIdentity:
    """Identity of one assembled candidate/correction-cycle state."""

    cycle_id: str
    cycle_index: int
    state_token: str

    def __post_init__(self) -> None:
        if not self.cycle_id.strip():
            raise ValueError("cycle_id must be non-empty")
        if self.cycle_index < 0:
            raise ValueError("cycle_index must be non-negative")
        if not self.state_token.strip():
            raise ValueError("state_token must be non-empty")


@dataclass(frozen=True)
class VectorFieldProvider:
    identity: CycleIdentity
    source_ref: str
    evaluator: VectorEvaluator

    def __post_init__(self) -> None:
        if not self.source_ref.strip():
            raise ValueError("source_ref must be non-empty")


@dataclass(frozen=True)
class ScalarFieldProvider:
    identity: CycleIdentity
    source_ref: str
    evaluator: ScalarEvaluator

    def __post_init__(self) -> None:
        if not self.source_ref.strip():
            raise ValueError("source_ref must be non-empty")


@dataclass(frozen=True)
class RestrictedForcingProvider:
    """Forcing provider whose restriction claim is externally traceable.

    Constructing this object is not independent validation of the restriction
    receipt. It prevents the defect API from accepting an anonymous free
    forcing or a caller-supplied residual in place of forcing.
    """

    identity: CycleIdentity
    source_ref: str
    restriction_receipt: str
    evaluator: VectorEvaluator

    def __post_init__(self) -> None:
        if not self.source_ref.strip():
            raise ValueError("source_ref must be non-empty")
        if not self.restriction_receipt.strip():
            raise ValueError("restriction_receipt must be non-empty")


@dataclass(frozen=True)
class DefectEvaluation:
    identity: CycleIdentity
    point_count: int
    viscosity: float
    spatial_step: float
    residual_values: np.ndarray
    vector_rms: float
    vector_max: float
    component_rms: np.ndarray

    def to_receipt(self) -> dict[str, object]:
        return {
            "identity": asdict(self.identity),
            "point_count": self.point_count,
            "viscosity": self.viscosity,
            "spatial_step": self.spatial_step,
            "vector_rms": self.vector_rms,
            "vector_max": self.vector_max,
            "component_rms": self.component_rms.tolist(),
            "residual_values": self.residual_values.tolist(),
        }


def _as_points(
    points: Sequence[Sequence[float]] | np.ndarray, *, label: str
) -> np.ndarray:
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 4 or arr.shape[0] == 0:
        raise ValueError(f"{label} must have shape (N, 4) with N >= 1")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{label} must contain only finite values")
    rows = [tuple(float(v) for v in row) for row in arr]
    if len(set(rows)) != len(rows):
        raise ValueError(f"{label} contains duplicate points")
    return arr


def _vector(
    provider: VectorFieldProvider | RestrictedForcingProvider, point: np.ndarray
) -> np.ndarray:
    value = np.asarray(provider.evaluator(*map(float, point)), dtype=float)
    if value.shape != (3,) or not np.all(np.isfinite(value)):
        raise ValueError(f"{provider.source_ref} must return a finite 3-vector")
    return value


def _scalar(provider: ScalarFieldProvider, point: np.ndarray) -> float:
    value = float(provider.evaluator(*map(float, point)))
    if not np.isfinite(value):
        raise ValueError(f"{provider.source_ref} must return a finite scalar")
    return value


def _shift(point: np.ndarray, axis: int, amount: float) -> np.ndarray:
    shifted = point.copy()
    shifted[axis] += amount
    return shifted


def _fd4_first_vector(
    provider: VectorFieldProvider, point: np.ndarray, axis: int, h: float
) -> np.ndarray:
    return (
        -_vector(provider, _shift(point, axis, 2.0 * h))
        + 8.0 * _vector(provider, _shift(point, axis, h))
        - 8.0 * _vector(provider, _shift(point, axis, -h))
        + _vector(provider, _shift(point, axis, -2.0 * h))
    ) / (12.0 * h)


def _fd4_second_vector(
    provider: VectorFieldProvider, point: np.ndarray, axis: int, h: float
) -> np.ndarray:
    return (
        -_vector(provider, _shift(point, axis, 2.0 * h))
        + 16.0 * _vector(provider, _shift(point, axis, h))
        - 30.0 * _vector(provider, point)
        + 16.0 * _vector(provider, _shift(point, axis, -h))
        - _vector(provider, _shift(point, axis, -2.0 * h))
    ) / (12.0 * h * h)


def _fd4_first_scalar(
    provider: ScalarFieldProvider, point: np.ndarray, axis: int, h: float
) -> float:
    return (
        -_scalar(provider, _shift(point, axis, 2.0 * h))
        + 8.0 * _scalar(provider, _shift(point, axis, h))
        - 8.0 * _scalar(provider, _shift(point, axis, -h))
        + _scalar(provider, _shift(point, axis, -2.0 * h))
    ) / (12.0 * h)


def _assert_same_cycle(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    forcing: RestrictedForcingProvider,
) -> CycleIdentity:
    identity = velocity.identity
    named = {
        "velocity_dt": velocity_dt.identity,
        "pressure": pressure.identity,
        "restricted_forcing": forcing.identity,
    }
    for name, other in named.items():
        if other != identity:
            raise ValueError(f"{name} does not share the velocity same-cycle identity")
    return identity


def evaluate_same_cycle_momentum_defect(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    points: Sequence[Sequence[float]] | np.ndarray,
    *,
    viscosity: float,
    spatial_step: float,
) -> DefectEvaluation:
    """Evaluate the raw same-cycle momentum defect from field providers only."""

    identity = _assert_same_cycle(
        velocity, velocity_dt, pressure, restricted_forcing
    )
    if not np.isfinite(viscosity) or viscosity < 0.0:
        raise ValueError("viscosity must be finite and non-negative")
    if not np.isfinite(spatial_step) or spatial_step <= 0.0:
        raise ValueError("spatial_step must be finite and positive")
    pts = _as_points(points, label="points")
    residuals: list[np.ndarray] = []
    for point in pts:
        u = _vector(velocity, point)
        u_t = _vector(velocity_dt, point)
        grad_u = np.column_stack(
            [
                _fd4_first_vector(velocity, point, axis, spatial_step)
                for axis in range(3)
            ]
        )
        lap_u = sum(
            (
                _fd4_second_vector(velocity, point, axis, spatial_step)
                for axis in range(3)
            ),
            start=np.zeros(3, dtype=float),
        )
        grad_p = np.asarray(
            [
                _fd4_first_scalar(pressure, point, axis, spatial_step)
                for axis in range(3)
            ],
            dtype=float,
        )
        forcing = _vector(restricted_forcing, point)
        convection = grad_u @ u
        residuals.append(
            u_t + convection + grad_p - viscosity * lap_u - forcing
        )
    values = np.vstack(residuals)
    point_norms = np.linalg.norm(values, axis=1)
    return DefectEvaluation(
        identity=identity,
        point_count=len(pts),
        viscosity=float(viscosity),
        spatial_step=float(spatial_step),
        residual_values=values,
        vector_rms=float(np.sqrt(np.mean(point_norms * point_norms))),
        vector_max=float(np.max(point_norms)),
        component_rms=np.sqrt(np.mean(values * values, axis=0)),
    )


def evaluate_disjoint_heldin_heldout(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    held_in_points: Sequence[Sequence[float]] | np.ndarray,
    held_out_points: Sequence[Sequence[float]] | np.ndarray,
    *,
    viscosity: float,
    spatial_step: float,
) -> dict[str, DefectEvaluation]:
    """Evaluate disjoint held-in/out sets without reusing a validation point."""

    held_in = _as_points(held_in_points, label="held_in_points")
    held_out = _as_points(held_out_points, label="held_out_points")
    held_in_rows = {tuple(float(v) for v in row) for row in held_in}
    held_out_rows = {tuple(float(v) for v in row) for row in held_out}
    if held_in_rows & held_out_rows:
        raise ValueError("held-in and held-out point sets must be disjoint")
    kwargs = {"viscosity": viscosity, "spatial_step": spatial_step}
    return {
        "held_in": evaluate_same_cycle_momentum_defect(
            velocity,
            velocity_dt,
            pressure,
            restricted_forcing,
            held_in,
            **kwargs,
        ),
        "held_out": evaluate_same_cycle_momentum_defect(
            velocity,
            velocity_dt,
            pressure,
            restricted_forcing,
            held_out,
            **kwargs,
        ),
    }


def truth_boundary() -> dict[str, object]:
    return {
        "momentum_defect_formula": "u_t + (u dot grad)u + grad(p) - nu*Delta(u) - f",
        "spatial_derivative_operator": "centered FD4",
        "actual_momentum_defect_formula_executable": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_target_allowed": False,
        "caller_supplied_pressure_gradient_allowed": False,
        "caller_supplied_gain_allowed": False,
        "same_cycle_identity_enforced": True,
        "restricted_forcing_provider_required": True,
        "heldin_heldout_disjointness_enforced": True,
        "restricted_forcing_semantics_independently_validated_here": False,
        "full_same_cycle_composite_requested_stress_materialized": False,
        "candidate_finite_head_mean_debt_materialized": False,
        "real_full_candidate_defect_consumed": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_momentum_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": 1.0e-3,
        "final_normalized_divergence_gate": 1.0e-5,
    }


def deterministic_receipt() -> dict[str, object]:
    identity = CycleIdentity("analytic-contract-regression-v1", 0, "fixed")
    velocity = VectorFieldProvider(
        identity, "analytic:u=(y,0,0)", lambda x, y, z, t: (y, 0.0, 0.0)
    )
    velocity_dt = VectorFieldProvider(
        identity, "analytic:u_t=0", lambda x, y, z, t: (0.0, 0.0, 0.0)
    )
    pressure = ScalarFieldProvider(
        identity, "analytic:p=2x+3z", lambda x, y, z, t: 2.0 * x + 3.0 * z
    )
    forcing = RestrictedForcingProvider(
        identity,
        "analytic:f=0",
        "regression-zero-forcing; independent of residual",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    split = evaluate_disjoint_heldin_heldout(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        held_in_points=[
            (-0.31, 0.27, 0.14, 0.37),
            (0.22, -0.19, -0.33, 0.63),
        ],
        held_out_points=[
            (0.11, 0.36, -0.21, 0.41),
            (-0.28, -0.17, 0.29, 0.59),
        ],
        viscosity=0.01,
        spatial_step=1.0e-3,
    )
    expected = np.asarray([2.0, 0.0, 3.0])
    max_error = max(
        float(np.max(np.abs(item.residual_values - expected)))
        for item in split.values()
    )
    return {
        "schema": "kokuno-a3-same-cycle-defect-contract-v1",
        "task": "KOKUNO-A3-SAME-CYCLE-DEFECT-CONTRACT-046",
        "provenance": {
            "parent_agent3_pr": 607,
            "parent_agent3_head": "a2fc36a44771cffb173e391d53fb04df220ab638",
            "purpose": "typed real-defect seam before full Agent-1 composite handoff exists",
        },
        "formula_regression": {
            "expected_residual": expected.tolist(),
            "expected_vector_rms": float(np.sqrt(13.0)),
            "maximum_absolute_error": max_error,
            "held_in": split["held_in"].to_receipt(),
            "held_out": split["held_out"].to_receipt(),
        },
        "truth_boundary": truth_boundary(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = deterministic_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
