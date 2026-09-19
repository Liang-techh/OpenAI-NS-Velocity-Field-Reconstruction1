"""Typed angular mean projection of the actual same-cycle momentum defect.

Kokuno Agent 3 needs the axisymmetric mean channels of the *real* candidate
defect before compact radial-stress reconstruction can act.  This module is a
narrow bridge from :mod:`kokuno_same_cycle_defect_contract`: callers provide
raw same-cycle ``u, u_t, p, f`` field providers, never a precomputed residual.
For each cylindrical ring we evaluate

    R = u_t + (u . grad)u + grad(p) - nu Delta(u) - f

at an equally spaced angular grid and form

    R_r     = e_r     . R,
    R_theta = e_theta . R,
    R_z     = e_z     . R,

then return their discrete angular means.  The theta/axial channels are the
typed input needed by the existing Agent-3 compact/radial stress machinery.

The equal-angle quadrature and this adapter are repository engineering
machinery.  They do not prove the Kokuno reconstruction, validate restricted
forcing semantics, or assess the final normalized Navier--Stokes residual.
"""

from __future__ import annotations

import argparse
import inspect
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
    evaluate_same_cycle_momentum_defect,
)


@dataclass(frozen=True)
class CylindricalRing:
    """One off-axis ring on which an angular mean is evaluated."""

    time: float
    radius: float
    axial_z: float
    angular_count: int = 32
    phase: float = 0.0

    def __post_init__(self) -> None:
        values = (self.time, self.radius, self.axial_z, self.phase)
        if not all(np.isfinite(value) for value in values):
            raise ValueError("ring parameters must be finite")
        if self.radius <= 0.0:
            raise ValueError("ring radius must be strictly positive")
        if self.angular_count < 8:
            raise ValueError("angular_count must be at least 8")


@dataclass(frozen=True)
class MeanDefectRingEvaluation:
    """Actual-defect angular means for one cylindrical ring."""

    identity: CycleIdentity
    ring: CylindricalRing
    viscosity: float
    spatial_step: float
    radial_mean: float
    theta_mean: float
    axial_mean: float
    radial_rms: float
    theta_rms: float
    axial_rms: float
    raw_vector_rms: float
    raw_vector_max: float

    def to_receipt(self) -> dict[str, object]:
        return {
            "identity": asdict(self.identity),
            "ring": asdict(self.ring),
            "viscosity": self.viscosity,
            "spatial_step": self.spatial_step,
            "radial_mean": self.radial_mean,
            "theta_mean": self.theta_mean,
            "axial_mean": self.axial_mean,
            "radial_rms": self.radial_rms,
            "theta_rms": self.theta_rms,
            "axial_rms": self.axial_rms,
            "raw_vector_rms": self.raw_vector_rms,
            "raw_vector_max": self.raw_vector_max,
        }


@dataclass(frozen=True)
class MeanDefectSplitEvaluation:
    """Disjoint held-in/out mean-defect projection for one cycle state."""

    identity: CycleIdentity
    held_in: tuple[MeanDefectRingEvaluation, ...]
    held_out: tuple[MeanDefectRingEvaluation, ...]

    def to_receipt(self) -> dict[str, object]:
        return {
            "identity": asdict(self.identity),
            "held_in": [item.to_receipt() for item in self.held_in],
            "held_out": [item.to_receipt() for item in self.held_out],
        }


def _angles(ring: CylindricalRing) -> np.ndarray:
    return ring.phase + 2.0 * np.pi * np.arange(ring.angular_count) / ring.angular_count


def _ring_points(ring: CylindricalRing) -> np.ndarray:
    angles = _angles(ring)
    x = ring.radius * np.cos(angles)
    y = ring.radius * np.sin(angles)
    z = np.full_like(angles, ring.axial_z, dtype=float)
    t = np.full_like(angles, ring.time, dtype=float)
    return np.column_stack((x, y, z, t))


def _point_keys(points: np.ndarray) -> set[tuple[float, float, float, float]]:
    # Generated held-in/out grids are deterministic.  Hex keys avoid accidental
    # decimal-string rounding while making exact reuse fail closed.
    return {
        tuple(float(value).hex() for value in row)  # type: ignore[arg-type]
        for row in points
    }


def _validate_ring_sequence(
    rings: Sequence[CylindricalRing], *, label: str
) -> tuple[CylindricalRing, ...]:
    values = tuple(rings)
    if not values:
        raise ValueError(f"{label} must contain at least one ring")
    point_keys: set[tuple[float, float, float, float]] = set()
    for ring in values:
        current = _point_keys(_ring_points(ring))
        if point_keys & current:
            raise ValueError(f"{label} reuses a physical quadrature point")
        point_keys.update(current)
    return values


def evaluate_actual_mean_defect_ring(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    ring: CylindricalRing,
    *,
    viscosity: float,
    spatial_step: float,
) -> MeanDefectRingEvaluation:
    """Project a raw-provider momentum defect to cylindrical angular means."""

    points = _ring_points(ring)
    defect = evaluate_same_cycle_momentum_defect(
        velocity,
        velocity_dt,
        pressure,
        restricted_forcing,
        points,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    angles = _angles(ring)
    cos_theta = np.cos(angles)
    sin_theta = np.sin(angles)
    values = defect.residual_values
    radial = cos_theta * values[:, 0] + sin_theta * values[:, 1]
    theta = -sin_theta * values[:, 0] + cos_theta * values[:, 1]
    axial = values[:, 2]

    return MeanDefectRingEvaluation(
        identity=defect.identity,
        ring=ring,
        viscosity=defect.viscosity,
        spatial_step=defect.spatial_step,
        radial_mean=float(np.mean(radial)),
        theta_mean=float(np.mean(theta)),
        axial_mean=float(np.mean(axial)),
        radial_rms=float(np.sqrt(np.mean(radial * radial))),
        theta_rms=float(np.sqrt(np.mean(theta * theta))),
        axial_rms=float(np.sqrt(np.mean(axial * axial))),
        raw_vector_rms=defect.vector_rms,
        raw_vector_max=defect.vector_max,
    )


def evaluate_disjoint_actual_mean_defect(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    held_in_rings: Sequence[CylindricalRing],
    held_out_rings: Sequence[CylindricalRing],
    *,
    viscosity: float,
    spatial_step: float,
) -> MeanDefectSplitEvaluation:
    """Evaluate typed mean channels on disjoint held-in/out cylindrical rings."""

    held_in = _validate_ring_sequence(held_in_rings, label="held_in_rings")
    held_out = _validate_ring_sequence(held_out_rings, label="held_out_rings")
    held_in_keys = set().union(*(_point_keys(_ring_points(ring)) for ring in held_in))
    held_out_keys = set().union(*(_point_keys(_ring_points(ring)) for ring in held_out))
    if held_in_keys & held_out_keys:
        raise ValueError("held-in and held-out ring quadrature points must be disjoint")

    kwargs = {"viscosity": viscosity, "spatial_step": spatial_step}
    projected_in = tuple(
        evaluate_actual_mean_defect_ring(
            velocity, velocity_dt, pressure, restricted_forcing, ring, **kwargs
        )
        for ring in held_in
    )
    projected_out = tuple(
        evaluate_actual_mean_defect_ring(
            velocity, velocity_dt, pressure, restricted_forcing, ring, **kwargs
        )
        for ring in held_out
    )
    identities = {item.identity for item in projected_in + projected_out}
    if len(identities) != 1:
        raise RuntimeError("mean-defect projection produced mixed cycle identities")
    return MeanDefectSplitEvaluation(
        identity=next(iter(identities)), held_in=projected_in, held_out=projected_out
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(evaluate_disjoint_actual_mean_defect)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "target",
        "gain",
        "normalized_score",
    }
    return {
        "source_defect_formula": "u_t + (u dot grad)u + grad(p) - nu*Delta(u) - f",
        "source_defect_operator": "kokuno_same_cycle_defect_contract centered FD4",
        "angular_projection": "equal-angle discrete mean of cylindrical components",
        "theta_projection_formula": "mean((-sin(theta), cos(theta), 0) dot R)",
        "axial_projection_formula": "mean((0, 0, 1) dot R)",
        "raw_same_cycle_providers_required": True,
        "caller_supplied_precomputed_defect_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "heldin_heldout_ring_points_disjoint": True,
        "mean_channels_materialized_from_actual_defect": True,
        "ready_as_typed_input_for_radial_stress_reconstruction": True,
        "restricted_forcing_semantics_independently_validated_here": False,
        "kokuno_theorem_machine_replayed_here": False,
        "full_same_cycle_composite_requested_stress_materialized": False,
        "candidate_finite_head_mean_debt_materialized": False,
        "real_full_candidate_defect_consumed": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
    }


def deterministic_receipt() -> dict[str, object]:
    """Nonzero analytic regression; forcing is zero and independent of R."""

    identity = CycleIdentity("analytic-rotating-mean-defect-v1", 0, "t0")
    velocity = VectorFieldProvider(
        identity,
        "analytic:u=t*(-y,x,0)",
        lambda x, y, z, t: (-t * y, t * x, 0.0),
    )
    velocity_dt = VectorFieldProvider(
        identity,
        "analytic:u_t=(-y,x,0)",
        lambda x, y, z, t: (-y, x, 0.0),
    )
    pressure = ScalarFieldProvider(
        identity, "analytic:p=0", lambda x, y, z, t: 0.0
    )
    forcing = RestrictedForcingProvider(
        identity,
        "analytic:f=0",
        "regression-zero-forcing; independent of residual",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    held_in_rings = (
        CylindricalRing(time=0.0, radius=0.20, axial_z=-0.17, angular_count=32),
        CylindricalRing(time=0.0, radius=0.34, axial_z=0.09, angular_count=32),
    )
    held_out_rings = (
        CylindricalRing(time=0.0, radius=0.27, axial_z=0.23, angular_count=32, phase=0.031),
        CylindricalRing(time=0.0, radius=0.41, axial_z=-0.29, angular_count=32, phase=0.047),
    )
    split = evaluate_disjoint_actual_mean_defect(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        held_in_rings,
        held_out_rings,
        viscosity=0.01,
        spatial_step=0.005,
    )
    all_items = split.held_in + split.held_out
    theta_errors = [abs(item.theta_mean - item.ring.radius) for item in all_items]
    radial_errors = [abs(item.radial_mean) for item in all_items]
    axial_errors = [abs(item.axial_mean) for item in all_items]
    return {
        "schema": "kokuno-a3-typed-mean-defect-projection-v1",
        "task": "KOKUNO-A3-TYPED-MEAN-DEFECT-PROJECTION-051",
        "provenance": {
            "parent_agent3_pr": 654,
            "parent_agent3_head": "b0c94639b3d4e83ff9c786e0c4f7b68a6bd6dae7",
            "structural_source": "Kokuno corrected reconstruction mean/radial correction architecture",
            "repository_role": "actual same-cycle defect -> angular mean channels before compact radial stress",
        },
        "analytic_regression": {
            "field": "u=t*(-y,x,0), u_t=(-y,x,0), p=0, f=0 evaluated at t=0",
            "expected_actual_defect": "R=(-y,x,0)",
            "expected_ring_channels": "mean(R_theta)=r, mean(R_r)=mean(R_z)=0",
            "maximum_theta_mean_error": float(max(theta_errors)),
            "maximum_radial_mean_error": float(max(radial_errors)),
            "maximum_axial_mean_error": float(max(axial_errors)),
            "split": split.to_receipt(),
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
