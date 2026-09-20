"""Strict-inner PA.10 plus oscillatory local transport seam for Kokuno Agent 2.

This module performs one bounded composition step for the currently executable
Agent-1 PA.10 *inner contraction-center* field and the frozen Agent-2 localized
complete-curl oscillatory field.  It exposes

    T_inner+osc = d_t u_inner+osc
                  + (u_inner+osc . grad) u_inner+osc
                  - nu Delta u_inner+osc,

with the repository physical viscosity fixed at ``nu=0.01`` and

    Delta u_inner+osc = Delta u_inner + Delta u_osc.

Agent 1 owns the strict-inner leading velocity/time derivative/self-advection/
Laplacian.  Agent 2 owns the oscillatory velocity and its derivative/nonlinear/
Laplacian seams.  This module only composes those typed quantities.  It does not
invent an outer/global leading continuation, pressure, restricted forcing,
Agent-3 mean/radial correction, or a complete Navier--Stokes residual.

Source/realization boundary
---------------------------
The corrected 2026-09-09 Kokuno reconstruction is structural provenance for the
upstream PA.10 inner formulas and localized complete-curl oscillatory structure.
The public-z pullback, repository-autonomous oscillatory choices, finite-
difference Laplacians/verifiers, and this transport composition are repository
realizations.  Nothing here is claimed paper-exact or globally source-complete.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
import math
from typing import Any, Protocol, runtime_checkable

import numpy as np

from .kokuno_public_inner_leading_oscillatory_full_advection import (
    evaluate_inner_leading_oscillatory_full_advection,
)
from .kokuno_public_inner_leading_oscillatory_time_derivative import (
    evaluate_inner_leading_oscillatory_time_derivative,
)
from .kokuno_public_oscillatory_spatial_laplacian import (
    evaluate_velocity_laplacian_fd6,
)
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-INNER-LEADING-OSCILLATORY-TRANSPORT-062"
SCHEMA = "kokuno-a2-inner-leading-oscillatory-transport-v1"
PARENT_AGENT2_PR = 830
PARENT_AGENT2_HEAD = "3479bed94ac6ef5c67547227f81a224ad43fb409"
AGENT1_LAPLACIAN_PR = 839
AGENT1_LAPLACIAN_HEAD = "e96ee90144976a992b62d76d4361a94eb16bd91e"
AGENT1_LAPLACIAN_BLOB = "74b0e185e8e0bbb08695d493e0a091c305a03006"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
VISCOSITY = 0.01
INNER_ADVECTION_SPATIAL_STEP = 1.0e-3
OSCILLATORY_ADVECTION_SPATIAL_STEP = 1.0e-3
OSCILLATORY_LAPLACIAN_SPATIAL_STEP = 1.0e-3
INDEPENDENT_FD4_STEPS = (4.0e-3, 2.0e-3, 1.0e-3)


@runtime_checkable
class InnerLeadingTransportBackend(Protocol):
    """Minimal strict-inner Agent-1 protocol; global leading semantics are absent."""

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...

    def self_advection(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...

    def velocity_laplacian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...


@dataclass(frozen=True)
class InnerLeadingOscillatoryTransportResult:
    inner_leading_velocity: np.ndarray
    oscillatory_velocity: np.ndarray
    inner_plus_oscillatory_velocity: np.ndarray
    inner_leading_velocity_dt: np.ndarray
    oscillatory_velocity_dt: np.ndarray
    inner_plus_oscillatory_velocity_dt: np.ndarray
    inner_plus_oscillatory_self_advection: np.ndarray
    inner_leading_laplacian: np.ndarray
    oscillatory_laplacian: np.ndarray
    inner_plus_oscillatory_laplacian: np.ndarray
    viscous_term: np.ndarray
    transport: np.ndarray
    viscosity: float


def _broadcast_xyzt(x: Any, y: Any, z: Any, t: Any) -> tuple[np.ndarray, ...]:
    arrays = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    if not all(np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("x, y, z and t must be finite and broadcastable")
    return tuple(arrays)


def _backend_vector(
    backend: InnerLeadingTransportBackend,
    method_name: str,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    method = getattr(backend, method_name, None)
    if not callable(method):
        raise TypeError(
            f"inner_leading_backend must expose callable {method_name}(x,y,z,t)"
        )
    value = np.asarray(method(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise ValueError(
            f"inner leading {method_name} must have shape {expected}, got {value.shape}"
        )
    if not np.all(np.isfinite(value)):
        raise ValueError(f"inner leading {method_name} must be finite")
    return value


def _require_exact_replay(left: np.ndarray, right: np.ndarray, label: str) -> None:
    if left.shape != right.shape or not np.array_equal(left, right):
        raise RuntimeError(f"{label} did not replay exactly across inherited seams")


def evaluate_inner_leading_oscillatory_transport(
    inner_leading_backend: InnerLeadingTransportBackend,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
) -> InnerLeadingOscillatoryTransportResult:
    """Compose the strict-inner local transport operator with fixed viscosity.

    The public interface deliberately exposes no viscosity or derivative-step
    knob.  Any Agent-1 stencil/domain failure propagates fail-closed.
    """
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)

    time_part = evaluate_inner_leading_oscillatory_time_derivative(
        inner_leading_backend, xb, yb, zb, tb
    )
    advection_part = evaluate_inner_leading_oscillatory_full_advection(
        inner_leading_backend,
        xb,
        yb,
        zb,
        tb,
        inner_leading_spatial_step=INNER_ADVECTION_SPATIAL_STEP,
        oscillatory_spatial_step=OSCILLATORY_ADVECTION_SPATIAL_STEP,
    )
    inner_laplacian = _backend_vector(
        inner_leading_backend, "velocity_laplacian", xb, yb, zb, tb
    )
    oscillatory_laplacian_result = evaluate_velocity_laplacian_fd6(
        xb,
        yb,
        zb,
        tb,
        spatial_step=OSCILLATORY_LAPLACIAN_SPATIAL_STEP,
    )
    oscillatory_laplacian = np.asarray(
        oscillatory_laplacian_result["laplacian"], dtype=float
    )
    oscillatory_velocity = np.asarray(
        oscillatory_laplacian_result["velocity"], dtype=float
    )

    _require_exact_replay(
        time_part.inner_leading_velocity,
        advection_part.inner_leading_velocity,
        "inner velocity time/advection",
    )
    _require_exact_replay(
        time_part.oscillatory_velocity,
        advection_part.oscillatory_velocity,
        "oscillatory velocity time/advection",
    )
    _require_exact_replay(
        time_part.oscillatory_velocity,
        oscillatory_velocity,
        "oscillatory velocity time/laplacian",
    )
    _require_exact_replay(
        time_part.inner_plus_oscillatory_velocity,
        advection_part.inner_plus_oscillatory_velocity,
        "composite velocity time/advection",
    )

    total_laplacian = inner_laplacian + oscillatory_laplacian
    viscous_term = -VISCOSITY * total_laplacian
    transport = (
        time_part.inner_plus_oscillatory_velocity_dt
        + advection_part.inner_plus_oscillatory_self_advection
        + viscous_term
    )
    expected = xb.shape + (3,)
    for name, value in (
        ("inner_laplacian", inner_laplacian),
        ("oscillatory_laplacian", oscillatory_laplacian),
        ("total_laplacian", total_laplacian),
        ("viscous_term", viscous_term),
        ("transport", transport),
    ):
        if value.shape != expected or not np.all(np.isfinite(value)):
            raise RuntimeError(f"{name} is non-finite or has the wrong shape")

    return InnerLeadingOscillatoryTransportResult(
        inner_leading_velocity=time_part.inner_leading_velocity,
        oscillatory_velocity=time_part.oscillatory_velocity,
        inner_plus_oscillatory_velocity=time_part.inner_plus_oscillatory_velocity,
        inner_leading_velocity_dt=time_part.inner_leading_velocity_dt,
        oscillatory_velocity_dt=time_part.oscillatory_velocity_dt,
        inner_plus_oscillatory_velocity_dt=time_part.inner_plus_oscillatory_velocity_dt,
        inner_plus_oscillatory_self_advection=
            advection_part.inner_plus_oscillatory_self_advection,
        inner_leading_laplacian=inner_laplacian,
        oscillatory_laplacian=oscillatory_laplacian,
        inner_plus_oscillatory_laplacian=total_laplacian,
        viscous_term=viscous_term,
        transport=transport,
        viscosity=VISCOSITY,
    )


def _validate_step(step: float) -> float:
    h = float(step)
    if not math.isfinite(h) or not 1.0e-5 <= h <= 5.0e-2:
        raise ValueError("step must be finite and lie in [1e-5,5e-2]")
    return h


def _total_public_velocity(
    backend: InnerLeadingTransportBackend,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    inner = _backend_vector(backend, "velocity", x, y, z, t)
    osc = np.asarray(velocity_osc(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if osc.shape != expected or not np.all(np.isfinite(osc)):
        raise RuntimeError("oscillatory velocity is non-finite or has the wrong shape")
    return inner + osc


def independent_fd4_inner_plus_oscillatory_transport(
    inner_leading_backend: InnerLeadingTransportBackend,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    """Implementation-distinct public-velocity-only FD4 transport verifier.

    This path never calls production ``velocity_dt``, ``self_advection`` or
    ``velocity_laplacian`` from Agent 1 and never calls A2 production time,
    advection or Laplacian evaluators.  It samples only the direct public sum.
    """
    h = _validate_step(step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)

    def value(xx: np.ndarray, yy: np.ndarray, zz: np.ndarray, tt: np.ndarray) -> np.ndarray:
        return _total_public_velocity(inner_leading_backend, xx, yy, zz, tt)

    dt = (
        value(xb, yb, zb, tb - 2.0 * h)
        - 8.0 * value(xb, yb, zb, tb - h)
        + 8.0 * value(xb, yb, zb, tb + h)
        - value(xb, yb, zb, tb + 2.0 * h)
    ) / (12.0 * h)

    u0 = value(xb, yb, zb, tb)
    speed = np.linalg.norm(u0, axis=-1)
    direction = np.zeros_like(u0)
    nonzero = speed > 1.0e-14
    direction[nonzero] = u0[nonzero] / speed[nonzero, None]
    line: dict[int, np.ndarray] = {}
    for offset in (-2, -1, 1, 2):
        displacement = offset * h * direction
        line[offset] = value(
            xb + displacement[..., 0],
            yb + displacement[..., 1],
            zb + displacement[..., 2],
            tb,
        )
    directional = (
        line[-2] - 8.0 * line[-1] + 8.0 * line[1] - line[2]
    ) / (12.0 * h)
    advection = speed[..., None] * directional

    laplacian = np.zeros_like(u0)
    base = (xb, yb, zb)
    for axis in range(3):
        shifted_values: dict[int, np.ndarray] = {}
        for offset in (-2, -1, 1, 2):
            shifted = [a for a in base]
            shifted[axis] = shifted[axis] + offset * h
            shifted_values[offset] = value(
                shifted[0], shifted[1], shifted[2], tb
            )
        laplacian += (
            -shifted_values[-2]
            + 16.0 * shifted_values[-1]
            - 30.0 * u0
            + 16.0 * shifted_values[1]
            - shifted_values[2]
        ) / (12.0 * h * h)

    out = dt + advection - VISCOSITY * laplacian
    if not np.all(np.isfinite(out)):
        raise RuntimeError("independent FD4 transport became non-finite")
    return out


def public_contract() -> dict[str, Any]:
    """Machine-readable scope/provenance boundary for downstream agents."""
    signature = inspect.signature(evaluate_inner_leading_oscillatory_transport)
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "alpha", "damping", "delta_y", "delta_a", "nu", "viscosity",
        "scientific_threshold", "spatial_step", "time_step",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "agent1_laplacian_pr": AGENT1_LAPLACIAN_PR,
        "agent1_laplacian_head": AGENT1_LAPLACIAN_HEAD,
        "agent1_laplacian_blob": AGENT1_LAPLACIAN_BLOB,
        "viscosity": VISCOSITY,
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        },
        "forbidden_public_inputs_present": sorted(forbidden.intersection(signature.parameters)),
        "truth_boundary": {
            "strict_inner_transport_executable": True,
            "viscosity_fixed_not_tunable": True,
            "inner_plus_oscillatory_laplacian_executable": True,
            "agent1_exact_head_ci_assumed_passed": False,
            "agent4_inner_field_independent_admission_assumed": False,
            "agent4_a2_oscillatory_independent_admission_assumed": False,
            "inner_leading_is_final_corrected_fixed_point": False,
            "global_leading_velocity_materialized": False,
            "outer_join_materialized": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
            "agent3_mean_radial_chain_reimplemented": False,
            "correction_velocity_included": False,
            "complete_kokuno_candidate_assembled": False,
            "complete_ns_residual": False,
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "paper_exact": False,
            "pde_validated": False,
        },
    }
