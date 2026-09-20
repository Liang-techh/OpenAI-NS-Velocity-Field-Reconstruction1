"""Strict-inner PA.10 plus oscillatory convective composition for Kokuno Agent 2.

This module closes one narrow nonlinear seam for the currently executable
*inner* Agent-1 PA.10 contraction-center field and the frozen Agent-2
oscillatory complete-curl field.  It exposes

    u_io = u_inner + u_osc

and

    (u_io . grad)u_io
      = (u_inner . grad)u_inner
      + (u_inner . grad)u_osc
      + (u_osc . grad)u_inner
      + (u_osc . grad)u_osc.

Agent 1 owns the analytic inner self-advection.  Agent 2 owns the frozen
oscillatory self-advection and the already-existing mixed cross-advection seam.
This module only performs the typed composition and does not invent a global
leading continuation, pressure, forcing, residual target, or Agent-3 mean/radial
correction.

Source/realization boundary
---------------------------
The corrected 2026-09-09 Kokuno reconstruction is provenance for the upstream
PA.10 inner formulas and localized complete-curl oscillatory structure.  The
public-z pullback, repository-autonomous oscillatory time law, Cartesian FD
operators in the inherited A2 diagnostics, and this nonlinear composition are
repository realizations.  Nothing here is paper-exact or a complete
Navier--Stokes residual.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
import math
from typing import Any, Protocol, runtime_checkable

import numpy as np

from .kokuno_public_inner_leading_oscillatory_cross_advection import (
    INNER_LEADING_FD4_STEP,
    OSCILLATORY_FD6_STEP,
    evaluate_inner_leading_oscillatory_cross_advection,
)
from .kokuno_public_oscillatory_self_advection import (
    evaluate_oscillatory_self_advection_fd6,
)
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-INNER-LEADING-OSCILLATORY-FULL-ADVECTION-061"
SCHEMA = "kokuno-a2-inner-leading-oscillatory-full-advection-v1"
PARENT_AGENT2_PR = 821
PARENT_AGENT2_HEAD = "110c7150df73114c482e9f1673ef3a928b95053e"
AGENT1_SELF_ADVECTION_PR = 829
AGENT1_SELF_ADVECTION_HEAD = "177190d770e0306e116674c5a2e1cbf977628760"
AGENT1_SELF_ADVECTION_BLOB = "b75f3b08ef3f1853550d8dfdef912980c15287c3"
AGENT2_CROSS_PR = 812
AGENT2_OSCILLATORY_SELF_PR = 779
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
DIRECTIONAL_FD4_STEPS = (4.0e-3, 2.0e-3, 1.0e-3)


@runtime_checkable
class InnerLeadingSelfAdvectionBackend(Protocol):
    """Minimal Agent-1-facing protocol; global leading semantics are not implied."""

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...

    def self_advection(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...


@dataclass(frozen=True)
class InnerLeadingOscillatoryFullAdvectionResult:
    inner_leading_velocity: np.ndarray
    oscillatory_velocity: np.ndarray
    inner_plus_oscillatory_velocity: np.ndarray
    inner_leading_self_advection: np.ndarray
    inner_advects_oscillation: np.ndarray
    oscillation_advects_inner: np.ndarray
    mixed_cross_advection: np.ndarray
    oscillatory_self_advection: np.ndarray
    inner_plus_oscillatory_self_advection: np.ndarray
    inner_leading_spatial_step: float
    oscillatory_spatial_step: float


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


def _validate_step(value: float, name: str) -> float:
    step = float(value)
    if not math.isfinite(step) or not 1.0e-5 <= step <= 5.0e-2:
        raise ValueError(f"{name} must be finite and lie in [1e-5,5e-2]")
    return step


def _backend_vector(
    backend: InnerLeadingSelfAdvectionBackend,
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
        raise RuntimeError(f"{label} did not replay exactly across inherited A2 seams")


def evaluate_inner_leading_oscillatory_full_advection(
    inner_leading_backend: InnerLeadingSelfAdvectionBackend,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    inner_leading_spatial_step: float = INNER_LEADING_FD4_STEP,
    oscillatory_spatial_step: float = OSCILLATORY_FD6_STEP,
) -> InnerLeadingOscillatoryFullAdvectionResult:
    """Return the strict-inner composite and its complete convective expansion.

    The pure inner term comes from the supplied Agent-1 backend.  The mixed and
    pure oscillatory terms are delegated to the already-frozen A2 seams.  Any
    Agent-1 inner-domain failure propagates; no outer continuation is invented.
    """
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    h_inner = _validate_step(inner_leading_spatial_step, "inner_leading_spatial_step")
    h_osc = _validate_step(oscillatory_spatial_step, "oscillatory_spatial_step")

    u_inner = _backend_vector(inner_leading_backend, "velocity", xb, yb, zb, tb)
    n_inner = _backend_vector(inner_leading_backend, "self_advection", xb, yb, zb, tb)

    cross = evaluate_inner_leading_oscillatory_cross_advection(
        inner_leading_backend,
        xb,
        yb,
        zb,
        tb,
        inner_leading_spatial_step=h_inner,
        oscillatory_spatial_step=h_osc,
    )
    osc = evaluate_oscillatory_self_advection_fd6(
        xb, yb, zb, tb, spatial_step=h_osc
    )
    u_osc = np.asarray(osc["velocity"], dtype=float)
    n_osc = np.asarray(osc["self_advection"], dtype=float)

    _require_exact_replay(u_inner, cross.inner_leading_velocity, "inner velocity")
    _require_exact_replay(u_osc, cross.oscillatory_velocity, "oscillatory velocity")

    total_velocity = u_inner + u_osc
    total_advection = n_inner + cross.cross_advection + n_osc
    if not np.all(np.isfinite(total_velocity)) or not np.all(np.isfinite(total_advection)):
        raise RuntimeError("inner-leading/oscillatory full advection became non-finite")

    return InnerLeadingOscillatoryFullAdvectionResult(
        inner_leading_velocity=u_inner,
        oscillatory_velocity=u_osc,
        inner_plus_oscillatory_velocity=total_velocity,
        inner_leading_self_advection=n_inner,
        inner_advects_oscillation=np.asarray(cross.inner_advects_oscillation, dtype=float),
        oscillation_advects_inner=np.asarray(cross.oscillation_advects_inner, dtype=float),
        mixed_cross_advection=np.asarray(cross.cross_advection, dtype=float),
        oscillatory_self_advection=n_osc,
        inner_plus_oscillatory_self_advection=total_advection,
        inner_leading_spatial_step=h_inner,
        oscillatory_spatial_step=h_osc,
    )


def directional_fd4_inner_plus_oscillatory_self_advection(
    inner_leading_backend: InnerLeadingSelfAdvectionBackend,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    directional_step: float,
) -> np.ndarray:
    """Independent public-velocity-only line derivative of the full sum.

    The direction is frozen to ``(u_inner+u_osc)/|u_inner+u_osc|`` at each base
    point.  This verifier never calls Agent-1 ``self_advection`` nor any A2
    Cartesian Jacobian/self-advection/cross-advection evaluator.
    """
    h = _validate_step(directional_step, "directional_step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)

    def total_velocity(xx: np.ndarray, yy: np.ndarray, zz: np.ndarray, tt: np.ndarray) -> np.ndarray:
        inner = _backend_vector(inner_leading_backend, "velocity", xx, yy, zz, tt)
        oscillatory = np.asarray(velocity_osc(xx, yy, zz, tt), dtype=float)
        expected = xx.shape + (3,)
        if oscillatory.shape != expected or not np.all(np.isfinite(oscillatory)):
            raise RuntimeError("velocity_osc directional evaluation failed")
        return inner + oscillatory

    base = total_velocity(xb, yb, zb, tb)
    speed = np.linalg.norm(base, axis=-1)
    if np.any(speed <= 1.0e-12):
        raise ValueError("directional audit requires nonzero composite velocity")
    direction = base / speed[..., None]

    values: dict[int, np.ndarray] = {}
    for offset in (-2, -1, 1, 2):
        displacement = offset * h * direction
        values[offset] = total_velocity(
            xb + displacement[..., 0],
            yb + displacement[..., 1],
            zb + displacement[..., 2],
            tb,
        )
    derivative = (
        values[-2] - 8.0 * values[-1] + 8.0 * values[1] - values[2]
    ) / (12.0 * h)
    out = speed[..., None] * derivative
    if not np.all(np.isfinite(out)):
        raise RuntimeError("directional composite self-advection became non-finite")
    return out


def public_contract() -> dict[str, Any]:
    """Machine-readable provenance and fail-closed scientific scope."""
    signature = inspect.signature(evaluate_inner_leading_oscillatory_full_advection)
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "alpha", "damping", "delta_y", "delta_a", "nu", "viscosity",
        "scientific_threshold",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "agent1_self_advection_pr": AGENT1_SELF_ADVECTION_PR,
        "agent1_self_advection_head": AGENT1_SELF_ADVECTION_HEAD,
        "agent1_self_advection_blob": AGENT1_SELF_ADVECTION_BLOB,
        "agent2_cross_pr": AGENT2_CROSS_PR,
        "agent2_oscillatory_self_pr": AGENT2_OSCILLATORY_SELF_PR,
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        },
        "forbidden_public_inputs_present": sorted(forbidden.intersection(signature.parameters)),
        "truth_boundary": {
            "inner_pa10_self_advection_consumable": True,
            "inner_plus_oscillatory_velocity_executable": True,
            "inner_plus_oscillatory_full_advection_executable": True,
            "production_reuses_agent1_analytic_inner_self_advection": True,
            "production_reuses_a2_complete_curl_oscillatory_field": True,
            "direct_verifier_uses_only_public_velocities": True,
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
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "paper_exact": False,
            "pde_validated": False,
        },
    }
