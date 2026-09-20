"""Typed inner-leading/oscillatory time-derivative seam for Kokuno Agent 2.

This module is intentionally inner-only.  Agent 1 PR #811 exposes an analytic
fixed-Cartesian-point time derivative for the executable PA.10 contraction-center
velocity inherited from #803.  That field is not the final corrected fixed point,
has no outer/global join, and has no matched global pressure.  Agent 2 therefore
consumes it only through a minimal ``velocity`` + ``velocity_dt`` backend protocol.

For an admitted call inside the Agent-1 inner domain this module exposes

    u_inner+osc = u_inner + u_osc,
    d_t u_inner+osc = d_t u_inner + d_t u_osc.

The oscillatory derivative is the already-frozen Agent-2 derivative of the
repository-autonomous covariance modulation.  No finite difference is used in
the production composition.  Any Agent-1 domain failure propagates; A2 does not
invent an outer continuation.

Source/realization boundary
---------------------------
The corrected 2026-09-09 Kokuno reconstruction is provenance for the upstream
PA.10 inner formulas and localized complete-curl oscillatory structure.  The
public-z pullback, repository-autonomous oscillatory time law, and this additive
composition are repository realizations.  Nothing here is paper-exact, a global
Kokuno field, or an independent Navier--Stokes validation.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
from typing import Any, Protocol, runtime_checkable

import numpy as np

from .kokuno_public_oscillatory_time_derivative import velocity_osc_dt
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-INNER-LEADING-OSCILLATORY-TIME-060"
SCHEMA = "kokuno-a2-inner-leading-oscillatory-time-v1"
PARENT_AGENT2_PR = 812
PARENT_AGENT2_HEAD = "39a55b6f80a4bf9dfa6e0d6bee7825ef0cf7e600"
AGENT1_PR = 811
AGENT1_HEAD = "6f16e837c3a95a6e54af5c2cb7725d094f667026"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
FD4_TIME_STEPS = (4.0e-3, 2.0e-3, 1.0e-3)


@runtime_checkable
class InnerLeadingVelocityTimeBackend(Protocol):
    """Minimal Agent-1-facing protocol; global leading semantics are not implied."""

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...


@dataclass(frozen=True)
class InnerLeadingOscillatoryTimeResult:
    inner_leading_velocity: np.ndarray
    oscillatory_velocity: np.ndarray
    inner_plus_oscillatory_velocity: np.ndarray
    inner_leading_velocity_dt: np.ndarray
    oscillatory_velocity_dt: np.ndarray
    inner_plus_oscillatory_velocity_dt: np.ndarray


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
    backend: InnerLeadingVelocityTimeBackend,
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


def evaluate_inner_leading_oscillatory_time_derivative(
    inner_leading_backend: InnerLeadingVelocityTimeBackend,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
) -> InnerLeadingOscillatoryTimeResult:
    """Return the strict-inner composite and its analytic/componentwise time derivative.

    The Agent-1 backend owns the PA.10 inner velocity and its time derivative.
    Agent 2 owns the unchanged oscillatory field and its frozen candidate time
    derivative.  This function only performs the typed additive composition.
    """
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    u_inner = _backend_vector(inner_leading_backend, "velocity", xb, yb, zb, tb)
    dt_inner = _backend_vector(inner_leading_backend, "velocity_dt", xb, yb, zb, tb)

    u_osc = np.asarray(velocity_osc(xb, yb, zb, tb), dtype=float)
    dt_osc = np.asarray(velocity_osc_dt(xb, yb, zb, tb), dtype=float)
    expected = xb.shape + (3,)
    if u_osc.shape != expected or dt_osc.shape != expected:
        raise RuntimeError("oscillatory velocity/time derivative returned an unexpected shape")
    if not np.all(np.isfinite(u_osc)) or not np.all(np.isfinite(dt_osc)):
        raise RuntimeError("oscillatory velocity/time derivative produced non-finite values")

    total = u_inner + u_osc
    dt_total = dt_inner + dt_osc
    if not np.all(np.isfinite(total)) or not np.all(np.isfinite(dt_total)):
        raise RuntimeError("inner-leading/oscillatory time composition became non-finite")

    return InnerLeadingOscillatoryTimeResult(
        inner_leading_velocity=u_inner,
        oscillatory_velocity=u_osc,
        inner_plus_oscillatory_velocity=total,
        inner_leading_velocity_dt=dt_inner,
        oscillatory_velocity_dt=dt_osc,
        inner_plus_oscillatory_velocity_dt=dt_total,
    )


def public_contract() -> dict[str, Any]:
    """Machine-readable truth/provenance boundary for downstream agents."""
    signature = inspect.signature(evaluate_inner_leading_oscillatory_time_derivative)
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "alpha", "damping", "delta_y", "delta_a", "nu", "viscosity",
        "scientific_threshold", "time_step",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "agent1_inner_velocity_dt_pr": AGENT1_PR,
        "agent1_inner_velocity_dt_head": AGENT1_HEAD,
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        },
        "forbidden_public_inputs_present": sorted(forbidden.intersection(signature.parameters)),
        "truth_boundary": {
            "inner_pa10_velocity_consumable": True,
            "inner_pa10_velocity_dt_consumable": True,
            "inner_plus_oscillatory_velocity_executable": True,
            "inner_plus_oscillatory_velocity_dt_executable": True,
            "production_composite_dt_uses_finite_difference": False,
            "agent1_exact_head_ci_assumed_passed": False,
            "agent4_inner_field_independent_admission_assumed": False,
            "inner_leading_is_final_corrected_fixed_point": False,
            "global_leading_velocity_materialized": False,
            "outer_join_materialized": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
            "agent3_mean_radial_chain_reimplemented": False,
            "complete_kokuno_candidate_assembled": False,
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "paper_exact": False,
            "pde_validated": False,
        },
    }
