"""Oscillatory transport increment on the strict-inner PA.10 background.

This Kokuno Agent-2 seam answers one narrow question without promoting it to a
Navier--Stokes residual: how much does the frozen localized complete-curl
oscillatory field change the pressure/forcing-free local transport of the
currently executable Agent-1 PA.10 *inner contraction-center* field?

With repository viscosity fixed at ``nu=0.01`` define

    T_inner = d_t u_inner + (u_inner . grad)u_inner - nu Delta u_inner,

    T_inner+osc = d_t(u_inner+u_osc)
                  + ((u_inner+u_osc) . grad)(u_inner+u_osc)
                  - nu Delta(u_inner+u_osc),

and expose the typed increment

    delta_T_osc = T_inner+osc - T_inner
                = d_t u_osc
                  + [(u_inner . grad)u_osc + (u_osc . grad)u_inner
                     + (u_osc . grad)u_osc]
                  - nu Delta u_osc.

Agent 1 owns the strict-inner leading field and its derivatives. Agent 2 owns
the frozen complete-curl oscillatory field and the inherited typed composition.
This module adds no amplitude/phase/carrier/support freedom, pressure, forcing,
residual target, correction velocity, or Agent-3 mean/radial operation.

The corrected Kokuno 2026-09-09 reconstruction remains structural provenance
for the upstream PA.10 and localized complete-curl organization. The public-z
pullback, repository-autonomous oscillatory choices, finite-difference
realizations, and this transport-delta diagnostic are repository realizations;
none is claimed paper-exact.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
import math
from typing import Any

import numpy as np

from .kokuno_public_inner_leading_oscillatory_transport import (
    AGENT1_LAPLACIAN_BLOB,
    AGENT1_LAPLACIAN_HEAD,
    AGENT1_LAPLACIAN_PR,
    INDEPENDENT_FD4_STEPS,
    InnerLeadingTransportBackend,
    VISCOSITY,
    evaluate_inner_leading_oscillatory_transport,
    independent_fd4_inner_plus_oscillatory_transport,
)

TASK = "KOKUNO-A2-INNER-LEADING-OSCILLATORY-TRANSPORT-DELTA-063"
SCHEMA = "kokuno-a2-inner-leading-oscillatory-transport-delta-v1"
PARENT_AGENT2_PR = 840
PARENT_AGENT2_HEAD = "016e3c152d7d11e4fedcd79b819df337f8942983"
PARENT_AGENT2_TRANSPORT_BLOB = "d6082ebd346ba8150321bf3d62f104c107a70300"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"


@dataclass(frozen=True)
class InnerLeadingOscillatoryTransportDeltaResult:
    inner_leading_transport: np.ndarray
    inner_plus_oscillatory_transport: np.ndarray
    oscillatory_transport_increment: np.ndarray
    oscillatory_time_increment: np.ndarray
    oscillatory_nonlinear_increment: np.ndarray
    oscillatory_viscous_increment: np.ndarray
    inner_leading_velocity: np.ndarray
    oscillatory_velocity: np.ndarray
    inner_plus_oscillatory_velocity: np.ndarray
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


def evaluate_inner_leading_oscillatory_transport_delta(
    inner_leading_backend: InnerLeadingTransportBackend,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
) -> InnerLeadingOscillatoryTransportDeltaResult:
    """Return the strict-inner transport before/after the frozen oscillation.

    This quantity intentionally excludes pressure and forcing and therefore is
    *not* a complete momentum residual. Any Agent-1 inner-domain/stencil failure
    propagates fail-closed.
    """
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    total = evaluate_inner_leading_oscillatory_transport(
        inner_leading_backend, xb, yb, zb, tb
    )

    inner_dt = _backend_vector(
        inner_leading_backend, "velocity_dt", xb, yb, zb, tb
    )
    inner_advection = _backend_vector(
        inner_leading_backend, "self_advection", xb, yb, zb, tb
    )
    inner_laplacian = _backend_vector(
        inner_leading_backend, "velocity_laplacian", xb, yb, zb, tb
    )

    if not np.array_equal(inner_dt, total.inner_leading_velocity_dt):
        raise RuntimeError("inner velocity_dt did not replay exactly across parent seam")
    if not np.array_equal(inner_laplacian, total.inner_leading_laplacian):
        raise RuntimeError("inner Laplacian did not replay exactly across parent seam")

    inner_transport = inner_dt + inner_advection - VISCOSITY * inner_laplacian
    nonlinear_increment = total.inner_plus_oscillatory_self_advection - inner_advection
    viscous_increment = -VISCOSITY * total.oscillatory_laplacian
    oscillatory_increment = (
        total.oscillatory_velocity_dt
        + nonlinear_increment
        + viscous_increment
    )

    delta_from_total = total.transport - inner_transport
    closure = float(np.max(np.abs(oscillatory_increment - delta_from_total)))
    if closure > 5.0e-13:
        raise RuntimeError(
            "oscillatory transport increment failed additive/decomposition closure"
        )

    expected = xb.shape + (3,)
    for name, value in (
        ("inner_transport", inner_transport),
        ("total_transport", total.transport),
        ("oscillatory_increment", oscillatory_increment),
        ("oscillatory_time_increment", total.oscillatory_velocity_dt),
        ("oscillatory_nonlinear_increment", nonlinear_increment),
        ("oscillatory_viscous_increment", viscous_increment),
    ):
        if value.shape != expected or not np.all(np.isfinite(value)):
            raise RuntimeError(f"{name} is non-finite or has the wrong shape")

    return InnerLeadingOscillatoryTransportDeltaResult(
        inner_leading_transport=inner_transport,
        inner_plus_oscillatory_transport=np.asarray(total.transport, dtype=float),
        oscillatory_transport_increment=oscillatory_increment,
        oscillatory_time_increment=np.asarray(total.oscillatory_velocity_dt, dtype=float),
        oscillatory_nonlinear_increment=nonlinear_increment,
        oscillatory_viscous_increment=viscous_increment,
        inner_leading_velocity=np.asarray(total.inner_leading_velocity, dtype=float),
        oscillatory_velocity=np.asarray(total.oscillatory_velocity, dtype=float),
        inner_plus_oscillatory_velocity=np.asarray(
            total.inner_plus_oscillatory_velocity, dtype=float
        ),
        viscosity=VISCOSITY,
    )


def _validate_step(step: float) -> float:
    h = float(step)
    if not math.isfinite(h) or not 1.0e-5 <= h <= 5.0e-2:
        raise ValueError("step must be finite and lie in [1e-5,5e-2]")
    return h


def independent_fd4_inner_transport(
    inner_leading_backend: InnerLeadingTransportBackend,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    """Public-inner-velocity-only FD4 reference for ``T_inner``.

    This path deliberately does not call Agent-1 production ``velocity_dt``,
    ``self_advection`` or ``velocity_laplacian``.
    """
    h = _validate_step(step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)

    def value(xx: np.ndarray, yy: np.ndarray, zz: np.ndarray, tt: np.ndarray) -> np.ndarray:
        return _backend_vector(inner_leading_backend, "velocity", xx, yy, zz, tt)

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
        raise RuntimeError("independent inner FD4 transport became non-finite")
    return out


def independent_fd4_oscillatory_transport_delta(
    inner_leading_backend: InnerLeadingTransportBackend,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    """Independent public-velocity-only FD4 before/after transport difference."""
    total = independent_fd4_inner_plus_oscillatory_transport(
        inner_leading_backend, x, y, z, t, step=step
    )
    inner = independent_fd4_inner_transport(
        inner_leading_backend, x, y, z, t, step=step
    )
    out = total - inner
    if not np.all(np.isfinite(out)):
        raise RuntimeError("independent oscillatory transport delta became non-finite")
    return out


def public_contract() -> dict[str, Any]:
    """Machine-readable scope/provenance boundary for downstream agents."""
    signature = inspect.signature(evaluate_inner_leading_oscillatory_transport_delta)
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
        "parent_agent2_transport_blob": PARENT_AGENT2_TRANSPORT_BLOB,
        "agent1_laplacian_pr": AGENT1_LAPLACIAN_PR,
        "agent1_laplacian_head": AGENT1_LAPLACIAN_HEAD,
        "agent1_laplacian_blob": AGENT1_LAPLACIAN_BLOB,
        "viscosity": VISCOSITY,
        "independent_fd4_steps": list(INDEPENDENT_FD4_STEPS),
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        },
        "forbidden_public_inputs_present": sorted(forbidden.intersection(signature.parameters)),
        "truth_boundary": {
            "strict_inner_before_after_transport_executable": True,
            "oscillatory_transport_increment_executable": True,
            "viscosity_fixed_not_tunable": True,
            "diagnostic_is_complete_ns_residual": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
            "correction_velocity_included": False,
            "inner_leading_is_final_corrected_fixed_point": False,
            "global_leading_velocity_materialized": False,
            "outer_join_materialized": False,
            "agent3_mean_radial_chain_reimplemented": False,
            "complete_kokuno_candidate_assembled": False,
            "heldout_ns_residual_assessed": False,
            "st006_same_protocol_comparison_valid": False,
            "residual_reduction_claimed": False,
            "paper_exact": False,
            "pde_validated": False,
        },
    }
