"""Typed inner-leading/oscillatory cross-advection seam for Kokuno Agent 2.

This module is deliberately narrower than a complete Kokuno candidate.  Agent 1
PR #803 is the first stack surface that exposes an executable *inner PA.10
contraction-center* Cartesian velocity.  It is not the corrected fixed point,
has no outer/global join, and has unresolved exact-head admission.  Agent 2
therefore consumes it only through a small ``velocity(x,y,z,t)`` backend
protocol and combines it with the already-frozen A2 oscillatory field.

For an admitted call inside the Agent-1 inner domain we expose

    u_inner_plus_osc = u_inner + u_osc

and the raw mixed convective contribution

    N_inner,osc = (u_inner . grad) u_osc + (u_osc . grad) u_inner.

The oscillatory Jacobian reuses A2's existing centered Cartesian FD6 diagnostic.
The inner-leading Jacobian is an independent centered Cartesian FD4 realization
of the supplied backend.  No pressure, forcing, residual target, Agent-3 mean,
or radial inverse is accepted here.

Source/realization boundary
---------------------------
The corrected Kokuno reconstruction is provenance for the upstream PA.10
coordinate/velocity formulas and localized complete-curl oscillatory structure.
The finite-difference derivatives and this mixed quadratic composition are
repository numerical realizations.  Nothing here is paper-exact or a full NS
residual.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
import math
from typing import Any, Protocol, runtime_checkable

import numpy as np

from .kokuno_public_oscillatory_vorticity_diagnostic import evaluate_vorticity_osc_fd6
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-INNER-LEADING-OSCILLATORY-CROSS-059"
SCHEMA = "kokuno-a2-inner-leading-oscillatory-cross-v1"
PARENT_AGENT2_PR = 804
PARENT_AGENT2_HEAD = "9bd1623f238a7827019a336aee435216fc38d393"
AGENT1_PR = 803
AGENT1_HEAD = "0666fff748d7c2659774b0aacba70042d045aab5"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
OSCILLATORY_FD6_STEP = 1.0e-3
INNER_LEADING_FD4_STEP = 1.0e-3
DIRECTIONAL_FD4_STEPS = (4.0e-3, 2.0e-3, 1.0e-3)


@runtime_checkable
class InnerLeadingVelocityBackend(Protocol):
    """Minimal Agent-1-facing protocol; global-field semantics are not implied."""

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...


@dataclass(frozen=True)
class InnerLeadingOscillatoryCrossResult:
    inner_leading_velocity: np.ndarray
    oscillatory_velocity: np.ndarray
    inner_plus_oscillatory_velocity: np.ndarray
    inner_leading_jacobian: np.ndarray
    oscillatory_jacobian: np.ndarray
    inner_advects_oscillation: np.ndarray
    oscillation_advects_inner: np.ndarray
    cross_advection: np.ndarray
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


def _backend_velocity(
    backend: InnerLeadingVelocityBackend,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    if not hasattr(backend, "velocity") or not callable(getattr(backend, "velocity")):
        raise TypeError("inner_leading_backend must expose callable velocity(x,y,z,t)")
    value = np.asarray(backend.velocity(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise ValueError(f"inner leading velocity must have shape {expected}, got {value.shape}")
    if not np.all(np.isfinite(value)):
        raise ValueError("inner leading velocity must be finite")
    return value


def _inner_leading_jacobian_fd4(
    backend: InnerLeadingVelocityBackend,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    spatial_step: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return public backend velocity and J[component,axis] by centered FD4."""
    h = _validate_step(spatial_step, "inner_leading_spatial_step")
    base = _backend_velocity(backend, x, y, z, t)
    jac = np.empty(x.shape + (3, 3), dtype=float)
    coords = [x, y, z]
    for axis in range(3):
        values: dict[int, np.ndarray] = {}
        for offset in (-2, -1, 1, 2):
            shifted = [c.copy() for c in coords]
            shifted[axis] = shifted[axis] + offset * h
            values[offset] = _backend_velocity(
                backend, shifted[0], shifted[1], shifted[2], t
            )
        jac[..., :, axis] = (
            values[-2] - 8.0 * values[-1] + 8.0 * values[1] - values[2]
        ) / (12.0 * h)
    if not np.all(np.isfinite(jac)):
        raise RuntimeError("inner leading FD4 Jacobian produced non-finite values")
    return base, jac


def _contract(jacobian: np.ndarray, advecting_velocity: np.ndarray) -> np.ndarray:
    return np.einsum("...ij,...j->...i", jacobian, advecting_velocity)


def evaluate_inner_leading_oscillatory_cross_advection(
    inner_leading_backend: InnerLeadingVelocityBackend,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    inner_leading_spatial_step: float = INNER_LEADING_FD4_STEP,
    oscillatory_spatial_step: float = OSCILLATORY_FD6_STEP,
) -> InnerLeadingOscillatoryCrossResult:
    """Evaluate the inner-center + oscillation sum and their mixed advection.

    Any Agent-1 inner-domain failure is allowed to propagate.  This is
    intentional: A2 does not invent an outer continuation for the inner field.
    """
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    h_lead = _validate_step(inner_leading_spatial_step, "inner_leading_spatial_step")
    h_osc = _validate_step(oscillatory_spatial_step, "oscillatory_spatial_step")

    u_lead, j_lead = _inner_leading_jacobian_fd4(
        inner_leading_backend, xb, yb, zb, tb, spatial_step=h_lead
    )
    osc = evaluate_vorticity_osc_fd6(xb, yb, zb, tb, spatial_step=h_osc)
    u_osc = np.asarray(osc["velocity"], dtype=float)
    j_osc = np.asarray(osc["velocity_gradient_fd6"], dtype=float)
    expected_velocity = xb.shape + (3,)
    expected_jacobian = xb.shape + (3, 3)
    if u_osc.shape != expected_velocity or j_osc.shape != expected_jacobian:
        raise RuntimeError("oscillatory diagnostic returned an unexpected shape")

    lead_advects_osc = _contract(j_osc, u_lead)
    osc_advects_lead = _contract(j_lead, u_osc)
    cross = lead_advects_osc + osc_advects_lead
    total = u_lead + u_osc
    if not all(
        np.all(np.isfinite(a))
        for a in (u_osc, j_osc, lead_advects_osc, osc_advects_lead, cross, total)
    ):
        raise RuntimeError("inner-leading/oscillatory composition produced non-finite values")

    return InnerLeadingOscillatoryCrossResult(
        inner_leading_velocity=u_lead,
        oscillatory_velocity=u_osc,
        inner_plus_oscillatory_velocity=total,
        inner_leading_jacobian=j_lead,
        oscillatory_jacobian=j_osc,
        inner_advects_oscillation=lead_advects_osc,
        oscillation_advects_inner=osc_advects_lead,
        cross_advection=cross,
        inner_leading_spatial_step=h_lead,
        oscillatory_spatial_step=h_osc,
    )


def directional_fd4_cross_advection(
    inner_leading_backend: InnerLeadingVelocityBackend,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    directional_step: float,
) -> np.ndarray:
    """Independent line-derivative verifier for the two mixed terms.

    It calls only the two public velocity providers.  It never calls either
    Cartesian Jacobian used by the production evaluator.
    """
    h = _validate_step(directional_step, "directional_step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    u_lead = _backend_velocity(inner_leading_backend, xb, yb, zb, tb)
    u_osc = np.asarray(velocity_osc(xb, yb, zb, tb), dtype=float)
    if u_osc.shape != xb.shape + (3,):
        raise RuntimeError("velocity_osc returned an unexpected shape")

    def line_derivative(
        evaluator: Any, advector: np.ndarray, *, label: str
    ) -> np.ndarray:
        speed = np.linalg.norm(advector, axis=-1)
        if np.any(speed <= 1.0e-12):
            raise ValueError(f"{label} directional audit requires nonzero advecting velocity")
        direction = advector / speed[..., None]
        samples: dict[int, np.ndarray] = {}
        for offset in (-2, -1, 1, 2):
            d = offset * h * direction
            samples[offset] = evaluator(
                xb + d[..., 0], yb + d[..., 1], zb + d[..., 2], tb
            )
        derivative = (
            samples[-2] - 8.0 * samples[-1] + 8.0 * samples[1] - samples[2]
        ) / (12.0 * h)
        return speed[..., None] * derivative

    def eval_osc(xx: np.ndarray, yy: np.ndarray, zz: np.ndarray, tt: np.ndarray) -> np.ndarray:
        out = np.asarray(velocity_osc(xx, yy, zz, tt), dtype=float)
        if out.shape != xx.shape + (3,) or not np.all(np.isfinite(out)):
            raise RuntimeError("velocity_osc directional evaluation failed")
        return out

    def eval_lead(xx: np.ndarray, yy: np.ndarray, zz: np.ndarray, tt: np.ndarray) -> np.ndarray:
        return _backend_velocity(inner_leading_backend, xx, yy, zz, tt)

    return line_derivative(eval_osc, u_lead, label="inner-leading") + line_derivative(
        eval_lead, u_osc, label="oscillatory"
    )


def public_contract() -> dict[str, Any]:
    """Machine-readable truth/provenance boundary for downstream agents."""
    signature = inspect.signature(evaluate_inner_leading_oscillatory_cross_advection)
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "alpha", "damping", "delta_y", "delta_a", "nu", "viscosity",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "agent1_inner_velocity_pr": AGENT1_PR,
        "agent1_inner_velocity_head": AGENT1_HEAD,
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        },
        "forbidden_public_inputs_present": sorted(forbidden.intersection(signature.parameters)),
        "truth_boundary": {
            "inner_pa10_contraction_center_velocity_consumable": True,
            "inner_plus_oscillatory_velocity_executable": True,
            "inner_oscillatory_cross_advection_executable": True,
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
