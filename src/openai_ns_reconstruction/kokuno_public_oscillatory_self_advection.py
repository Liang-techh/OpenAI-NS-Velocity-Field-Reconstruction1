"""Nonlinear self-advection diagnostic for the frozen Kokuno Agent-2 oscillation.

This module adds no oscillatory degree of freedom.  It consumes the already
frozen public velocity ``velocity_osc`` and the existing FD6 spatial-gradient
diagnostic, then exposes the quadratic term

    (u_osc . grad) u_osc.

That term is needed later when a real leading field is available and the
staged Navier--Stokes momentum contribution of ``u_lead + u_osc`` is assembled.
It is intentionally kept separate from Agent-3 mean/radial correction logic:
this module neither averages the nonlinear term into a mean defect nor builds a
correction from it.

Provenance boundary
-------------------
The corrected Kokuno reconstruction is structural provenance for the upstream
localized complete-curl oscillatory field.  The Cartesian FD6 gradient and the
quadratic composition in this module are executable repository diagnostics,
not hidden paper data and not a paper-exact Navier--Stokes residual.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_public_oscillatory_vorticity_diagnostic import (
    evaluate_vorticity_osc_fd6,
)
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-PUBLIC-OSCILLATORY-SELF-ADVECTION-055"
SCHEMA = "kokuno-a2-public-oscillatory-self-advection-v1"
PARENT_AGENT2_PR = 769
PARENT_AGENT2_HEAD = "f298116969c4096cccc38d037656f59cab509ed1"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
DIRECTIONAL_FD4_STEPS = (0.008, 0.004, 0.002)
GRADIENT_FD6_STEP = 0.001


def _broadcast_xyzt(x: Any, y: Any, z: Any, t: Any) -> tuple[np.ndarray, ...]:
    arrays = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    if not all(np.all(np.isfinite(array)) for array in arrays):
        raise ValueError("x, y, z and t must be finite and broadcastable")
    return tuple(arrays)


def _validate_step(step: float) -> float:
    h = float(step)
    if not math.isfinite(h) or not (1.0e-5 <= h <= 5.0e-2):
        raise ValueError("spatial_step must be finite and lie in [1e-5, 5e-2]")
    return h


def _self_advection_from_jacobian(
    velocity: np.ndarray, jacobian: np.ndarray
) -> np.ndarray:
    """Return ``J @ u`` for ``J[component, derivative_axis]``."""
    u = np.asarray(velocity, dtype=float)
    j = np.asarray(jacobian, dtype=float)
    if u.shape[-1:] != (3,):
        raise ValueError("velocity must end in component dimension 3")
    if j.shape != u.shape[:-1] + (3, 3):
        raise ValueError("jacobian must have shape velocity.shape[:-1] + (3,3)")
    if not np.all(np.isfinite(u)) or not np.all(np.isfinite(j)):
        raise ValueError("velocity and jacobian must be finite")
    return np.einsum("...ij,...j->...i", j, u)


def evaluate_oscillatory_self_advection_fd6(
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float = GRADIENT_FD6_STEP,
) -> dict[str, np.ndarray | float | str]:
    """Evaluate frozen ``u_osc`` and ``(u_osc . grad)u_osc``.

    The gradient comes from the existing centered Cartesian FD6 diagnostic.
    No pressure, forcing, residual, correction target, gain, damping, or
    Agent-3 mean/radial state is accepted by this API.
    """
    h = _validate_step(spatial_step)
    diagnostic = evaluate_vorticity_osc_fd6(x, y, z, t, spatial_step=h)
    velocity = np.asarray(diagnostic["velocity"], dtype=float)
    gradient = np.asarray(diagnostic["velocity_gradient_fd6"], dtype=float)
    self_advection = _self_advection_from_jacobian(velocity, gradient)
    if not np.all(np.isfinite(self_advection)):
        raise RuntimeError("oscillatory self-advection produced non-finite values")
    return {
        "velocity": velocity,
        "velocity_gradient": gradient,
        "self_advection": self_advection,
        "divergence": np.asarray(diagnostic["divergence_fd6"], dtype=float),
        "vorticity": np.asarray(diagnostic["vorticity_fd6"], dtype=float),
        "spatial_step": h,
        "gradient_operator": "centered_cartesian_fd6",
        "self_advection_contract": "einsum('...ij,...j->...i', grad_u, u)",
    }


def _directional_fd4_self_advection(
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    directional_step: float,
) -> np.ndarray:
    """Independent frozen-direction FD4 check of ``(u.grad)u``.

    At each sample point, the direction is frozen to ``n=u/|u|`` and the
    directional derivative of the public velocity is evaluated along the line
    ``x+s n``.  Multiplying by ``|u|`` gives ``(u.grad)u`` without using the
    Cartesian gradient implementation under audit.
    """
    h = _validate_step(directional_step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    u0 = np.asarray(velocity_osc(xb, yb, zb, tb), dtype=float)
    if u0.shape != xb.shape + (3,):
        raise RuntimeError("velocity_osc returned unexpected shape")
    speed = np.linalg.norm(u0, axis=-1)
    if np.any(speed <= 1.0e-12):
        raise ValueError("directional audit requires nonzero oscillatory velocity")
    direction = u0 / speed[..., None]

    values: dict[int, np.ndarray] = {}
    for offset in (-2, -1, 1, 2):
        displacement = offset * h * direction
        values[offset] = np.asarray(
            velocity_osc(
                xb + displacement[..., 0],
                yb + displacement[..., 1],
                zb + displacement[..., 2],
                tb,
            ),
            dtype=float,
        )
    directional_derivative = (
        values[-2] - 8.0 * values[-1] + 8.0 * values[1] - values[2]
    ) / (12.0 * h)
    return speed[..., None] * directional_derivative


def _vector_rms(value: np.ndarray) -> float:
    a = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(np.sum(a * a, axis=-1))))


def _relative_rms(reference: np.ndarray, value: np.ndarray) -> float:
    ref = np.asarray(reference, dtype=float)
    val = np.asarray(value, dtype=float)
    return _vector_rms(val - ref) / max(_vector_rms(ref), np.finfo(float).tiny)


def _relative_sampled_max(reference: np.ndarray, value: np.ndarray) -> float:
    ref = np.asarray(reference, dtype=float)
    val = np.asarray(value, dtype=float)
    diff = np.linalg.norm(val - ref, axis=-1)
    denom = max(float(np.max(np.linalg.norm(ref, axis=-1))), np.finfo(float).tiny)
    return float(np.max(diff) / denom)


def _fresh_offgrid_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fresh support-interior points, frozen before Actions execution."""
    radius = np.asarray(
        [0.43, 0.51, 0.60, 0.69, 0.78, 0.87,
         0.47, 0.56, 0.65, 0.74, 0.83, 0.92,
         0.50, 0.59, 0.68, 0.77, 0.86, 0.95],
        dtype=float,
    )
    angle = np.asarray(
        [0.23, 0.87, 1.51, 2.18, 2.79, -2.61,
         -1.97, -1.31, -0.62, 0.41, 1.09, 1.74,
         2.43, -2.88, -2.19, -1.48, -0.79, 0.66],
        dtype=float,
    )
    z = np.asarray(
        [-1.08, -0.73, -0.39, -0.05, 0.31, 0.67,
         1.02, -0.94, -0.58, -0.21, 0.17, 0.54,
         0.89, 1.16, -1.14, -0.81, 0.46, 0.06],
        dtype=float,
    )
    t = np.asarray(([0.42] * 6) + ([0.50] * 6) + ([0.58] * 6), dtype=float)
    return radius * np.cos(angle), radius * np.sin(angle), z, t


def _manufactured_algebra_check() -> dict[str, float]:
    x = np.asarray([0.2, -0.3, 0.5, -0.7], dtype=float)
    y = np.asarray([-0.4, 0.6, 0.1, -0.2], dtype=float)
    z = np.asarray([0.8, -0.5, 0.3, 0.4], dtype=float)
    u = np.stack((x + y, y + z, z + x), axis=-1)
    jacobian = np.broadcast_to(
        np.asarray(
            [
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 1.0],
                [1.0, 0.0, 1.0],
            ],
            dtype=float,
        ),
        u.shape[:-1] + (3, 3),
    )
    actual = _self_advection_from_jacobian(u, jacobian)
    expected = np.stack(
        (u[..., 0] + u[..., 1], u[..., 1] + u[..., 2], u[..., 0] + u[..., 2]),
        axis=-1,
    )
    return {"absolute_max_error": float(np.max(np.abs(actual - expected)))}


def build_receipt() -> dict[str, Any]:
    x, y, z, t = _fresh_offgrid_cloud()
    production = evaluate_oscillatory_self_advection_fd6(
        x, y, z, t, spatial_step=GRADIENT_FD6_STEP
    )
    reference = np.asarray(production["self_advection"], dtype=float)
    directional = [
        _directional_fd4_self_advection(
            x, y, z, t, directional_step=step
        )
        for step in DIRECTIONAL_FD4_STEPS
    ]
    errors = [_relative_rms(reference, value) for value in directional]
    max_errors = [_relative_sampled_max(reference, value) for value in directional]
    refinement = [
        errors[0] / max(errors[1], np.finfo(float).tiny),
        errors[1] / max(errors[2], np.finfo(float).tiny),
    ]

    exterior_x = np.asarray([2.20, -2.24, 0.0, 0.0, 1.70, -1.72], dtype=float)
    exterior_y = np.asarray([0.0, 0.0, 2.22, -2.26, 1.70, -1.72], dtype=float)
    exterior_z = np.asarray([0.0, 0.0, 0.0, 0.0, 2.20, -2.23], dtype=float)
    exterior_t = np.full(exterior_x.shape, 0.50, dtype=float)
    exterior = evaluate_oscillatory_self_advection_fd6(
        exterior_x,
        exterior_y,
        exterior_z,
        exterior_t,
        spatial_step=GRADIENT_FD6_STEP,
    )
    support_exterior_absolute_max = float(
        max(
            np.max(np.abs(np.asarray(exterior["velocity"], dtype=float))),
            np.max(np.abs(np.asarray(exterior["self_advection"], dtype=float))),
        )
    )

    manufactured = _manufactured_algebra_check()
    failed: list[str] = []
    if _vector_rms(reference) < 1.0e-8:
        failed.append("self_advection_nontrivial")
    if refinement[0] < 8.0 or refinement[1] < 8.0:
        failed.append("directional_fd4_refinement")
    if errors[-1] > 1.0e-3:
        failed.append("finest_directional_relative_rms")
    if max_errors[-1] > 2.0e-3:
        failed.append("finest_directional_relative_sampled_max")
    if manufactured["absolute_max_error"] > 1.0e-14:
        failed.append("manufactured_algebra")
    if support_exterior_absolute_max > 1.0e-12:
        failed.append("support_exterior")

    signature = inspect.signature(evaluate_oscillatory_self_advection_fd6)
    forbidden = {
        "residual",
        "defect",
        "mean",
        "stress",
        "pressure",
        "forcing",
        "target",
        "gain",
        "alpha",
        "damping",
        "delta_a",
        "delta_y",
    }
    if forbidden.intersection(signature.parameters):
        failed.append("forbidden_public_inputs")

    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
            "structural_scope": "localized oscillatory complete-curl field only",
        },
        "sample_count": int(x.size),
        "sample_times": [0.42, 0.50, 0.58],
        "gradient_fd6_step": GRADIENT_FD6_STEP,
        "directional_fd4_steps": list(DIRECTIONAL_FD4_STEPS),
        "self_advection_rms": _vector_rms(reference),
        "directional_relative_rms": errors,
        "directional_relative_sampled_max": max_errors,
        "directional_refinement_ratios": refinement,
        "manufactured_algebra": manufactured,
        "support_exterior_absolute_max": support_exterior_absolute_max,
        "failed_guards": failed,
        "truth_boundary": {
            "oscillatory_velocity_candidate_changed": False,
            "amplitude_phase_carrier_support_retuned": False,
            "self_advection_executable": True,
            "self_advection_is_ns_residual": False,
            "agent3_mean_radial_chain_reimplemented": False,
            "agent3_mean_defect_contract_emitted": False,
            "real_global_leading_field_bound": False,
            "full_composite_velocity_available": False,
            "heldout_ns_momentum_residual_assessed": False,
            "residual_reduction_claimed": False,
            "paper_exact": False,
            "pde_validated": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = build_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    if receipt["failed_guards"]:
        raise SystemExit(
            "oscillatory self-advection guards failed: "
            + ", ".join(receipt["failed_guards"])
        )
    print("oscillatory_self_advection=PASS")


if __name__ == "__main__":
    main()
