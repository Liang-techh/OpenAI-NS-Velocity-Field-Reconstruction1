"""Independent vector-identity cross-check for the frozen Kokuno A2 Laplacian.

This module does not create or retune an oscillatory velocity candidate.  It
uses the already-frozen public ``velocity_osc(x,y,z,t)`` and compares two
numerically distinct routes to the Cartesian vector Laplacian,

    Delta u = grad(div u) - curl(curl u).

The right-hand identity is reconstructed with nested centered FD4 first
derivatives, whereas the comparison Laplacian comes from Agent-2 #708's
centered FD6 pure-second-derivative diagnostic.  This is an internal numerical
consistency check for future viscous-term assembly, not a Navier--Stokes
residual and not a replacement for independent Agent-4 validation.

The corrected Kokuno 2026-09-09 reconstruction supplies structural provenance
for the upstream localized complete-curl field.  The public-z realization and
both finite-difference routes used here are repository diagnostics, not hidden
paper-exact data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_oscillatory_spatial_laplacian import (
    evaluate_velocity_laplacian_fd6,
)
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-LAPLACIAN-VECTOR-IDENTITY-048"
SCHEMA = "kokuno-a2-laplacian-vector-identity-v1"
PARENT_AGENT2_PR = 708
PARENT_AGENT2_HEAD = "9f8bae37c4d3b2559bee6db050655702e4c058e9"
VELOCITY_PR = 561
VELOCITY_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"

# Centered fourth-order first derivative:
# [f(-2h)-8f(-h)+8f(h)-f(2h)]/(12h).
_FD4_D1_OFFSETS = (-2, -1, 1, 2)
_FD4_D1_COEFFICIENTS = (1.0, -8.0, 8.0, -1.0)
_FD4_D1_DENOMINATOR = 12.0

VectorProvider = Callable[[Any, Any, Any, Any], np.ndarray]
ScalarProvider = Callable[[Any, Any, Any, Any], np.ndarray]


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


def _validate_step(step: float, *, name: str) -> float:
    h = float(step)
    if not math.isfinite(h) or not (1.0e-5 <= h <= 5.0e-2):
        raise ValueError(f"{name} must be finite and lie in [1e-5, 5e-2]")
    return h


def _vector_value(
    provider: VectorProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise RuntimeError(
            f"vector provider returned shape {value.shape}, expected {expected}"
        )
    if not np.all(np.isfinite(value)):
        raise RuntimeError("vector provider returned non-finite values")
    return value


def _scalar_value(
    provider: ScalarProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    if value.shape != x.shape:
        raise RuntimeError(
            f"scalar provider returned shape {value.shape}, expected {x.shape}"
        )
    if not np.all(np.isfinite(value)):
        raise RuntimeError("scalar provider returned non-finite values")
    return value


def _fd4_vector_gradient(
    provider: VectorProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> np.ndarray:
    """Return ``[..., component, derivative_axis]`` with centered FD4."""
    h = _validate_step(spatial_step, name="spatial_step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = (xb, yb, zb)
    derivatives: list[np.ndarray] = []
    for axis in range(3):
        accum = np.zeros(xb.shape + (3,), dtype=float)
        for offset, coefficient in zip(
            _FD4_D1_OFFSETS, _FD4_D1_COEFFICIENTS, strict=True
        ):
            shifted = [a for a in base]
            shifted[axis] = shifted[axis] + offset * h
            accum += coefficient * _vector_value(
                provider, shifted[0], shifted[1], shifted[2], tb
            )
        derivatives.append(accum / (_FD4_D1_DENOMINATOR * h))
    return np.stack(derivatives, axis=-1)


def _fd4_scalar_gradient(
    provider: ScalarProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> np.ndarray:
    """Return ``[..., derivative_axis]`` with centered FD4."""
    h = _validate_step(spatial_step, name="spatial_step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = (xb, yb, zb)
    derivatives: list[np.ndarray] = []
    for axis in range(3):
        accum = np.zeros(xb.shape, dtype=float)
        for offset, coefficient in zip(
            _FD4_D1_OFFSETS, _FD4_D1_COEFFICIENTS, strict=True
        ):
            shifted = [a for a in base]
            shifted[axis] = shifted[axis] + offset * h
            accum += coefficient * _scalar_value(
                provider, shifted[0], shifted[1], shifted[2], tb
            )
        derivatives.append(accum / (_FD4_D1_DENOMINATOR * h))
    return np.stack(derivatives, axis=-1)


def _curl_from_gradient(gradient: np.ndarray) -> np.ndarray:
    grad = np.asarray(gradient, dtype=float)
    if grad.shape[-2:] != (3, 3):
        raise ValueError("gradient must end in shape (3,3)")
    return np.stack(
        (
            grad[..., 2, 1] - grad[..., 1, 2],
            grad[..., 0, 2] - grad[..., 2, 0],
            grad[..., 1, 0] - grad[..., 0, 1],
        ),
        axis=-1,
    )


def _fd4_curl(
    provider: VectorProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> np.ndarray:
    return _curl_from_gradient(
        _fd4_vector_gradient(provider, x, y, z, t, spatial_step=spatial_step)
    )


def _fd4_divergence(
    provider: VectorProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> np.ndarray:
    grad = _fd4_vector_gradient(
        provider, x, y, z, t, spatial_step=spatial_step
    )
    return np.trace(grad, axis1=-2, axis2=-1)


def _fd4_laplacian_vector_identity(
    provider: VectorProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``grad(div u)``, ``curl(curl u)`` and their difference.

    Both outer operators use centered FD4.  The inner divergence/curl are
    recomputed at every shifted outer-stencil point, deliberately avoiding the
    pure-second-derivative implementation from #708.
    """
    h = _validate_step(spatial_step, name="spatial_step")

    def divergence_provider(xp: Any, yp: Any, zp: Any, tp: Any) -> np.ndarray:
        return _fd4_divergence(
            provider, xp, yp, zp, tp, spatial_step=h
        )

    def curl_provider(xp: Any, yp: Any, zp: Any, tp: Any) -> np.ndarray:
        return _fd4_curl(provider, xp, yp, zp, tp, spatial_step=h)

    grad_div = _fd4_scalar_gradient(
        divergence_provider, x, y, z, t, spatial_step=h
    )
    curl_curl = _fd4_curl(
        curl_provider, x, y, z, t, spatial_step=h
    )
    identity_laplacian = grad_div - curl_curl
    return grad_div, curl_curl, identity_laplacian


def evaluate_laplacian_vector_identity(
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    first_derivative_step: float = 0.0045,
    laplacian_step: float = 0.0045,
) -> dict[str, np.ndarray | float | str]:
    """Cross-check the frozen velocity Laplacian through two numerical routes."""
    h1 = _validate_step(first_derivative_step, name="first_derivative_step")
    h2 = _validate_step(laplacian_step, name="laplacian_step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    velocity = _vector_value(velocity_osc, xb, yb, zb, tb)
    grad_div, curl_curl, identity_laplacian = _fd4_laplacian_vector_identity(
        velocity_osc, xb, yb, zb, tb, spatial_step=h1
    )
    fd6 = evaluate_velocity_laplacian_fd6(
        xb, yb, zb, tb, spatial_step=h2
    )
    laplacian_fd6 = np.asarray(fd6["laplacian"], dtype=float)
    return {
        "velocity": velocity,
        "grad_div_fd4": grad_div,
        "curl_curl_fd4": curl_curl,
        "identity_laplacian_fd4": identity_laplacian,
        "minus_curl_curl_fd4": -curl_curl,
        "laplacian_fd6": laplacian_fd6,
        "identity_minus_fd6": identity_laplacian - laplacian_fd6,
        "first_derivative_step": h1,
        "laplacian_step": h2,
        "identity_operator": "nested_centered_cartesian_fd4_first_derivatives",
        "comparison_operator": "centered_cartesian_fd6_pure_second_derivatives",
    }


def _vector_rms(value: np.ndarray) -> float:
    a = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(np.sum(a * a, axis=-1))))


def _offgrid_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fresh support-interior cloud, distinct from #669/#694/#708 receipts."""
    radius = np.asarray(
        [0.43, 0.51, 0.60, 0.69, 0.78, 0.87, 0.96, 1.05, 1.14,
         0.47, 0.56, 0.65, 0.74, 0.83, 0.92, 1.01, 1.10, 1.19],
        dtype=float,
    )
    angle = np.asarray(
        [0.29, 0.98, 1.67, 2.36, 3.02, -2.48, -1.71, -0.97, -0.26,
         0.52, 1.25, 1.95, 2.69, -2.80, -2.04, -1.30, -0.57, 0.79],
        dtype=float,
    )
    z = np.asarray(
        [-1.16, -0.79, -0.44, -0.10, 0.27, 0.60, 0.95, 1.18, -1.02,
         -0.66, -0.29, 0.12, 0.46, 0.81, 1.09, -1.25, 0.69, 0.02],
        dtype=float,
    )
    t = np.asarray(([0.37] * 6) + ([0.50] * 6) + ([0.63] * 6), dtype=float)
    return radius * np.cos(angle), radius * np.sin(angle), z, t


def _manufactured_provider(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    del tb
    return np.stack(
        (
            xb**4 + xb * yb * zb + 2.0 * yb**2,
            yb**4 + xb**2 * zb + 3.0 * xb * zb**2,
            zb**4 + xb * yb**2 + 2.0 * xb**2,
        ),
        axis=-1,
    )


def _manufactured_laplacian(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    del tb
    return np.stack(
        (
            12.0 * xb**2 + 4.0,
            12.0 * yb**2 + 2.0 * zb + 6.0 * xb,
            12.0 * zb**2 + 2.0 * xb + 4.0,
        ),
        axis=-1,
    )


def build_receipt() -> dict[str, Any]:
    x, y, z, t = _offgrid_cloud()
    steps = (0.018, 0.009, 0.0045)
    comparison_step = 0.0045
    levels = [
        evaluate_laplacian_vector_identity(
            x,
            y,
            z,
            t,
            first_derivative_step=h,
            laplacian_step=comparison_step,
        )
        for h in steps
    ]
    identity = [
        np.asarray(level["identity_laplacian_fd4"], dtype=float)
        for level in levels
    ]
    diff_rms = [
        _vector_rms(identity[0] - identity[1]),
        _vector_rms(identity[1] - identity[2]),
    ]
    refinement = diff_rms[0] / max(diff_rms[1], np.finfo(float).tiny)
    finest = levels[-1]
    finest_identity = np.asarray(finest["identity_laplacian_fd4"], dtype=float)
    finest_fd6 = np.asarray(finest["laplacian_fd6"], dtype=float)
    finest_difference = finest_identity - finest_fd6
    fd6_rms = _vector_rms(finest_fd6)
    fd6_sampled_max = float(np.max(np.linalg.norm(finest_fd6, axis=-1)))
    identity_rms = _vector_rms(finest_identity)
    relative_rms = _vector_rms(finest_difference) / max(
        fd6_rms, np.finfo(float).tiny
    )
    relative_sampled_max = float(
        np.max(np.linalg.norm(finest_difference, axis=-1))
    ) / max(fd6_sampled_max, np.finfo(float).tiny)
    grad_div_rms = _vector_rms(np.asarray(finest["grad_div_fd4"], dtype=float))
    grad_div_relative_rms = grad_div_rms / max(
        identity_rms, np.finfo(float).tiny
    )

    mx = np.asarray([-0.41, 0.23, 0.52, -0.31], dtype=float)
    my = np.asarray([0.19, -0.37, 0.28, 0.44], dtype=float)
    mz = np.asarray([0.33, -0.21, -0.46, 0.17], dtype=float)
    mt = np.asarray([0.37, 0.45, 0.55, 0.63], dtype=float)
    _, _, manufactured = _fd4_laplacian_vector_identity(
        _manufactured_provider, mx, my, mz, mt, spatial_step=0.013
    )
    manufactured_exact = _manufactured_laplacian(mx, my, mz, mt)
    manufactured_abs_max = float(np.max(np.abs(manufactured - manufactured_exact)))
    manufactured_scale = max(float(np.max(np.abs(manufactured_exact))), 1.0)
    manufactured_rel_max = manufactured_abs_max / manufactured_scale

    # Each nested FD4 route reaches at most 4*h from the center.  These points
    # stay safely outside the registered compact support for the coarse h.
    exterior = evaluate_laplacian_vector_identity(
        np.asarray([0.0, 1.80, 0.65]),
        np.asarray([0.0, 0.0, 0.0]),
        np.asarray([0.0, 0.0, 2.35]),
        np.asarray([0.50, 0.50, 0.50]),
        first_derivative_step=steps[0],
        laplacian_step=steps[0],
    )
    exterior_abs_max = max(
        float(np.max(np.abs(np.asarray(exterior[key], dtype=float))))
        for key in (
            "velocity",
            "grad_div_fd4",
            "curl_curl_fd4",
            "identity_laplacian_fd4",
            "laplacian_fd6",
            "identity_minus_fd6",
        )
    )

    guards = {
        "identity_laplacian_nontrivial_rms_min": 1.0e-8,
        "nested_fd4_successive_difference_refinement_ratio_min": 8.0,
        "finest_identity_vs_fd6_relative_rms_max": 2.0e-2,
        "finest_identity_vs_fd6_relative_sampled_max_max": 5.0e-2,
        "manufactured_relative_max_error_max": 1.0e-8,
        "support_exterior_absolute_max": 1.0e-12,
    }
    failed: list[str] = []
    if identity_rms < guards["identity_laplacian_nontrivial_rms_min"]:
        failed.append("identity_laplacian_nontrivial_rms")
    if refinement < guards["nested_fd4_successive_difference_refinement_ratio_min"]:
        failed.append("nested_fd4_successive_difference_refinement_ratio")
    if relative_rms > guards["finest_identity_vs_fd6_relative_rms_max"]:
        failed.append("finest_identity_vs_fd6_relative_rms")
    if (
        relative_sampled_max
        > guards["finest_identity_vs_fd6_relative_sampled_max_max"]
    ):
        failed.append("finest_identity_vs_fd6_relative_sampled_max")
    if manufactured_rel_max > guards["manufactured_relative_max_error_max"]:
        failed.append("manufactured_relative_max_error")
    if exterior_abs_max > guards["support_exterior_absolute_max"]:
        failed.append("support_exterior_absolute_max")

    identity_payload = {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "velocity_pr": VELOCITY_PR,
        "velocity_head": VELOCITY_HEAD,
        "identity_operator": "nested_centered_cartesian_fd4_first_derivatives",
        "comparison_operator": "centered_cartesian_fd6_pure_second_derivatives",
        "first_derivative_steps": list(steps),
        "comparison_laplacian_step": comparison_step,
    }
    identity_sha256 = hashlib.sha256(
        json.dumps(identity_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    return {
        **identity_payload,
        "identity_sha256": identity_sha256,
        "sample_count": int(t.size),
        "verification_times": sorted({float(v) for v in t}),
        "nested_fd4_ladder": {
            "steps": list(steps),
            "identity_laplacian_rms": [_vector_rms(value) for value in identity],
            "successive_difference_rms": diff_rms,
            "successive_difference_refinement_ratio": refinement,
        },
        "finest_identity_laplacian_rms": identity_rms,
        "finest_fd6_laplacian_rms": fd6_rms,
        "finest_identity_vs_fd6_relative_rms": relative_rms,
        "finest_identity_vs_fd6_relative_sampled_max": relative_sampled_max,
        "finest_grad_div_relative_rms": grad_div_relative_rms,
        "manufactured_polynomial": {
            "step": 0.013,
            "absolute_max_error": manufactured_abs_max,
            "relative_max_error": manufactured_rel_max,
            "degree_bound": 4,
        },
        "support_exterior_absolute_max": exterior_abs_max,
        "guards": guards,
        "failed_guards": failed,
        "self_diagnostic_passed": failed == [],
        "provenance": {
            "source_structure": (
                "inherits corrected-Kokuno localized oscillatory complete-curl "
                "structure from the frozen upstream public velocity"
            ),
            "public_z_pullback": "repository_public_coordinate_realization",
            "nested_fd4_operator": "repository_numerical_diagnostic",
            "fd6_laplacian_operator": "repository_numerical_diagnostic_from_pr_708",
        },
        "truth_boundary": {
            "velocity_candidate_changed": False,
            "oscillatory_coefficients_retuned": False,
            "source_formula_changed": False,
            "laplacian_vector_identity_self_diagnostic_only": True,
            "viscous_term_formula_or_pressure_fit_changed": False,
            "independent_agent4_validation_replaced": False,
            "independent_agent4_vector_potential_audit_required": True,
            "agent3_momentum_defect_or_correction_duplicated": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = build_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
