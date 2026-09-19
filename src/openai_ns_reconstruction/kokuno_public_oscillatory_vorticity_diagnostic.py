"""FD6 spatial-derivative diagnostics for the frozen Kokuno Agent-2 field.

This module does **not** create or retune an oscillatory velocity candidate.  It
numerically differentiates the already-frozen public provider from Agent-2
#561 and exposes a batchable velocity-gradient / vorticity diagnostic for
morphology and downstream debugging.

Provenance boundary
-------------------
Kokuno's corrected 2026-09-09 reconstruction motivates the localized
vector-potential / complete-curl structure used by the upstream field.  The
centered Cartesian FD6 operator and the morphology observables in this module
are repository diagnostics, not Kokuno source formulas and not paper-exact
quantities.  They are intentionally separate from Agent-4 independent
validation and from Agent-3 momentum-defect/correction machinery.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-PUBLIC-VORTICITY-DIAGNOSTIC-045"
SCHEMA = "kokuno-a2-public-oscillatory-vorticity-diagnostic-v1"
PARENT_AGENT2_PR = 643
PARENT_AGENT2_HEAD = "dfea7e7297d07154adbb59c30e74c81410cf0d58"
VELOCITY_PR = 561
VELOCITY_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"

# A centered sixth-order first derivative.  Offsets are ordered explicitly so
# the implementation is easy to audit and is not delegated to numpy.gradient.
_FD6_OFFSETS = (-3, -2, -1, 1, 2, 3)
_FD6_COEFFICIENTS = (-1.0, 9.0, -45.0, 45.0, -9.0, 1.0)
_FD6_DENOMINATOR = 60.0

Provider = Callable[[Any, Any, Any, Any], np.ndarray]


def _broadcast_xyz_t(x: Any, y: Any, z: Any, t: Any) -> tuple[np.ndarray, ...]:
    arrays = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    if not all(np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("x, y, z and t must be finite and broadcastable")
    return tuple(arrays)


def _validate_step(spatial_step: float) -> float:
    h = float(spatial_step)
    if not math.isfinite(h) or not (1.0e-5 <= h <= 5.0e-2):
        raise ValueError("spatial_step must be finite and lie in [1e-5, 5e-2]")
    return h


def _provider_value(
    provider: Provider, x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray
) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise RuntimeError(
            f"velocity provider returned shape {value.shape}, expected {expected}"
        )
    if not np.all(np.isfinite(value)):
        raise RuntimeError("velocity provider returned non-finite values")
    return value


def _fd6_velocity_gradient(
    provider: Provider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> np.ndarray:
    """Return ``[..., component, derivative_axis]`` using centered FD6."""
    h = _validate_step(spatial_step)
    xb, yb, zb, tb = _broadcast_xyz_t(x, y, z, t)
    base = (xb, yb, zb)
    derivatives: list[np.ndarray] = []
    for axis in range(3):
        accum = np.zeros(xb.shape + (3,), dtype=float)
        for offset, coeff in zip(_FD6_OFFSETS, _FD6_COEFFICIENTS, strict=True):
            shifted = [a for a in base]
            shifted[axis] = shifted[axis] + offset * h
            accum += coeff * _provider_value(
                provider, shifted[0], shifted[1], shifted[2], tb
            )
        derivatives.append(accum / (_FD6_DENOMINATOR * h))
    return np.stack(derivatives, axis=-1)


def evaluate_vorticity_osc_fd6(
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float = 0.0045,
) -> dict[str, np.ndarray | float | str]:
    """Evaluate frozen velocity, FD6 gradient, vorticity and divergence.

    The result is a numerical diagnostic only.  In particular, the returned
    vorticity is *not* an additional source formula and must not replace an
    independent validator when assessing the Navier--Stokes residual.
    """
    h = _validate_step(spatial_step)
    xb, yb, zb, tb = _broadcast_xyz_t(x, y, z, t)
    u = _provider_value(velocity_osc, xb, yb, zb, tb)
    grad = _fd6_velocity_gradient(
        velocity_osc, xb, yb, zb, tb, spatial_step=h
    )
    omega = np.stack(
        (
            grad[..., 2, 1] - grad[..., 1, 2],
            grad[..., 0, 2] - grad[..., 2, 0],
            grad[..., 1, 0] - grad[..., 0, 1],
        ),
        axis=-1,
    )
    divergence = np.trace(grad, axis1=-2, axis2=-1)
    return {
        "velocity": u,
        "velocity_gradient_fd6": grad,
        "vorticity_fd6": omega,
        "divergence_fd6": divergence,
        "spatial_step": h,
        "operator": "centered_cartesian_fd6",
    }


def _rms(value: np.ndarray) -> float:
    a = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(a * a)))


def _offgrid_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Fixed before receipt execution.  Radii and |z| stay away from the compact
    # support edges by more than the largest 3h stencil reach.
    radius = np.asarray(
        [0.36, 0.44, 0.53, 0.61, 0.69, 0.78, 0.86, 0.95, 1.03,
         1.11, 0.40, 0.49, 0.58, 0.73, 0.82, 0.91, 1.00, 1.16],
        dtype=float,
    )
    angle = np.asarray(
        [0.17, 0.91, 1.63, 2.31, 2.87, -2.44, -1.78, -1.02, -0.41,
         0.36, 1.14, 1.92, 2.56, -2.91, -2.12, -1.37, -0.72, 0.63],
        dtype=float,
    )
    z = np.asarray(
        [-1.18, -0.83, -0.47, -0.12, 0.24, 0.61, 0.98, 1.27, -1.05,
         -0.69, -0.31, 0.08, 0.43, 0.76, 1.12, -1.31, 0.55, -0.02],
        dtype=float,
    )
    t = np.asarray(([0.34] * 6) + ([0.50] * 6) + ([0.66] * 6), dtype=float)
    return radius * np.cos(angle), radius * np.sin(angle), z, t


def _weighted_morphology(
    x: np.ndarray, z: np.ndarray, vector: np.ndarray
) -> dict[str, float]:
    weight = np.sum(np.asarray(vector, dtype=float) ** 2, axis=-1)
    total = float(np.sum(weight))
    if not math.isfinite(total) or total <= 0.0:
        raise RuntimeError("morphology weight is not positive and finite")
    r = np.abs(np.asarray(x, dtype=float))
    zz = np.asarray(z, dtype=float)
    r_mean = float(np.sum(weight * r) / total)
    z_mean = float(np.sum(weight * zz) / total)
    r_var = float(np.sum(weight * (r - r_mean) ** 2) / total)
    z_var = float(np.sum(weight * (zz - z_mean) ** 2) / total)
    aspect = math.sqrt(z_var / max(r_var, np.finfo(float).tiny))
    central = (r <= 0.70) & (np.abs(zz) <= 0.80)
    return {
        "weight_sum": total,
        "radial_centroid": r_mean,
        "axial_centroid": z_mean,
        "radial_rms_width": math.sqrt(max(r_var, 0.0)),
        "axial_rms_width": math.sqrt(max(z_var, 0.0)),
        "axial_to_radial_aspect": aspect,
        "central_weight_fraction": float(np.sum(weight[central]) / total),
    }


def _plane_morphology(resolution: int) -> dict[str, Any]:
    n = int(resolution)
    if n < 9 or n % 2 == 0:
        raise ValueError("morphology resolution must be an odd integer >=9")
    x1 = np.linspace(-1.25, 1.25, n)
    z1 = np.linspace(-1.70, 1.70, n)
    xx, zz = np.meshgrid(x1, z1, indexing="ij")
    yy = np.zeros_like(xx)
    tt = np.full_like(xx, 0.50)
    diagnostic = evaluate_vorticity_osc_fd6(
        xx, yy, zz, tt, spatial_step=0.0045
    )
    u = np.asarray(diagnostic["velocity"], dtype=float)
    omega = np.asarray(diagnostic["vorticity_fd6"], dtype=float)
    return {
        "resolution": n,
        "plane": "public y=0 meridional plane at t=0.50",
        "x_window": [-1.25, 1.25],
        "z_window": [-1.70, 1.70],
        "vorticity_derivative_step": 0.0045,
        "velocity_rms": _rms(u),
        "velocity_max_abs": float(np.max(np.abs(u))),
        "vorticity_rms": _rms(omega),
        "vorticity_max_abs": float(np.max(np.abs(omega))),
        "velocity_morphology": _weighted_morphology(xx, zz, u),
        "vorticity_morphology": _weighted_morphology(xx, zz, omega),
    }


def build_receipt() -> dict[str, Any]:
    """Build the preregistered three-level derivative/morphology receipt."""
    x, y, z, t = _offgrid_cloud()
    steps = (0.018, 0.009, 0.0045)
    levels = [
        evaluate_vorticity_osc_fd6(x, y, z, t, spatial_step=h) for h in steps
    ]
    omega = [np.asarray(level["vorticity_fd6"], dtype=float) for level in levels]
    divergence = [
        np.asarray(level["divergence_fd6"], dtype=float) for level in levels
    ]
    omega_rms = [_rms(value) for value in omega]
    diff_rms = [_rms(omega[0] - omega[1]), _rms(omega[1] - omega[2])]
    refinement_ratio = diff_rms[0] / max(diff_rms[1], np.finfo(float).tiny)
    finest_relative_divergence_rms = _rms(divergence[-1]) / max(
        omega_rms[-1], np.finfo(float).tiny
    )
    finest_relative_divergence_max = float(np.max(np.abs(divergence[-1]))) / max(
        float(np.max(np.abs(omega[-1]))), np.finfo(float).tiny
    )

    # Exterior points are farther than 3*h_coarse from the registered support,
    # so the complete FD6 stencil must remain exactly in a zero region.
    exterior = evaluate_vorticity_osc_fd6(
        np.asarray([0.0, 1.55, 0.60]),
        np.asarray([0.0, 0.0, 0.0]),
        np.asarray([0.0, 0.0, 2.15]),
        np.asarray([0.50, 0.50, 0.50]),
        spatial_step=steps[0],
    )
    exterior_abs_max = max(
        float(np.max(np.abs(np.asarray(exterior["velocity"], dtype=float)))),
        float(np.max(np.abs(np.asarray(exterior["vorticity_fd6"], dtype=float)))),
        float(np.max(np.abs(np.asarray(exterior["divergence_fd6"], dtype=float)))),
    )

    morphology = [_plane_morphology(n) for n in (17, 25, 33)]
    finite_morphology = all(
        math.isfinite(float(level[key]))
        for level in morphology
        for key in (
            "velocity_rms",
            "velocity_max_abs",
            "vorticity_rms",
            "vorticity_max_abs",
        )
    ) and all(
        math.isfinite(float(level[group][key]))
        for level in morphology
        for group in ("velocity_morphology", "vorticity_morphology")
        for key in (
            "radial_centroid",
            "axial_centroid",
            "radial_rms_width",
            "axial_rms_width",
            "axial_to_radial_aspect",
            "central_weight_fraction",
        )
    )

    # Scientific guards were selected before exact-head CI.  They characterize
    # this A2 self-diagnostic only; they are not final PDE acceptance gates.
    guards = {
        "vorticity_nontrivial_rms_min": 1.0e-8,
        "fd6_vorticity_refinement_ratio_min": 20.0,
        "finest_relative_divergence_rms_max": 5.0e-5,
        "finest_relative_divergence_max_max": 1.0e-4,
        "support_exterior_absolute_max": 1.0e-12,
        "morphology_finite_required": True,
    }
    failed: list[str] = []
    if not (omega_rms[-1] >= guards["vorticity_nontrivial_rms_min"]):
        failed.append("vorticity_nontrivial_rms")
    if not (refinement_ratio >= guards["fd6_vorticity_refinement_ratio_min"]):
        failed.append("fd6_vorticity_refinement_ratio")
    if not (
        finest_relative_divergence_rms
        <= guards["finest_relative_divergence_rms_max"]
    ):
        failed.append("finest_relative_divergence_rms")
    if not (
        finest_relative_divergence_max
        <= guards["finest_relative_divergence_max_max"]
    ):
        failed.append("finest_relative_divergence_max")
    if not (exterior_abs_max <= guards["support_exterior_absolute_max"]):
        failed.append("support_exterior_absolute_max")
    if not finite_morphology:
        failed.append("morphology_finite")

    identity_payload = {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "velocity_pr": VELOCITY_PR,
        "velocity_head": VELOCITY_HEAD,
        "provider": "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc",
        "operator": "centered_cartesian_fd6",
        "steps": list(steps),
        "morphology_resolutions": [17, 25, 33],
    }
    identity_sha256 = hashlib.sha256(
        json.dumps(identity_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    return {
        **identity_payload,
        "identity_sha256": identity_sha256,
        "derivative_resolution_ladder": {
            "steps": list(steps),
            "vorticity_rms": omega_rms,
            "successive_difference_rms": diff_rms,
            "successive_difference_refinement_ratio": float(refinement_ratio),
            "divergence_rms": [_rms(value) for value in divergence],
            "finest_relative_divergence_rms": float(
                finest_relative_divergence_rms
            ),
            "finest_relative_divergence_max": float(
                finest_relative_divergence_max
            ),
        },
        "support_exterior_absolute_max": float(exterior_abs_max),
        "morphology_resolution_ladder": morphology,
        "guards": guards,
        "failed_guards": failed,
        "self_diagnostic_passed": not failed,
        "truth_boundary": {
            "velocity_candidate_changed": False,
            "oscillatory_coefficients_retuned": False,
            "source_formula_changed": False,
            "numerical_fd6_diagnostic_only": True,
            "independent_agent4_validation_replaced": False,
            "agent3_momentum_defect_or_correction_duplicated": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = build_receipt()
    raw = json.dumps(receipt, indent=2, sort_keys=True)
    if args.output is None:
        print(raw)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(raw + "\n", encoding="utf-8")
        print(args.output)
        print("identity_sha256=", receipt["identity_sha256"])
        print("failed_guards=", receipt["failed_guards"])


if __name__ == "__main__":
    _main()
