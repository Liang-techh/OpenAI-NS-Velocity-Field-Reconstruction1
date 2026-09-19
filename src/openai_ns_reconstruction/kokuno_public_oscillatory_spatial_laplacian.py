"""FD6 spatial second-derivative/Laplacian diagnostic for frozen Kokuno A2 velocity.

This module does not create or retune an oscillatory candidate. It consumes the
already-frozen public ``velocity_osc(x,y,z,t)`` and exposes a batch-safe centered
Cartesian sixth-order second-derivative diagnostic. The resulting Laplacian is
useful for later viscous-term assembly, but it is not itself a Navier--Stokes
residual and it does not replace independent Agent-4 validation.

The corrected Kokuno reconstruction is structural provenance for the upstream
localized complete-curl field. The public-z realization and the finite-
difference operator here are repository choices, not paper-exact hidden data.
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

TASK = "KOKUNO-A2-PUBLIC-SPATIAL-LAPLACIAN-047"
SCHEMA = "kokuno-a2-public-spatial-laplacian-v1"
PARENT_AGENT2_PR = 694
PARENT_AGENT2_HEAD = "9c4d1feb514fe44ca16793130f7bd1acfda95e26"
VELOCITY_PR = 561
VELOCITY_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"

# Seven-point centered sixth-order second derivative.
_FD6_D2_OFFSETS = (-3, -2, -1, 0, 1, 2, 3)
_FD6_D2_COEFFICIENTS = (2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0)
_FD6_D2_DENOMINATOR = 180.0

Provider = Callable[[Any, Any, Any, Any], np.ndarray]


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


def _validate_spatial_step(step: float) -> float:
    h = float(step)
    if not math.isfinite(h) or not (1.0e-5 <= h <= 5.0e-2):
        raise ValueError("spatial_step must be finite and lie in [1e-5, 5e-2]")
    return h


def _provider_value(
    provider: Provider,
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


def _fd6_second_derivatives(
    provider: Provider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> np.ndarray:
    """Return ``[..., component, axis]`` pure second derivatives."""
    h = _validate_spatial_step(spatial_step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = (xb, yb, zb)
    derivatives: list[np.ndarray] = []
    for axis in range(3):
        accum = np.zeros(xb.shape + (3,), dtype=float)
        for offset, coefficient in zip(
            _FD6_D2_OFFSETS, _FD6_D2_COEFFICIENTS, strict=True
        ):
            shifted = [a for a in base]
            shifted[axis] = shifted[axis] + offset * h
            accum += coefficient * _provider_value(
                provider, shifted[0], shifted[1], shifted[2], tb
            )
        derivatives.append(accum / (_FD6_D2_DENOMINATOR * h * h))
    return np.stack(derivatives, axis=-1)


def evaluate_velocity_laplacian_fd6(
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float = 0.0045,
) -> dict[str, np.ndarray | float | str]:
    """Evaluate frozen velocity, pure second derivatives and Cartesian Laplacian."""
    h = _validate_spatial_step(spatial_step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    velocity = _provider_value(velocity_osc, xb, yb, zb, tb)
    second = _fd6_second_derivatives(
        velocity_osc, xb, yb, zb, tb, spatial_step=h
    )
    laplacian = np.sum(second, axis=-1)
    return {
        "velocity": velocity,
        "second_derivatives": second,
        "laplacian": laplacian,
        "spatial_step": h,
        "spatial_operator": "centered_cartesian_fd6_second_derivative",
    }


def _vector_rms(value: np.ndarray) -> float:
    a = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(np.sum(a * a, axis=-1))))


def _offgrid_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fresh support-interior cloud, distinct from #669/#694 receipts."""
    radius = np.asarray(
        [0.41, 0.49, 0.58, 0.67, 0.76, 0.85, 0.94, 1.03, 1.12,
         0.45, 0.54, 0.63, 0.72, 0.81, 0.90, 0.99, 1.08, 1.17],
        dtype=float,
    )
    angle = np.asarray(
        [0.17, 0.91, 1.58, 2.22, 2.87, -2.61, -1.92, -1.19, -0.43,
         0.36, 1.14, 1.88, 2.56, -2.94, -2.17, -1.42, -0.69, 0.68],
        dtype=float,
    )
    z = np.asarray(
        [-1.19, -0.83, -0.51, -0.16, 0.21, 0.57, 0.92, 1.20, -1.05,
         -0.70, -0.33, 0.08, 0.43, 0.77, 1.12, -1.30, 0.63, -0.01],
        dtype=float,
    )
    t = np.asarray(([0.36] * 6) + ([0.50] * 6) + ([0.64] * 6), dtype=float)
    return radius * np.cos(angle), radius * np.sin(angle), z, t


def _manufactured_provider(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    del tb
    return np.stack(
        (
            xb**6 + 2.0 * yb**4 + 3.0 * zb**2,
            0.5 * xb**4 - yb**6 + zb**4,
            2.0 * xb**2 + 3.0 * yb**2 + 4.0 * zb**6,
        ),
        axis=-1,
    )


def _manufactured_laplacian(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    del tb
    return np.stack(
        (
            30.0 * xb**4 + 24.0 * yb**2 + 6.0,
            6.0 * xb**2 - 30.0 * yb**4 + 12.0 * zb**2,
            10.0 + 120.0 * zb**4,
        ),
        axis=-1,
    )


def build_receipt() -> dict[str, Any]:
    x, y, z, t = _offgrid_cloud()
    steps = (0.018, 0.009, 0.0045)
    levels = [
        evaluate_velocity_laplacian_fd6(x, y, z, t, spatial_step=h)
        for h in steps
    ]
    lap = [np.asarray(level["laplacian"], dtype=float) for level in levels]
    diff_rms = [_vector_rms(lap[0] - lap[1]), _vector_rms(lap[1] - lap[2])]
    refinement = diff_rms[0] / max(diff_rms[1], np.finfo(float).tiny)
    finest_rms = _vector_rms(lap[-1])
    finest_max = float(np.max(np.linalg.norm(lap[-1], axis=-1)))

    mx = np.asarray([-0.41, 0.23, 0.52, -0.31], dtype=float)
    my = np.asarray([0.19, -0.37, 0.28, 0.44], dtype=float)
    mz = np.asarray([0.33, -0.21, -0.46, 0.17], dtype=float)
    mt = np.asarray([0.36, 0.44, 0.56, 0.64], dtype=float)
    manufactured = np.sum(
        _fd6_second_derivatives(
            _manufactured_provider, mx, my, mz, mt, spatial_step=0.013
        ),
        axis=-1,
    )
    manufactured_exact = _manufactured_laplacian(mx, my, mz, mt)
    manufactured_abs_max = float(
        np.max(np.abs(manufactured - manufactured_exact))
    )
    manufactured_scale = max(float(np.max(np.abs(manufactured_exact))), 1.0)
    manufactured_rel_max = manufactured_abs_max / manufactured_scale

    exterior = evaluate_velocity_laplacian_fd6(
        np.asarray([0.0, 1.75, 0.65]),
        np.asarray([0.0, 0.0, 0.0]),
        np.asarray([0.0, 0.0, 2.30]),
        np.asarray([0.50, 0.50, 0.50]),
        spatial_step=steps[0],
    )
    exterior_abs_max = max(
        float(np.max(np.abs(np.asarray(exterior[key], dtype=float))))
        for key in ("velocity", "second_derivatives", "laplacian")
    )

    guards = {
        "laplacian_nontrivial_rms_min": 1.0e-8,
        "fd6_successive_difference_refinement_ratio_min": 16.0,
        "manufactured_relative_max_error_max": 5.0e-8,
        "support_exterior_absolute_max": 1.0e-12,
    }
    failed: list[str] = []
    if finest_rms < guards["laplacian_nontrivial_rms_min"]:
        failed.append("laplacian_nontrivial_rms")
    if refinement < guards["fd6_successive_difference_refinement_ratio_min"]:
        failed.append("fd6_successive_difference_refinement_ratio")
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
        "spatial_operator": "centered_cartesian_fd6_second_derivative",
        "spatial_steps": list(steps),
    }
    identity_sha256 = hashlib.sha256(
        json.dumps(identity_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    return {
        **identity_payload,
        "identity_sha256": identity_sha256,
        "sample_count": int(t.size),
        "verification_times": sorted({float(v) for v in t}),
        "laplacian_ladder": {
            "steps": list(steps),
            "rms": [_vector_rms(value) for value in lap],
            "sampled_max": [
                float(np.max(np.linalg.norm(value, axis=-1))) for value in lap
            ],
            "successive_difference_rms": diff_rms,
            "successive_difference_refinement_ratio": refinement,
        },
        "finest_laplacian_rms": finest_rms,
        "finest_laplacian_sampled_max": finest_max,
        "manufactured_polynomial": {
            "step": 0.013,
            "absolute_max_error": manufactured_abs_max,
            "relative_max_error": manufactured_rel_max,
            "degree_bound": 6,
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
            "fd6_second_derivative_operator": "repository_numerical_diagnostic",
        },
        "truth_boundary": {
            "velocity_candidate_changed": False,
            "oscillatory_coefficients_retuned": False,
            "source_formula_changed": False,
            "numerical_spatial_laplacian_diagnostic_only": True,
            "viscous_term_ready_for_future_composite_diagnostics": True,
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
