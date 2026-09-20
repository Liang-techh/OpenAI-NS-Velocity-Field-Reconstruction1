"""Independent vector-identity audit for the A2 signed-amplitude correction Laplacian.

This module cross-checks the already-materialized divergence-free correction
``delta u = curl(delta A)`` through two numerically distinct Cartesian routes:

    Delta(delta u) = grad(div(delta u)) - curl(curl(delta u)).

The vector-identity side uses nested centered FD4 first derivatives, while the
comparison side uses Agent-2 #752's centered FD6 pure-second-derivative
Laplacian. The corrected Kokuno reconstruction supplies structural provenance
for the upstream localized complete-curl field; the compact radial spline,
public-z pullback, autonomous correction time lift, and both finite-difference
operators are repository realizations rather than paper-exact hidden data.

No Agent-3 defect/mean/radial inverse, pressure, forcing, residual, target,
gain, damping choice, or scientific threshold enters the public evaluator.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_signed_amplitude_complete_curl import (
    SignedRadialAmplitudeCorrection,
    correction_velocity,
    profile_from_delta_a,
)
from .kokuno_public_signed_amplitude_spatial_laplacian import (
    evaluate_correction_velocity_laplacian_fd6,
)

TASK = "KOKUNO-A2-SIGNED-AMPLITUDE-LAPLACIAN-VECTOR-IDENTITY-053"
SCHEMA = "kokuno-a2-signed-amplitude-laplacian-vector-identity-v1"
PARENT_AGENT2_PR = 752
PARENT_AGENT2_HEAD = "5ed44296743f0d4166f6ffb80838fe8adce7cfbf"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
FD4_STEPS = (0.004, 0.002, 0.001)
FD6_COMPARISON_STEP = 0.001

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
    if not all(np.all(np.isfinite(array)) for array in arrays):
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
    h = _validate_step(spatial_step, name="first_derivative_step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = (xb, yb, zb)
    derivatives: list[np.ndarray] = []
    for axis in range(3):
        accum = np.zeros(xb.shape + (3,), dtype=float)
        for offset, coefficient in zip(
            _FD4_D1_OFFSETS, _FD4_D1_COEFFICIENTS, strict=True
        ):
            shifted = [array for array in base]
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
    h = _validate_step(spatial_step, name="first_derivative_step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = (xb, yb, zb)
    derivatives: list[np.ndarray] = []
    for axis in range(3):
        accum = np.zeros(xb.shape, dtype=float)
        for offset, coefficient in zip(
            _FD4_D1_OFFSETS, _FD4_D1_COEFFICIENTS, strict=True
        ):
            shifted = [array for array in base]
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
    gradient = _fd4_vector_gradient(
        provider, x, y, z, t, spatial_step=spatial_step
    )
    return np.trace(gradient, axis1=-2, axis2=-1)


def _fd4_laplacian_vector_identity(
    provider: VectorProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    h = _validate_step(spatial_step, name="first_derivative_step")

    def divergence_provider(xp: Any, yp: Any, zp: Any, tp: Any) -> np.ndarray:
        return _fd4_divergence(provider, xp, yp, zp, tp, spatial_step=h)

    def curl_provider(xp: Any, yp: Any, zp: Any, tp: Any) -> np.ndarray:
        return _fd4_curl(provider, xp, yp, zp, tp, spatial_step=h)

    grad_div = _fd4_scalar_gradient(
        divergence_provider, x, y, z, t, spatial_step=h
    )
    curl_curl = _fd4_curl(curl_provider, x, y, z, t, spatial_step=h)
    return grad_div, curl_curl, grad_div - curl_curl


def evaluate_correction_laplacian_vector_identity(
    correction: SignedRadialAmplitudeCorrection,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    first_derivative_step: float = 0.001,
    laplacian_step: float = 0.001,
) -> dict[str, np.ndarray | float | str]:
    """Cross-check ``Delta(delta u)`` through nested FD4 and pure-second FD6."""
    h1 = _validate_step(first_derivative_step, name="first_derivative_step")
    h2 = _validate_step(laplacian_step, name="laplacian_step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)

    def provider(xx: Any, yy: Any, zz: Any, tt: Any) -> np.ndarray:
        return correction_velocity(correction, xx, yy, zz, tt)

    velocity = _vector_value(provider, xb, yb, zb, tb)
    grad_div, curl_curl, identity_laplacian = _fd4_laplacian_vector_identity(
        provider, xb, yb, zb, tb, spatial_step=h1
    )
    fd6 = evaluate_correction_velocity_laplacian_fd6(
        correction, xb, yb, zb, tb, spatial_step=h2
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


def _vector_rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def _verification_profile() -> SignedRadialAmplitudeCorrection:
    radii = np.linspace(0.33, 1.17, 23)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(math.pi * s) ** 8
    delta_a = np.stack(
        (
            0.0087 * bump * (1.0 + 0.08 * np.cos(2.0 * math.pi * s)),
            -0.0061 * bump * (1.0 - 0.06 * np.sin(2.0 * math.pi * s)),
        ),
        axis=-1,
    )
    delta_a[[0, -1], :] = 0.0
    return profile_from_delta_a(
        radii,
        delta_a,
        reference_time=0.50,
        producer_kind="analytic-regression-not-agent3-real-candidate",
        provenance=(
            "fresh compact signed radial mechanics profile for correction vector-identity "
            "Laplacian cross-check only"
        ),
        source_mean_amplitude_differential_certified=False,
    )


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r0 = 0.33
    dr = (1.17 - 0.33) / 22.0
    cell_ids = (2, 5, 8, 11, 14, 17)
    radii = np.asarray([r0 + (index + 0.5) * dr for index in cell_ids], dtype=float)
    z_block = np.asarray((-0.57, -0.34, -0.11, 0.15, 0.38, 0.58), dtype=float)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    count = 0
    for time in (0.40, 0.50, 0.60):
        for radius, z_value in zip(radii, z_block, strict=True):
            angle = 0.271 + (count + 1) * golden
            rows.append(
                (
                    radius * math.cos(angle),
                    radius * math.sin(angle),
                    float(z_value),
                    time,
                )
            )
            count += 1
    array = np.asarray(rows, dtype=float)
    return array[:, 0], array[:, 1], array[:, 2], array[:, 3]


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


def truth_boundary() -> dict[str, bool]:
    return {
        "correction_laplacian_vector_identity_executable": True,
        "fd6_parent_laplacian_cross_checked_by_independent_route": True,
        "vector_potential_first_complete_curl_reused": True,
        "candidate_amplitude_phase_support_retuned": False,
        "agent3_mean_radial_chain_reimplemented": False,
        "source_agent2_complete_curl_certified_for_full_candidate": False,
        "independent_agent4_correction_vector_potential_audit_required": True,
        "real_agent3_delta_a_bound": False,
        "real_full_candidate_correction_cycle_run": False,
        "heldout_ns_momentum_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def verification_receipt() -> dict[str, Any]:
    correction = _verification_profile()
    x, y, z, t = _verification_cloud()
    levels = [
        evaluate_correction_laplacian_vector_identity(
            correction,
            x,
            y,
            z,
            t,
            first_derivative_step=step,
            laplacian_step=FD6_COMPARISON_STEP,
        )
        for step in FD4_STEPS
    ]
    identity = [
        np.asarray(level["identity_laplacian_fd4"], dtype=float)
        for level in levels
    ]
    successive = [
        _vector_rms(identity[0] - identity[1]),
        _vector_rms(identity[1] - identity[2]),
    ]
    refinement = successive[0] / max(successive[1], np.finfo(float).tiny)

    finest = levels[-1]
    identity_finest = np.asarray(finest["identity_laplacian_fd4"], dtype=float)
    fd6_finest = np.asarray(finest["laplacian_fd6"], dtype=float)
    difference = identity_finest - fd6_finest
    fd6_rms = _vector_rms(fd6_finest)
    fd6_max = float(np.max(np.linalg.norm(fd6_finest, axis=-1)))
    relative_rms = _vector_rms(difference) / max(fd6_rms, np.finfo(float).tiny)
    relative_max = float(np.max(np.linalg.norm(difference, axis=-1))) / max(
        fd6_max, np.finfo(float).tiny
    )
    grad_div_ratio = _vector_rms(np.asarray(finest["grad_div_fd4"])) / max(
        _vector_rms(identity_finest), np.finfo(float).tiny
    )

    mx = np.asarray((0.37, -0.41, 0.52, -0.28), dtype=float)
    my = np.asarray((-0.26, 0.33, 0.19, -0.47), dtype=float)
    mz = np.asarray((0.21, -0.18, 0.31, -0.36), dtype=float)
    mt = np.full(mx.shape, 0.50)
    manufactured_numeric = _fd4_laplacian_vector_identity(
        _manufactured_provider,
        mx,
        my,
        mz,
        mt,
        spatial_step=FD4_STEPS[-1],
    )[2]
    manufactured_exact = _manufactured_laplacian(mx, my, mz, mt)
    manufactured_error = manufactured_numeric - manufactured_exact
    manufactured_scale = float(np.max(np.linalg.norm(manufactured_exact, axis=-1)))
    manufactured_relative_max = float(
        np.max(np.linalg.norm(manufactured_error, axis=-1))
    ) / max(manufactured_scale, np.finfo(float).tiny)

    exterior_x = np.asarray((1.42, -1.47, 1.51, -1.56), dtype=float)
    exterior_y = np.asarray((0.17, -0.21, -0.16, 0.19), dtype=float)
    exterior_z = np.asarray((-0.41, -0.12, 0.18, 0.43), dtype=float)
    exterior_t = np.asarray((0.40, 0.50, 0.50, 0.60), dtype=float)
    exterior = evaluate_correction_laplacian_vector_identity(
        correction,
        exterior_x,
        exterior_y,
        exterior_z,
        exterior_t,
        first_derivative_step=FD4_STEPS[0],
        laplacian_step=FD6_COMPARISON_STEP,
    )
    support_exterior_absolute_max = max(
        float(np.max(np.abs(np.asarray(exterior[key], dtype=float))))
        for key in (
            "velocity",
            "grad_div_fd4",
            "curl_curl_fd4",
            "identity_laplacian_fd4",
            "laplacian_fd6",
        )
    )

    guards = {
        "identity_laplacian_nontrivial_rms": _vector_rms(identity_finest) >= 1.0e-8,
        "successive_difference_refinement": refinement >= 8.0,
        "identity_vs_fd6_relative_rms": relative_rms <= 2.0e-2,
        "identity_vs_fd6_relative_max": relative_max <= 5.0e-2,
        "manufactured_relative_max_error": manufactured_relative_max <= 1.0e-8,
        "support_exterior": support_exterior_absolute_max <= 1.0e-12,
    }
    failed = [name for name, passed in guards.items() if not passed]
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "source_corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
        "source_corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        "sample_count": int(x.size),
        "fd4_steps": list(FD4_STEPS),
        "fd6_comparison_step": FD6_COMPARISON_STEP,
        "successive_difference_rms": successive,
        "successive_difference_refinement_ratio": refinement,
        "finest_identity_laplacian_rms": _vector_rms(identity_finest),
        "finest_fd6_laplacian_rms": fd6_rms,
        "identity_vs_fd6_relative_rms": relative_rms,
        "identity_vs_fd6_relative_sampled_max": relative_max,
        "grad_div_over_identity_laplacian_rms": grad_div_ratio,
        "manufactured_polynomial": {
            "relative_max_error": manufactured_relative_max,
            "operator": "nested_centered_cartesian_fd4_first_derivatives",
        },
        "support_exterior_absolute_max": support_exterior_absolute_max,
        "guards": guards,
        "failed_guards": failed,
        "public_signature": str(inspect.signature(evaluate_correction_laplacian_vector_identity)),
        "truth_boundary": truth_boundary(),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    receipt = verification_receipt()
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    else:
        print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0 if not receipt["failed_guards"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
