"""Correction-induced local transport increment for the frozen A2 curl family.

This module adds no new oscillatory or correction degree of freedom.  It combines
already-exposed Agent-2 quantities for an already-materialized
``SignedRadialAmplitudeCorrection`` into

    delta T_A2 = d_t(delta u)
               + (u_osc . grad) delta u
               + (delta u . grad) u_osc
               + (delta u . grad) delta u
               - nu Delta(delta u),

with the repository physical viscosity frozen at ``nu=0.01``.

Equivalently, in the absence of a leading field, pressure, and forcing, this is
the correction-induced change in the local transport operator
``d_t u + (u.grad)u - nu Delta u`` when moving from ``u_osc`` to
``u_osc + delta u``.  It is deliberately *not* called a Navier--Stokes
residual: leading/correction cross terms, pressure gradient, restricted forcing,
and the still-missing global leading field are outside this seam.

The corrected Kokuno 2026-09-09 reconstruction is structural provenance only
for the upstream localized oscillatory / signed complete-curl construction.
The public-z pullback, compact radial spline, autonomous time lift, Cartesian
finite differences, viscosity binding, and the transport composition below are
repository realizations, not paper-exact hidden data.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_correction_self_advection import (
    evaluate_correction_self_advection_fd4,
)
from .kokuno_public_oscillatory_correction_cross_advection import (
    evaluate_oscillatory_correction_cross_advection,
)
from .kokuno_public_signed_amplitude_complete_curl import (
    SignedRadialAmplitudeCorrection,
    correction_velocity,
    correction_velocity_dt,
    profile_from_delta_a,
)
from .kokuno_public_signed_amplitude_spatial_laplacian import (
    evaluate_correction_velocity_laplacian_fd6,
)
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-CORRECTION-TRANSPORT-INCREMENT-058"
SCHEMA = "kokuno-a2-correction-transport-increment-v1"
PARENT_AGENT2_PR = 795
PARENT_AGENT2_HEAD = "af99abe67feb45e29614972ab791f227ea004aa6"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
VISCOSITY = 0.01
OSCILLATORY_FD6_STEP = 0.001
CORRECTION_FD4_STEP = 0.001
CORRECTION_LAPLACIAN_FD6_STEP = 0.001
INDEPENDENT_FD4_STEPS = (0.004, 0.002, 0.001)

VectorProvider = Callable[[Any, Any, Any, Any], np.ndarray]


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


def _validate_step(step: float, name: str) -> float:
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
        raise RuntimeError(f"vector provider returned {value.shape}, expected {expected}")
    if not np.all(np.isfinite(value)):
        raise RuntimeError("vector provider returned non-finite values")
    return value


def evaluate_correction_transport_increment(
    correction: SignedRadialAmplitudeCorrection,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    oscillatory_spatial_step: float = OSCILLATORY_FD6_STEP,
    correction_spatial_step: float = CORRECTION_FD4_STEP,
    correction_laplacian_step: float = CORRECTION_LAPLACIAN_FD6_STEP,
) -> dict[str, Any]:
    """Evaluate the A2 correction-induced local transport increment.

    The caller supplies only an already-materialized correction witness and
    derivative-resolution settings.  Viscosity is frozen at the repository
    physical value ``0.01`` and is not a caller-tunable amplitude/residual knob.
    """
    ho = _validate_step(oscillatory_spatial_step, "oscillatory_spatial_step")
    hc = _validate_step(correction_spatial_step, "correction_spatial_step")
    hl = _validate_step(correction_laplacian_step, "correction_laplacian_step")

    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    dt = np.asarray(correction_velocity_dt(correction, xb, yb, zb, tb), dtype=float)
    cross = evaluate_oscillatory_correction_cross_advection(
        correction,
        xb,
        yb,
        zb,
        tb,
        oscillatory_spatial_step=ho,
        correction_spatial_step=hc,
    )
    self_term = evaluate_correction_self_advection_fd4(
        correction, xb, yb, zb, tb, spatial_step=hc
    )
    lap = evaluate_correction_velocity_laplacian_fd6(
        correction, xb, yb, zb, tb, spatial_step=hl
    )

    cross_value = np.asarray(cross["cross_advection"], dtype=float)
    self_value = np.asarray(self_term["correction_self_advection"], dtype=float)
    lap_value = np.asarray(lap["laplacian"], dtype=float)
    viscous_value = -VISCOSITY * lap_value
    increment = dt + cross_value + self_value + viscous_value

    expected = xb.shape + (3,)
    for name, value in (
        ("correction_velocity_dt", dt),
        ("cross_advection", cross_value),
        ("correction_self_advection", self_value),
        ("correction_laplacian", lap_value),
        ("transport_increment", increment),
    ):
        if value.shape != expected or not np.all(np.isfinite(value)):
            raise RuntimeError(f"{name} is non-finite or has the wrong shape")

    correction_values = (
        np.asarray(cross["correction_velocity"], dtype=float),
        np.asarray(self_term["correction_velocity"], dtype=float),
        np.asarray(lap["velocity"], dtype=float),
    )
    if not all(np.array_equal(correction_values[0], value) for value in correction_values[1:]):
        raise RuntimeError("upstream A2 correction velocity seams do not replay exactly")

    return {
        "correction_velocity": correction_values[0],
        "correction_velocity_dt": dt,
        "oscillation_advects_correction": np.asarray(
            cross["oscillation_advects_correction"], dtype=float
        ),
        "correction_advects_oscillation": np.asarray(
            cross["correction_advects_oscillation"], dtype=float
        ),
        "cross_advection": cross_value,
        "correction_self_advection": self_value,
        "correction_laplacian": lap_value,
        "correction_viscous_term": viscous_value,
        "correction_transport_increment": increment,
        "viscosity": VISCOSITY,
        "oscillatory_spatial_step": ho,
        "correction_spatial_step": hc,
        "correction_laplacian_step": hl,
        "handoff_contract": {
            "consumer_lane": "Kokuno Agent 3 mean/radial bookkeeping",
            "quantity_kind": "A2 correction-induced local transport increment",
            "complete_ns_defect": False,
            "includes_correction_time_derivative": True,
            "includes_oscillation_correction_cross_advection": True,
            "includes_correction_self_advection": True,
            "includes_correction_viscosity": True,
            "includes_leading_cross_terms": False,
            "includes_pressure_gradient": False,
            "includes_restricted_forcing": False,
        },
    }


def _fd4_time_derivative(
    provider: VectorProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    h = _validate_step(step, "independent_time_step")
    return (
        _vector_value(provider, x, y, z, t - 2.0 * h)
        - 8.0 * _vector_value(provider, x, y, z, t - h)
        + 8.0 * _vector_value(provider, x, y, z, t + h)
        - _vector_value(provider, x, y, z, t + 2.0 * h)
    ) / (12.0 * h)


def _fd4_laplacian(
    provider: VectorProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    h = _validate_step(step, "independent_spatial_step")
    base = (x, y, z)
    center = _vector_value(provider, x, y, z, t)
    lap = np.zeros_like(center)
    for axis in range(3):
        values: dict[int, np.ndarray] = {}
        for offset in (-2, -1, 1, 2):
            shifted = [array for array in base]
            shifted[axis] = shifted[axis] + offset * h
            values[offset] = _vector_value(
                provider, shifted[0], shifted[1], shifted[2], t
            )
        lap += (
            -values[-2]
            + 16.0 * values[-1]
            - 30.0 * center
            + 16.0 * values[1]
            - values[2]
        ) / (12.0 * h * h)
    return lap


def _directional_fd4_self_advection(
    provider: VectorProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    h = _validate_step(step, "independent_directional_step")
    velocity = _vector_value(provider, x, y, z, t)
    speed = np.linalg.norm(velocity, axis=-1)
    direction = np.zeros_like(velocity)
    nonzero = speed > 1.0e-14
    direction[nonzero] = velocity[nonzero] / speed[nonzero, None]
    values: dict[int, np.ndarray] = {}
    for offset in (-2, -1, 1, 2):
        displacement = offset * h * direction
        values[offset] = _vector_value(
            provider,
            x + displacement[..., 0],
            y + displacement[..., 1],
            z + displacement[..., 2],
            t,
        )
    derivative = (
        values[-2] - 8.0 * values[-1] + 8.0 * values[1] - values[2]
    ) / (12.0 * h)
    return speed[..., None] * derivative


def _independent_transport(
    provider: VectorProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    return (
        _fd4_time_derivative(provider, x, y, z, t, step)
        + _directional_fd4_self_advection(provider, x, y, z, t, step)
        - VISCOSITY * _fd4_laplacian(provider, x, y, z, t, step)
    )


def _independent_direct_difference(
    correction: SignedRadialAmplitudeCorrection,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)

    def osc_provider(xx: Any, yy: Any, zz: Any, tt: Any) -> np.ndarray:
        return np.asarray(velocity_osc(xx, yy, zz, tt), dtype=float)

    def total_provider(xx: Any, yy: Any, zz: Any, tt: Any) -> np.ndarray:
        return np.asarray(velocity_osc(xx, yy, zz, tt), dtype=float) + np.asarray(
            correction_velocity(correction, xx, yy, zz, tt), dtype=float
        )

    return _independent_transport(total_provider, xb, yb, zb, tb, step) - _independent_transport(
        osc_provider, xb, yb, zb, tb, step
    )


def _vector_rms(value: np.ndarray) -> float:
    array = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def _relative_rms(reference: np.ndarray, value: np.ndarray) -> float:
    return _vector_rms(np.asarray(value) - np.asarray(reference)) / max(
        _vector_rms(reference), np.finfo(float).tiny
    )


def _relative_sampled_max(reference: np.ndarray, value: np.ndarray) -> float:
    ref = np.asarray(reference, dtype=float)
    val = np.asarray(value, dtype=float)
    error = np.linalg.norm(val - ref, axis=-1)
    scale = max(float(np.max(np.linalg.norm(ref, axis=-1))), np.finfo(float).tiny)
    return float(np.max(error) / scale)


def _verification_profile() -> SignedRadialAmplitudeCorrection:
    radii = np.linspace(0.32, 1.18, 31)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(math.pi * s) ** 10
    delta_a = np.stack(
        (
            0.0077 * bump * (1.0 + 0.05 * np.cos(2.0 * math.pi * s)),
            -0.0054 * bump * (1.0 - 0.04 * np.sin(2.0 * math.pi * s)),
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
            "fresh compact signed radial mechanics profile for A2 correction transport "
            "increment verification only; no Agent-3 defect, mean, stress, inverse, "
            "damping, pressure, forcing, or residual input"
        ),
        source_mean_amplitude_differential_certified=False,
    )


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r0 = 0.32
    dr = (1.18 - 0.32) / 30.0
    cell_ids = (3, 8, 12, 17, 22, 26)
    radii = np.asarray([r0 + (index + 0.5) * dr for index in cell_ids], dtype=float)
    z_block = np.asarray((-0.54, -0.31, -0.08, 0.17, 0.35, 0.53), dtype=float)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    count = 0
    for time in (0.45, 0.50, 0.55):
        for radius, z_value in zip(radii, z_block, strict=True):
            theta = 0.317 + (count + 1) * golden
            rows.append(
                (
                    radius * math.cos(theta),
                    radius * math.sin(theta),
                    float(z_value),
                    time,
                )
            )
            count += 1
    array = np.asarray(rows, dtype=float)
    return array[:, 0], array[:, 1], array[:, 2], array[:, 3]


def build_receipt() -> dict[str, Any]:
    correction = _verification_profile()
    x, y, z, t = _verification_cloud()
    production = evaluate_correction_transport_increment(
        correction,
        x,
        y,
        z,
        t,
        oscillatory_spatial_step=OSCILLATORY_FD6_STEP,
        correction_spatial_step=CORRECTION_FD4_STEP,
        correction_laplacian_step=CORRECTION_LAPLACIAN_FD6_STEP,
    )
    reference = np.asarray(production["correction_transport_increment"], dtype=float)

    direct = [
        _independent_direct_difference(correction, x, y, z, t, step=step)
        for step in INDEPENDENT_FD4_STEPS
    ]
    successive = [
        _vector_rms(direct[0] - direct[1]),
        _vector_rms(direct[1] - direct[2]),
    ]
    refinement = successive[0] / max(successive[1], np.finfo(float).tiny)
    relative_rms = _relative_rms(reference, direct[-1])
    relative_max = _relative_sampled_max(reference, direct[-1])

    explicit_sum = (
        np.asarray(production["correction_velocity_dt"], dtype=float)
        + np.asarray(production["cross_advection"], dtype=float)
        + np.asarray(production["correction_self_advection"], dtype=float)
        + np.asarray(production["correction_viscous_term"], dtype=float)
    )
    assembly_max = float(np.max(np.abs(reference - explicit_sum)))

    exterior = evaluate_correction_transport_increment(
        correction,
        np.asarray((2.20, -2.22, 0.0, 1.72), dtype=float),
        np.asarray((0.0, 0.0, 2.21, 1.72), dtype=float),
        np.asarray((0.0, 0.0, 0.0, 2.20), dtype=float),
        np.full(4, 0.50, dtype=float),
    )
    exterior_abs_max = max(
        float(np.max(np.abs(np.asarray(exterior[key], dtype=float))))
        for key in (
            "correction_velocity",
            "correction_velocity_dt",
            "cross_advection",
            "correction_self_advection",
            "correction_laplacian",
            "correction_transport_increment",
        )
    )

    signature = inspect.signature(evaluate_correction_transport_increment)
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "delta_y",
        "delta_a", "normalized_score", "scientific_threshold", "viscosity", "nu",
    }
    failed: list[str] = []
    if _vector_rms(reference) < 1.0e-10:
        failed.append("correction_transport_increment_nontrivial")
    if refinement < 6.0:
        failed.append("independent_fd4_refinement")
    if relative_rms > 1.0e-2:
        failed.append("independent_finest_relative_rms")
    if relative_max > 2.0e-2:
        failed.append("independent_finest_relative_sampled_max")
    if assembly_max > 1.0e-13:
        failed.append("termwise_assembly")
    if exterior_abs_max > 1.0e-12:
        failed.append("support_exterior")
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
            "structural_scope": "localized oscillatory and signed complete-curl fields",
        },
        "repository_physical_contract": {
            "viscosity": VISCOSITY,
            "viscosity_caller_tunable": False,
            "quantity": (
                "correction-induced local transport increment relative to frozen u_osc; "
                "not complete NS residual"
            ),
        },
        "sample_count": int(t.size),
        "sample_times": [0.45, 0.50, 0.55],
        "oscillatory_fd6_step": OSCILLATORY_FD6_STEP,
        "correction_fd4_step": CORRECTION_FD4_STEP,
        "correction_laplacian_fd6_step": CORRECTION_LAPLACIAN_FD6_STEP,
        "independent_fd4_steps": list(INDEPENDENT_FD4_STEPS),
        "correction_transport_increment_rms": _vector_rms(reference),
        "independent_successive_difference_rms": successive,
        "independent_refinement_ratio": float(refinement),
        "independent_finest_relative_rms": float(relative_rms),
        "independent_finest_relative_sampled_max": float(relative_max),
        "termwise_assembly_absolute_max": assembly_max,
        "support_exterior_absolute_max": exterior_abs_max,
        "failed_guards": failed,
        "handoff_contract": production["handoff_contract"],
        "truth_boundary": {
            "oscillatory_or_correction_candidate_changed": False,
            "amplitude_phase_carrier_support_retuned": False,
            "vector_potential_complete_curl_reused": True,
            "correction_transport_increment_executable": True,
            "agent3_raw_transport_handoff_available": True,
            "agent3_mean_radial_chain_reimplemented": False,
            "real_agent3_delta_a_bound": False,
            "real_global_leading_field_bound": False,
            "leading_cross_terms_included": False,
            "pressure_gradient_included": False,
            "restricted_forcing_included": False,
            "complete_ns_defect": False,
            "full_composite_velocity_available": False,
            "heldout_ns_momentum_residual_assessed": False,
            "residual_reduction_claimed": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "pde_validated": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    receipt = build_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0 if not receipt["failed_guards"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
