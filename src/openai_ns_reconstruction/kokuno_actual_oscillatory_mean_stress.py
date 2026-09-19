"""Materialize the admitted oscillatory candidate's real mean-stress burden.

Kokuno Agent 3 owns the mean-correction / compact radial-stress lane.  This
module consumes the *actual frozen Agent-2 public field* rather than accepting a
caller-supplied residual array.  It applies a fixed raw Navier--Stokes operator
component

    R_osc = d_t w + (w . grad) w - nu Delta w

with no fitted pressure and no fitted forcing, ring-averages its physical
cylindrical components, and reconstructs the compact theta/e=2 and z/e=1
radial stresses with the corrected-reader moment-complement inverse

    M_e F = integral r^e F dr,
    P_e F = F - b_e M_e F,       integral r^e b_e dr = 1,
    sigma_e(r) = -r^(-e) integral_0^r s^e P_e(s) ds.

The production time derivative is Agent-2 #579's exact derivative of the frozen
repository-autonomous time law.  Spatial first/second derivatives are evaluated
with one fixed centered FD4 stencil in *public Cartesian coordinates*.

This is deliberately narrower than the eventual same-cycle composite defect:
it materializes the real oscillatory self-defect / quadratic-mean burden of the
independently admitted Agent-2 #561 field.  Agent-1 leading cross terms, matched
pressure, restricted forcing, finite-head weighting, signed covariance inverse,
and any public correction velocity remain downstream.  In particular, this
module cannot accept a hand-filled/surrogate defect as success evidence.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_public_oscillatory_time_derivative import velocity_osc_dt
from .kokuno_public_z_pullback_velocity import default_field, velocity_osc

TASK = "KOKUNO-A3-ACTUAL-OSCILLATORY-MEAN-STRESS-039"
SCHEMA = "kokuno-a3-actual-oscillatory-mean-stress-v1"

# Exact upstream identities already independently admitted by A3 #571.
ADMITTED_AGENT2_PR = 561
ADMITTED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
INDEPENDENT_AGENT4_PR = 563
INDEPENDENT_AGENT4_HEAD = "9b0f86012c53fa8e32a19f766dbc150931870425"
INDEPENDENT_AGENT4_WORKFLOW = 35422203621
PARENT_AGENT2_DERIVATIVE_PR = 579
PARENT_AGENT2_DERIVATIVE_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"

DEFAULT_NU = 0.01
DEFAULT_FD4_STEP = 0.005
DEFAULT_TIME = 0.50
DEFAULT_Z = 0.08
DEFAULT_RADIAL_COUNT = 49
DEFAULT_ANGULAR_COUNT = 48


def _finite_scalar(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _vector_rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    if array.ndim < 1 or array.shape[-1] != 3:
        raise ValueError("vector array must end in Cartesian/cylindrical dimension 3")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def _scalar_rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def _cumulative_trapezoid(values: np.ndarray, coordinates: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    coordinates = np.asarray(coordinates, dtype=float)
    if values.ndim != 1 or coordinates.ndim != 1 or values.shape != coordinates.shape:
        raise ValueError("cumulative trapezoid requires matching one-dimensional arrays")
    if coordinates.size < 3 or np.any(np.diff(coordinates) <= 0.0):
        raise ValueError("coordinates must be strictly increasing with at least three nodes")
    out = np.zeros_like(values)
    increments = 0.5 * (values[1:] + values[:-1]) * np.diff(coordinates)
    out[1:] = np.cumsum(increments)
    return out


def _compact_cos8_bump(radii: np.ndarray, center: float, halfwidth: float) -> np.ndarray:
    radii = np.asarray(radii, dtype=float)
    s = (radii - center) / halfwidth
    out = np.zeros_like(s)
    mask = np.abs(s) < 1.0
    if np.any(mask):
        out[mask] = np.cos(0.5 * math.pi * s[mask]) ** 8
    return out


def _compact_radial_stress(
    radii: np.ndarray,
    source: np.ndarray,
    *,
    exponent: int,
    bump_center: float,
    bump_halfwidth: float,
) -> dict[str, Any]:
    """Apply the discrete moment-complement form of Kokuno's compact inverse."""
    r = np.asarray(radii, dtype=float)
    f = np.asarray(source, dtype=float)
    if r.ndim != 1 or f.shape != r.shape or r.size < 9:
        raise ValueError("radial stress requires matching one-dimensional arrays with >=9 nodes")
    if exponent not in (1, 2):
        raise ValueError("only the source axial e=1 and angular e=2 channels are supported")
    if np.any(r <= 0.0) or not np.all(np.isfinite(f)):
        raise ValueError("radii must be positive and source finite")

    weight = r ** exponent
    moment = float(_cumulative_trapezoid(weight * f, r)[-1])
    bump_raw = _compact_cos8_bump(r, bump_center, bump_halfwidth)
    bump_weighted_integral = float(_cumulative_trapezoid(weight * bump_raw, r)[-1])
    if not math.isfinite(bump_weighted_integral) or bump_weighted_integral <= 0.0:
        raise RuntimeError("compact bump lost positive weighted normalization")
    bump = bump_raw / bump_weighted_integral
    complement = f - bump * moment
    complement_moment = float(_cumulative_trapezoid(weight * complement, r)[-1])
    primitive = _cumulative_trapezoid(weight * complement, r)
    stress = -primitive / weight

    # This is the exact algebraic RHS of (d_r+e/r)sigma_e=-F+b_e M_e.
    reconstructed_force = -f + bump * moment
    return {
        "exponent": exponent,
        "weighted_moment": moment,
        "bump_weighted_integral": float(
            _cumulative_trapezoid(weight * bump, r)[-1]
        ),
        "moment_complement_weighted_moment": complement_moment,
        "stress": stress,
        "stress_rms": _scalar_rms(stress),
        "stress_max_abs": float(np.max(np.abs(stress))),
        "stress_inner_edge": float(stress[0]),
        "stress_outer_edge": float(stress[-1]),
        "reconstructed_force_rms": _scalar_rms(reconstructed_force),
    }


def _fd4_spatial_operator(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    step: float,
    nu: float,
) -> dict[str, np.ndarray]:
    """Evaluate the raw oscillatory self-operator using public value calls only."""
    h = _finite_scalar(step, "step")
    viscosity = _finite_scalar(nu, "nu")
    if h <= 0.0 or viscosity < 0.0:
        raise ValueError("step must be positive and nu nonnegative")

    coords = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]
    if not (coords[0].shape == coords[1].shape == coords[2].shape):
        raise ValueError("x/y/z must have the same shape")
    tt = np.broadcast_to(np.asarray(t, dtype=float), coords[0].shape)
    if not all(np.all(np.isfinite(c)) for c in (*coords, tt)):
        raise ValueError("spacetime coordinates must be finite")

    base = np.asarray(velocity_osc(coords[0], coords[1], coords[2], tt), dtype=float)
    dt = np.asarray(velocity_osc_dt(coords[0], coords[1], coords[2], tt), dtype=float)
    if base.shape != coords[0].shape + (3,) or dt.shape != base.shape:
        raise RuntimeError("public velocity/time-derivative shape contract changed")

    first: list[np.ndarray] = []
    second: list[np.ndarray] = []
    for axis in range(3):
        shifted: dict[int, np.ndarray] = {}
        for multiple in (-2, -1, 1, 2):
            args = [c.copy() for c in coords]
            args[axis] = args[axis] + multiple * h
            shifted[multiple] = np.asarray(
                velocity_osc(args[0], args[1], args[2], tt), dtype=float
            )
        first.append(
            (
                -shifted[2]
                + 8.0 * shifted[1]
                - 8.0 * shifted[-1]
                + shifted[-2]
            )
            / (12.0 * h)
        )
        second.append(
            (
                -shifted[2]
                + 16.0 * shifted[1]
                - 30.0 * base
                + 16.0 * shifted[-1]
                - shifted[-2]
            )
            / (12.0 * h * h)
        )

    gradient = np.stack(first, axis=-2)  # (..., spatial derivative, velocity component)
    laplacian = second[0] + second[1] + second[2]
    convective = np.einsum("...j,...jk->...k", base, gradient)
    linear = dt - viscosity * laplacian
    residual = linear + convective
    return {
        "velocity": base,
        "time_derivative": dt,
        "gradient": gradient,
        "laplacian": laplacian,
        "linear": linear,
        "quadratic": convective,
        "raw_self_residual": residual,
    }


def _cartesian_to_cylindrical(vectors: np.ndarray, theta: np.ndarray) -> np.ndarray:
    vectors = np.asarray(vectors, dtype=float)
    theta = np.asarray(theta, dtype=float)
    c = np.cos(theta)
    s = np.sin(theta)
    radial = c * vectors[..., 0] + s * vectors[..., 1]
    angular = -s * vectors[..., 0] + c * vectors[..., 1]
    return np.stack((radial, angular, vectors[..., 2]), axis=-1)


def materialize_actual_oscillatory_mean_stress(
    *,
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    time: float = DEFAULT_TIME,
    z: float = DEFAULT_Z,
    fd4_step: float = DEFAULT_FD4_STEP,
    nu: float = DEFAULT_NU,
) -> dict[str, Any]:
    """Materialize the real A2 oscillatory self-defect mean and compact stresses.

    There is intentionally no ``defect`` or ``residual`` argument: the operator
    is evaluated directly on the admitted public candidate.  The returned
    stresses are therefore candidate-generated data, not caller-supplied target
    data.  They remain only one component of the future full composite defect.

    The registered field is exactly zero at the two radial support boundaries.
    We therefore evaluate the public provider only on strict interior radial
    nodes and insert those two known zero-extension boundary values before the
    compact radial integral.  This avoids asking upstream complete-curl code to
    normalize a zero tangent at a support endpoint while preserving the compact
    support/moment contract.
    """
    if int(radial_count) != radial_count or radial_count < 17:
        raise ValueError("radial_count must be an integer >=17")
    if int(angular_count) != angular_count or angular_count < 16:
        raise ValueError("angular_count must be an integer >=16")
    radial_count = int(radial_count)
    angular_count = int(angular_count)
    time_value = _finite_scalar(time, "time")
    z_value = _finite_scalar(z, "z")
    h = _finite_scalar(fd4_step, "fd4_step")
    if h <= 0.0:
        raise ValueError("fd4_step must be positive")

    field = default_field()
    if not (field.time_min <= time_value <= field.time_max):
        raise ValueError("time is outside the admitted candidate interval")
    radial_min = field.radial_center - field.radial_halfwidth
    radial_max = field.radial_center + field.radial_halfwidth
    if radial_min <= 0.0:
        raise RuntimeError("admitted field no longer has an axis-safe radial annulus")
    if abs(z_value - field.axial_center) >= field.axial_halfwidth:
        raise ValueError("z must lie strictly inside the admitted axial support")

    radii = np.linspace(radial_min, radial_max, radial_count)
    sample_radii = radii[1:-1]
    boundary_margin = min(
        float(sample_radii[0] - radial_min),
        float(radial_max - sample_radii[-1]),
    )
    if 2.0 * h >= boundary_margin:
        raise ValueError(
            "FD4 stencil must remain strictly inside the registered radial support; "
            "increase radial boundary margin or decrease fd4_step"
        )

    # Mid-cell angles avoid privileging the Cartesian axes while retaining an exact uniform ring.
    angles = 2.0 * math.pi * (np.arange(angular_count, dtype=float) + 0.5) / angular_count
    rr, tt = np.meshgrid(sample_radii, angles, indexing="ij")
    x = (rr * np.cos(tt)).reshape(-1)
    y = (rr * np.sin(tt)).reshape(-1)
    zz = np.full_like(x, z_value)
    times = np.full_like(x, time_value)

    operator = _fd4_spatial_operator(x, y, zz, times, step=h, nu=nu)
    interior_count = radial_count - 2
    cylindrical: dict[str, np.ndarray] = {}
    for name in ("linear", "quadratic", "raw_self_residual"):
        cyl = _cartesian_to_cylindrical(operator[name], tt.reshape(-1))
        cylindrical[name] = cyl.reshape(interior_count, angular_count, 3)

    ring_mean_interior = {
        name: np.mean(values, axis=1) for name, values in cylindrical.items()
    }
    closure = ring_mean_interior["raw_self_residual"] - (
        ring_mean_interior["linear"] + ring_mean_interior["quadratic"]
    )

    # The public field has a registered smooth zero extension outside the annulus,
    # so the exact endpoint values of all three operator components are inserted
    # as zero rather than obtained by floating-point evaluation on the boundary.
    ring_mean_full: dict[str, np.ndarray] = {}
    for name, interior in ring_mean_interior.items():
        full = np.zeros((radial_count, 3), dtype=float)
        full[1:-1] = interior
        ring_mean_full[name] = full

    # requestedStress ordering follows the corrected formal construction:
    # theta residual -> physicalBarSigma_2, axial residual -> physicalBarSigma_1.
    theta_stress = _compact_radial_stress(
        radii,
        ring_mean_full["raw_self_residual"][:, 1],
        exponent=2,
        bump_center=field.radial_center,
        bump_halfwidth=field.radial_halfwidth,
    )
    axial_stress = _compact_radial_stress(
        radii,
        ring_mean_full["raw_self_residual"][:, 2],
        exponent=1,
        bump_center=field.radial_center,
        bump_halfwidth=field.radial_halfwidth,
    )

    requested_stress = np.stack((theta_stress["stress"], axial_stress["stress"]), axis=-1)
    raw_mean = ring_mean_interior["raw_self_residual"]
    quadratic_mean = ring_mean_interior["quadratic"]
    linear_mean = ring_mean_interior["linear"]
    raw_mean_rms = _vector_rms(raw_mean)
    quadratic_mean_rms = _vector_rms(quadratic_mean)
    linear_mean_rms = _vector_rms(linear_mean)

    signature = inspect.signature(materialize_actual_oscillatory_mean_stress)
    prohibited_names = {"defect", "residual", "target", "stress"}
    caller_supplied_target_parameters = [
        name for name in signature.parameters if name.lower() in prohibited_names
    ]

    return {
        "task": TASK,
        "schema": SCHEMA,
        "provenance": {
            "kokuno_public_structure": "corrected 2026-09-09 reconstruction; compact moment-complement inverse",
            "admitted_agent2_pr": ADMITTED_AGENT2_PR,
            "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
            "independent_agent4_pr": INDEPENDENT_AGENT4_PR,
            "independent_agent4_head": INDEPENDENT_AGENT4_HEAD,
            "independent_agent4_workflow": INDEPENDENT_AGENT4_WORKFLOW,
            "parent_agent2_derivative_pr": PARENT_AGENT2_DERIVATIVE_PR,
            "parent_agent2_derivative_head": PARENT_AGENT2_DERIVATIVE_HEAD,
            "velocity_provider": "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc",
            "time_derivative_provider": "openai_ns_reconstruction.kokuno_public_oscillatory_time_derivative:velocity_osc_dt",
        },
        "state": {
            "time": time_value,
            "z": z_value,
            "nu": float(nu),
            "fd4_step": h,
            "radial_count": radial_count,
            "radial_evaluation_count": interior_count,
            "angular_count": angular_count,
            "radial_support": [float(radial_min), float(radial_max)],
            "radial_boundary_margin": boundary_margin,
            "boundary_values_inserted_from_registered_zero_extension": True,
        },
        "mean_operator": {
            "kind": "physical_uniform_angle_mean_of_raw_oscillatory_self_operator",
            "raw_self_operator": "dt(w)+(w.grad)w-nu*Delta(w)",
            "pressure_used": False,
            "forcing_used": False,
            "raw_mean_vector_rms": raw_mean_rms,
            "raw_mean_sampled_max": float(np.max(np.linalg.norm(raw_mean, axis=-1))),
            "quadratic_mean_vector_rms": quadratic_mean_rms,
            "linear_mean_vector_rms": linear_mean_rms,
            "linear_to_raw_mean_rms_ratio": linear_mean_rms / max(raw_mean_rms, 1.0e-300),
            "quadratic_to_raw_mean_rms_ratio": quadratic_mean_rms / max(raw_mean_rms, 1.0e-300),
            "theta_raw_mean_rms": _scalar_rms(raw_mean[:, 1]),
            "axial_raw_mean_rms": _scalar_rms(raw_mean[:, 2]),
            "radial_raw_mean_rms": _scalar_rms(raw_mean[:, 0]),
            "mean_decomposition_closure_max_abs": float(np.max(np.abs(closure))),
        },
        "requested_stress_component": {
            "ordering": ["theta_e2", "axial_e1"],
            "requested_stress_rms": float(np.sqrt(np.mean(np.sum(requested_stress**2, axis=-1)))),
            "requested_stress_max_abs": float(np.max(np.abs(requested_stress))),
            "theta_e2": {key: value for key, value in theta_stress.items() if key != "stress"},
            "axial_e1": {key: value for key, value in axial_stress.items() if key != "stress"},
            "theta_stress_values": theta_stress["stress"].tolist(),
            "axial_stress_values": axial_stress["stress"].tolist(),
            "radial_nodes": radii.tolist(),
        },
        "anti_surrogate_contract": {
            "caller_supplied_target_parameters": caller_supplied_target_parameters,
            "surrogate_defect_used": False,
            "actual_public_velocity_evaluated_internally": True,
        },
        "truth_boundary": {
            "admitted_agent2_velocity_consumed": True,
            "exact_candidate_time_derivative_consumed": True,
            "real_oscillatory_self_defect_component_consumed": True,
            "real_oscillatory_quadratic_mean_component_measured": True,
            "oscillatory_requested_stress_component_materialized": True,
            "full_same_cycle_composite_requested_stress_materialized": False,
            "agent1_leading_cross_terms_included": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
            "autonomous_finite_head_factor_applied": False,
            "candidate_finite_head_mean_debt_materialized": False,
            "signed_mean_inverse_input_ready": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_rerun_allowed": False,
            "finite_correction_cycle_run": False,
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--radial-count", type=int, default=DEFAULT_RADIAL_COUNT)
    parser.add_argument("--angular-count", type=int, default=DEFAULT_ANGULAR_COUNT)
    parser.add_argument("--time", type=float, default=DEFAULT_TIME)
    parser.add_argument("--z", type=float, default=DEFAULT_Z)
    parser.add_argument("--fd4-step", type=float, default=DEFAULT_FD4_STEP)
    parser.add_argument("--nu", type=float, default=DEFAULT_NU)
    args = parser.parse_args()
    receipt = materialize_actual_oscillatory_mean_stress(
        radial_count=args.radial_count,
        angular_count=args.angular_count,
        time=args.time,
        z=args.z,
        fd4_step=args.fd4_step,
        nu=args.nu,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
