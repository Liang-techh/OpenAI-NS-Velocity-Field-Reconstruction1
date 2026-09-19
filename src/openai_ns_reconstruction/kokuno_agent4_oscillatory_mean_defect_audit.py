"""Independent Agent-4 audit of Agent-3's oscillatory mean-defect input.

The audit consumes only the already-admitted public black-box interfaces

    velocity_osc(x, y, z, t)
    velocity_osc_dt(x, y, z, t)

and independently reconstructs the raw oscillatory self-operator

    Q[w] = d_t w + (w . grad) w - nu Delta w

with a centered Cartesian FD6 spatial operator. It never calls Agent-3's FD4
self-defect or compact radial-stress helpers. Physical-angle means use
Gauss--Legendre quadrature rather than Agent-3's uniform midpoint ring.

This validates the component-level numerical input entering the correction
lane. It is not a full Navier--Stokes residual because global leading/cross
terms, matched pressure, restricted forcing and a public correction velocity do
not yet coexist in one candidate. Final project gates remain unchanged:
normalized momentum max/L2 <= 1e-3 and divergence max/L2 <= 1e-5.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
from numpy.polynomial.legendre import leggauss

from .kokuno_public_oscillatory_time_derivative import velocity_osc_dt
from .kokuno_public_z_pullback_velocity import default_field, velocity_osc

TASK = "KOKUNO-A4-OSCILLATORY-MEAN-DEFECT-INDEPENDENT-AUDIT-043"
SCHEMA = "kokuno-agent4-oscillatory-mean-defect-audit-v1"

AUDITED_AGENT3_PR = 581
AUDITED_AGENT3_HEAD = "c95dab915ae25e689e7c82bf1b8614e205f9274c"
AUDITED_AGENT3_WORKFLOW = 35428028030
AUDITED_AGENT3_ARTIFACT_ID = 10579891653
AUDITED_AGENT3_ARTIFACT_DIGEST = (
    "sha256:269a6948fde38c5ec1510c44917ee857dfaaa3f01fac5dbb3b9f6519d62a5cf5"
)
PARENT_AGENT3_MEAN_DEBT_PR = 586
PARENT_AGENT3_MEAN_DEBT_HEAD = "ca4da20c90ccd7a5046f79fb2e7965ec05b6086f"
ADMITTED_AGENT2_PR = 561
ADMITTED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
AUDITED_AGENT2_DERIVATIVE_PR = 579
AUDITED_AGENT2_DERIVATIVE_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"
INDEPENDENT_AGENT4_VELOCITY_PR = 563
INDEPENDENT_AGENT4_VELOCITY_HEAD = "9b0f86012c53fa8e32a19f766dbc150931870425"
INDEPENDENT_AGENT4_TIME_DERIVATIVE_PR = 580
INDEPENDENT_AGENT4_TIME_DERIVATIVE_HEAD = "6f7edb4d66dcb162a8510a40e4227f48f3c57e28"

RADIAL_COUNT = 25
ANGULAR_ORDER = 32
TIME = 0.50
Z = 0.08
NU = 0.01
FD6_STEPS = (0.01, 0.005, 0.0025)

AGENT3_RECEIPT_TARGET = {
    "raw_mean_vector_rms": 96194924.96271932,
    "raw_mean_sampled_max": 225706907.52261892,
    "radial_raw_mean_rms": 86751074.56430492,
    "theta_raw_mean_rms": 31210509.501432016,
    "axial_raw_mean_rms": 27452117.353359643,
    "quadratic_mean_vector_rms": 96194924.96273328,
    "linear_mean_vector_rms": 0.00033125340454790233,
}

GUARDS = {
    "finest_receipt_raw_mean_rms_relative": 5.0e-3,
    "finest_receipt_raw_mean_max_relative": 1.0e-2,
    "finest_receipt_component_rms_relative": 1.0e-2,
    "finest_receipt_quadratic_rms_relative": 5.0e-3,
    "finest_fd6_profile_change": 2.0e-4,
    "minimum_fd6_profile_refinement_ratio": 8.0,
    "minimum_quadratic_to_raw_ratio": 0.99,
    "maximum_linear_to_raw_ratio": 1.0e-7,
    "maximum_decomposition_closure_relative": 1.0e-12,
    "minimum_quadratic_sign_flip_response": 1.0,
}

VelocityEvaluator = Callable[[Any, Any, Any, Any], np.ndarray]


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.ndim < 1 or arr.shape[-1] != 3:
        raise ValueError("vector array must end in dimension 3")
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def _scalar_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(arr * arr)))


def _relative_vector_rms(delta: np.ndarray, reference: np.ndarray) -> float:
    return _vector_rms(delta) / max(_vector_rms(reference), 1.0e-300)


def _relative_scalar(value: float, reference: float) -> float:
    return abs(float(value) - float(reference)) / max(abs(float(reference)), 1.0e-300)


def _fd6_spatial_operator(
    velocity: VelocityEvaluator,
    velocity_dt: VelocityEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    step: float,
    nu: float,
    base: np.ndarray | None = None,
    dt: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Independent centered-FD6 Cartesian self-operator using public values."""
    h = float(step)
    viscosity = float(nu)
    if not math.isfinite(h) or h <= 0.0:
        raise ValueError("step must be positive and finite")
    if not math.isfinite(viscosity) or viscosity < 0.0:
        raise ValueError("nu must be finite and nonnegative")

    coords = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]
    if not (coords[0].shape == coords[1].shape == coords[2].shape):
        raise ValueError("x/y/z must share shape")
    tt = np.broadcast_to(np.asarray(t, dtype=float), coords[0].shape)
    if not all(np.all(np.isfinite(v)) for v in (*coords, tt)):
        raise ValueError("coordinates/time must be finite")

    if base is None:
        base_arr = np.asarray(velocity(coords[0], coords[1], coords[2], tt), dtype=float)
    else:
        base_arr = np.asarray(base, dtype=float)
    if dt is None:
        dt_arr = np.asarray(velocity_dt(coords[0], coords[1], coords[2], tt), dtype=float)
    else:
        dt_arr = np.asarray(dt, dtype=float)
    expected = coords[0].shape + (3,)
    if base_arr.shape != expected or dt_arr.shape != expected:
        raise RuntimeError("public velocity/time-derivative shape contract changed")

    first: list[np.ndarray] = []
    second: list[np.ndarray] = []
    for axis in range(3):
        shifted: dict[int, np.ndarray] = {}
        for multiple in (-3, -2, -1, 1, 2, 3):
            args = [v.copy() for v in coords]
            args[axis] = args[axis] + multiple * h
            shifted[multiple] = np.asarray(
                velocity(args[0], args[1], args[2], tt), dtype=float
            )
        d1 = (
            -shifted[-3]
            + 9.0 * shifted[-2]
            - 45.0 * shifted[-1]
            + 45.0 * shifted[1]
            - 9.0 * shifted[2]
            + shifted[3]
        ) / (60.0 * h)
        d2 = (
            2.0 * shifted[-3]
            - 27.0 * shifted[-2]
            + 270.0 * shifted[-1]
            - 490.0 * base_arr
            + 270.0 * shifted[1]
            - 27.0 * shifted[2]
            + 2.0 * shifted[3]
        ) / (180.0 * h * h)
        first.append(d1)
        second.append(d2)

    gradient = np.stack(first, axis=-2)
    laplacian = second[0] + second[1] + second[2]
    quadratic = np.einsum("...j,...jk->...k", base_arr, gradient)
    linear = dt_arr - viscosity * laplacian
    raw = linear + quadratic
    return {
        "velocity": base_arr,
        "time_derivative": dt_arr,
        "gradient": gradient,
        "laplacian": laplacian,
        "linear": linear,
        "quadratic": quadratic,
        "raw_self_residual": raw,
    }


def _cartesian_to_cylindrical(vectors: np.ndarray, theta: np.ndarray) -> np.ndarray:
    values = np.asarray(vectors, dtype=float)
    angle = np.asarray(theta, dtype=float)
    c = np.cos(angle)
    s = np.sin(angle)
    radial = c * values[..., 0] + s * values[..., 1]
    angular = -s * values[..., 0] + c * values[..., 1]
    return np.stack((radial, angular, values[..., 2]), axis=-1)


def _ring_cloud() -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    field = default_field()
    radial_min = float(field.radial_center - field.radial_halfwidth)
    radial_max = float(field.radial_center + field.radial_halfwidth)
    radii = np.linspace(radial_min, radial_max, RADIAL_COUNT)
    sample_radii = radii[1:-1]

    nodes, weights = leggauss(ANGULAR_ORDER)
    theta = math.pi * nodes
    mean_weights = 0.5 * weights

    rr, aa = np.meshgrid(sample_radii, theta, indexing="ij")
    x = (rr * np.cos(aa)).reshape(-1)
    y = (rr * np.sin(aa)).reshape(-1)
    z = np.full_like(x, Z)
    t = np.full_like(x, TIME)
    return sample_radii, theta, mean_weights, x, y, z, t


def _weighted_ring_mean(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != (RADIAL_COUNT - 2, ANGULAR_ORDER, 3):
        raise ValueError("unexpected ring tensor shape")
    w = np.asarray(weights, dtype=float)
    if w.shape != (ANGULAR_ORDER,):
        raise ValueError("unexpected angular weights")
    return np.einsum("j,rjc->rc", w, arr)


def _metric_row(
    operator: dict[str, np.ndarray], theta: np.ndarray, weights: np.ndarray, step: float
) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    n_r = RADIAL_COUNT - 2
    means: dict[str, np.ndarray] = {}
    angle_flat = np.tile(theta, n_r)
    for key in ("linear", "quadratic", "raw_self_residual"):
        converted = _cartesian_to_cylindrical(operator[key], angle_flat).reshape(
            n_r, ANGULAR_ORDER, 3
        )
        means[key] = _weighted_ring_mean(converted, weights)

    closure = means["raw_self_residual"] - (means["linear"] + means["quadratic"])
    raw = means["raw_self_residual"]
    quadratic = means["quadratic"]
    linear = means["linear"]
    row = {
        "step": float(step),
        "raw_mean_vector_rms": _vector_rms(raw),
        "raw_mean_sampled_max": float(np.max(np.linalg.norm(raw, axis=-1), initial=0.0)),
        "radial_raw_mean_rms": _scalar_rms(raw[:, 0]),
        "theta_raw_mean_rms": _scalar_rms(raw[:, 1]),
        "axial_raw_mean_rms": _scalar_rms(raw[:, 2]),
        "quadratic_mean_vector_rms": _vector_rms(quadratic),
        "linear_mean_vector_rms": _vector_rms(linear),
        "quadratic_to_raw_mean_rms_ratio": _vector_rms(quadratic)
        / max(_vector_rms(raw), 1.0e-300),
        "linear_to_raw_mean_rms_ratio": _vector_rms(linear)
        / max(_vector_rms(raw), 1.0e-300),
        "decomposition_closure_max_abs": float(np.max(np.abs(closure), initial=0.0)),
        "decomposition_closure_relative": _vector_rms(closure)
        / max(_vector_rms(raw), 1.0e-300),
    }
    return row, means


def run_audit() -> dict[str, Any]:
    sample_radii, theta, weights, x, y, z, t = _ring_cloud()
    base = np.asarray(velocity_osc(x, y, z, t), dtype=float)
    dt = np.asarray(velocity_osc_dt(x, y, z, t), dtype=float)

    rows: list[dict[str, float]] = []
    profiles: list[dict[str, np.ndarray]] = []
    for step in FD6_STEPS:
        operator = _fd6_spatial_operator(
            velocity_osc,
            velocity_osc_dt,
            x,
            y,
            z,
            t,
            step=step,
            nu=NU,
            base=base,
            dt=dt,
        )
        row, means = _metric_row(operator, theta, weights, step)
        rows.append(row)
        profiles.append(means)

    profile_changes = [
        _relative_vector_rms(
            profiles[i]["raw_self_residual"] - profiles[i + 1]["raw_self_residual"],
            profiles[i + 1]["raw_self_residual"],
        )
        for i in range(len(profiles) - 1)
    ]
    refinement_ratio = profile_changes[0] / max(profile_changes[1], 1.0e-300)

    finest = rows[-1]
    receipt_relative = {
        name: _relative_scalar(finest[name], reference)
        for name, reference in AGENT3_RECEIPT_TARGET.items()
    }

    finest_means = profiles[-1]
    raw = finest_means["raw_self_residual"]
    raw_norm = np.linalg.norm(raw, axis=-1)
    max_index = int(np.argmax(raw_norm))
    profile_concentration = {
        "max_vector_norm_radius": float(sample_radii[max_index]),
        "max_vector_norm": float(raw_norm[max_index]),
        "max_vector_components_cylindrical": raw[max_index].tolist(),
        "component_rms": {
            "radial": _scalar_rms(raw[:, 0]),
            "theta": _scalar_rms(raw[:, 1]),
            "axial": _scalar_rms(raw[:, 2]),
        },
    }
    sign_flip = finest_means["linear"] - finest_means["quadratic"]
    sign_flip_response = _relative_vector_rms(sign_flip - raw, raw)

    failed: list[str] = []
    if receipt_relative["raw_mean_vector_rms"] > GUARDS["finest_receipt_raw_mean_rms_relative"]:
        failed.append("finest_receipt_raw_mean_rms_relative")
    if receipt_relative["raw_mean_sampled_max"] > GUARDS["finest_receipt_raw_mean_max_relative"]:
        failed.append("finest_receipt_raw_mean_max_relative")
    component_names = ("radial_raw_mean_rms", "theta_raw_mean_rms", "axial_raw_mean_rms")
    if max(receipt_relative[name] for name in component_names) > GUARDS[
        "finest_receipt_component_rms_relative"
    ]:
        failed.append("finest_receipt_component_rms_relative")
    if receipt_relative["quadratic_mean_vector_rms"] > GUARDS[
        "finest_receipt_quadratic_rms_relative"
    ]:
        failed.append("finest_receipt_quadratic_rms_relative")
    if profile_changes[-1] > GUARDS["finest_fd6_profile_change"]:
        failed.append("finest_fd6_profile_change")
    if refinement_ratio < GUARDS["minimum_fd6_profile_refinement_ratio"]:
        failed.append("minimum_fd6_profile_refinement_ratio")
    if finest["quadratic_to_raw_mean_rms_ratio"] < GUARDS["minimum_quadratic_to_raw_ratio"]:
        failed.append("minimum_quadratic_to_raw_ratio")
    if finest["linear_to_raw_mean_rms_ratio"] > GUARDS["maximum_linear_to_raw_ratio"]:
        failed.append("maximum_linear_to_raw_ratio")
    if finest["decomposition_closure_relative"] > GUARDS[
        "maximum_decomposition_closure_relative"
    ]:
        failed.append("maximum_decomposition_closure_relative")
    if sign_flip_response < GUARDS["minimum_quadratic_sign_flip_response"]:
        failed.append("minimum_quadratic_sign_flip_response")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "provenance": {
            "audited_agent3_pr": AUDITED_AGENT3_PR,
            "audited_agent3_head": AUDITED_AGENT3_HEAD,
            "audited_agent3_workflow": AUDITED_AGENT3_WORKFLOW,
            "audited_agent3_artifact_id": AUDITED_AGENT3_ARTIFACT_ID,
            "audited_agent3_artifact_digest": AUDITED_AGENT3_ARTIFACT_DIGEST,
            "parent_agent3_mean_debt_pr": PARENT_AGENT3_MEAN_DEBT_PR,
            "parent_agent3_mean_debt_head": PARENT_AGENT3_MEAN_DEBT_HEAD,
            "admitted_agent2_pr": ADMITTED_AGENT2_PR,
            "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
            "agent2_derivative_pr": AUDITED_AGENT2_DERIVATIVE_PR,
            "agent2_derivative_head": AUDITED_AGENT2_DERIVATIVE_HEAD,
            "independent_agent4_velocity_pr": INDEPENDENT_AGENT4_VELOCITY_PR,
            "independent_agent4_velocity_head": INDEPENDENT_AGENT4_VELOCITY_HEAD,
            "independent_agent4_time_derivative_pr": INDEPENDENT_AGENT4_TIME_DERIVATIVE_PR,
            "independent_agent4_time_derivative_head": INDEPENDENT_AGENT4_TIME_DERIVATIVE_HEAD,
        },
        "state": {
            "radial_count": RADIAL_COUNT,
            "radial_evaluation_count": RADIAL_COUNT - 2,
            "angular_quadrature": "gauss_legendre",
            "angular_order": ANGULAR_ORDER,
            "time": TIME,
            "z": Z,
            "nu": NU,
            "fd6_steps": list(FD6_STEPS),
            "radial_nodes_evaluated": sample_radii.tolist(),
        },
        "agent3_receipt_target": AGENT3_RECEIPT_TARGET,
        "guards_frozen_before_actions": GUARDS,
        "independent_fd6_rows": rows,
        "raw_mean_profile_relative_changes": profile_changes,
        "raw_mean_profile_refinement_ratio": refinement_ratio,
        "finest_vs_agent3_receipt_relative": receipt_relative,
        "quadratic_sign_flip_mutation_relative_rms": sign_flip_response,
        "finest_profile_concentration": profile_concentration,
        "failed_guards": failed,
        "oscillatory_mean_defect_input_preflight_passed": len(failed) == 0,
        "independence": {
            "public_black_box_values_only": True,
            "agent3_fd4_operator_reused": False,
            "agent3_uniform_midpoint_ring_reused": False,
            "agent3_compact_radial_stress_reused": False,
            "training_tensor_or_loss_read": False,
            "pressure_fitted": False,
            "forcing_fitted": False,
            "spatial_operator": "centered_cartesian_fd6",
            "angular_operator": "gauss_legendre_physical_angle_mean",
        },
        "truth_boundary": {
            "oscillatory_component_mean_defect_input_only": True,
            "full_same_cycle_composite_requested_stress_materialized": False,
            "candidate_finite_head_mean_debt_materialized": False,
            "signed_mean_inverse_input_ready": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
            "pressure_assessed": False,
            "restricted_forcing_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "residual_reduction_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "final_momentum_gate": 1.0e-3,
            "final_divergence_gate": 1.0e-5,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
