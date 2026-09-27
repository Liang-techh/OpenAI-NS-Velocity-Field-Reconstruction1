"""Bounded endpoint acceleration projection for the refined 264-control wave.

The balanced refined candidate is replayed at ``delta_k = 1e-6`` on the
44,400-point cache grid.  A new correction has velocity

    0.5 * (t - t0)^2 * V a

and pressure

    (t - t0) * P b,

so both corrections vanish at the reference time.  Degree-2 compact
divergence-free Fourier bases for modes 0 through 4 provide a 324-real-column
linearized endpoint response.  A ridge fit is compared with one actual
finite-difference endpoint replay; the quadratic convection term is reported
as the nonlinear mismatch.  No endpoint moments, cones, shape constraints, or
trajectory acceptance are claimed for the new correction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from broad_meridional_constrained import load_saved_field  # noqa: E402
from full_wave_tangent import LocalPotentialField  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from supported_fourier_analytic_jets import basis_jets  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402
from wave_momentum_projection import _metric  # noqa: E402


ROOT_CACHE = ROOT / "refined_wave_momentum_cache.npz"
BALANCED_PATH = ROOT / "balanced_refined_tangent.json"
FROZEN_PATH = ROOT / "full_wave_frozen_cache.json"
OUTPUT_PATH = ROOT / "endpoint_acceleration_projection.json"

DELTA_K = 1.0e-6
RIDGE = 1.0e-6
DEGREE = 2
MODES = (0, 1, 2, 3, 4)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(report: dict) -> None:
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _packed_complex(value):
    value = np.asarray(value, dtype=complex)
    return np.stack((value.real, value.imag), axis=-1).tolist()


def _unpack_wave(value):
    value = np.asarray(value, dtype=float)
    if value.shape != (27, 2):
        raise ValueError(f"Expected 27 complex wave coefficients, got {value.shape}")
    return value[:, 0] + 1j * value[:, 1]


def _unpack_mixed_affine_control(control):
    """Decode 264 controls: mode0 degree3, mode1 degree2, mode2 degree3."""

    control = np.asarray(control, dtype=float)
    if control.shape != (264,):
        raise ValueError(f"Expected 264 affine controls, got {control.shape}")
    result = {}
    cursor = 0
    for mode, degree in ((0, 3), (1, 2), (2, 3)):
        q = (degree + 1) ** 2
        velocity_count = 3 * q
        block_count = velocity_count + q if mode == 0 else 2 * (velocity_count + q)
        block = control[cursor:cursor + block_count]
        cursor += block_count
        if mode == 0:
            result[mode] = (
                block[:velocity_count].astype(complex),
                block[velocity_count:].astype(complex),
                degree,
            )
        else:
            velocity = block[:2 * velocity_count:2] + 1j * block[1:2 * velocity_count:2]
            pressure_start = 2 * velocity_count
            pressure = (
                block[pressure_start:pressure_start + 2 * q:2]
                + 1j * block[pressure_start + 1:pressure_start + 2 * q:2]
            )
            result[mode] = (velocity, pressure, degree)
    if cursor != len(control):
        raise ValueError(f"Consumed {cursor} of {len(control)} affine controls")
    return result


def _unpack_acceleration_control(control):
    """Decode 324 controls: degree2 velocity acceleration and pressure slope."""

    control = np.asarray(control, dtype=float)
    expected = 36 + 4 * 72
    if control.shape != (expected,):
        raise ValueError(f"Expected {expected} acceleration controls, got {control.shape}")
    result = {}
    cursor = 0
    for mode in MODES:
        q = (DEGREE + 1) ** 2
        velocity_count = 3 * q
        block_count = velocity_count + q if mode == 0 else 2 * (velocity_count + q)
        block = control[cursor:cursor + block_count]
        cursor += block_count
        if mode == 0:
            result[mode] = (
                block[:velocity_count].astype(complex),
                block[velocity_count:].astype(complex),
            )
        else:
            velocity = block[:2 * velocity_count:2] + 1j * block[1:2 * velocity_count:2]
            pressure_start = 2 * velocity_count
            pressure = (
                block[pressure_start:pressure_start + 2 * q:2]
                + 1j * block[pressure_start + 1:pressure_start + 2 * q:2]
            )
            result[mode] = (velocity, pressure)
    return result


class MixedAffineField:
    """Reference mean plus current wave and balanced affine controls."""

    def __init__(self, mean, center, widths, carrier, wave, control, tau0):
        self.mean = mean
        self.center = tuple(float(v) for v in center)
        self.widths = tuple(float(v) for v in widths)
        self.carrier = np.asarray(carrier, dtype=float)
        self.wave = np.asarray(wave, dtype=complex)
        self.affine = _unpack_mixed_affine_control(control)
        self.tau0 = float(tau0)
        self.nu = float(mean.nu)

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        tau = float(np.asarray(tau).ravel()[0])
        velocity, pressure = self.mean.fields(points, tau)
        velocity = velocity.copy()
        pressure = pressure.copy()
        physical_delta = self.tau0 - tau
        for mode in (0, 1, 2):
            derivative, pressure_coeff, degree = self.affine[mode]
            initial = self.wave if mode == 1 else np.zeros_like(derivative)
            carrier = mode * self.carrier
            basis_velocity, basis_pressure, _ = basis_data(
                points, self.center, self.widths, mode, degree, carrier
            )
            coefficient = initial + physical_delta * derivative
            velocity += np.einsum("niq,q->ni", basis_velocity, coefficient).real
            pressure += np.einsum("nq,q->n", basis_pressure, pressure_coeff).real
        return velocity, pressure


class AccelerationCorrectionField:
    """Add zero-at-reference quadratic velocity and linear pressure terms."""

    def __init__(self, base, center, widths, carrier, acceleration, pressure_slope, tau0):
        self.base = base
        self.center = tuple(float(v) for v in center)
        self.widths = tuple(float(v) for v in widths)
        self.carrier = np.asarray(carrier, dtype=float)
        self.acceleration = acceleration
        self.pressure_slope = pressure_slope
        self.tau0 = float(tau0)
        self.nu = float(base.nu)

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        tau = float(np.asarray(tau).ravel()[0])
        velocity, pressure = self.base.fields(points, tau)
        velocity = velocity.copy()
        pressure = pressure.copy()
        physical_delta = self.tau0 - tau
        for mode in MODES:
            acceleration = self.acceleration[mode]
            pressure_slope = self.pressure_slope[mode]
            basis_velocity, basis_pressure, _ = basis_data(
                points, self.center, self.widths, mode, DEGREE, mode * self.carrier
            )
            velocity += 0.5 * physical_delta**2 * np.einsum(
                "niq,q->ni", basis_velocity, acceleration
            ).real
            pressure += physical_delta * np.einsum(
                "nq,q->n", basis_pressure, pressure_slope
            ).real
        return velocity, pressure


def _metric_difference(a, b, weights):
    return _metric(np.asarray(a) - np.asarray(b), weights)


def _max_detail(value, points, weights):
    magnitudes = np.linalg.norm(value, axis=1)
    index = int(np.argmax(magnitudes))
    return {
        "index": index,
        "point": np.asarray(points[index], dtype=float).tolist(),
        "weight": float(weights[index]),
        "vector": np.asarray(value[index], dtype=float).tolist(),
        "norm": float(magnitudes[index]),
    }


def _acceleration_design(points, center, widths, carrier, dt, base_velocity, base_gradient, nu):
    """Build the linearized endpoint residual response for modes 0 through 4."""

    point_count = len(points)
    column_count = 36 + 4 * 72
    design = np.empty((3 * point_count, column_count), dtype=float)
    layout = []
    cursor = 0
    for mode in MODES:
        values, gradients, diffusion, _, pressure_gradient = basis_jets(
            points, center, widths, mode, DEGREE, mode * np.asarray(carrier), nu
        )
        q = (DEGREE + 1) ** 2
        velocity_count = 3 * q
        real_imag = (("real", 1.0),) if mode == 0 else (("real", 1.0), ("imag", -1.0))
        for kind, tensor, count in (
            ("velocity", values, velocity_count),
            ("pressure", pressure_gradient, q),
        ):
            for index in range(count):
                for component, sign in real_imag:
                    tensor_column = tensor[:, :, index] * sign
                    if kind == "velocity":
                        gradient_column = gradients[:, :, :, index] * sign
                        diffusion_column = diffusion[:, :, index] * sign
                        response = (
                            dt * tensor_column
                            + 0.5 * dt**2 * diffusion_column
                            + 0.5 * dt**2 * (
                                np.einsum("nij,nj->ni", gradient_column, base_velocity)
                                + np.einsum("nij,nj->ni", base_gradient, tensor_column)
                            )
                        )
                    else:
                        response = dt * tensor_column
                    # Re(V*c) gives Re(V) for a real control and -Im(V)
                    # for an imaginary control. The sign was applied above;
                    # select the matching component rather than duplicating Re.
                    design[:, cursor] = getattr(response, component).reshape(-1)
                    layout.append({
                        "mode": mode,
                        "kind": kind,
                        "component": component,
                        "index": index,
                        "coefficient_units": "velocity acceleration" if kind == "velocity" else "pressure slope",
                    })
                    cursor += 1
    if cursor != column_count:
        raise ValueError(f"Acceleration layout consumed {cursor} of {column_count}")
    return design, layout


def _ridge_fit(design, residual, weights, ridge=RIDGE):
    row_weight = np.repeat(np.sqrt(np.asarray(weights, dtype=float)), 3)
    weighted_design = design * row_weight[:, None]
    weighted_rhs = -np.asarray(residual, dtype=float).reshape(-1) * row_weight
    scales = np.maximum(np.linalg.norm(weighted_design, axis=0), 1.0e-300)
    normalized = weighted_design / scales[None, :]
    singular = np.linalg.svd(normalized, compute_uv=False)
    gram = normalized.T @ normalized + float(ridge) * np.eye(normalized.shape[1])
    rhs = normalized.T @ weighted_rhs
    scaled = np.linalg.solve(gram, rhs)
    coefficients = scaled / scales
    predicted = residual + (design @ coefficients).reshape(-1, 3)
    return coefficients, predicted, {
        "ridge": float(ridge),
        "column_count": int(design.shape[1]),
        "singular_values": singular.tolist(),
        "rank_rcond_1e-10": int(np.sum(singular > 1.0e-10 * singular[0])),
        "condition_number": float(singular[0] / max(singular[-1], 1.0e-300)),
        "column_scales_min": float(np.min(scales)),
        "column_scales_max": float(np.max(scales)),
    }


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    balanced_raw = BALANCED_PATH.read_bytes()
    frozen_raw = FROZEN_PATH.read_bytes()
    balanced = json.loads(balanced_raw)
    frozen = json.loads(frozen_raw)
    if balanced.get("status") != "completed":
        raise ValueError("balanced_refined_tangent.json is not complete")
    if not balanced.get("assembled_feasible", False):
        raise ValueError("Balanced source does not report an assembled feasible candidate")
    with np.load(ROOT_CACHE, allow_pickle=False) as loaded:
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
    if points.shape != (44400, 3):
        raise ValueError(f"Expected refined cache points (44400,3), got {points.shape}")
    base_control = np.asarray(balanced["selected"]["tangent_coefficients"], dtype=float)
    if base_control.shape != (264,):
        raise ValueError(f"Expected balanced 264 controls, got {base_control.shape}")
    wave = _unpack_wave(balanced["selected"]["coefficients_original"])
    geometry = frozen["inputs"]["wave"]
    center = tuple(float(v) for v in geometry["center"])
    widths = tuple(float(v) for v in geometry["widths"])
    carrier = np.asarray(geometry["carrier"], dtype=float)
    tau0 = float(frozen["inputs"]["mean"]["tau"])
    k0 = float(frozen["inputs"]["mean"]["k"])
    tau1 = 0.5 * 2.0 ** (-(k0 + DELTA_K))
    dt = tau0 - tau1
    hspace = float(frozen["timesteps"]["hspace"])
    htime = float(frozen["timesteps"]["htime"])
    nu = float(frozen["inputs"]["mean"]["nu"])
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "constraints_maintained": False,
        "wave_integrated": False,
        "scope": (
            "Bounded actual endpoint acceleration projection. The balanced 264-control "
            "affine candidate is replayed at delta_k=1e-6; a separate degree-2 mode "
            "0..4 correction has zero velocity and pressure at t0. A ridge linearized "
            "endpoint fit is followed by one actual finite-difference replay. Endpoint "
            "moments, cones, shape, trajectory, and recursion are diagnostic only."
        ),
        "sources": {
            "balanced": {"path": BALANCED_PATH.name, "sha256": _sha256(BALANCED_PATH)},
            "frozen_geometry": {"path": FROZEN_PATH.name, "sha256": _sha256(FROZEN_PATH)},
            "refined_cache": {"path": ROOT_CACHE.name, "sha256": _sha256(ROOT_CACHE)},
        },
        "inputs": {
            "point_count": int(len(points)),
            "delta_k": DELTA_K,
            "k0": k0,
            "tau0": tau0,
            "tau1": tau1,
            "physical_time_increment": dt,
            "hspace": hspace,
            "htime": htime,
            "viscosity": nu,
            "degree": DEGREE,
            "modes": list(MODES),
            "acceleration_control_count": 324,
            "ridge": RIDGE,
            "base_control_count": 264,
            "initial_wave_coefficients": _packed_complex(wave),
            "base_tangent_coefficients": base_control.tolist(),
            "acceleration_layout": "mode0 [27 velocity acceleration, 9 pressure slope] real; modes1..4 each interleaved velocity [Re,-Im] then pressure [Re,-Im]",
        },
        "source_constraints": {
            "reference_assembled_feasible": bool(balanced.get("assembled_feasible", False)),
            "reference_constraints_scope": balanced.get("constraints_maintained_scope"),
            "reference_selected_moment_error": balanced["selected"].get("assembled_moment_max_abs"),
            "reference_selected_cone_min_margin": balanced["selected"].get("assembled_cone_min_margin"),
            "reference_selected_endpoint_margin": balanced["selected"].get("endpoint_margin"),
            "new_acceleration_constraints_enforced": False,
            "new_acceleration_initial_state_unchanged": True,
        },
    }
    _save(report)
    print(json.dumps({"stage": "sources_loaded", "point_count": len(points), "delta_t": dt}), flush=True)

    mean, mean_report = load_saved_field()
    install_in_field(mean)
    base_field = MixedAffineField(mean, center, widths, carrier, wave, base_control, tau0)
    print(json.dumps({"stage": "base_endpoint_replay_started", "point_count": len(points)}), flush=True)
    base_endpoint_jets = jets(base_field, points, tau1, hspace, htime)
    base_endpoint_residual = momentum(base_endpoint_jets)
    base_velocity = base_endpoint_jets[0]
    base_gradient = base_endpoint_jets[1]
    report["base_endpoint"] = {
        "method": "actual MixedAffineField plus affine_momentum.jets five-point FD",
        "metric": _metric(base_endpoint_residual, weights),
        "max_residual": _max_detail(base_endpoint_residual, points, weights),
        "mean_loader_status": mean_report.get("status", "loaded"),
        "reference_affine_candidate_metric": balanced["selected"]["training_momentum"],
    }
    report["status"] = "base_endpoint_replayed"
    _save(report)
    print(json.dumps({"stage": "base_endpoint_replayed", "l2": report["base_endpoint"]["metric"]["momentum_volume_L2"]}), flush=True)

    print(json.dumps({"stage": "acceleration_design_started", "point_count": len(points)}), flush=True)
    design, layout = _acceleration_design(
        points, center, widths, carrier, dt, base_velocity, base_gradient, nu
    )
    report["acceleration_design"] = {
        "shape": list(design.shape),
        "layout": layout,
        "response_definition": (
            "dt*V + 0.5*dt^2*(-nu*lapV + (gradV)u_endpoint + (gradU_endpoint)V) "
            "for velocity acceleration; dt*gradP for pressure slope"
        ),
    }
    report["status"] = "acceleration_design_assembled"
    _save(report)
    print(json.dumps({"stage": "acceleration_design_assembled", "shape": list(design.shape)}), flush=True)

    coefficients, linear_predicted, fit_data = _ridge_fit(
        design, base_endpoint_residual, weights, RIDGE
    )
    report["ridge_fit"] = fit_data
    report["ridge_fit"]["linearized_endpoint_metric"] = _metric(linear_predicted, weights)
    report["ridge_fit"]["linearized_correction_metric"] = _metric(
        linear_predicted - base_endpoint_residual, weights
    )
    report["ridge_fit"]["coefficient_norm"] = float(np.linalg.norm(coefficients))
    report["ridge_fit"]["coefficient_max_abs"] = float(np.max(np.abs(coefficients)))
    report["ridge_fit"]["coefficients"] = coefficients.tolist()
    report["ridge_fit"]["layout"] = layout
    report["status"] = "ridge_fit_complete"
    _save(report)
    print(json.dumps({
        "stage": "ridge_fit_complete",
        "linear_l2": report["ridge_fit"]["linearized_endpoint_metric"]["momentum_volume_L2"],
        "coefficient_norm": report["ridge_fit"]["coefficient_norm"],
    }), flush=True)

    acceleration = _unpack_acceleration_control(coefficients)
    acceleration_velocity = {mode: value[0] for mode, value in acceleration.items()}
    acceleration_pressure = {mode: value[1] for mode, value in acceleration.items()}
    corrected_field = AccelerationCorrectionField(
        base_field, center, widths, carrier,
        acceleration_velocity, acceleration_pressure, tau0,
    )
    # Cheap exact construction check at the reference time: both correction
    # fields are zero algebraically, before the endpoint replay begins.
    sample = points[: min(8, len(points))]
    base_t0 = base_field.fields(sample, tau0)
    corrected_t0 = corrected_field.fields(sample, tau0)
    report["initial_invariance_check"] = {
        "sample_count": int(len(sample)),
        "velocity_max_abs": float(np.max(np.abs(corrected_t0[0] - base_t0[0]))),
        "pressure_max_abs": float(np.max(np.abs(corrected_t0[1] - base_t0[1]))),
        "statement": "correction is exactly zero at tau0 by its delta^2 and delta factors",
    }
    print(json.dumps({"stage": "corrected_endpoint_replay_started", "point_count": len(points)}), flush=True)
    corrected_endpoint_jets = jets(corrected_field, points, tau1, hspace, htime)
    corrected_endpoint_residual = momentum(corrected_endpoint_jets)
    nonlinear_difference = corrected_endpoint_residual - linear_predicted
    report["corrected_endpoint"] = {
        "method": "actual AccelerationCorrectionField plus affine_momentum.jets five-point FD",
        "metric": _metric(corrected_endpoint_residual, weights),
        "max_residual": _max_detail(corrected_endpoint_residual, points, weights),
        "linear_prediction_metric": _metric(linear_predicted, weights),
        "nonlinear_mismatch_metric": _metric(nonlinear_difference, weights),
        "nonlinear_mismatch_relative_to_linear_correction": float(
            _metric(nonlinear_difference, weights)["momentum_volume_L2"]
            / max(_metric(linear_predicted - base_endpoint_residual, weights)["momentum_volume_L2"], 1.0e-300)
        ),
        "quadratic_convection_scope": "difference between actual endpoint replay and linearized response; includes q·grad(q) and FD roundoff/truncation",
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "base_l2": report["base_endpoint"]["metric"]["momentum_volume_L2"],
        "linear_l2": report["ridge_fit"]["linearized_endpoint_metric"]["momentum_volume_L2"],
        "actual_l2": report["corrected_endpoint"]["metric"]["momentum_volume_L2"],
        "nonlinear_l2": report["corrected_endpoint"]["nonlinear_mismatch_metric"]["momentum_volume_L2"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
