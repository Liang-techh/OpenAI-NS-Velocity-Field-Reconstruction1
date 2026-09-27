"""Endpoint shape oracle for the frozen 324-control acceleration layout.

The 7776-point enriched endpoint cache already contains the endpoint velocity
and gradient of the balanced 264-control parent.  This module reconstructs
that parent directly from the cache, then adds the analytic degree-2
mode-0..4 acceleration velocity/gradient response.  Pressure-slope columns
are zero for velocity and shape observables, as required by the acceleration
field construction.

The oracle is a groundwork diagnostic for a later constrained acceleration
fit.  It performs no mean finite-difference replay, optimization, PDE check,
trajectory integration, or adoption decision.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from enriched_endpoint_shape_cache import EnrichedEndpointShape  # noqa: E402
from endpoint_acceleration_projection import _unpack_acceleration_control  # noqa: E402
from supported_fourier_analytic_jets import basis_jets  # noqa: E402
from wave_endpoint_shape_cache import _observables_and_jacobian  # noqa: E402


DEFAULT_CACHE = ROOT / "enriched_endpoint_shape_cache.npz"
DEFAULT_CACHE_REPORT = ROOT / "enriched_endpoint_shape_cache.json"
DEFAULT_BALANCED = ROOT / "balanced_refined_tangent.json"
DEFAULT_ACCELERATION_REPORT = ROOT / "endpoint_acceleration_projection.json"
DEFAULT_ACTUAL_SHAPE_REPORT = ROOT / "endpoint_acceleration_shape.json"
DEFAULT_GEOMETRY = ROOT / "full_wave_frozen_cache.json"
DEFAULT_OUTPUT = ROOT / "acceleration_shape_oracle.json"

DEGREE = 2
MODES = (0, 1, 2, 3, 4)
CONTROL_COUNT = 324


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _acceleration_layout():
    """Return the exact interleaved layout consumed by _unpack_acceleration_control."""

    layout = []
    q = (DEGREE + 1) ** 2
    velocity_count = 3 * q
    for mode in MODES:
        components = ("real",) if mode == 0 else ("real", "imag")
        for kind, count, units in (
            ("velocity", velocity_count, "velocity acceleration"),
            ("pressure", q, "pressure slope"),
        ):
            for index in range(count):
                for component in components:
                    layout.append(
                        {
                            "mode": int(mode),
                            "kind": kind,
                            "component": component,
                            "index": int(index),
                            "coefficient_units": units,
                        }
                    )
    if len(layout) != CONTROL_COUNT:
        raise ValueError(f"Acceleration layout has {len(layout)} columns, expected {CONTROL_COUNT}")
    return layout


def _analytic_acceleration_response(points, geometry, dt, nu):
    """Build endpoint velocity/J response for acceleration controls.

    For velocity acceleration columns the endpoint factor is
    ``0.5 * dt**2``.  Pressure slope columns have exactly zero velocity and
    gradient response, because pressure does not enter endpoint shape.
    """

    points = np.asarray(points, dtype=float)
    center = geometry["center"]
    widths = geometry["widths"]
    carrier = np.asarray(geometry["carrier"], dtype=float)
    factor = 0.5 * float(dt) ** 2
    velocity_blocks = []
    gradient_blocks = []
    for mode in MODES:
        values, gradients, _, _, _ = basis_jets(
            points,
            center,
            widths,
            mode,
            DEGREE,
            mode * carrier,
            nu,
        )
        q = (DEGREE + 1) ** 2
        velocity_count = 3 * q
        if mode == 0:
            components = (("real", 1.0),)
        else:
            components = (("real", 1.0), ("imag", -1.0))
        for index in range(velocity_count):
            for _component, sign in components:
                velocity_blocks.append(factor * sign * values[:, :, index].real)
                gradient_blocks.append(factor * sign * gradients[:, :, :, index].real)
        # Pressure-slope columns have no velocity or gradient response.
        for _index in range(q):
            for _component, _sign in components:
                velocity_blocks.append(np.zeros((len(points), 3), dtype=float))
                gradient_blocks.append(np.zeros((len(points), 3, 3), dtype=float))
    velocity = np.stack(velocity_blocks, axis=-1)
    gradient = np.stack(gradient_blocks, axis=-1)
    if velocity.shape != (len(points), 3, CONTROL_COUNT):
        raise ValueError(f"Unexpected acceleration velocity response shape {velocity.shape}")
    if gradient.shape != (len(points), 3, 3, CONTROL_COUNT):
        raise ValueError(f"Unexpected acceleration gradient response shape {gradient.shape}")
    return velocity, gradient


class AccelerationShapeOracle:
    """Evaluate endpoint shape observables for acceleration controls.

    ``evaluate(control)`` accepts the 324-real-column layout used by
    ``endpoint_acceleration_projection._unpack_acceleration_control`` and
    returns ``(values, jacobian)`` with shapes ``(3,)`` and ``(3,324)``.
    The state at zero acceleration control is the balanced parent endpoint
    reconstructed from ``enriched_endpoint_shape_cache.npz``.
    """

    def __init__(
        self,
        cache_path=DEFAULT_CACHE,
        balanced_path=DEFAULT_BALANCED,
        geometry_path=DEFAULT_GEOMETRY,
    ):
        cache_path = Path(cache_path)
        balanced_path = Path(balanced_path)
        geometry_path = Path(geometry_path)
        enriched = EnrichedEndpointShape(cache_path)
        balanced = json.loads(balanced_path.read_text(encoding="utf-8"))
        geometry = json.loads(geometry_path.read_text(encoding="utf-8"))["inputs"]["wave"]
        parent_control = np.asarray(balanced["selected"]["tangent_coefficients"], dtype=float)
        if parent_control.shape != (264,):
            raise ValueError(f"Expected balanced parent control shape (264,), got {parent_control.shape}")
        if not np.array_equal(np.asarray(enriched.coefficients_original), np.asarray(balanced["selected"]["coefficients_original"], dtype=float)):
            raise ValueError("Endpoint cache and balanced parent use different locked waves")
        if enriched.point_count != 7776 or enriched.control_count != 264:
            raise ValueError("Unexpected enriched endpoint cache dimensions")

        self.cache_path = cache_path
        self.balanced_path = balanced_path
        self.geometry_path = geometry_path
        self.points = np.asarray(enriched.points, dtype=float)
        self.weights = np.asarray(enriched.weights, dtype=float)
        self.point_count = int(len(self.points))
        self.reference = np.asarray(enriched.reference, dtype=float)
        self.coefficients_original = np.asarray(enriched.coefficients_original, dtype=float)
        self.parent_control = parent_control
        self.parent_velocity = enriched.velocity_offset + np.einsum(
            "ncq,q->nc", enriched.velocity_control, parent_control
        )
        self.parent_gradient = enriched.gradient_offset + np.einsum(
            "ncdq,q->ncd", enriched.gradient_control, parent_control
        )
        self.tau0 = float(enriched.tau0)
        self.tau_endpoint = float(enriched.tau_endpoint)
        self.dt = float(enriched.tau0 - enriched.tau_endpoint)
        self.k0 = float(enriched.k0)
        self.k_endpoint = float(enriched.k_endpoint)
        self.nu = float(json.loads(geometry_path.read_text(encoding="utf-8"))["inputs"]["mean"]["nu"])
        self.geometry = geometry
        self.control_count = CONTROL_COUNT
        self.layout = _acceleration_layout()
        self.velocity_control, self.gradient_control = _analytic_acceleration_response(
            self.points, geometry, self.dt, self.nu
        )
        self.parent_values, _, self.parent_detail = _observables_and_jacobian(
            self.parent_velocity, self.parent_gradient, self.weights, self.points
        )

    def evaluate(self, control):
        """Return endpoint shape values and the exact 3 by 324 Jacobian."""

        control = np.asarray(control, dtype=float)
        if control.shape != (CONTROL_COUNT,):
            raise ValueError(f"control must have shape ({CONTROL_COUNT},), got {control.shape}")
        velocity = self.parent_velocity + np.einsum(
            "ncq,q->nc", self.velocity_control, control
        )
        gradient = self.parent_gradient + np.einsum(
            "ncdq,q->ncd", self.gradient_control, control
        )
        values, jacobian, detail = _observables_and_jacobian(
            velocity,
            gradient,
            self.weights,
            self.points,
            self.velocity_control,
            self.gradient_control,
        )
        return values, jacobian


# A concise alias for downstream fit code.
AccelerationShape = AccelerationShapeOracle


def _fractional_change(values, reference):
    values = np.asarray(values, dtype=float)
    reference = np.asarray(reference, dtype=float)
    signs = np.array((-1.0, 1.0, np.sign(reference[2])), dtype=float)
    return signs * (values - reference) / np.maximum(np.abs(reference), 1.0e-30)


def run(output_path=DEFAULT_OUTPUT):
    started = time.perf_counter()
    output_path = Path(output_path)
    acceleration_report = json.loads(DEFAULT_ACCELERATION_REPORT.read_text(encoding="utf-8"))
    actual_shape_report = json.loads(DEFAULT_ACTUAL_SHAPE_REPORT.read_text(encoding="utf-8"))
    if actual_shape_report.get("status") != "completed":
        raise ValueError("Frozen actual acceleration shape report is not complete")
    oracle = AccelerationShapeOracle()
    ridge_control = np.asarray(acceleration_report["ridge_fit"]["coefficients"], dtype=float)
    if ridge_control.shape != (CONTROL_COUNT,):
        raise ValueError("Frozen acceleration projection does not contain 324 coefficients")

    zero_control = np.zeros(CONTROL_COUNT, dtype=float)
    parent_values, _, parent_detail = _observables_and_jacobian(
        oracle.parent_velocity, oracle.parent_gradient, oracle.weights, oracle.points
    )
    corrected_values, corrected_jacobian = oracle.evaluate(ridge_control)

    # The grouped actual replay is read as frozen evidence only.  It uses the
    # existing five-point FD velocity gradient, while this oracle uses the
    # requested analytic acceleration gradient.  Keeping both rows explicit
    # prevents an optimizer from silently mixing the two precision conventions.
    actual_reference_row = actual_shape_report["rows"]["balanced_reference"]
    actual_corrected_row = actual_shape_report["rows"]["corrected_endpoint"]
    actual_reference_values = np.array(
        (
            actual_reference_row["enstrophy_radial_rms"],
            actual_reference_row["enstrophy_aspect_ratio"],
            actual_reference_row["enstrophy_weighted_angular_speed"],
        ),
        dtype=float,
    )
    actual_corrected_values = np.array(
        (
            actual_corrected_row["enstrophy_radial_rms"],
            actual_corrected_row["enstrophy_aspect_ratio"],
            actual_corrected_row["enstrophy_weighted_angular_speed"],
        ),
        dtype=float,
    )

    # Use a scale-aware direction around the frozen ridge point so the tiny
    # dt^2 acceleration response is resolved by centered finite differences.
    direction = np.linspace(-1.0, 1.0, CONTROL_COUNT) * np.maximum(np.abs(ridge_control), 1.0)
    epsilon = 1.0e-6
    plus_values, _ = oracle.evaluate(ridge_control + epsilon * direction)
    minus_values, _ = oracle.evaluate(ridge_control - epsilon * direction)
    finite_difference = (plus_values - minus_values) / (2.0 * epsilon)
    analytic_directional = corrected_jacobian @ direction
    error = finite_difference - analytic_directional
    relative_error = float(np.max(np.abs(error)) / max(np.max(np.abs(finite_difference)), 1.0e-30))

    decoded = _unpack_acceleration_control(zero_control)
    pressure_columns = [
        index for index, entry in enumerate(oracle.layout) if entry["kind"] == "pressure"
    ]
    pressure_velocity_norm = float(np.max(np.abs(oracle.velocity_control[:, :, pressure_columns])))
    pressure_gradient_norm = float(np.max(np.abs(oracle.gradient_control[:, :, :, pressure_columns])))
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": (
            "Endpoint acceleration shape oracle only. The balanced parent endpoint "
            "u,J is reconstructed from the frozen enriched endpoint cache; degree-2 "
            "mode-0..4 acceleration V,J responses are analytic and pressure-slope "
            "columns are zero. No mean replay, acceleration fit, trajectory, PDE, "
            "or adoption claim is made."
        ),
        "sources": {
            "enriched_endpoint_cache": {"path": oracle.cache_path.name, "sha256": _sha(oracle.cache_path)},
            "enriched_endpoint_cache_report": {"path": DEFAULT_CACHE_REPORT.name, "sha256": _sha(DEFAULT_CACHE_REPORT)},
            "balanced_parent": {"path": oracle.balanced_path.name, "sha256": _sha(oracle.balanced_path)},
            "acceleration_projection": {"path": DEFAULT_ACCELERATION_REPORT.name, "sha256": _sha(DEFAULT_ACCELERATION_REPORT)},
            "actual_shape_replay": {"path": DEFAULT_ACTUAL_SHAPE_REPORT.name, "sha256": _sha(DEFAULT_ACTUAL_SHAPE_REPORT)},
            "frozen_geometry": {"path": oracle.geometry_path.name, "sha256": _sha(oracle.geometry_path)},
        },
        "inputs": {
            "point_count": int(oracle.point_count),
            "control_count": CONTROL_COUNT,
            "degree": DEGREE,
            "modes": list(MODES),
            "k0": oracle.k0,
            "k_endpoint": oracle.k_endpoint,
            "tau0": oracle.tau0,
            "tau_endpoint": oracle.tau_endpoint,
            "physical_time_increment": oracle.dt,
            "viscosity": oracle.nu,
            "layout": oracle.layout,
            "layout_source": "endpoint_acceleration_projection._unpack_acceleration_control",
            "pressure_column_count": int(len(pressure_columns)),
            "velocity_response_factor": float(0.5 * oracle.dt**2),
            "pressure_response_factor": float(oracle.dt),
        },
        "parent_reconstruction": {
            "method": "enriched_endpoint_shape_cache velocity_offset/gradient_offset plus balanced 264-control response",
            "parent_control_norm": float(np.linalg.norm(oracle.parent_control)),
            "parent_values": parent_values.tolist(),
            "parent_detail": {
                "energy": float(parent_detail["energy"]),
                "zmean": float(parent_detail["zmean"]),
                "radial_rms": float(parent_detail["radial_rms"]),
                "axial_rms": float(parent_detail["axial_rms"]),
                "spin": float(parent_detail["spin"]),
            },
            "full_mean_replay_performed": False,
        },
        "reference_endpoint": {
            "observables": oracle.reference.tolist(),
            "signed_fractional_parent_change": _fractional_change(parent_values, oracle.reference).tolist(),
        },
        "corrected_endpoint": {
            "control_source": DEFAULT_ACCELERATION_REPORT.name,
            "control_norm": float(np.linalg.norm(ridge_control)),
            "observables": corrected_values.tolist(),
            "signed_fractional_change": _fractional_change(corrected_values, oracle.reference).tolist(),
            "jacobian_shape": list(corrected_jacobian.shape),
        },
        "analytic_vs_frozen_actual_shape": {
            "frozen_source": DEFAULT_ACTUAL_SHAPE_REPORT.name,
            "frozen_method": "grouped actual acceleration shape replay with five-point FD velocity gradient",
            "analytic_oracle_reference": oracle.reference.tolist(),
            "analytic_oracle_corrected": corrected_values.tolist(),
            "frozen_actual_reference": actual_reference_values.tolist(),
            "frozen_actual_corrected": actual_corrected_values.tolist(),
            "analytic_raw_fractional_change": ((corrected_values - oracle.reference) / np.maximum(np.abs(oracle.reference), 1.0e-30)).tolist(),
            "frozen_actual_raw_fractional_change": ((actual_corrected_values - actual_reference_values) / np.maximum(np.abs(actual_reference_values), 1.0e-30)).tolist(),
            "analytic_minus_frozen_corrected": (corrected_values - actual_corrected_values).tolist(),
            "precision_note": (
                "The analytic response is the optimizer-facing oracle. The frozen "
                "actual row is retained for calibration; its five-point FD gradient "
                "changes the aspect response at this endpoint."
            ),
        },
        "response_checks": {
            "pressure_columns": pressure_columns,
            "pressure_velocity_max_abs": pressure_velocity_norm,
            "pressure_gradient_max_abs": pressure_gradient_norm,
            "pressure_columns_exactly_zero": bool(pressure_velocity_norm == 0.0 and pressure_gradient_norm == 0.0),
            "layout_decode_zero_control_modes": sorted(int(mode) for mode in decoded),
            "velocity_response_shape": list(oracle.velocity_control.shape),
            "gradient_response_shape": list(oracle.gradient_control.shape),
        },
        "directional_gradient_check": {
            "base_control_source": DEFAULT_ACCELERATION_REPORT.name,
            "epsilon": epsilon,
            "direction_norm": float(np.linalg.norm(direction)),
            "analytic_directional": analytic_directional.tolist(),
            "finite_difference_directional": finite_difference.tolist(),
            "max_abs_error": float(np.max(np.abs(error))),
            "relative_error": relative_error,
        },
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "point_count": oracle.point_count,
        "parent_values": parent_values.tolist(),
        "corrected_values": corrected_values.tolist(),
        "directional_relative_error": relative_error,
        "pressure_columns_zero": report["response_checks"]["pressure_columns_exactly_zero"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output)
