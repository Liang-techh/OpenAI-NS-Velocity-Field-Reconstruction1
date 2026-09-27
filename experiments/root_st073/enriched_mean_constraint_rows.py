"""Degree-3 mode-0 constraint rows for the frozen mean/wave geometry.

The degree-3 mode-0 family has 16 scalar potential functions, hence 48 real
velocity-derivative controls and 16 pressure controls.  This module assembles
the integrated moment rows, sampled cone rows, and fixed-cylinder shape rows
without solving a fit.  The existing degree-2 rows are embedded by an
explicit multi-index map and compared numerically before the report is
frozen.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from broad_shear_dynamic_control import load_saved_field as load_dynamic
from full_wave_tangent import _basis_columns
from grouped_joined_field import install_in_field
from joined_field import coordinates
from supported_fourier_analytic_jets import basis_jets
from supported_fourier_basis import basis_data
from vortex_state_observables import fixed_cylinder
from wave_mean_cone_projection import _cone_geometry
from wave_shape_tangent_rows import (
    _build_field,
    _curl,
    _metric_rate,
    _state_with_fd_jets,
)


ROOT = Path(__file__).resolve().parent
MEAN_PATH = ROOT / "broad_meridional_constrained.json"
DYNAMIC_PATH = ROOT / "broad_shear_dynamic_control.json"
CANDIDATE_PATH = ROOT / "wave_moment_cone_tangent.json"
SNAPSHOT_PATH = ROOT / "full_wave_frozen_cache.json"
MOMENT_REFERENCE_PATH = ROOT / "wave_dynamics_moment_compatibility.json"
CONE_REFERENCE_PATH = ROOT / "wave_mean_cone_projection.json"
SHAPE_REFERENCE_PATH = ROOT / "wave_shape_tangent_rows.json"
OUTPUT_PATH = ROOT / "enriched_mean_constraint_rows.json"

DEGREE = 3
OLD_DEGREE = 2
Q = (DEGREE + 1) ** 2
OLD_Q = (OLD_DEGREE + 1) ** 2
MODE0_CONTROL_COUNT = 4 * Q
OLD_MODE0_CONTROL_COUNT = 4 * OLD_Q
MOMENT_ORDER = 96
MOMENT_ZFACTOR = 0.002
CONE_ORDER = 64
NU = 0.01


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(report: dict, path: Path) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _degree2_to_degree3_map():
    """Map old mode-0 columns to degree-3 columns by radial/axial powers."""
    mapping = []
    for component in range(3):
        for old_index in range(OLD_Q):
            radial_power, axial_power = divmod(old_index, OLD_DEGREE + 1)
            mapping.append(component * Q + radial_power * (DEGREE + 1) + axial_power)
    for old_index in range(OLD_Q):
        radial_power, axial_power = divmod(old_index, OLD_DEGREE + 1)
        mapping.append(3 * Q + radial_power * (DEGREE + 1) + axial_power)
    return np.asarray(mapping, dtype=int)


def _row_layout():
    layout = []
    for index in range(3 * Q):
        radial_power, axial_power = divmod(index % Q, DEGREE + 1)
        layout.append(
            {
                "mode": 0,
                "kind": "velocity",
                "component": "real",
                "index": index,
                "component_block": index // Q,
                "radial_power": radial_power,
                "axial_power": axial_power,
            }
        )
    for index in range(Q):
        radial_power, axial_power = divmod(index, DEGREE + 1)
        layout.append(
            {
                "mode": 0,
                "kind": "pressure",
                "component": "real",
                "index": index,
                "radial_power": radial_power,
                "axial_power": axial_power,
            }
        )
    return layout


def _radial_integrals(inner, radius_outer, z, tau, breaks, order, center, widths):
    """Mode-0 radial integrals, generalized to the degree-3 basis."""
    x_max = float(inner.p.X_max)
    q_value = float(
        coordinates(0.0, z / np.sqrt(inner.nu), tau, inner.h)["q"]
    )
    radius_inner = np.sqrt(2.0 * inner.nu * q_value * x_max)
    edges = [0.0, radius_inner, radius_outer]
    edges.extend(
        radius_inner + (radius_outer - radius_inner) * float(value)
        for value in breaks
    )
    edges = sorted(set(np.clip(edges, 0.0, radius_outer)))
    nodes, node_weights = leggauss(int(order))
    radii = np.concatenate(
        [
            0.5 * (lo + hi) + 0.5 * (hi - lo) * nodes
            for lo, hi in zip(edges[:-1], edges[1:])
        ]
    )
    weights = np.concatenate(
        [
            0.5 * (hi - lo) * node_weights
            for lo, hi in zip(edges[:-1], edges[1:])
        ]
    )
    points = np.column_stack(
        (radii, np.zeros_like(radii), np.full_like(radii, z))
    )
    velocity, pressure, _ = basis_data(
        points, center, widths, 0, DEGREE, (0.0, 0.0)
    )
    velocity = velocity.real
    pressure = pressure.real
    return {
        "velocity_theta": np.einsum(
            "n,n,nq->q", weights, radii**2, velocity[:, 1, :]
        ),
        "velocity_axial": np.einsum(
            "n,n,nq->q", weights, radii, velocity[:, 2, :]
        ),
        "pressure": np.einsum("n,n,nq->q", weights, radii, pressure),
    }


def _mode0_moment_rows(inner, tau, breaks, center, widths):
    """Return order-96 integrated moment rows with 64 degree-3 controls."""
    htime = 1.0e-4 * tau
    response = np.zeros((4, MODE0_CONTROL_COUNT), dtype=float)
    metadata = []
    for row_index, eta in enumerate((-0.2, 0.2)):
        outer = inner.from_similarity(
            inner.p.X_max * np.array([16.0**2]), [eta], tau
        )[0]
        radius_outer, z = float(outer[0]), float(outer[2])
        zp = inner.from_similarity(
            inner.p.X_max * np.array([1.0]), [eta + 0.01], tau
        )[0, 2]
        zm = inner.from_similarity(
            inner.p.X_max * np.array([1.0]), [eta - 0.01], tau
        )[0, 2]
        hz = MOMENT_ZFACTOR * abs(zp - zm) / 0.02
        tau_offsets = (-2, -1, 1, 2)
        theta_states = []
        axial_states = []
        for offset in tau_offsets:
            tau_j = tau + offset * htime
            state = _radial_integrals(
                inner,
                radius_outer,
                z,
                tau_j,
                breaks,
                MOMENT_ORDER,
                center,
                widths,
            )
            coefficient = -offset * htime
            theta_states.append(coefficient * state["velocity_theta"])
            axial_states.append(coefficient * state["velocity_axial"])
        dtau_theta = (
            theta_states[0]
            - 8.0 * theta_states[1]
            + 8.0 * theta_states[2]
            - theta_states[3]
        ) / (12.0 * htime)
        dtau_axial = (
            axial_states[0]
            - 8.0 * axial_states[1]
            + 8.0 * axial_states[2]
            - axial_states[3]
        ) / (12.0 * htime)
        z_states = []
        for offset in (-2, -1, 1, 2):
            state = _radial_integrals(
                inner,
                radius_outer,
                z + offset * hz,
                tau,
                breaks,
                MOMENT_ORDER,
                center,
                widths,
            )
            z_states.append(state["pressure"])
        dz_pressure = (
            z_states[0]
            - 8.0 * z_states[1]
            + 8.0 * z_states[2]
            - z_states[3]
        ) / (12.0 * hz)
        response[row_index * 2, : 3 * Q] = dtau_theta / radius_outer**2
        response[row_index * 2 + 1, : 3 * Q] = dtau_axial / radius_outer
        response[row_index * 2 + 1, 3 * Q :] = -dz_pressure / radius_outer
        metadata.append(
            {
                "eta": float(eta),
                "outer_radius": radius_outer,
                "z": z,
                "z_step": hz,
                "time_step": htime,
                "row_order": ["theta", "axial"],
            }
        )
    return response, metadata


def _mode0_cone_rows(saved, inner, k, breaks, center, widths):
    points, panels, H = _cone_geometry(saved, inner, k, breaks)
    target_response = np.zeros((len(panels), 2, MODE0_CONTROL_COUNT), dtype=float)
    for index, (sl, radii, weights, radius) in enumerate(panels):
        velocity, _, pressure_gradient = basis_data(
            points[sl], center, widths, 0, DEGREE, (0.0, 0.0)
        )
        response = np.concatenate((velocity.real, pressure_gradient.real), axis=2)
        target_response[index, 0] = -np.einsum(
            "n,n,nq->q", weights, radii**2 / radius**2, response[:, 1, :]
        )
        target_response[index, 1] = -np.einsum(
            "n,n,nq->q", weights, radii / radius, response[:, 2, :]
        )
    transformed = np.einsum("rab,rbi->rai", H, target_response)
    rows = transformed.reshape(len(panels) * 3, MODE0_CONTROL_COUNT)
    return rows, target_response, H, points, panels


def _mode0_shape_rows(candidate, snapshot, geometry, shape_reference):
    """Build fixed-cylinder shape rows using analytic degree-3 mode-0 jets."""
    field, base, dynamic_report, _, _, _ = _build_field(candidate, snapshot)
    k0 = float(dynamic_report["k"])
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    points, weights, domain = fixed_cylinder(base.inner, k0, angles=12)
    domain = dict(domain)
    domain["point_count"] = int(len(points))
    velocity, omega, density, angular_speed, reference = _state_with_fd_jets(
        field, points, weights, tau0
    )
    reference.update(
        {"density": density, "omega": omega, "angular_speed": angular_speed}
    )
    values, gradients, _, _, _ = basis_jets(
        points,
        geometry["center"],
        geometry["widths"],
        0,
        DEGREE,
        (0.0, 0.0),
        NU,
    )
    curls = _curl(gradients)
    blocks = []
    for index in range(3 * Q):
        blocks.append(
            _metric_rate(
                values[:, :, index].real,
                curls[:, :, index].real,
                points,
                weights,
                reference,
                np.log(2.0) * tau0,
            )
        )
    blocks.extend(np.zeros(3, dtype=float) for _ in range(Q))
    rows = np.asarray(blocks, dtype=float).T
    return rows, domain, reference


def _mapped_error(new_rows, old_rows, mapping):
    mapped = np.asarray(new_rows)[:, mapping]
    old = np.asarray(old_rows)
    return float(np.max(np.abs(mapped - old))), mapped


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    saved = json.loads(MEAN_PATH.read_text(encoding="utf-8"))
    moment_reference = json.loads(MOMENT_REFERENCE_PATH.read_text(encoding="utf-8"))
    cone_reference = json.loads(CONE_REFERENCE_PATH.read_text(encoding="utf-8"))
    shape_reference = json.loads(SHAPE_REFERENCE_PATH.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE_PATH.read_text(encoding="utf-8"))
    snapshot = json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))
    dynamic, dynamic_report = load_dynamic(DYNAMIC_PATH)
    install_in_field(dynamic)
    inner = dynamic.inner
    k = float(saved["k"])
    tau = float(saved["tau"])
    breaks = [float(value) for value in saved["quadrature"]["radial_split_breaks"]]
    geometry = cone_reference["wave"]
    center = geometry["center"]
    widths = geometry["widths"]

    report = {
        "status": "running",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "degree": DEGREE,
        "mode0_control_count": MODE0_CONTROL_COUNT,
        "reference_k": k,
        "tau": tau,
        "nu": float(saved["nu"]),
        "physical_time": -tau,
        "source_hashes": {
            "candidate": _sha256(CANDIDATE_PATH),
            "mean": _sha256(MEAN_PATH),
            "dynamic": _sha256(DYNAMIC_PATH),
            "moment_reference": _sha256(MOMENT_REFERENCE_PATH),
            "cone_reference": _sha256(CONE_REFERENCE_PATH),
            "shape_reference": _sha256(SHAPE_REFERENCE_PATH),
        },
        "scope": (
            "Degree-3 mode-0 instantaneous linear constraint rows at the saved "
            "geometry. No optimization, trajectory, PDE, or scale-recursion claim."
        ),
    }
    _save(report, output_path)

    moment_rows, moment_geometry = _mode0_moment_rows(
        inner, tau, breaks, center, widths
    )
    _save({**report, "stage": "moment_rows_assembled"}, output_path)
    cone_rows, cone_targets, cone_H, cone_points, cone_panels = _mode0_cone_rows(
        saved, inner, k, breaks, center, widths
    )
    _save({**report, "stage": "cone_rows_assembled"}, output_path)
    shape_rows, shape_domain, shape_reference_state = _mode0_shape_rows(
        candidate, snapshot, geometry, shape_reference
    )

    mapping = _degree2_to_degree3_map()
    old_moment_rows = np.asarray(
        moment_reference["reusable_moment_linearization"]["moment_rows"],
        dtype=float,
    )
    old_cone_rows = np.asarray(cone_reference["cone_control_rows"], dtype=float)
    old_shape_rows = np.asarray(shape_reference["shape_control_rows"], dtype=float)
    moment_error, mapped_moment = _mapped_error(
        moment_rows, old_moment_rows[:, :OLD_MODE0_CONTROL_COUNT], mapping
    )
    cone_error, mapped_cone = _mapped_error(cone_rows, old_cone_rows[:, :OLD_MODE0_CONTROL_COUNT], mapping)
    shape_error, mapped_shape = _mapped_error(shape_rows, old_shape_rows[:, :OLD_MODE0_CONTROL_COUNT], mapping)

    old_tangent = np.asarray(
        candidate["selected"]["tangent_coefficients"], dtype=float
    )[:OLD_MODE0_CONTROL_COUNT]
    expanded_tangent = np.zeros(MODE0_CONTROL_COUNT, dtype=float)
    expanded_tangent[mapping] = old_tangent
    selected_response_errors = {
        "moment": float(
            np.max(
                np.abs(
                    moment_rows @ expanded_tangent
                    - old_moment_rows[:, :OLD_MODE0_CONTROL_COUNT] @ old_tangent
                )
            )
        ),
        "cone": float(
            np.max(np.abs(cone_rows @ expanded_tangent - old_cone_rows[:, :OLD_MODE0_CONTROL_COUNT] @ old_tangent))
        ),
        "shape": float(
            np.max(np.abs(shape_rows @ expanded_tangent - old_shape_rows[:, :OLD_MODE0_CONTROL_COUNT] @ old_tangent))
        ),
    }

    mlin = moment_reference["reusable_moment_linearization"]
    report.update(
        {
            "geometry": {
                "radial_split_breaks": breaks,
                "moment_order": MOMENT_ORDER,
                "moment_zfactor": MOMENT_ZFACTOR,
                "cone_order": CONE_ORDER,
                "cone_location_count": len(cone_panels),
                "cone_inequality_count": int(len(cone_rows)),
                "cone_grouped_point_count": int(len(cone_points)),
                "shape_domain": shape_domain,
                "shape_point_count": int(shape_domain["point_count"]),
            },
            "wave": {
                "mode": geometry["mode"],
                "degree": geometry["degree"],
                "center": center,
                "widths": widths,
                "carrier": geometry["carrier"],
            },
            "row_layout": _row_layout(),
            "mode0_layout": _row_layout(),
            "degree2_to_degree3_column_map": mapping.tolist(),
            "moment_control_rows": moment_rows.tolist(),
            "moment_rows": moment_rows.tolist(),
            "cone_control_rows": cone_rows.tolist(),
            "shape_control_rows": shape_rows.tolist(),
            "baseline_moments": mlin["baseline_moments"],
            "cone_baseline": cone_reference["cone_baseline"],
            "cone_lower": cone_reference["cone_lower"],
            "shape_baseline": shape_reference["shape_baseline"],
            "wave_moment_forms": mlin["wave_moment_forms"],
            "cone_wave_forms": cone_reference["cone_wave_forms"],
            "validation": {
                "embedded_degree2_moment_max_abs": moment_error,
                "embedded_degree2_cone_max_abs": cone_error,
                "embedded_degree2_shape_max_abs": shape_error,
                "embedded_degree2_selected_response_max_abs": selected_response_errors,
                "embedded_degree2_moment_rows": mapped_moment.tolist(),
                "embedded_degree2_cone_rows": mapped_cone.tolist(),
                "embedded_degree2_shape_rows": mapped_shape.tolist(),
                "shape_reference_observables": shape_reference["reference_observables"],
                "shape_reference_state_replay": {
                    key: float(value)
                    for key, value in shape_reference_state.items()
                    if key not in {"density", "omega", "angular_speed"}
                },
                "shape_reference_observable_errors": {
                    name: float(
                        shape_reference_state[name]
                        - shape_reference["reference_observables_detail"][name]
                    )
                    for name in (
                        "enstrophy_radial_rms",
                        "enstrophy_aspect_ratio",
                        "enstrophy_weighted_angular_speed",
                    )
                },
            },
            "baseline_sources": {
                "moment": MOMENT_REFERENCE_PATH.name,
                "cone": CONE_REFERENCE_PATH.name,
                "shape": SHAPE_REFERENCE_PATH.name,
                "note": "Baseline vectors and wave quadratic forms are copied unchanged from frozen degree-2 reports.",
            },
        }
    )
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report, output_path)
    print(
        json.dumps(
            {
                "status": report["status"],
                "moment_rows": list(moment_rows.shape),
                "cone_rows": list(cone_rows.shape),
                "shape_rows": list(shape_rows.shape),
                "embedded_degree2_errors": {
                    "moment": moment_error,
                    "cone": cone_error,
                    "shape": shape_error,
                },
                "elapsed_seconds": report["elapsed_seconds"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
