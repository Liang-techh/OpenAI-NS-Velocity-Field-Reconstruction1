"""Constraint and endpoint-shape rows for the localized degree-2 patch.

The local patch uses the same mode-0 moment and cone row definitions as the
frozen mean compatibility machinery, evaluated at a new compact support.  Its
36 mode-0 columns can affect the instantaneous mean moments and cones.  The
144 nonzero-angular columns have zero instantaneous mean rows.  Endpoint shape
uses all 180 local velocity columns with the cached physical time increment;
pressure columns are exactly zero for endpoint velocity and gradient.

This module assembles rows only.  It performs no optimization, mean replay,
finite-time evolution, PDE check, or acceptance decision.
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

from broad_shear_dynamic_control import load_saved_field as load_dynamic_field  # noqa: E402
from enriched_endpoint_shape_cache import EnrichedEndpointShape  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from supported_fourier_analytic_jets import basis_jets  # noqa: E402
from wave_dynamics_mean_compatibility import _mode0_moment_rows  # noqa: E402
from wave_endpoint_shape_cache import _observables_and_jacobian  # noqa: E402
from wave_mean_cone_projection import _cone_geometry, _mode0_rows  # noqa: E402


MEAN_PATH = ROOT / "broad_meridional_constrained.json"
DYNAMIC_PATH = ROOT / "broad_shear_dynamic_control.json"
MOMENT_REFERENCE_PATH = ROOT / "wave_dynamics_mean_compatibility.json"
CONE_REFERENCE_PATH = ROOT / "wave_mean_cone_projection.json"
ENRICHED_ROWS_PATH = ROOT / "enriched_mean_constraint_rows.json"
ENDPOINT_CACHE_PATH = ROOT / "enriched_endpoint_shape_cache.npz"
ENDPOINT_CACHE_REPORT_PATH = ROOT / "enriched_endpoint_shape_cache.json"
CANDIDATE_PATH = ROOT / "enriched_mean_endpoint_tangent.json"
GEOMETRY_PATH = ROOT / "full_wave_frozen_cache.json"
OUTPUT_PATH = ROOT / "localized_constraint_rows.json"

PATCH_CENTER = np.asarray((9.4e-4, -5.49e-5), dtype=float)
PATCH_WIDTHS = np.asarray((2.5e-4, 2.0e-4), dtype=float)
DEGREE = 2
MODES = (0, 1, 2)
NU = 0.01


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _layout():
    q = (DEGREE + 1) ** 2
    result = []
    for mode, components in ((0, ("real",)), (1, ("real", "imag")), (2, ("real", "imag"))):
        for kind, count in (("velocity", 3 * q), ("pressure", q)):
            for index in range(count):
                for component in components:
                    result.append({
                        "mode": mode,
                        "kind": kind,
                        "component": component,
                        "index": index,
                    })
    if len(result) != 180:
        raise ValueError(f"Unexpected localized layout length {len(result)}")
    return result


def _local_endpoint_response(points, center, widths, carrier, dtau):
    """Return analytic endpoint velocity/gradient response for 180 controls."""

    return _local_endpoint_response_with_carrier(points, center, widths, carrier, dtau)


def _local_endpoint_response_with_carrier(points, center, widths, carrier, dtau):
    points = np.asarray(points, dtype=float)
    carriers = {
        0: np.zeros(2, dtype=float),
        1: np.asarray(carrier, dtype=float),
        2: 2.0 * np.asarray(carrier, dtype=float),
    }
    blocks_v = []
    blocks_j = []
    layout = []
    precision = []
    q = (DEGREE + 1) ** 2
    for mode in MODES:
        values, gradients, _, _, _ = basis_jets(
            points, center, widths, mode, DEGREE, carriers[mode], NU
        )
        if mode == 0:
            for index in range(3 * q):
                blocks_v.append(float(dtau) * values[:, :, index].real)
                blocks_j.append(float(dtau) * gradients[:, :, :, index].real)
                layout.append((mode, "velocity", "real", index))
                precision.append("analytic_basis_jets_times_dtau")
            for index in range(q):
                blocks_v.append(np.zeros((len(points), 3), dtype=float))
                blocks_j.append(np.zeros((len(points), 3, 3), dtype=float))
                layout.append((mode, "pressure", "real", index))
                precision.append("zero_velocity_pressure")
        else:
            for index in range(3 * q):
                blocks_v.append(float(dtau) * values[:, :, index].real)
                blocks_j.append(float(dtau) * gradients[:, :, :, index].real)
                layout.append((mode, "velocity", "real", index))
                precision.append("analytic_basis_jets_times_dtau")
                blocks_v.append(-float(dtau) * values[:, :, index].imag)
                blocks_j.append(-float(dtau) * gradients[:, :, :, index].imag)
                layout.append((mode, "velocity", "imag", index))
                precision.append("analytic_basis_jets_times_dtau")
            for index in range(q):
                blocks_v.append(np.zeros((len(points), 3), dtype=float))
                blocks_j.append(np.zeros((len(points), 3, 3), dtype=float))
                layout.append((mode, "pressure", "real", index))
                precision.append("zero_velocity_pressure")
                blocks_v.append(np.zeros((len(points), 3), dtype=float))
                blocks_j.append(np.zeros((len(points), 3, 3), dtype=float))
                layout.append((mode, "pressure", "imag", index))
                precision.append("zero_velocity_pressure")
    velocity = np.stack(blocks_v, axis=-1)
    gradient = np.stack(blocks_j, axis=-1)
    if velocity.shape != (len(points), 3, 180):
        raise ValueError(f"Unexpected local endpoint velocity shape {velocity.shape}")
    if gradient.shape != (len(points), 3, 3, 180):
        raise ValueError(f"Unexpected local endpoint gradient shape {gradient.shape}")
    return velocity, gradient, layout, precision


def _support_report(geometry):
    original_center = np.asarray(geometry["center"], dtype=float)
    original_widths = np.asarray(geometry["widths"], dtype=float)
    original_lower = original_center - original_widths
    original_upper = original_center + original_widths
    patch_lower = PATCH_CENTER - PATCH_WIDTHS
    patch_upper = PATCH_CENTER + PATCH_WIDTHS
    lower_margin = patch_lower - original_lower
    upper_margin = original_upper - patch_upper
    strict_inside = bool(np.all(lower_margin > 0.0) and np.all(upper_margin > 0.0))
    return {
        "coordinates": ["radius", "z"],
        "patch_center": PATCH_CENTER.tolist(),
        "patch_widths": PATCH_WIDTHS.tolist(),
        "patch_bounds": {"lower": patch_lower.tolist(), "upper": patch_upper.tolist()},
        "original_wave_center": original_center.tolist(),
        "original_wave_widths": original_widths.tolist(),
        "original_wave_bounds": {"lower": original_lower.tolist(), "upper": original_upper.tolist()},
        "strictly_inside_original_support": strict_inside,
        "patch_outside_original_support": bool(not strict_inside),
        "added_support_outside_original_support": bool(not strict_inside),
        "lower_containment_margin": lower_margin.tolist(),
        "upper_containment_margin": upper_margin.tolist(),
    }


def _shape_fd_check(parent_velocity, parent_gradient, local_velocity, local_gradient,
                    weights, points, jacobian):
    direction = np.linspace(-1.0, 1.0, local_velocity.shape[-1], dtype=float)
    direction /= np.linalg.norm(direction)
    fd_step = 1.0e6
    delta_velocity = np.einsum("ncq,q->nc", local_velocity, direction)
    delta_gradient = np.einsum("ncdq,q->ncd", local_gradient, direction)
    plus_values, _, _ = _observables_and_jacobian(
        parent_velocity + fd_step * delta_velocity,
        parent_gradient + fd_step * delta_gradient,
        weights,
        points,
    )
    minus_values, _, _ = _observables_and_jacobian(
        parent_velocity - fd_step * delta_velocity,
        parent_gradient - fd_step * delta_gradient,
        weights,
        points,
    )
    finite_difference = (plus_values - minus_values) / (2.0 * fd_step)
    analytic = jacobian @ direction
    error = finite_difference - analytic
    return {
        "direction_norm": 1.0,
        "step": fd_step,
        "analytic": analytic.tolist(),
        "finite_difference": finite_difference.tolist(),
        "max_abs_error": float(np.max(np.abs(error))),
        "relative_error": float(np.max(np.abs(error)) / max(np.max(np.abs(finite_difference)), 1.0e-30)),
    }


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    mean_report = json.loads(MEAN_PATH.read_text(encoding="utf-8"))
    moment_reference = json.loads(MOMENT_REFERENCE_PATH.read_text(encoding="utf-8"))
    cone_reference = json.loads(CONE_REFERENCE_PATH.read_text(encoding="utf-8"))
    enriched_rows = json.loads(ENRICHED_ROWS_PATH.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE_PATH.read_text(encoding="utf-8"))
    geometry = json.loads(GEOMETRY_PATH.read_text(encoding="utf-8"))["inputs"]["wave"]
    if candidate.get("status") != "completed" or candidate.get("control_count") != 264:
        raise ValueError("Expected completed enriched_mean_endpoint_tangent 264-control candidate")
    candidate_control = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    if candidate_control.shape != (264,):
        raise ValueError(f"Unexpected selected candidate control shape {candidate_control.shape}")

    # Assemble the moment response using the same degree-2 reusable definition
    # used by the frozen mean compatibility rows, but at the local patch support.
    dynamic, dynamic_report = load_dynamic_field(DYNAMIC_PATH)
    grouped_replacements = int(install_in_field(dynamic))
    tau = float(mean_report["tau"])
    k = float(mean_report["k"])
    breaks = [float(value) for value in mean_report["quadrature"]["radial_split_breaks"]]
    moment_rows, moment_geometry = _mode0_moment_rows(
        dynamic.inner, tau, breaks, PATCH_CENTER, PATCH_WIDTHS
    )

    # Reuse the frozen grouped cone geometry and H transforms.  _mode0_rows
    # returns the 81-row full tangent layout with only mode-0 columns filled.
    saved = mean_report
    cone_points, cone_panels, cone_H = _cone_geometry(saved, dynamic.inner, k, breaks)
    cone_target_response, cone_rows = _mode0_rows(
        cone_points, cone_panels, cone_H, PATCH_CENTER, PATCH_WIDTHS
    )
    if moment_rows.shape != (4, 180) or cone_rows.shape != (81, 180):
        raise ValueError(f"Unexpected local row shapes {moment_rows.shape} {cone_rows.shape}")

    # Baseline vectors and initial-wave quadratic forms are copied from the
    # existing frozen reports; the local patch changes only linear rows.
    moment_linearization = moment_reference["reusable_moment_linearization"]
    baseline_moments = np.asarray(moment_linearization["baseline_moments"], dtype=float)
    wave_moment_forms = np.asarray(moment_linearization["wave_moment_forms"], dtype=float)
    cone_baseline = np.asarray(cone_reference["cone_baseline"], dtype=float)
    cone_wave_forms = np.asarray(cone_reference["cone_wave_forms"], dtype=float)
    cone_lower = np.asarray(cone_reference["cone_lower"], dtype=float)
    if baseline_moments.shape != (4,) or wave_moment_forms.shape != (4, 54, 54):
        raise ValueError("Unexpected frozen moment baseline/form shapes")
    if cone_baseline.shape != (81,) or cone_wave_forms.shape != (81, 54, 54) or cone_lower.shape != (81,):
        raise ValueError("Unexpected frozen cone baseline/form shapes")

    endpoint = EnrichedEndpointShape(ENDPOINT_CACHE_PATH)
    if endpoint.point_count != 7776:
        raise ValueError(f"Expected 7776 endpoint points, got {endpoint.point_count}")
    parent_velocity = endpoint.velocity_offset + np.einsum(
        "ncq,q->nc", endpoint.velocity_control, candidate_control
    )
    parent_gradient = endpoint.gradient_offset + np.einsum(
        "ncdq,q->ncd", endpoint.gradient_control, candidate_control
    )
    parent_values, _, parent_detail = _observables_and_jacobian(
        parent_velocity, parent_gradient, endpoint.weights, endpoint.points
    )
    endpoint_values_check, endpoint_full_jacobian = endpoint.evaluate(candidate_control)
    if not np.allclose(parent_values, endpoint_values_check, rtol=0.0, atol=1.0e-13):
        raise ValueError("Parent endpoint reconstruction disagrees with EnrichedEndpointShape.evaluate")
    local_velocity, local_gradient, endpoint_layout, endpoint_precision = (
        _local_endpoint_response_with_carrier(
            endpoint.points,
            PATCH_CENTER,
            PATCH_WIDTHS,
            np.asarray(geometry["carrier"], dtype=float),
            endpoint.dtau,
        )
    )
    endpoint_shape_values, endpoint_shape_jacobian, endpoint_detail = _observables_and_jacobian(
        parent_velocity,
        parent_gradient,
        endpoint.weights,
        endpoint.points,
        local_velocity,
        local_gradient,
    )
    shape_fd = _shape_fd_check(
        parent_velocity,
        parent_gradient,
        local_velocity,
        local_gradient,
        endpoint.weights,
        endpoint.points,
        endpoint_shape_jacobian,
    )

    mode0_tail_moment = float(np.max(np.abs(moment_rows[:, 36:])))
    mode0_tail_cone = float(np.max(np.abs(cone_rows[:, 36:])))
    pressure_indices = [index for index, entry in enumerate(endpoint_layout) if entry[1] == "pressure"]
    pressure_velocity_max = float(np.max(np.abs(local_velocity[:, :, pressure_indices])))
    pressure_gradient_max = float(np.max(np.abs(local_gradient[:, :, :, pressure_indices])))
    source_paths = {
        "mean": MEAN_PATH,
        "dynamic": DYNAMIC_PATH,
        "moment_reference": MOMENT_REFERENCE_PATH,
        "cone_reference": CONE_REFERENCE_PATH,
        "enriched_constraint_rows": ENRICHED_ROWS_PATH,
        "endpoint_cache": ENDPOINT_CACHE_PATH,
        "endpoint_cache_report": ENDPOINT_CACHE_REPORT_PATH,
        "candidate": CANDIDATE_PATH,
        "geometry": GEOMETRY_PATH,
    }
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": (
            "Localized degree-2 constraint rows and endpoint shape response only. "
            "The 264-control enriched_mean_endpoint_tangent candidate is the parent "
            "state. New moment/cone rows are linear responses; frozen initial-wave "
            "baselines and quadratic forms are copied unchanged. No optimization, "
            "mean replay, finite-time, PDE, or adoption claim."
        ),
        "source": CANDIDATE_PATH.name,
        "source_sha256": _sha256(CANDIDATE_PATH),
        "source_hashes": {name: {"path": str(path.name), "sha256": _sha256(path)} for name, path in source_paths.items()},
        "inputs": {
            "patch_center": PATCH_CENTER.tolist(),
            "patch_widths": PATCH_WIDTHS.tolist(),
            "degree": DEGREE,
            "modes": list(MODES),
            "tau": tau,
            "k": k,
            "nu": NU,
            "physical_time": -tau,
            "grouped_backend_replacements": grouped_replacements,
            "candidate_control_count": 264,
            "local_control_count": 180,
            "candidate_control_norm": float(np.linalg.norm(candidate_control)),
            "endpoint_point_count": int(endpoint.point_count),
            "endpoint_dtau": float(endpoint.dtau),
            "endpoint_carrier": np.asarray(geometry["carrier"], dtype=float).tolist(),
            "endpoint_parent_control_source": CANDIDATE_PATH.name,
        },
        "support_containment": _support_report(geometry),
        "row_layout": {
            "local_layout": _layout(),
            "mode0_columns": 36,
            "mode1_columns": 72,
            "mode2_columns": 72,
            "moment_rows": "4x180; only columns 0:36 are nonzero",
            "cone_rows": "81x180; only columns 0:36 are nonzero",
            "endpoint_response": "velocity/gradient response for all 180 local controls; pressure columns zero",
        },
        "moment_rows": moment_rows.tolist(),
        "cone_rows": cone_rows.tolist(),
        "moment_geometry": moment_geometry,
        "cone_geometry": {
            "location_count": int(len(cone_panels)),
            "inequality_count": int(len(cone_rows)),
            "grouped_point_count": int(len(cone_points)),
            "cone_order": 64,
            "H_shape": list(np.asarray(cone_H).shape),
        },
        "row_checks": {
            "moment_shape": list(moment_rows.shape),
            "cone_shape": list(cone_rows.shape),
            "moment_angular_tail_max_abs": mode0_tail_moment,
            "cone_angular_tail_max_abs": mode0_tail_cone,
            "mean_zero_angular_columns": bool(mode0_tail_moment == 0.0 and mode0_tail_cone == 0.0),
            "moment_baseline_shape": list(baseline_moments.shape),
            "moment_wave_form_shape": list(wave_moment_forms.shape),
            "cone_baseline_shape": list(cone_baseline.shape),
            "cone_wave_form_shape": list(cone_wave_forms.shape),
            "cone_lower_shape": list(cone_lower.shape),
            "frozen_initial_wave_baseline_unchanged": True,
            "baseline_sources": {
                "moment_baseline": MOMENT_REFERENCE_PATH.name,
                "moment_wave_forms": MOMENT_REFERENCE_PATH.name,
                "cone_baseline": CONE_REFERENCE_PATH.name,
                "cone_wave_forms": CONE_REFERENCE_PATH.name,
                "cone_lower": CONE_REFERENCE_PATH.name,
            },
        },
        "frozen_baselines": {
            "moment_baseline": baseline_moments.tolist(),
            "moment_wave_forms_shape": list(wave_moment_forms.shape),
            "cone_baseline": cone_baseline.tolist(),
            "cone_lower": cone_lower.tolist(),
            "cone_wave_forms_shape": list(cone_wave_forms.shape),
        },
        "endpoint_shape": {
            "parent_observables": parent_values.tolist(),
            "parent_detail": {
                "energy": float(parent_detail["energy"]),
                "zmean": float(parent_detail["zmean"]),
                "radial_rms": float(parent_detail["radial_rms"]),
                "axial_rms": float(parent_detail["axial_rms"]),
                "spin": float(parent_detail["spin"]),
            },
            "parent_reconstruction_max_abs_vs_cache_api": float(np.max(np.abs(parent_values - endpoint_values_check))),
            "response_velocity_shape": list(local_velocity.shape),
            "response_gradient_shape": list(local_gradient.shape),
            "response_velocity_max_abs": float(np.max(np.abs(local_velocity))),
            "response_gradient_max_abs": float(np.max(np.abs(local_gradient))),
            "pressure_column_count": int(len(pressure_indices)),
            "pressure_velocity_max_abs": pressure_velocity_max,
            "pressure_gradient_max_abs": pressure_gradient_max,
            "pressure_response_exactly_zero": bool(pressure_velocity_max == 0.0 and pressure_gradient_max == 0.0),
            "shape_values_at_parent": endpoint_shape_values.tolist(),
            "shape_jacobian": endpoint_shape_jacobian.tolist(),
            "shape_jacobian_shape": list(endpoint_shape_jacobian.shape),
            "endpoint_full_parent_jacobian_shape": list(endpoint_full_jacobian.shape),
            "layout": [
                {"mode": int(mode), "kind": kind, "component": component, "index": int(index), "precision": endpoint_precision[i]}
                for i, (mode, kind, component, index) in enumerate(endpoint_layout)
            ],
            "directional_jacobian_check": shape_fd,
        },
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "moment_rows": list(moment_rows.shape),
        "cone_rows": list(cone_rows.shape),
        "moment_angular_tail_max_abs": mode0_tail_moment,
        "cone_angular_tail_max_abs": mode0_tail_cone,
        "endpoint_shape_jacobian": list(endpoint_shape_jacobian.shape),
        "endpoint_jacobian_relative_error": shape_fd["relative_error"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
