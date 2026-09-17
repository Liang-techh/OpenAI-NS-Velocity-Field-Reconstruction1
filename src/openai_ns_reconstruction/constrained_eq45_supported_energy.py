"""Independent kinetic-energy audit for the support-connected Eq. (4.5) child.

The physical-space support transform changes the public ``[u,v,w]`` field in
its exterior collar. Energy evidence from the untapered parent therefore cannot
be inherited. This module consumes only the public ``at_points(points, time)``
output of :class:`Eq45SupportedVelocityCandidate` and applies independent
tensor-product Gauss--Legendre quadrature.

The audit reads the preregistered CR001 energy definition, quadrature ladder and
thresholds from ``configs/constraints.json``. It reports those gates exactly as
registered; it never relaxes them. Passing or failing this energy audit is
independent of PDE validity and visual correspondence.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.polynomial.legendre import leggauss

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


SCHEMA = "eq45_supported_energy_audit_v1"


def _as_box_half_width(config: Mapping[str, Any]) -> float:
    box = np.asarray(config["domain"]["evaluation_box"], dtype=float)
    if box.shape != (3, 2) or not np.all(np.isfinite(box)):
        raise ValueError("evaluation_box must be finite with shape (3,2)")
    if not np.allclose(box[:, 0], -box[:, 1], rtol=0.0, atol=0.0):
        raise ValueError("this audit requires a symmetric Cartesian evaluation box")
    half_widths = box[:, 1]
    if not np.allclose(half_widths, half_widths[0], rtol=0.0, atol=0.0):
        raise ValueError("this audit requires one common cube half-width")
    if half_widths[0] <= 0.0:
        raise ValueError("evaluation box half-width must be positive")
    return float(half_widths[0])


def _validate_orders(values: Sequence[int]) -> tuple[int, ...]:
    orders = tuple(int(v) for v in values)
    if len(orders) < 3 or any(v < 2 for v in orders):
        raise ValueError("at least three quadrature orders >=2 are required")
    if any(float(raw) != int(raw) for raw in values):
        raise ValueError("quadrature orders must be integers")
    if any(b <= a for a, b in zip(orders, orders[1:])):
        raise ValueError("quadrature orders must be strictly increasing")
    return orders


def _validate_times(values: Sequence[float], start: float, end: float) -> tuple[float, ...]:
    times = tuple(float(v) for v in values)
    if not times or not np.all(np.isfinite(times)):
        raise ValueError("times must be nonempty and finite")
    if any(t < start or t > end for t in times):
        raise ValueError("audit time lies outside the candidate interval")
    return times


def gauss_energy_and_components(
    field,
    *,
    time: float,
    half_width: float,
    order: int,
) -> dict[str, float]:
    """Integrate total/radial/swirl/axial kinetic energy on one cube.

    Velocity is sampled only through ``field.at_points``. Cylindrical
    decomposition is reconstructed independently from Cartesian output.
    """
    if half_width <= 0.0 or not np.isfinite(half_width):
        raise ValueError("half_width must be positive and finite")
    if int(order) != order or int(order) < 2:
        raise ValueError("order must be an integer >=2")

    nodes, weights = leggauss(int(order))
    coord = half_width * nodes
    weight = half_width * weights
    xx, yy = np.meshgrid(coord, coord, indexing="ij")
    ww = np.multiply.outer(weight, weight)
    flat_xy = np.column_stack((xx.ravel(), yy.ravel()))
    flat_w = ww.ravel()

    totals = np.zeros(4, dtype=float)  # total, radial, swirl, axial
    for z, wz in zip(coord, weight):
        points = np.empty((flat_xy.shape[0], 3), dtype=float)
        points[:, :2] = flat_xy
        points[:, 2] = z
        velocity = np.asarray(field.at_points(points, float(time)), dtype=float)
        if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
            raise ValueError("public velocity must return finite values with shape (n,3)")

        x = points[:, 0]
        y = points[:, 1]
        r = np.hypot(x, y)
        u = velocity[:, 0]
        v = velocity[:, 1]
        w = velocity[:, 2]
        radial = np.divide(x * u + y * v, r, out=np.zeros_like(r), where=r > 0.0)
        swirl = np.divide(-y * u + x * v, r, out=np.zeros_like(r), where=r > 0.0)

        density = np.column_stack(
            (
                u * u + v * v + w * w,
                radial * radial,
                swirl * swirl,
                w * w,
            )
        )
        totals += 0.5 * float(wz) * np.sum(flat_w[:, None] * density, axis=0)

    if not np.all(np.isfinite(totals)) or totals[0] <= 0.0:
        raise RuntimeError("energy quadrature returned a nonfinite or inactive field")
    closure = abs(float(totals[0] - np.sum(totals[1:]))) / float(totals[0])
    if closure > 5e-12:
        raise RuntimeError("cylindrical component energies do not close to total energy")
    return {
        "total": float(totals[0]),
        "radial": float(totals[1]),
        "swirl": float(totals[2]),
        "axial": float(totals[3]),
        "component_closure_relative_error": float(closure),
    }


def audit_supported_eq45_energy(
    child: Eq45SupportedVelocityCandidate,
    *,
    constraints_path: str | Path,
    times: Sequence[float] = (0.25, 0.5, 0.75),
) -> dict[str, Any]:
    """Revalidate CR001 energy gates after the physical support transform."""
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("child must be Eq45SupportedVelocityCandidate")

    config = json.loads(Path(constraints_path).read_text(encoding="utf-8"))
    if config.get("experiment_id") != "compact_axisymmetric_window_v1":
        raise ValueError("unexpected constrained experiment identity")
    if config.get("domain", {}).get("physical") != "R^3":
        raise ValueError("CR001 physical domain drifted")
    if config.get("domain", {}).get("support") != "r < 2 and abs(z) < 2":
        raise ValueError("CR001 support contract drifted")

    half_width = _as_box_half_width(config)
    taper = child.taper
    if taper.radial_support != half_width or taper.axial_half_height != half_width:
        raise ValueError("supported-child taper does not match the CR001 support box")

    validation = config["validation"]
    orders = _validate_orders(validation["quadrature_orders_per_axis"])
    audit_times = _validate_times(times, child.time_start, child.time_end)
    nontriviality = config["nontriviality"]
    reference_time = float(nontriviality["reference_time"])
    if reference_time not in audit_times:
        raise ValueError("audit times must include the preregistered reference time")

    reference_energy = float(nontriviality["reference_energy"])
    reference_tolerance = float(nontriviality["reference_energy_abs_tolerance"])
    min_energy = float(nontriviality["minimum_energy_each_validation_time"])
    max_energy = float(nontriviality["maximum_energy_each_validation_time"])
    quadrature_tolerance = float(validation["thresholds"]["energy_quadrature_relative_change"])
    if min(reference_energy, reference_tolerance, min_energy, max_energy, quadrature_tolerance) <= 0.0:
        raise ValueError("energy contract contains a nonpositive value")

    rows: list[dict[str, Any]] = []
    for time in audit_times:
        per_order = []
        for order in orders:
            energy = gauss_energy_and_components(
                child,
                time=time,
                half_width=half_width,
                order=order,
            )
            per_order.append({"order": int(order), **energy})

        finest = float(per_order[-1]["total"])
        medium = float(per_order[-2]["total"])
        fine_change = abs(finest - medium) / finest

        # Same-order parent comparison quantifies what the support transform
        # changed without treating parent finite-window energy as an R^3 result.
        compare_order = int(orders[1])
        parent_energy = gauss_energy_and_components(
            child.parent,
            time=time,
            half_width=half_width,
            order=compare_order,
        )
        child_compare = next(row for row in per_order if row["order"] == compare_order)
        parent_total = float(parent_energy["total"])
        child_total_compare = float(child_compare["total"])
        transform_fraction = (child_total_compare - parent_total) / parent_total

        rows.append(
            {
                "time": float(time),
                "quadrature": per_order,
                "finest_total_energy": finest,
                "medium_to_finest_relative_change": float(fine_change),
                "quadrature_gate_pass": bool(fine_change <= quadrature_tolerance),
                "energy_range_gate_pass": bool(min_energy <= finest <= max_energy),
                "parent_comparison_order": compare_order,
                "parent_finite_window_total_energy": parent_total,
                "child_same_order_total_energy": child_total_compare,
                "support_transform_energy_change_fraction": float(transform_fraction),
                "child_same_order_component_fractions": {
                    name: float(child_compare[name]) / child_total_compare
                    for name in ("radial", "swirl", "axial")
                },
                "parent_same_order_component_fractions": {
                    name: float(parent_energy[name]) / parent_total
                    for name in ("radial", "swirl", "axial")
                },
            }
        )

    reference_row = next(row for row in rows if row["time"] == reference_time)
    reference_observed = float(reference_row["finest_total_energy"])
    reference_abs_error = abs(reference_observed - reference_energy)

    return {
        "schema": SCHEMA,
        "task_id": "CR007-EQ45-SUPPORTED-ENERGY-REVALIDATION-017",
        "candidate_schema": child.to_dict()["schema"],
        "candidate_sha256": child.sha256,
        "parent_sha256": child.parent_sha256,
        "public_velocity_interface": "Eq45SupportedVelocityCandidate.at_points",
        "kinetic_energy_definition": nontriviality["kinetic_energy_definition"],
        "physical_domain": config["domain"]["physical"],
        "support_contract": config["domain"]["support"],
        "integration_cube_half_width": half_width,
        "quadrature_orders_per_axis": list(orders),
        "quadrature_relative_change_threshold": quadrature_tolerance,
        "reference_time": reference_time,
        "reference_energy": reference_energy,
        "reference_energy_abs_tolerance": reference_tolerance,
        "reference_observed_finest_energy": reference_observed,
        "reference_energy_abs_error": float(reference_abs_error),
        "reference_energy_gate_pass": bool(reference_abs_error <= reference_tolerance),
        "all_quadrature_gates_pass": bool(all(row["quadrature_gate_pass"] for row in rows)),
        "all_energy_range_gates_pass": bool(all(row["energy_range_gate_pass"] for row in rows)),
        "rows": rows,
        "interpretation": {
            "child_energy_is_post_support_transform_evidence": True,
            "parent_comparison_is_finite_window_context_only": True,
            "support_transform_changes_velocity": True,
            "energy_gate_does_not_imply_pde_validity": True,
            "energy_gate_does_not_imply_visual_correspondence": True,
        },
        "truth_boundary": {
            "velocity_changed_by_this_audit": False,
            "physical_support_connection_implemented": True,
            "physical_support_validated": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
