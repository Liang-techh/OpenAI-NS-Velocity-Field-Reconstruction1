"""Whole fixed-domain momentum diagnostic for the unified global candidate.

The loader assembles the frozen global mean, disjoint collar correction, and
endpoint acceleration correction.  This script integrates the complete
Cartesian finite-difference momentum residual over one fixed bounding
cylinder at the reference state and at ``k0 + 1e-6``.  The cylinder includes
the axis only through Gauss nodes (no singular point is sampled), and the
finite-difference stencil is allowed to cross the support boundary.

The result is a bounded sampled diagnostic.  It does not assert a continuum
estimate, a whole-space theorem, a trajectory, or scale-recursion acceptance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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

from global_acceleration_candidate import load  # noqa: E402
from global_collar_tangent import BOXES  # noqa: E402
from global_support_energy import panel_nodes, support_bounds  # noqa: E402
from affine_momentum import jets, momentum  # noqa: E402


OUTPUT_PATH = ROOT / "global_full_momentum.json"
SNAPSHOT_PATH = ROOT / "full_wave_frozen_cache.json"
GLOBAL_SOURCE_PATH = ROOT / "global_axial_extension.py"
SUPPORT_SOURCE_PATH = ROOT / "global_support_energy.py"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixed_edges(localized, geometry, tau0, tau_endpoint, hspace):
    """Return one fixed partition covering all supports and FD margin."""

    bounds0 = support_bounds(localized, geometry, tau0)
    bounds1 = support_bounds(localized, geometry, tau_endpoint)
    collar_radius = max(float(box[1]) for box in BOXES)
    collar_z_lower = min(float(box[2]) for box in BOXES)
    collar_z_upper = max(float(box[3]) for box in BOXES)
    fd_margin = 2.0 * float(hspace)
    radius = max(float(bounds0["radius"]), float(bounds1["radius"]), collar_radius) + fd_margin
    z_lower = min(float(bounds0["z_lower"]), float(bounds1["z_lower"]), collar_z_lower) - fd_margin
    z_upper = max(float(bounds0["z_upper"]), float(bounds1["z_upper"]), collar_z_upper) + fd_margin
    center = np.asarray(geometry["center"], dtype=float)
    widths = np.asarray(geometry["widths"], dtype=float)
    radial_edges = [
        0.0,
        float(center[0] - widths[0]),
        float(center[0] + widths[0]),
        float(localized.radial_radii(tau0)[0]),
        float(localized.radial_radii(tau_endpoint)[0]),
        float(localized.radial_radii(tau0)[1]),
        float(localized.radial_radii(tau_endpoint)[1]),
        radius,
    ]
    axial_edges = [
        z_lower,
        float(center[1] - widths[1]),
        float(center[1] + widths[1]),
        0.0,
        z_upper,
    ]
    # Include every registered collar boundary so the wave, collar, and rest
    # masks do not straddle a requested support edge inside a panel.
    for r0, r1, z0, z1 in BOXES:
        radial_edges.extend((float(r0), float(r1)))
        axial_edges.extend((float(z0), float(z1)))
    radial_edges = sorted(set(value for value in radial_edges if 0.0 <= value <= radius))
    axial_edges = sorted(set(value for value in axial_edges if z_lower <= value <= z_upper))
    if radial_edges[0] != 0.0 or radial_edges[-1] != radius:
        raise ValueError("Radial partition does not span the fixed cylinder")
    if axial_edges[0] != z_lower or axial_edges[-1] != z_upper:
        raise ValueError("Axial partition does not span the fixed cylinder")
    return {
        "bounds_reference": bounds0,
        "bounds_endpoint": bounds1,
        "raw_collar_bounds": {
            "radius_upper": collar_radius,
            "z_lower": collar_z_lower,
            "z_upper": collar_z_upper,
        },
        "finite_difference_stencil_margin": fd_margin,
        "radius": radius,
        "z_lower": z_lower,
        "z_upper": z_upper,
        "radial_edges": radial_edges,
        "axial_edges": axial_edges,
    }


def _volume_grid(edges, order=4, angles=20):
    radius, radial_weights = panel_nodes(edges["radial_edges"], order)
    z, axial_weights = panel_nodes(edges["axial_edges"], order)
    theta = 2.0 * np.pi * np.arange(int(angles), dtype=float) / float(angles)
    rr, zz, aa = np.meshgrid(radius, z, theta, indexing="ij")
    points = np.column_stack(
        (
            rr.ravel() * np.cos(aa.ravel()),
            rr.ravel() * np.sin(aa.ravel()),
            zz.ravel(),
        )
    )
    weights = (
        radial_weights[:, None, None]
        * axial_weights[None, :, None]
        * np.ones((1, 1, int(angles)), dtype=float)
        * rr
        * (2.0 * np.pi / float(angles))
    ).ravel()
    return points, weights


def _partition_masks(points, geometry):
    center = np.asarray(geometry["center"], dtype=float)
    widths = np.asarray(geometry["widths"], dtype=float)
    radius = np.hypot(points[:, 0], points[:, 1])
    wave = (
        (np.abs(radius - center[0]) < widths[0])
        & (np.abs(points[:, 2] - center[1]) < widths[1])
    )
    collar = np.zeros(len(points), dtype=bool)
    for r0, r1, z0, z1 in BOXES:
        collar |= (
            (radius > float(r0))
            & (radius < float(r1))
            & (points[:, 2] > float(z0))
            & (points[:, 2] < float(z1))
        )
    if np.any(wave & collar):
        raise ValueError("Wave and collar quadrature partitions overlap")
    rest = ~(wave | collar)
    return {"wave_patch": wave, "collar_union": collar, "rest": rest}


def _metrics(residual, points, weights, mask):
    selected = np.asarray(mask, dtype=bool)
    if not np.any(selected):
        return {
            "point_count": 0,
            "physical_volume": 0.0,
            "momentum_max": 0.0,
            "momentum_max_point": None,
            "momentum_volume_L2_squared": 0.0,
            "momentum_volume_L2": 0.0,
        }
    norms = np.linalg.norm(residual[selected], axis=1)
    local_weights = weights[selected]
    local_points = points[selected]
    index = int(np.argmax(norms))
    squared = float(local_weights @ np.sum(residual[selected] ** 2, axis=1))
    return {
        "point_count": int(np.sum(selected)),
        "physical_volume": float(np.sum(local_weights)),
        "momentum_max": float(norms[index]),
        "momentum_max_point": local_points[index].tolist(),
        "momentum_volume_L2_squared": squared,
        "momentum_volume_L2": float(math.sqrt(max(squared, 0.0))),
    }


def _state_metrics(residual, points, weights, masks):
    norms = np.linalg.norm(residual, axis=1)
    index = int(np.argmax(norms))
    total_squared = float(weights @ np.sum(residual**2, axis=1))
    partitions = {name: _metrics(residual, points, weights, mask) for name, mask in masks.items()}
    partition_squared = float(
        sum(row["momentum_volume_L2_squared"] for row in partitions.values())
    )
    return {
        "point_count": int(len(points)),
        "physical_volume": float(np.sum(weights)),
        "momentum_max": float(norms[index]),
        "momentum_max_point": points[index].tolist(),
        "momentum_volume_L2_squared": total_squared,
        "momentum_volume_L2": float(math.sqrt(max(total_squared, 0.0))),
        "partition_squared_L2_sum": partition_squared,
        "partition_squared_L2_absolute_difference": float(abs(total_squared - partition_squared)),
        "partitions": partitions,
    }


def _outside_probe(full, edges, tau, hspace):
    points = np.asarray(
        [
            [edges["radius"] + 3.0 * hspace, 0.0, 0.0],
            [0.0, 0.0, edges["z_upper"] + 3.0 * hspace],
            [0.0, 0.0, edges["z_lower"] - 3.0 * hspace],
        ],
        dtype=float,
    )
    velocity, pressure = full.fields(points, tau)
    return {
        "points": points.tolist(),
        "velocity_max_abs": float(np.max(np.abs(velocity))),
        "pressure_max_abs": float(np.max(np.abs(pressure))),
        "zero_outside_probe": bool(
            np.max(np.abs(velocity)) <= 1.0e-12
            and np.max(np.abs(pressure)) <= 1.0e-12
        ),
    }


def run(output_path=OUTPUT_PATH, order=4, angles=20):
    started = time.perf_counter()
    output_path = Path(output_path)
    full, global_base, localized, snapshot, reports, hashes = load()
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    k0 = float(snapshot["inputs"]["mean"]["k"])
    delta_k = 1.0e-6
    k_endpoint = k0 + delta_k
    tau_endpoint = tau0 * 2.0 ** (-delta_k)
    geometry = snapshot["inputs"]["wave"]
    edges = _fixed_edges(
        localized,
        geometry,
        tau0,
        tau_endpoint,
        float(snapshot["timesteps"]["hspace"]),
    )
    points, weights = _volume_grid(edges, order=order, angles=angles)
    masks = _partition_masks(points, geometry)
    prior_report = None
    if output_path.exists():
        try:
            candidate_prior = json.loads(output_path.read_text(encoding="utf-8"))
            if candidate_prior.get("status") == "completed":
                prior_report = candidate_prior
        except (OSError, json.JSONDecodeError):
            prior_report = None
    report = {
        "status": "running",
        "accepted": False,
        "pde_validated": False,
        "finite_time_evolution_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Complete Cartesian momentum residual over one fixed bounding cylinder "
            "covering the reference and endpoint supports. Gauss nodes avoid the "
            "axis; FD stencils are evaluated across support boundaries. This is a "
            "bounded sampled diagnostic with no whole-R3, continuum, trajectory, "
            "or scale-recursion acceptance claim."
        ),
        "sources": {
            "global_acceleration_candidate": {
                "path": "global_acceleration_candidate.py",
                "sha256": _sha256(ROOT / "global_acceleration_candidate.py"),
            },
            "global_source": {
                "path": GLOBAL_SOURCE_PATH.name,
                "sha256": _sha256(GLOBAL_SOURCE_PATH),
            },
            "support_bounds_source": {
                "path": SUPPORT_SOURCE_PATH.name,
                "sha256": _sha256(SUPPORT_SOURCE_PATH),
            },
            "balanced": {"sha256": hashes["balanced"]},
            "acceleration": {"sha256": hashes["acceleration"]},
            "collar": {"sha256": hashes["collar"]},
            "frozen_cache": {
                "path": SNAPSHOT_PATH.name,
                "sha256": _sha256(SNAPSHOT_PATH),
            },
        },
        "states": {
            "reference": {"k": k0, "tau": tau0},
            "endpoint": {"k": k_endpoint, "tau": tau_endpoint, "delta_k": delta_k},
            "physical_time": "t=-tau",
        },
        "quadrature": {
            "radial_order_per_panel": int(order),
            "axial_order_per_panel": int(order),
            "angle_order": int(angles),
            "point_count": int(len(points)),
            "physical_volume": float(np.sum(weights)),
            "radial_edges": edges["radial_edges"],
            "axial_edges": edges["axial_edges"],
            "bounds_reference": edges["bounds_reference"],
            "bounds_endpoint": edges["bounds_endpoint"],
            "raw_collar_bounds": edges["raw_collar_bounds"],
            "finite_difference_stencil_margin": edges["finite_difference_stencil_margin"],
            "axis_handling": "Gauss nodes lie strictly inside the first radial panel; r=0 is not sampled.",
            "jacobian": "2*pi*r dr dz with evenly spaced theta weights summing to 2*pi",
            "partition_point_counts": {name: int(np.sum(mask)) for name, mask in masks.items()},
            "partition_overlap_count": int(np.sum(masks["wave_patch"] & masks["collar_union"])),
            "partition_coverage_count": int(np.sum(masks["wave_patch"] | masks["collar_union"] | masks["rest"])),
            "quadrature_converged": False,
        },
        "support_enclosure": {
            "constructed_components": (
                "The global axial/radial cutoff localizes the mean; the compact "
                "wave and acceleration bases use the frozen rectangular bump; "
                "the collar correction is supported only on the three BOXES."
            ),
            "fixed_cylinder_includes_reference_and_endpoint": True,
            "fd_stencils_cross_support_boundary": True,
            "exterior_zero_claim": (
                "Pointwise zero exterior follows from these compact support factors "
                "for the constructed field on the registered time slab; the probe "
                "checks are numerical evidence, not a continuum theorem."
            ),
        },
        "prior_domain_replay": (
            {
                "status": "invalid_support_enclosure",
                "reason": (
                    "The first order-4 grid used global_support_energy bounds only; "
                    "those bounds omitted the collar extension to r=0.0105 and |z|=0.001."
                ),
                "report": prior_report,
            }
            if prior_report is not None
            else None
        ),
        "rows": [],
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    for label, k, tau in (("reference", k0, tau0), ("endpoint", k_endpoint, tau_endpoint)):
        residual = momentum(jets(full, points, tau, float(snapshot["timesteps"]["hspace"]), float(snapshot["timesteps"]["htime"])))
        row = {
            "label": label,
            "k": float(k),
            "tau": float(tau),
            "delta_k": float(k - k0),
            "finite_difference_steps": {
                "hspace": float(snapshot["timesteps"]["hspace"]),
                "htime": float(snapshot["timesteps"]["htime"]),
            },
            "metrics": _state_metrics(residual, points, weights, masks),
            "outside_support_probe": _outside_probe(
                full, edges, tau, float(snapshot["timesteps"]["hspace"])
            ),
        }
        report["rows"].append(row)
        report["status"] = f"{label}_evaluated"
        output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(
            json.dumps(
                {
                    "label": label,
                    "point_count": len(points),
                    "momentum_max": row["metrics"]["momentum_max"],
                    "momentum_volume_L2": row["metrics"]["momentum_volume_L2"],
                }
            ),
            flush=True,
        )

    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": report["status"],
                "output": str(output_path),
                "point_count": len(points),
                "elapsed_seconds": report["elapsed_seconds"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--order", type=int, default=4)
    parser.add_argument("--angles", type=int, default=20)
    args = parser.parse_args()
    run(args.output, order=args.order, angles=args.angles)
