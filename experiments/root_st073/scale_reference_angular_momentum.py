"""Global angular-momentum diagnostic for the prescribed scale family.

For a compact divergence-free field with finite angular momentum, the pure
anisotropic global map gives

    M_z(s) = s**(3/2 - 2*h) M_z(1).

The angular integral removes every nonzero Fourier mode, so only the
axisymmetric swirl component contributes.  A smooth unforced Navier--Stokes
field conserves this quantity; the displayed scale family therefore requires
the corresponding net torque unless the reference angular momentum vanishes
or an outer compensating field is supplied.  This file measures the global
quantity for the current compact construction.  It makes no PDE or recursive
acceptance claim.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from global_collar_tangent import BOXES  # noqa: E402
from global_support_energy import panel_nodes, support_bounds  # noqa: E402
from scale_reference_candidate import load_reference  # noqa: E402


OUTPUT_PATH = ROOT / "scale_reference_angular_momentum.json"
REFERENCE_PATH = ROOT / "scale_reference_trust_nonlinear_fit.json"
FROZEN_PATH = ROOT / "full_wave_frozen_cache.json"
QUADRATURE_ORDERS = (4, 8)
ANGLE_COUNT = 12


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _edges(reference: Any) -> dict[str, Any]:
    tau0 = float(reference.tau0)
    geometry = reference.snapshot["inputs"]["wave"]
    bounds = support_bounds(reference.localized, geometry, tau0, extra_boxes=BOXES)
    radial = [0.0, float(bounds["radius"])]
    axial = [float(bounds["z_lower"]), float(bounds["z_upper"]), 0.0]
    center = np.asarray(geometry["center"], dtype=float)
    widths = np.asarray(geometry["widths"], dtype=float)
    radial.extend((float(center[0] - widths[0]), float(center[0] + widths[0])))
    axial.extend((float(center[1] - widths[1]), float(center[1] + widths[1])))
    for patch_center, patch_widths in zip(reference.center, reference.widths):
        radial.extend((float(patch_center[0] - patch_widths[0]), float(patch_center[0] + patch_widths[0])))
        axial.extend((float(patch_center[1] - patch_widths[1]), float(patch_center[1] + patch_widths[1])))
    radial.extend(float(value) for value in reference.localized.radial_radii(tau0))
    for r0, r1, z0, z1 in BOXES:
        radial.extend((float(r0), float(r1)))
        axial.extend((float(z0), float(z1)))
    radial_edges = sorted(set(value for value in radial if 0.0 <= value <= bounds["radius"]))
    axial_edges = sorted(set(value for value in axial if bounds["z_lower"] <= value <= bounds["z_upper"]))
    if radial_edges[0] != 0.0 or radial_edges[-1] != float(bounds["radius"]):
        raise ValueError("Radial partition does not span the compact support")
    if axial_edges[0] != float(bounds["z_lower"]) or axial_edges[-1] != float(bounds["z_upper"]):
        raise ValueError("Axial partition does not span the compact support")
    return {
        "radial_edges": radial_edges,
        "axial_edges": axial_edges,
        "bounds": bounds,
    }


def _volume_grid(edges: dict[str, Any], order: int, angles: int = ANGLE_COUNT) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    radius, radial_weights = panel_nodes(edges["radial_edges"], order)
    z, axial_weights = panel_nodes(edges["axial_edges"], order)
    theta = 2.0 * np.pi * np.arange(int(angles), dtype=float) / float(angles)
    rr, zz, aa = np.meshgrid(radius, z, theta, indexing="ij")
    points = np.column_stack((
        rr.ravel() * np.cos(aa.ravel()),
        rr.ravel() * np.sin(aa.ravel()),
        zz.ravel(),
    ))
    weights = (
        radial_weights[:, None, None]
        * axial_weights[None, :, None]
        * np.ones((1, 1, int(angles)), dtype=float)
        * rr
        * (2.0 * np.pi / float(angles))
    ).ravel()
    return points, weights, {
        "radius": radius,
        "radial_weights": radial_weights,
        "z": z,
        "axial_weights": axial_weights,
        "angles": theta,
        "angle_count": int(angles),
    }


def _angular_integrand(velocity: np.ndarray, points: np.ndarray) -> np.ndarray:
    theta = np.arctan2(points[:, 1], points[:, 0])
    u_theta = -np.sin(theta) * velocity[:, 0] + np.cos(theta) * velocity[:, 1]
    radius = np.hypot(points[:, 0], points[:, 1])
    return radius * u_theta


def _integral(values: np.ndarray, weights: np.ndarray) -> float:
    return float(
        np.asarray(weights, dtype=float).reshape(-1)
        @ np.asarray(values, dtype=float).reshape(-1)
    )


def _one_grid(reference: Any, base_reference: Any, edges: dict[str, Any], order: int) -> dict[str, Any]:
    points, weights, metadata = _volume_grid(edges, order)
    full_velocity = np.asarray(reference.velocity(points), dtype=float)
    base_velocity = np.asarray(base_reference.velocity(points), dtype=float)
    correction_velocity = full_velocity - base_velocity
    if full_velocity.shape != (len(points), 3) or not np.all(np.isfinite(full_velocity)):
        raise ValueError("Reference velocity has an invalid shape or non-finite values")
    full_integrand = _angular_integrand(full_velocity, points)
    base_integrand = _angular_integrand(base_velocity, points)
    correction_integrand = _angular_integrand(correction_velocity, points)
    # The grid is ordered (r, z, theta), so angular means can be integrated
    # without relying on cancellation between separate panels.
    nr, nz, na = len(metadata["radius"]), len(metadata["z"]), metadata["angle_count"]
    full_mean = full_integrand.reshape(nr, nz, na).mean(axis=2)
    base_mean = base_integrand.reshape(nr, nz, na).mean(axis=2)
    correction_mean = correction_integrand.reshape(nr, nz, na).mean(axis=2)
    cylindrical_factor = (
        metadata["radial_weights"][:, None]
        * metadata["axial_weights"][None, :]
        * (2.0 * np.pi)
        * metadata["radius"][:, None] ** 1
    )
    full_mean_mz = _integral(full_mean, cylindrical_factor)
    base_mean_mz = _integral(base_mean, cylindrical_factor)
    correction_mean_mz = _integral(correction_mean, cylindrical_factor)
    full_mz = _integral(full_integrand, weights)
    base_mz = _integral(base_integrand, weights)
    correction_mz = _integral(correction_integrand, weights)
    return {
        "panel_order_per_interval": int(order),
        "angle_count": int(na),
        "point_count": int(len(points)),
        "physical_volume": float(np.sum(weights)),
        "M_z": full_mz,
        "M_z_base_without_new_velocity": base_mz,
        "M_z_new_velocity_correction": correction_mz,
        "M_z_axisymmetric_angular_mean": full_mean_mz,
        "M_z_base_axisymmetric_angular_mean": base_mean_mz,
        "M_z_correction_axisymmetric_angular_mean": correction_mean_mz,
        "nonaxisymmetric_cancellation_error": full_mz - full_mean_mz,
        "base_nonaxisymmetric_cancellation_error": base_mz - base_mean_mz,
        "correction_nonaxisymmetric_cancellation_error": correction_mz - correction_mean_mz,
        "velocity_max_abs": float(np.max(np.abs(full_velocity))),
        "correction_velocity_max_abs": float(np.max(np.abs(correction_velocity))),
    }


def run(output_path: Path = OUTPUT_PATH) -> dict[str, Any]:
    started = time.perf_counter()
    output_path = Path(output_path)
    reference = load_reference(step=REFERENCE_PATH)
    edges = _edges(reference)
    rows = []
    for order in QUADRATURE_ORDERS:
        row = _one_grid(reference, reference.base_reference, edges, order)
        rows.append(row)
        print(json.dumps({"order": order, "point_count": row["point_count"], "M_z": row["M_z"]}), flush=True)
    finest = rows[-1]
    coarse = rows[-2]
    mz0 = float(finest["M_z_axisymmetric_angular_mean"])
    h = float(reference.localized.inner.h)
    tau0 = float(reference.tau0)
    beta = 1.5 - 2.0 * h
    scale_rows = []
    for scale in (1.0, 0.5, 0.25):
        scale_rows.append({
            "scale": scale,
            "M_z_predicted": float(scale ** beta * mz0),
            "required_net_torque": float(-beta * mz0 / tau0 * scale ** (0.5 - 2.0 * h)),
        })
    report: dict[str, Any] = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": (
            "Global compact angular-momentum quadrature for the prescribed reference "
            "velocity and its pure anisotropic scale law. This diagnoses the net-torque "
            "condition; it does not assert PDE, trajectory, or recursive-scale acceptance."
        ),
        "sources": {
            "reference_report": {"path": REFERENCE_PATH.name, "sha256": _sha256(REFERENCE_PATH)},
            "reference_loader": {"path": "scale_reference_candidate.py", "sha256": _sha256(ROOT / "scale_reference_candidate.py")},
            "frozen_geometry": {"path": FROZEN_PATH.name, "sha256": _sha256(FROZEN_PATH)},
            "support_quadrature": {"path": "global_support_energy.py", "sha256": _sha256(ROOT / "global_support_energy.py")},
            "collar_geometry": {"path": "global_collar_tangent.py", "sha256": _sha256(ROOT / "global_collar_tangent.py")},
        },
        "reference": {
            "tau0": tau0,
            "h": h,
            "viscosity": float(reference.nu),
            "source_semantics": reference.control_semantics,
        },
        "support": {
            "bounds": edges["bounds"],
            "radial_edges": edges["radial_edges"],
            "axial_edges": edges["axial_edges"],
            "collar_boxes": [list(box) for box in BOXES],
            "finite_compact_support_used": True,
        },
        "quadrature": {
            "orders": list(QUADRATURE_ORDERS),
            "angle_count": ANGLE_COUNT,
            "jacobian": "r dr dz dtheta",
            "integrand": "r u_theta, so cylindrical volume weighting gives r^2 u_theta dr dz dtheta",
            "rows": rows,
            "coarse_to_fine_abs_difference": float(finest["M_z"] - coarse["M_z"]),
            "coarse_to_fine_relative_difference": float(
                abs(finest["M_z"] - coarse["M_z"]) / max(abs(finest["M_z"]), 1.0e-300)
            ),
            "angular_modes_cancellation_checked": True,
        },
        "angular_momentum": {
            "M_z_reference": float(finest["M_z"]),
            "M_z_axisymmetric_reference": mz0,
            "M_z_base_without_new_velocity": float(finest["M_z_base_without_new_velocity"]),
            "M_z_new_velocity_correction": float(finest["M_z_new_velocity_correction"]),
            "base_plus_correction_reconstruction_error": float(
                finest["M_z"] - finest["M_z_base_without_new_velocity"] - finest["M_z_new_velocity_correction"]
            ),
            "relative_new_correction": float(
                finest["M_z_new_velocity_correction"] / max(abs(finest["M_z"]), 1.0e-300)
            ),
            "nonaxisymmetric_cancellation_error": float(finest["nonaxisymmetric_cancellation_error"]),
        },
        "pure_global_scale_law": {
            "exponent": beta,
            "formula": "M_z(s) = s^(3/2 - 2h) M_z(1)",
            "rows": scale_rows,
            "torque_formula": "dM_z/dt = -(3/2 - 2h) M_z(1) / tau0 * s^(1/2 - 2h)",
            "interpretation": (
                "A nonzero measured M_z requires net torque for this prescribed global "
                "map. A moment-neutral reference or an outer angular-momentum compensation "
                "can remove that requirement; this diagnostic does not choose either." 
            ),
        },
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "M_z": mz0, "required_torque_at_s1": scale_rows[0]["required_net_torque"], "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
