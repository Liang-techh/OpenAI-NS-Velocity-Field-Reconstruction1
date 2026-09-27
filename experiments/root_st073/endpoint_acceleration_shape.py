"""Fixed-cylinder shape replay for the frozen endpoint acceleration candidate.

This is a velocity-only diagnostic.  It reconstructs the balanced affine
field from ``endpoint_acceleration_projection.json`` and adds the frozen
zero-at-reference acceleration correction.  Spatial derivatives use the
same five-point stencil as ``vortex_state_observables.observe`` but are
batched over all stencil points.  No full momentum replay is performed.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from broad_meridional_constrained import load_saved_field
from endpoint_acceleration_projection import (
    AccelerationCorrectionField,
    MixedAffineField,
    _unpack_acceleration_control,
    _unpack_wave,
    _unpack_mixed_affine_control,
)
from grouped_joined_field import install_in_field
from vortex_state_observables import fixed_cylinder, observe


ROOT = Path(__file__).resolve().parent
PROJECTION = ROOT / "endpoint_acceleration_projection.json"
FROZEN = ROOT / "full_wave_frozen_cache.json"
OUTPUT = ROOT / "endpoint_acceleration_shape.json"


OBSERVABLE_KEYS = (
    "sampled_vorticity_max",
    "cylinder_enstrophy",
    "enstrophy_radial_rms",
    "enstrophy_axial_rms",
    "enstrophy_aspect_ratio",
    "enstrophy_weighted_angular_speed",
    "sampled_swirl_max",
    "sampled_circulation_max",
    "cylinder_kinetic_energy",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metric_delta(reference, value):
    """Return signed and relative changes for scalar shape observables."""

    result = {}
    for key in OBSERVABLE_KEYS:
        a = float(reference[key])
        b = float(value[key])
        result[key] = {
            "absolute": b - a,
            "relative_to_reference": (b - a) / max(abs(a), 1.0e-300),
        }
    return result


def _batched_observe(field, points, weights, k):
    """Equivalent to ``vortex_state_observables.observe`` with one batch.

    The only change is evaluating all 13 spatial stencil point sets in one
    ``fields`` call, which keeps this bounded shape replay from repeatedly
    entering the joined-field point loop.  The formulas and stencil are the
    existing diagnostic's formulas.
    """

    points = np.asarray(points, dtype=float)
    weights = np.asarray(weights, dtype=float)
    tau = 0.5 * 2.0 ** (-float(k))
    h = 5.0e-4 * np.sqrt(float(field.nu) * tau)
    offsets = [np.zeros(3, dtype=float)]
    for axis in np.eye(3):
        offsets.extend([-2.0 * h * axis, -h * axis, h * axis, 2.0 * h * axis])
    stacked = np.concatenate([points + offset for offset in offsets], axis=0)
    stacked_velocity = np.asarray(field.fields(stacked, tau)[0], dtype=float)
    count = len(points)
    u = stacked_velocity[:count]
    gradient = np.empty((count, 3, 3), dtype=float)
    cursor = count
    for axis_index in range(3):
        um2 = stacked_velocity[cursor:cursor + count]
        um = stacked_velocity[cursor + count:cursor + 2 * count]
        up = stacked_velocity[cursor + 2 * count:cursor + 3 * count]
        up2 = stacked_velocity[cursor + 3 * count:cursor + 4 * count]
        gradient[:, :, axis_index] = (um2 - 8.0 * um + 8.0 * up - up2) / (12.0 * h)
        cursor += 4 * count

    omega = np.column_stack((
        gradient[:, 2, 1] - gradient[:, 1, 2],
        gradient[:, 0, 2] - gradient[:, 2, 0],
        gradient[:, 1, 0] - gradient[:, 0, 1],
    ))
    density = np.sum(omega * omega, axis=1)
    enstrophy = float(weights @ density)
    if not np.isfinite(enstrophy) or enstrophy <= 0.0:
        raise ValueError("Nonpositive or nonfinite sampled enstrophy")
    radius = np.linalg.norm(points[:, :2], axis=1)
    swirl = (-points[:, 1] * u[:, 0] + points[:, 0] * u[:, 1]) / radius
    zmean = float(weights @ (density * points[:, 2]) / enstrophy)
    radial_rms = float(np.sqrt(weights @ (density * radius * radius) / enstrophy))
    axial_rms = float(np.sqrt(
        weights @ (density * (points[:, 2] - zmean) ** 2) / enstrophy
    ))
    return {
        "k": float(k),
        "tau": float(tau),
        "point_count": int(count),
        "sampled_vorticity_max": float(np.sqrt(density.max())),
        "cylinder_enstrophy": enstrophy,
        "enstrophy_radial_rms": radial_rms,
        "enstrophy_axial_rms": axial_rms,
        "enstrophy_aspect_ratio": float(axial_rms / radial_rms),
        "enstrophy_weighted_angular_speed": float(weights @ (density * swirl / radius) / enstrophy),
        "sampled_swirl_max": float(np.max(np.abs(swirl))),
        "sampled_circulation_max": float(np.max(np.abs(2.0 * np.pi * radius * swirl))),
        "cylinder_kinetic_energy": float(0.5 * (weights @ np.sum(u * u, axis=1))),
    }


def _build_fields(projection, frozen):
    geometry = frozen["inputs"]["wave"]
    center = tuple(float(value) for value in geometry["center"])
    widths = tuple(float(value) for value in geometry["widths"])
    carrier = np.asarray(geometry["carrier"], dtype=float)
    tau0 = float(frozen["inputs"]["mean"]["tau"])
    base_control = np.asarray(projection["inputs"]["base_tangent_coefficients"], dtype=float)
    if base_control.shape != (264,):
        raise ValueError(f"Expected 264 base tangent controls, got {base_control.shape}")
    wave = _unpack_wave(projection["inputs"]["initial_wave_coefficients"])
    mean, mean_report = load_saved_field()
    install_in_field(mean)
    base = MixedAffineField(mean, center, widths, carrier, wave, base_control, tau0)
    acceleration_coefficients = np.asarray(projection["ridge_fit"]["coefficients"], dtype=float)
    acceleration = _unpack_acceleration_control(acceleration_coefficients)
    acceleration_velocity = {mode: values[0] for mode, values in acceleration.items()}
    acceleration_pressure = {mode: values[1] for mode, values in acceleration.items()}
    corrected = AccelerationCorrectionField(
        base, center, widths, carrier,
        acceleration_velocity, acceleration_pressure, tau0,
    )
    return base, corrected, mean, mean_report, geometry, tau0


def run(output_path=OUTPUT):
    started = time.perf_counter()
    output_path = Path(output_path)
    projection_raw = PROJECTION.read_bytes()
    frozen_raw = FROZEN.read_bytes()
    projection = json.loads(projection_raw)
    frozen = json.loads(frozen_raw)
    if projection.get("status") != "completed":
        raise ValueError("endpoint_acceleration_projection.json is not completed")
    if projection.get("accepted", True):
        raise ValueError("Expected the diagnostic projection to remain unaccepted")
    base, corrected, mean, mean_report, geometry, tau0 = _build_fields(projection, frozen)
    k0 = float(frozen["inputs"]["mean"]["k"])
    delta_k = float(projection["inputs"]["delta_k"])
    k1 = k0 + delta_k
    points, weights, domain = fixed_cylinder(mean.inner, k0, radial_order=6, axial_order=12, angles=12)
    report = {
        "status": "running",
        "accepted": False,
        "pde_validated": False,
        "constraints_maintained": False,
        "source": {
            "projection": {"path": PROJECTION.name, "sha256": _sha256(PROJECTION)},
            "frozen_geometry": {"path": FROZEN.name, "sha256": _sha256(FROZEN)},
            "mean_loader_status": mean_report.get("status", "loaded"),
        },
        "scope": (
            "Fixed physical-cylinder shape diagnostic for the frozen endpoint "
            "acceleration projection. The balanced affine field and its zero-at-reference "
            "degree-2 mode-0..4 correction are compared at t0 and k0+delta_k. "
            "No full momentum replay, endpoint constraint solve, PDE, trajectory, or recursion acceptance."
        ),
        "inputs": {
            "reference_k": k0,
            "endpoint_k": k1,
            "delta_k": delta_k,
            "reference_tau": tau0,
            "endpoint_tau": float(0.5 * 2.0 ** (-k1)),
            "point_count": int(len(points)),
            "geometry": geometry,
            "domain": domain,
            "observable_method": "batched five-point Cartesian velocity FD, equivalent to vortex_state_observables.observe",
        },
        "rows": {},
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "shape_replay_started", "point_count": len(points), "rows": 4}), flush=True)

    # The reference correction is exactly zero by construction.  Evaluate the
    # balanced reference directly, then record the corrected reference as an
    # exact identity rather than paying for a duplicate FD stencil.
    print(json.dumps({"stage": "balanced_reference_started"}), flush=True)
    balanced_reference = _batched_observe(base, points, weights, k0)
    report["rows"]["balanced_reference"] = balanced_reference
    report["rows"]["corrected_reference"] = dict(balanced_reference)
    report["rows"]["corrected_reference"]["identity_from_balanced"] = True
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "balanced_reference_complete", "enstrophy": balanced_reference["cylinder_enstrophy"]}), flush=True)

    print(json.dumps({"stage": "balanced_endpoint_started"}), flush=True)
    balanced_endpoint = _batched_observe(base, points, weights, k1)
    report["rows"]["balanced_endpoint"] = balanced_endpoint
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "balanced_endpoint_complete", "enstrophy": balanced_endpoint["cylinder_enstrophy"]}), flush=True)

    print(json.dumps({"stage": "corrected_endpoint_started"}), flush=True)
    corrected_endpoint = _batched_observe(corrected, points, weights, k1)
    report["rows"]["corrected_endpoint"] = corrected_endpoint
    report["correction_endpoint_delta"] = _metric_delta(balanced_endpoint, corrected_endpoint)
    report["balanced_endpoint_from_reference"] = _metric_delta(balanced_reference, balanced_endpoint)
    report["corrected_endpoint_from_reference"] = _metric_delta(balanced_reference, corrected_endpoint)
    report["reference_correction_identity"] = {
        "velocity_max_abs": 0.0,
        "pressure_max_abs": 0.0,
        "statement": "AccelerationCorrectionField carries 0.5*(tau0-tau)^2 velocity and (tau0-tau) pressure factors."
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "stage": "completed",
        "balanced_reference_enstrophy": balanced_reference["cylinder_enstrophy"],
        "balanced_endpoint_enstrophy": balanced_endpoint["cylinder_enstrophy"],
        "corrected_endpoint_enstrophy": corrected_endpoint["cylinder_enstrophy"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
