"""Renderer-independent morphology fingerprints for continuous ST054 velocity fields.

This module samples only ``velocity(x,y,z,t)`` and reports autonomous cylindrical
kinematic diagnostics.  It does not compare pixels, infer hidden OpenAI data,
retune a candidate, or perform Navier--Stokes acceptance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Protocol, Sequence

import numpy as np

from .st054_snapshot import available_st054_models, load_st054

SCHEMA = "st054_cylindrical_morphology_fingerprint_v1"
DEFAULT_RADII = (0.4, 0.7, 1.0, 1.3)
DEFAULT_Z_MAGNITUDES = (0.3, 0.8)
DEFAULT_TIMES = (0.25, 0.5, 0.75)
DEFAULT_AZIMUTHS = 16
CORE_RADII = tuple(float(value) for value in np.linspace(0.2, 1.6, 15))
SWIRL_SPEED_FLOOR = 1.0e-12

EXTERNAL_METHOD_SCREEN = (
    {
        "repo": "pyvista/pyvista",
        "screened_commit": "f749a1b0a10a5a4c3c5ca3eedbc8f7860c6f9d9e",
        "license": "MIT",
        "classification": "screened_not_adopted",
        "scope": "streamline and derivative/filter APIs for vector-field visualization",
        "difference": (
            "PyVista can consume the VTK handoff already merged in PR #765, but this "
            "increment needs only direct callable-to-cylindrical statistics.  Adding a "
            "VTK/PyVista dependency here would duplicate the delivery path and make the "
            "diagnostic less implementation-direct."
        ),
    },
)

TRUTH_BOUNDARY = {
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


class VelocityField(Protocol):
    def velocity(self, x: Any, y: Any, z: Any, t: float) -> np.ndarray:
        """Return Cartesian velocity with final component axis of length three."""


def _finite_positive(values: Sequence[float], name: str) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if not result or not all(np.isfinite(result)) or not all(value > 0.0 for value in result):
        raise ValueError(f"{name} must contain finite positive values")
    return result


def _finite_times(values: Sequence[float]) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if not result or not all(np.isfinite(result)):
        raise ValueError("times must contain finite values")
    return result


def _azimuth_grid(count: int) -> np.ndarray:
    count = int(count)
    if count < 8 or count % 2:
        raise ValueError("azimuth_count must be an even integer >= 8")
    return 2.0 * np.pi * np.arange(count, dtype=float) / count


def _evaluate_ring(
    field: VelocityField,
    *,
    radius: float,
    z: float,
    time: float,
    azimuth_count: int,
) -> dict[str, float | int]:
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("radius must be finite and positive")
    if not np.isfinite(z) or z == 0.0:
        raise ValueError("z must be finite and nonzero for the signed axial-stretch diagnostic")
    if not np.isfinite(time):
        raise ValueError("time must be finite")

    phi = _azimuth_grid(azimuth_count)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    x = radius * cos_phi
    y = radius * sin_phi
    zz = np.full_like(x, float(z))
    velocity = np.asarray(field.velocity(x, y, zz, float(time)), dtype=float)
    if velocity.shape != (phi.size, 3) or not np.isfinite(velocity).all():
        raise RuntimeError("velocity callable returned malformed/non-finite ring data")

    u_r = velocity[:, 0] * cos_phi + velocity[:, 1] * sin_phi
    u_theta = -velocity[:, 0] * sin_phi + velocity[:, 1] * cos_phi
    u_z = velocity[:, 2]
    speed = np.linalg.norm(velocity, axis=1)
    circulation_speed = np.abs(u_theta)
    angular_rate = u_theta / radius
    away_from_midplane = np.sign(z) * u_z
    swirl_valid = circulation_speed > SWIRL_SPEED_FLOOR
    inward_swirl = (u_r < 0.0) & swirl_valid
    inward_to_circulation = np.zeros_like(u_r)
    inward_to_circulation[swirl_valid] = -u_r[swirl_valid] / circulation_speed[swirl_valid]

    return {
        "radius": float(radius),
        "z": float(z),
        "time": float(time),
        "azimuth_count": int(phi.size),
        "mean_u_r": float(np.mean(u_r)),
        "std_u_r": float(np.std(u_r)),
        "mean_u_theta": float(np.mean(u_theta)),
        "std_u_theta": float(np.std(u_theta)),
        "mean_abs_u_theta": float(np.mean(circulation_speed)),
        "mean_angular_rate": float(np.mean(angular_rate)),
        "std_angular_rate": float(np.std(angular_rate)),
        "mean_away_from_midplane_u_z": float(np.mean(away_from_midplane)),
        "std_away_from_midplane_u_z": float(np.std(away_from_midplane)),
        "mean_speed": float(np.mean(speed)),
        "max_speed": float(np.max(speed)),
        "inward_fraction": float(np.mean(u_r < 0.0)),
        "swirl_nonzero_fraction": float(np.mean(swirl_valid)),
        "inward_swirl_fraction": float(np.mean(inward_swirl)),
        "mean_inward_to_circulation_ratio": (
            float(np.mean(inward_to_circulation[swirl_valid])) if np.any(swirl_valid) else 0.0
        ),
        "away_from_midplane_fraction": float(np.mean(away_from_midplane > 0.0)),
    }


def _core_speed_proxy(
    field: VelocityField,
    *,
    time: float,
    radii: Sequence[float],
    azimuth_count: int,
) -> dict[str, Any]:
    radii_tuple = _finite_positive(radii, "core_radii")
    phi = _azimuth_grid(azimuth_count)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    ring_speed: list[float] = []
    for radius in radii_tuple:
        x = radius * cos_phi
        y = radius * sin_phi
        z = np.zeros_like(x)
        velocity = np.asarray(field.velocity(x, y, z, float(time)), dtype=float)
        if velocity.shape != (phi.size, 3) or not np.isfinite(velocity).all():
            raise RuntimeError("velocity callable returned malformed/non-finite core data")
        ring_speed.append(float(np.mean(np.linalg.norm(velocity, axis=1))))

    radius_array = np.asarray(radii_tuple, dtype=float)
    speed_array = np.asarray(ring_speed, dtype=float)
    weight = np.maximum(speed_array, 0.0)
    weight_sum = float(np.sum(weight))
    weighted_radius = (
        float(np.sqrt(np.sum(weight * radius_array * radius_array) / weight_sum))
        if weight_sum > 0.0
        else None
    )
    peak_index = int(np.argmax(speed_array))
    return {
        "time": float(time),
        "radii": list(radii_tuple),
        "mean_speed_by_radius": ring_speed,
        "speed_weighted_rms_radius": weighted_radius,
        "peak_mean_speed": float(speed_array[peak_index]),
        "peak_speed_radius": float(radius_array[peak_index]),
        "proxy_note": (
            "Autonomous z=0 speed-concentration proxy only; no OpenAI numerical core radius "
            "or frame-to-physical-time target is asserted."
        ),
    }


def validate_fingerprint_receipt(receipt: dict[str, Any]) -> None:
    """Fail closed if an autonomous diagnostic is laundered into a source/PDE claim."""
    if receipt.get("schema") != SCHEMA:
        raise ValueError("unexpected cylindrical morphology schema")
    if receipt.get("cylindrical_morphology_diagnostic_ready") is not True:
        raise ValueError("diagnostic readiness must be explicitly true")
    protocol = receipt.get("protocol")
    if not isinstance(protocol, dict):
        raise ValueError("missing morphology protocol")
    for key in ("source_numeric_targets_used", "renderer_or_camera_used", "pixel_loss_used"):
        if protocol.get(key) is not False:
            raise ValueError(f"{key} must remain false")
    for key, expected in TRUTH_BOUNDARY.items():
        if receipt.get(key) is not expected:
            raise ValueError(f"truth boundary promoted: {key}")
    screen = receipt.get("external_method_screen")
    if not isinstance(screen, list) or len(screen) != 1:
        raise ValueError("external method screen drift")
    expected_screen = EXTERNAL_METHOD_SCREEN[0]
    for key in ("repo", "screened_commit", "license", "classification", "scope", "difference"):
        if screen[0].get(key) != expected_screen[key]:
            raise ValueError(f"external method screen drift: {key}")
    measurements = receipt.get("measurements")
    if not isinstance(measurements, dict):
        raise ValueError("missing morphology measurements")
    digest = hashlib.sha256(
        json.dumps(measurements, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    if receipt.get("measurement_sha256") != digest:
        raise ValueError("morphology measurement digest mismatch")


def fingerprint_velocity_field(
    field: VelocityField,
    *,
    model_id: str = "custom",
    radii: Sequence[float] = DEFAULT_RADII,
    z_magnitudes: Sequence[float] = DEFAULT_Z_MAGNITUDES,
    times: Sequence[float] = DEFAULT_TIMES,
    azimuth_count: int = DEFAULT_AZIMUTHS,
    core_radii: Sequence[float] = CORE_RADII,
) -> dict[str, Any]:
    """Build a deterministic, renderer-independent cylindrical morphology receipt.

    The reported quantities are autonomous repository diagnostics.  They are useful
    for later candidate-to-candidate comparisons, but none is an OpenAI-supplied
    numerical target and none is a PDE acceptance statistic.
    """
    radii_tuple = _finite_positive(radii, "radii")
    z_tuple = _finite_positive(z_magnitudes, "z_magnitudes")
    times_tuple = _finite_times(times)
    core_radii_tuple = _finite_positive(core_radii, "core_radii")
    _azimuth_grid(azimuth_count)

    rows: list[dict[str, float | int]] = []
    radial_variation: list[dict[str, Any]] = []
    for time in times_tuple:
        for z_abs in z_tuple:
            for z in (-z_abs, z_abs):
                group: list[dict[str, float | int]] = []
                for radius in radii_tuple:
                    row = _evaluate_ring(
                        field,
                        radius=radius,
                        z=z,
                        time=time,
                        azimuth_count=azimuth_count,
                    )
                    rows.append(row)
                    group.append(row)
                circulation = np.asarray([float(row["mean_abs_u_theta"]) for row in group])
                angular = np.asarray([float(row["mean_angular_rate"]) for row in group])
                radial_variation.append(
                    {
                        "time": float(time),
                        "z": float(z),
                        "circulation_speed_span_across_radii": float(np.ptp(circulation)),
                        "angular_rate_span_across_radii": float(np.ptp(angular)),
                        "all_rings_inward": bool(all(float(row["mean_u_r"]) < 0.0 for row in group)),
                        "all_rings_swirl_nonzero": bool(
                            all(float(row["swirl_nonzero_fraction"]) == 1.0 for row in group)
                        ),
                    }
                )

    core = [
        _core_speed_proxy(
            field,
            time=time,
            radii=core_radii_tuple,
            azimuth_count=azimuth_count,
        )
        for time in times_tuple
    ]
    first_radius = core[0]["speed_weighted_rms_radius"]
    last_radius = core[-1]["speed_weighted_rms_radius"]
    time_evolution = {
        "first_time": float(times_tuple[0]),
        "last_time": float(times_tuple[-1]),
        "speed_weighted_radius_decreased": (
            bool(last_radius < first_radius)
            if first_radius is not None and last_radius is not None
            else None
        ),
        "peak_mean_speed_increased": bool(core[-1]["peak_mean_speed"] > core[0]["peak_mean_speed"]),
        "interpretation": (
            "Autonomous proxy trend only.  It is not a source-correspondence verdict and "
            "does not recover OpenAI camera, scale, or frame-time registration."
        ),
    }

    measurement_payload = {
        "ring_rows": rows,
        "radial_variation": radial_variation,
        "core_speed_proxy": core,
        "time_evolution_proxy": time_evolution,
    }
    measurement_sha256 = hashlib.sha256(
        json.dumps(measurement_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "model_id": str(model_id),
        "source_commit": getattr(field, "source_commit", None),
        "source_mat_sha256": getattr(field, "mat_sha256", None),
        "protocol": {
            "radii": list(radii_tuple),
            "z_magnitudes": list(z_tuple),
            "times": list(times_tuple),
            "azimuth_count": int(azimuth_count),
            "core_radii": list(core_radii_tuple),
            "coordinate_conversion": "u_r=u_x*cos(phi)+u_y*sin(phi); u_theta=-u_x*sin(phi)+u_y*cos(phi)",
            "spiral_proxy": "-u_r/abs(u_theta) where abs(u_theta)>1e-12",
            "source_numeric_targets_used": False,
            "renderer_or_camera_used": False,
            "pixel_loss_used": False,
        },
        "measurements": measurement_payload,
        "measurement_sha256": measurement_sha256,
        "external_method_screen": list(EXTERNAL_METHOD_SCREEN),
        "cylindrical_morphology_diagnostic_ready": True,
        **TRUTH_BOUNDARY,
        "scope": "autonomous_renderer_independent_velocity_kinematics_not_source_or_pde_acceptance",
    }
    validate_fingerprint_receipt(receipt)
    return receipt


def fingerprint_st054(model_id: str) -> dict[str, Any]:
    field = load_st054(model_id)
    return fingerprint_velocity_field(field, model_id=model_id)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=available_st054_models(), default="ST054-Q2")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = fingerprint_st054(args.model)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"model_id": args.model, "measurement_sha256": receipt["measurement_sha256"]}))


if __name__ == "__main__":
    main()
