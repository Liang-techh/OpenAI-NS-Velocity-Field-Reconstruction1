"""Truth-bounded center-plane radial circulation profile for retained ST006.

This module is a visualization diagnostic only. It samples the published
``velocity(x,y,z,t)->[u,v,w]`` interface on fixed physical rings; it does not
fit the candidate, a camera, a hidden frame time, pressure, or forcing.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import simpson

from . import CANDIDATE_SHA256, load_best

SCHEMA = "st006_centerplane_radial_circulation_v1"
REFERENCE_TIMES = np.array([0.25, 0.50, 0.75], dtype=float)
RADII = np.linspace(0.05, 1.50, 31, dtype=float)
N_THETA = 64
Z_PLANE = 0.0
PROVENANCE = "retained_ST006_public_velocity_api_centerplane_fixed_rings"

_FALSE_TRUTH_FLAGS = {
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
    "hidden_frame_time_inferred": False,
    "camera_fitted": False,
    "visual_acceptance_threshold_selected": False,
    "used_for_pde_acceptance": False,
}


def _jsonable_float_array(values: np.ndarray) -> list[float]:
    return [float(x) for x in np.asarray(values, dtype=float).tolist()]


def _report_digest(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def _sample_centerplane_profile(field: Any) -> dict[str, Any]:
    if getattr(field, "sha256", None) != CANDIDATE_SHA256:
        raise ValueError("radial-circulation audit is bound to retained ST006 identity")

    theta = 2.0 * np.pi * np.arange(N_THETA, dtype=float) / N_THETA
    c = np.cos(theta)
    s = np.sin(theta)
    rr = np.repeat(RADII, N_THETA)
    cc = np.tile(c, RADII.size)
    ss = np.tile(s, RADII.size)
    points = np.column_stack((rr * cc, rr * ss, np.full(rr.shape, Z_PLANE)))

    per_time: list[dict[str, Any]] = []
    for t in REFERENCE_TIMES:
        uvw = np.asarray(field.at_points(points, float(t)), dtype=float)
        if uvw.shape != points.shape:
            raise ValueError("velocity output must have shape (n,3)")
        if not np.isfinite(uvw).all():
            raise ValueError("velocity output must be finite")

        u = uvw[:, 0].reshape(RADII.size, N_THETA)
        v = uvw[:, 1].reshape(RADII.size, N_THETA)
        w = uvw[:, 2].reshape(RADII.size, N_THETA)
        C = c[None, :]
        S = s[None, :]
        ur = u * C + v * S
        utheta = -u * S + v * C
        speed = np.sqrt(u * u + v * v + w * w)

        mean_ur = np.mean(ur, axis=1)
        mean_utheta = np.mean(utheta, axis=1)
        mean_abs_utheta = np.mean(np.abs(utheta), axis=1)
        mean_uz = np.mean(w, axis=1)
        mean_speed = np.mean(speed, axis=1)
        circulation = 2.0 * np.pi * RADII * mean_utheta

        component_means = np.stack((mean_ur, mean_utheta, mean_uz), axis=1)
        cyl = np.stack((ur, utheta, w), axis=2)
        ring_deviation = np.sqrt(np.mean((cyl - component_means[:, None, :]) ** 2, axis=(1, 2)))
        sampled_velocity_rms = float(np.sqrt(np.mean(uvw * uvw)))
        if not np.isfinite(sampled_velocity_rms) or sampled_velocity_rms <= 0.0:
            raise ValueError("exact-zero/nonfinite sampled velocity is not a valid visualization audit")

        swirl_integral = float(simpson(mean_abs_utheta, x=RADII))
        swirl_roundoff_floor = float(
            128.0
            * np.finfo(float).eps
            * (RADII[-1] - RADII[0])
            * np.max(mean_speed)
        )
        if not np.isfinite(swirl_integral) or swirl_integral <= swirl_roundoff_floor:
            raise ValueError("center-plane sampled swirl is zero/nonfinite at floating-point scale")
        swirl_weighted_radius = float(
            simpson(RADII * mean_abs_utheta, x=RADII) / swirl_integral
        )

        peak_swirl_idx = int(np.argmax(mean_abs_utheta))
        peak_speed_idx = int(np.argmax(mean_speed))
        relative_range = float(
            (np.max(mean_abs_utheta) - np.min(mean_abs_utheta)) / np.max(mean_abs_utheta)
        )

        per_time.append(
            {
                "time": float(t),
                "mean_radial_velocity": _jsonable_float_array(mean_ur),
                "mean_circulating_velocity": _jsonable_float_array(mean_utheta),
                "mean_abs_circulating_speed": _jsonable_float_array(mean_abs_utheta),
                "mean_axial_velocity": _jsonable_float_array(mean_uz),
                "mean_speed": _jsonable_float_array(mean_speed),
                "circulation": _jsonable_float_array(circulation),
                "ring_component_deviation_rms": _jsonable_float_array(ring_deviation),
                "sampled_velocity_rms": sampled_velocity_rms,
                "max_ring_component_deviation_over_velocity_rms": float(
                    np.max(ring_deviation) / sampled_velocity_rms
                ),
                "peak_abs_circulating_speed": float(mean_abs_utheta[peak_swirl_idx]),
                "peak_abs_circulating_speed_radius": float(RADII[peak_swirl_idx]),
                "swirl_weighted_radius": swirl_weighted_radius,
                "circulating_speed_relative_radial_range": relative_range,
                "peak_mean_speed": float(mean_speed[peak_speed_idx]),
                "peak_mean_speed_radius": float(RADII[peak_speed_idx]),
                "inward_mean_radial_ring_count": int(np.count_nonzero(mean_ur < 0.0)),
                "outward_mean_radial_ring_count": int(np.count_nonzero(mean_ur > 0.0)),
                "zero_sign_mean_radial_ring_count": int(np.count_nonzero(mean_ur == 0.0)),
            }
        )

    first = per_time[0]
    last = per_time[-1]
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "candidate": "ST006",
        "candidate_sha256": CANDIDATE_SHA256,
        "provenance": PROVENANCE,
        "measurement_contract": {
            "plane": "z=0",
            "times": _jsonable_float_array(REFERENCE_TIMES),
            "radii": _jsonable_float_array(RADII),
            "theta_samples_per_ring": N_THETA,
            "theta_rule": "uniform endpoint-excluded [0,2*pi)",
            "radial_integration": "scipy.integrate.simpson on the fixed radius grid",
            "coordinate_frame": "registered physical Cartesian frame; no rotation/registration",
            "target_values_or_visual_thresholds": None,
        },
        "per_time": per_time,
        "endpoint_deltas": {
            "peak_abs_circulating_speed": float(
                last["peak_abs_circulating_speed"] - first["peak_abs_circulating_speed"]
            ),
            "peak_abs_circulating_speed_radius": float(
                last["peak_abs_circulating_speed_radius"]
                - first["peak_abs_circulating_speed_radius"]
            ),
            "swirl_weighted_radius": float(
                last["swirl_weighted_radius"] - first["swirl_weighted_radius"]
            ),
            "peak_mean_speed": float(last["peak_mean_speed"] - first["peak_mean_speed"]),
            "peak_mean_speed_radius": float(
                last["peak_mean_speed_radius"] - first["peak_mean_speed_radius"]
            ),
        },
        "interpretation_boundary": {
            "candidate_side_only": True,
            "may_support_public_observable_review": [
                "circulating speed depends on radius",
                "central-region contraction/speed trend, if independently paired to public evidence",
            ],
            "does_not_infer": [
                "OpenAI numerical velocity",
                "OpenAI physical radius scale",
                "OpenAI hidden frame times",
                "OpenAI camera pose",
                "OpenAI streamline seeds",
            ],
            **_FALSE_TRUTH_FLAGS,
        },
    }
    payload["report_sha256"] = _report_digest(payload)
    return payload


def measure_retained_st006_radial_circulation() -> dict[str, Any]:
    """Measure the frozen retained ST006 field on the preregistered ring contract."""
    return _sample_centerplane_profile(load_best())


def write_report(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing report: {path}")
    report = measure_retained_st006_radial_circulation()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    report = (
        write_report(args.output)
        if args.output is not None
        else measure_retained_st006_radial_circulation()
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
