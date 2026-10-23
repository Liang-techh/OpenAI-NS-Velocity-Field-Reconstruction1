from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass(frozen=True)
class StreamlineCurvatureReport:
    s_fraction: np.ndarray
    points: np.ndarray
    curvature: np.ndarray
    arclength: float
    mean_curvature: float
    rms_curvature: float
    max_curvature: float
    total_turning: float
    endpoint_chord: float
    tortuosity: float
    metadata: Mapping[str, object]


def measure_streamline_curvature_profile(
    points: np.ndarray,
    *,
    sample_count: int = 129,
    provenance: str,
) -> StreamlineCurvatureReport:
    """Measure fixed-geometry curvature of one already-computed 3-D streamline.

    This is a visualization diagnostic only. It does not integrate a streamline,
    fit a camera, alter velocity, smooth data, or define a visual pass threshold.
    """
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 3 or p.shape[0] < 5:
        raise ValueError("points must have shape (n,3) with n>=5")
    if not np.all(np.isfinite(p)):
        raise ValueError("points must be finite")
    if not isinstance(sample_count, int) or sample_count < 33 or sample_count % 2 == 0:
        raise ValueError("sample_count must be an odd integer >=33")
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("nonempty provenance is required")

    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    if np.any(seg <= 0.0):
        raise ValueError("consecutive duplicate/zero-length samples are not allowed")
    s = np.concatenate(([0.0], np.cumsum(seg)))
    parameter_length = float(s[-1])
    chord = float(np.linalg.norm(p[-1] - p[0]))
    if not np.isfinite(parameter_length) or parameter_length <= 0.0 or chord <= 0.0:
        raise ValueError("streamline must have positive length and endpoint chord")

    q = np.linspace(0.0, parameter_length, sample_count)
    xyz = np.empty((sample_count, 3), dtype=float)
    d1 = np.empty_like(xyz)
    d2 = np.empty_like(xyz)
    for j in range(3):
        cs = CubicSpline(s, p[:, j], bc_type="natural", extrapolate=False)
        xyz[:, j] = cs(q)
        d1[:, j] = cs(q, 1)
        d2[:, j] = cs(q, 2)

    speed = np.linalg.norm(d1, axis=1)
    if np.any(~np.isfinite(speed)) or np.min(speed) <= 1e-12:
        raise ValueError("interpolated curve has a degenerate tangent")
    cross = np.cross(d1, d2)
    kappa = np.linalg.norm(cross, axis=1) / speed**3
    if not np.all(np.isfinite(kappa)):
        raise ValueError("nonfinite curvature")

    dq = np.diff(q)
    ds_weight = speed
    arclength = float(np.sum(0.5 * (ds_weight[1:] + ds_weight[:-1]) * dq))
    if not np.isfinite(arclength) or arclength <= 0.0:
        raise ValueError("interpolated curve has nonpositive arclength")
    turning_integrand = kappa * ds_weight
    total_turning = float(
        np.sum(0.5 * (turning_integrand[1:] + turning_integrand[:-1]) * dq)
    )
    kappa_sq_integrand = kappa * kappa * ds_weight
    kappa_sq_integral = float(
        np.sum(0.5 * (kappa_sq_integrand[1:] + kappa_sq_integrand[:-1]) * dq)
    )
    rms = float(np.sqrt(kappa_sq_integral / arclength))
    mean = total_turning / arclength
    max_k = float(np.max(kappa))

    s_fraction = q / parameter_length
    for arr in (s_fraction, xyz, kappa):
        arr.setflags(write=False)

    metadata = {
        "claim_scope": "visualization_streamline_geometry_only",
        "interpolation": "scipy.interpolate.CubicSpline natural exact interpolation",
        "smoothing_applied": False,
        "camera_or_registration_fit": False,
        "velocity_changed": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "provenance": provenance,
    }
    return StreamlineCurvatureReport(
        s_fraction=s_fraction,
        points=xyz,
        curvature=kappa,
        arclength=arclength,
        mean_curvature=mean,
        rms_curvature=rms,
        max_curvature=max_k,
        total_turning=total_turning,
        endpoint_chord=chord,
        tortuosity=arclength / chord,
        metadata=metadata,
    )
