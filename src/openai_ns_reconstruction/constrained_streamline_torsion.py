from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.interpolate import make_interp_spline


@dataclass(frozen=True)
class StreamlineTorsionReport:
    s_fraction: np.ndarray
    points: np.ndarray
    torsion: np.ndarray
    valid_mask: np.ndarray
    arclength: float
    valid_arclength_fraction: float
    mean_signed_torsion: float
    mean_abs_torsion: float
    rms_torsion: float
    max_abs_torsion: float
    integrated_signed_torsion: float
    integrated_abs_torsion: float
    handedness_coherence: float
    metadata: Mapping[str, object]


def measure_streamline_torsion_profile(
    points: np.ndarray,
    *,
    sample_count: int = 257,
    provenance: str,
) -> StreamlineTorsionReport:
    """Measure intrinsic 3-D torsion of one already-computed streamline.

    The supplied polyline is interpolated exactly with a quintic B-spline in
    cumulative-polyline-length parameter. No smoothing, registration, camera
    fitting, path reversal, or velocity modification is performed. Torsion is
    a visualization-geometry diagnostic only and is never PDE evidence.
    """
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 3 or p.shape[0] < 6:
        raise ValueError("points must have shape (n,3) with n>=6")
    if not np.all(np.isfinite(p)):
        raise ValueError("points must be finite")
    if not isinstance(sample_count, int) or sample_count < 65 or sample_count % 2 == 0:
        raise ValueError("sample_count must be an odd integer >=65")
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("nonempty provenance is required")

    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    if np.any(~np.isfinite(seg)) or np.any(seg <= 0.0):
        raise ValueError("consecutive duplicate/zero-length samples are not allowed")
    s = np.concatenate(([0.0], np.cumsum(seg)))
    parameter_length = float(s[-1])
    if not np.isfinite(parameter_length) or parameter_length <= 0.0:
        raise ValueError("streamline must have positive length")

    q = np.linspace(0.0, parameter_length, sample_count)
    spline = make_interp_spline(s, p, k=5, axis=0, check_finite=True)
    xyz = np.asarray(spline(q), dtype=float)
    d1 = np.asarray(spline(q, 1), dtype=float)
    d2 = np.asarray(spline(q, 2), dtype=float)
    d3 = np.asarray(spline(q, 3), dtype=float)
    if not all(np.all(np.isfinite(arr)) for arr in (xyz, d1, d2, d3)):
        raise ValueError("nonfinite spline derivative")

    speed = np.linalg.norm(d1, axis=1)
    if np.any(speed <= 0.0):
        raise ValueError("interpolated curve has a degenerate tangent")
    cross = np.cross(d1, d2)
    denom = np.einsum("ij,ij->i", cross, cross)
    denom_scale = float(np.max(denom))
    if not np.isfinite(denom_scale) or denom_scale <= 0.0:
        raise ValueError("torsion is undefined for a zero-curvature curve")

    numerical_floor = np.finfo(float).eps * max(1.0, denom_scale) * 1024.0
    valid = denom > numerical_floor
    if not np.any(valid):
        raise ValueError("torsion is numerically undefined on the sampled curve")

    numerator = np.einsum("ij,ij->i", cross, d3)
    torsion = np.zeros(sample_count, dtype=float)
    torsion[valid] = numerator[valid] / denom[valid]
    if not np.all(np.isfinite(torsion[valid])):
        raise ValueError("nonfinite torsion")

    def trapz(values: np.ndarray) -> float:
        return float(np.sum(0.5 * (values[1:] + values[:-1]) * np.diff(q)))

    # Integrate against physical spline arclength. Invalid near-zero-curvature
    # samples contribute no torsion and are reported separately, not imputed.
    valid_weight = valid.astype(float) * speed
    total_arclength = trapz(speed)
    valid_arclength = trapz(valid_weight)
    if not np.isfinite(total_arclength) or total_arclength <= 0.0:
        raise ValueError("interpolated curve has nonpositive arclength")
    if not np.isfinite(valid_arclength) or valid_arclength <= 0.0:
        raise ValueError("torsion has no valid arclength support")

    signed_integrand = torsion * valid_weight
    abs_integrand = np.abs(torsion) * valid_weight
    square_integrand = torsion * torsion * valid_weight
    integrated_signed = trapz(signed_integrand)
    integrated_abs = trapz(abs_integrand)
    rms = float(np.sqrt(trapz(square_integrand) / valid_arclength))
    mean_signed = integrated_signed / valid_arclength
    mean_abs = integrated_abs / valid_arclength
    max_abs = float(np.max(np.abs(torsion[valid])))
    coherence = 0.0 if integrated_abs == 0.0 else abs(integrated_signed) / integrated_abs

    s_fraction = q / parameter_length
    for arr in (s_fraction, xyz, torsion, valid):
        arr.setflags(write=False)

    metadata = {
        "claim_scope": "visualization_streamline_geometry_only",
        "interpolation": "scipy.interpolate.make_interp_spline exact quintic B-spline",
        "smoothing_applied": False,
        "parameterization": "cumulative input-polyline length",
        "torsion_formula": "dot(cross(r1,r2),r3)/norm(cross(r1,r2))^2",
        "near_zero_curvature_policy": "exclude only machine-scale denominator samples and report valid_arclength_fraction",
        "camera_or_registration_fit": False,
        "path_reversal_or_reflection_fit": False,
        "velocity_changed": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "provenance": provenance,
    }
    return StreamlineTorsionReport(
        s_fraction=s_fraction,
        points=xyz,
        torsion=torsion,
        valid_mask=valid,
        arclength=total_arclength,
        valid_arclength_fraction=valid_arclength / total_arclength,
        mean_signed_torsion=mean_signed,
        mean_abs_torsion=mean_abs,
        rms_torsion=rms,
        max_abs_torsion=max_abs,
        integrated_signed_torsion=integrated_signed,
        integrated_abs_torsion=integrated_abs,
        handedness_coherence=coherence,
        metadata=metadata,
    )
