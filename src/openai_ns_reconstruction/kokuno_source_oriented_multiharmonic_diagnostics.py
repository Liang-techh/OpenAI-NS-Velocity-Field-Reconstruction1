"""Three-resolution covariance diagnostics for oriented Kokuno complete curls.

K2-OSC-108 adds one bounded repository-autonomous azimuthal frame orientation
around the provider-driven, spatial/pulse-support-certified complete-curl family.
This increment does not change that velocity field.  It consumes only the public
Cartesian ``velocity(x,y,z,t)`` surface and checks, with an implementation-
distinct centered FD2 operator, that the oriented field behaves consistently
with rigid z-axis covariance at three frozen source-chart-equivalent resolutions.

For corresponding parent/lab probes it records

* Cartesian divergence of the oriented and alpha=0 fields;
* numerical divergence invariance under the rigid frame rotation;
* full Cartesian gradient covariance ``G_alpha = Q G_0 Q^T``;
* vorticity-vector covariance ``omega_alpha = Q omega_0``;
* off-grid public-velocity covariance; and
* exact-zero behavior in the common certified axis core and outside every
  certified radial support envelope.

The derivative ladder is fixed at ``0.02 / 0.01 / 0.005`` in source-chart units
and mapped to physical Cartesian steps using the corrected Q scaling.  No
caller supplies derivative steps, scientific gates, residuals, pressure,
forcing, gains, optimizers, or tolerances.

This is numerical covariance/derivative evidence only.  It is not recovery of
the source-exact Kokuno frame vectors/orientation and it does not assess a
complete Navier--Stokes residual.
"""
from __future__ import annotations

import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from . import kokuno_source_azimuthal_frame_oriented_multiharmonic_velocity as _parent_module
from .kokuno_source_azimuthal_frame_oriented_multiharmonic_velocity import (
    KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity,
    azimuthal_frame_orientation_contract,
)
from .kokuno_source_physicalized_localized_harmonic import source_chart_to_physical_rzt

TASK = "K2-OSC-109"
SCHEMA = "kokuno-a2-oriented-multiharmonic-three-resolution-diagnostics-v1"
PARENT_A2_PR = 1170
PARENT_A2_HEAD = "59d3e6d64fcdfeaf145cbee67e68b8a25e89489c"
PARENT_A2_SOURCE_BLOB = "ca8acae2c3a53b79528cb55f45502443f6ffb2b6"

CHART_FD_STEPS = (2.0e-2, 1.0e-2, 5.0e-3)

# These frozen points lie well inside the manufactured provider used by the
# exact-head receipt.  For arbitrary external providers they are observations,
# not an assertion that every term must be active at every probe.
_R = np.asarray((0.30, 0.36, 0.43, 0.51, 0.59, 0.68), dtype=float)
_THETA = np.asarray((-1.08, -0.53, 0.08, 0.69, 1.31, 1.88), dtype=float)
_Z = np.asarray((-0.18, -0.11, -0.03, 0.05, 0.13, 0.22), dtype=float)
_T = np.asarray((0.31, 0.38, 0.45, 0.52, 0.59, 0.66), dtype=float)
_OFFGRID_SHIFT = {"R": 0.017, "theta": 0.043, "Z": -0.019, "T": 0.011}


def _git_blob_sha1(path: str | Path) -> str:
    data = Path(path).read_bytes()
    prefix = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(prefix + data).hexdigest()


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _parent_identity() -> None:
    path = getattr(_parent_module, "__file__", None)
    if not path or _git_blob_sha1(path) != PARENT_A2_SOURCE_BLOB:
        raise RuntimeError("exact K2-OSC-108 parent source blob drifted")
    contract = azimuthal_frame_orientation_contract()
    if contract.get("bounded_azimuthal_frame_orientation_materialized") is not True:
        raise RuntimeError("K2-OSC-108 orientation prerequisite is unavailable")
    if contract.get("self_contained_velocity_xyzt_provider") is not False:
        raise RuntimeError("K2-OSC-108 source-provider truth boundary drifted")


def _validate_field(
    field: Any,
) -> KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity:
    _parent_identity()
    if not isinstance(field, KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity):
        raise ValueError(
            "field must be KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity"
        )
    semantic = str(getattr(field, "semantic_sha256", ""))
    if len(semantic) != 64 or any(c not in "0123456789abcdef" for c in semantic):
        raise ValueError("field semantic identity must be a lowercase 64-hex digest")
    return field


def _rotation_z(alpha: float) -> np.ndarray:
    c = math.cos(float(alpha))
    s = math.sin(float(alpha))
    return np.asarray(((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0)), dtype=float)


def _rotate_points(x: np.ndarray, y: np.ndarray, alpha: float) -> tuple[np.ndarray, np.ndarray]:
    c = math.cos(float(alpha))
    s = math.sin(float(alpha))
    return c * x - s * y, s * x + c * y


def _parent_probe_coordinates(field, *, offgrid: bool):
    R = _R.copy()
    theta = _THETA.copy()
    Z = _Z.copy()
    T = _T.copy()
    if offgrid:
        R += _OFFGRID_SHIFT["R"]
        theta += _OFFGRID_SHIFT["theta"]
        Z += _OFFGRID_SHIFT["Z"]
        T += _OFFGRID_SHIFT["T"]
    scaling = field.parent.scaling
    r, z, t = source_chart_to_physical_rzt(R, Z, T, ell=scaling.ell, h=scaling.h)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, z, t, R, theta, Z, T


def _velocity(field, x, y, z, t) -> np.ndarray:
    value = np.asarray(field.velocity(x, y, z, t), dtype=float)
    expected = np.broadcast(np.asarray(x), np.asarray(y), np.asarray(z), np.asarray(t)).shape + (3,)
    if value.shape != expected:
        raise RuntimeError(f"public velocity returned shape {value.shape}, expected {expected}")
    if not np.all(np.isfinite(value)):
        raise RuntimeError("public velocity returned nonfinite values")
    return value


def _fd_gradient(field, x, y, z, t, *, chart_step: float) -> np.ndarray:
    scaling = field.parent.scaling
    h_xy = float(chart_step) * float(scaling.radial_scale)
    h_z = float(chart_step) * float(scaling.axial_scale)
    if not all(math.isfinite(v) and v > 0.0 for v in (h_xy, h_z)):
        raise RuntimeError("physical finite-difference step is not finite and positive")
    for coordinate, step in ((x, h_xy), (y, h_xy), (z, h_z)):
        if (
            np.any(coordinate + step == coordinate)
            or np.any(coordinate - step == coordinate)
            or np.any(coordinate + step == coordinate - step)
        ):
            raise RuntimeError("finite-difference perturbation collapsed in binary64")

    du_dx = (_velocity(field, x + h_xy, y, z, t) - _velocity(field, x - h_xy, y, z, t)) / (2.0 * h_xy)
    du_dy = (_velocity(field, x, y + h_xy, z, t) - _velocity(field, x, y - h_xy, z, t)) / (2.0 * h_xy)
    du_dz = (_velocity(field, x, y, z + h_z, t) - _velocity(field, x, y, z - h_z, t)) / (2.0 * h_z)
    return np.stack((du_dx, du_dy, du_dz), axis=-1)


def _div_vort(gradient: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if gradient.ndim != 3 or gradient.shape[1:] != (3, 3):
        raise RuntimeError("Cartesian gradient must have shape (n,3,3)")
    div = gradient[:, 0, 0] + gradient[:, 1, 1] + gradient[:, 2, 2]
    vort = np.stack(
        (
            gradient[:, 2, 1] - gradient[:, 1, 2],
            gradient[:, 0, 2] - gradient[:, 2, 0],
            gradient[:, 1, 0] - gradient[:, 0, 1],
        ),
        axis=-1,
    )
    return div, vort


def _rms(array: np.ndarray) -> float:
    arr = np.asarray(array, dtype=float)
    return float(np.sqrt(np.mean(arr * arr)))


def _vector_rms(array: np.ndarray) -> float:
    arr = np.asarray(array, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def _tensor_rms(array: np.ndarray) -> float:
    arr = np.asarray(array, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=(-2, -1)))))


def _resolution_record(g0: np.ndarray, ga: np.ndarray, q: np.ndarray, chart_step: float, field) -> dict[str, Any]:
    d0, w0 = _div_vort(g0)
    da, wa = _div_vort(ga)
    g_expected = np.einsum("ij,njk,kl->nil", q, g0, q.T)
    w_expected = np.einsum("ij,nj->ni", q, w0)
    grad_error = ga - g_expected
    vort_error = wa - w_expected
    div_error = da - d0

    grad_scale = max(_tensor_rms(g_expected), np.finfo(float).tiny)
    vort_scale = max(_vector_rms(w_expected), np.finfo(float).tiny)
    div_scale = max(_tensor_rms(g_expected), _tensor_rms(ga), np.finfo(float).tiny)
    vort_mag = np.linalg.norm(wa, axis=-1)
    w_rms = _vector_rms(wa)
    axial = _rms(wa[:, 2])
    scaling = field.parent.scaling
    return {
        "chart_step": float(chart_step),
        "physical_xy_step": float(chart_step * scaling.radial_scale),
        "physical_z_step": float(chart_step * scaling.axial_scale),
        "reference_divergence_max_abs": float(np.max(np.abs(d0))),
        "reference_divergence_rms": _rms(d0),
        "oriented_divergence_max_abs": float(np.max(np.abs(da))),
        "oriented_divergence_rms": _rms(da),
        "divergence_covariance_rms": _rms(div_error),
        "divergence_covariance_normalized_rms": _rms(div_error) / div_scale,
        "gradient_covariance_relative_rms": _tensor_rms(grad_error) / grad_scale,
        "vorticity_covariance_relative_rms": _vector_rms(vort_error) / vort_scale,
        "oriented_vorticity_rms": w_rms,
        "oriented_vorticity_max": float(np.max(vort_mag)),
        "oriented_axial_vorticity_rms_fraction": axial / max(w_rms, np.finfo(float).tiny),
    }


def _support_observation(field) -> dict[str, Any]:
    terms = field.parent.terms
    if not terms:
        raise RuntimeError("oriented parent must contain at least one harmonic term")
    min_core = min(float(term.spatial.harmonic.provider.axis_zero_radius_chart) for term in terms)
    max_outer = max(float(term.spatial.radial_support_max_chart) for term in terms)
    scaling = field.parent.scaling

    R_axis = np.asarray((0.0, 0.5 * min_core), dtype=float)
    z_axis_chart = np.zeros_like(R_axis)
    T_axis = np.full_like(R_axis, 0.5)
    r_axis, z_axis, t_axis = source_chart_to_physical_rzt(
        R_axis, z_axis_chart, T_axis, ell=scaling.ell, h=scaling.h
    )
    x_axis, y_axis = _rotate_points(r_axis, np.zeros_like(r_axis), field.frame_azimuth)
    axis_u = _velocity(field, x_axis, y_axis, z_axis, t_axis)

    R_outer = np.asarray((max_outer + 0.08, max_outer + 0.17), dtype=float)
    z_outer_chart = np.asarray((0.0, 0.11), dtype=float)
    T_outer = np.asarray((0.47, 0.61), dtype=float)
    r_outer, z_outer, t_outer = source_chart_to_physical_rzt(
        R_outer, z_outer_chart, T_outer, ell=scaling.ell, h=scaling.h
    )
    theta_outer = np.asarray((0.37, -1.12), dtype=float)
    xp = r_outer * np.cos(theta_outer)
    yp = r_outer * np.sin(theta_outer)
    x_outer, y_outer = _rotate_points(xp, yp, field.frame_azimuth)
    outer_u = _velocity(field, x_outer, y_outer, z_outer, t_outer)
    return {
        "common_axis_core_exact_zero": bool(np.array_equal(axis_u, np.zeros_like(axis_u))),
        "axis_probe_R": R_axis.tolist(),
        "outer_spatial_support_exact_zero": bool(np.array_equal(outer_u, np.zeros_like(outer_u))),
        "outer_probe_R": R_outer.tolist(),
        "outer_support_basis": "strictly outside max provider-certified radial_support_max_chart",
    }


def audit_oriented_multiharmonic_field(
    field: KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity,
) -> dict[str, Any]:
    """Run the frozen public-velocity rigid-covariance diagnostic.

    The only scientific input is the already-constructed field.  No derivative
    step, gate, pressure, forcing, residual, optimizer, or tolerance is caller
    configurable.
    """
    field = _validate_field(field)
    alpha = field.frame_azimuth
    q = _rotation_z(alpha)
    reference = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        field.parent, frame_azimuth=0.0
    )

    xp, yp, z, t, R, theta, Z, T = _parent_probe_coordinates(field, offgrid=False)
    x, y = _rotate_points(xp, yp, alpha)
    u0 = _velocity(reference, xp, yp, z, t)
    ua = _velocity(field, x, y, z, t)
    u_expected = np.einsum("ij,nj->ni", q, u0)
    base_velocity_covariance = ua - u_expected

    resolutions: list[dict[str, Any]] = []
    for step in CHART_FD_STEPS:
        g0 = _fd_gradient(reference, xp, yp, z, t, chart_step=step)
        ga = _fd_gradient(field, x, y, z, t, chart_step=step)
        resolutions.append(_resolution_record(g0, ga, q, step, field))

    xpo, ypo, zo, to, Ro, thetao, Zo, To = _parent_probe_coordinates(field, offgrid=True)
    xo, yo = _rotate_points(xpo, ypo, alpha)
    u0o = _velocity(reference, xpo, ypo, zo, to)
    uao = _velocity(field, xo, yo, zo, to)
    uao_expected = np.einsum("ij,nj->ni", q, u0o)
    offgrid_error = uao - uao_expected

    support = _support_observation(field)
    if not support["common_axis_core_exact_zero"]:
        raise RuntimeError("provider-certified common axis core is not exact zero")
    if not support["outer_spatial_support_exact_zero"]:
        raise RuntimeError("provider-certified radial exterior is not exact zero")

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "parent": {
            "pr": PARENT_A2_PR,
            "head": PARENT_A2_HEAD,
            "source_blob_sha1": PARENT_A2_SOURCE_BLOB,
        },
        "field_semantic_sha256": field.semantic_sha256,
        "reference_alpha0_semantic_sha256": reference.semantic_sha256,
        "frame_azimuth": float(alpha),
        "protocol": {
            "chart_fd_steps": list(CHART_FD_STEPS),
            "probe_count": int(R.size),
            "correspondence": "lab probes are Q(alpha) times frozen parent-frame probes",
            "derivatives_consume_public_velocity_only": True,
            "scientific_acceptance_gate": None,
            "offgrid_shift": dict(_OFFGRID_SHIFT),
        },
        "probe_chart_parent_frame": {
            "R": R.tolist(),
            "theta": theta.tolist(),
            "Z": Z.tolist(),
            "T": T.tolist(),
        },
        "velocity_covariance_observation": {
            "base_probe_relative_rms": _vector_rms(base_velocity_covariance)
            / max(_vector_rms(u_expected), np.finfo(float).tiny),
            "base_probe_max_abs": float(np.max(np.abs(base_velocity_covariance))),
            "reference_speed_rms": _vector_rms(u0),
            "oriented_speed_rms": _vector_rms(ua),
            "nonzero_oriented_probe_count": int(np.count_nonzero(np.linalg.norm(ua, axis=-1) > 0.0)),
        },
        "resolutions": resolutions,
        "offgrid_observation": {
            "finite": bool(np.all(np.isfinite(uao))),
            "velocity_covariance_relative_rms": _vector_rms(offgrid_error)
            / max(_vector_rms(uao_expected), np.finfo(float).tiny),
            "velocity_covariance_max_abs": float(np.max(np.abs(offgrid_error))),
            "chart_R": Ro.tolist(),
            "chart_theta": thetao.tolist(),
            "chart_Z": Zo.tolist(),
            "chart_T": To.tolist(),
        },
        "support_observation": support,
        "truth_boundary": oriented_multiharmonic_diagnostics_contract(),
    }
    report["diagnostic_sha256"] = _sha256(report)
    return report


def oriented_multiharmonic_diagnostics_contract() -> dict[str, Any]:
    parent = azimuthal_frame_orientation_contract()
    if parent["bounded_azimuthal_frame_orientation_materialized"] is not True:
        raise RuntimeError("K2-OSC-108 orientation prerequisite missing")
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_a2_pr": PARENT_A2_PR,
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_a2_source_blob": PARENT_A2_SOURCE_BLOB,
        "oriented_public_velocity_three_resolution_diagnostic_materialized": True,
        "public_velocity_only_cartesian_derivatives_materialized": True,
        "rigid_rotation_divergence_covariance_numerically_assessed": True,
        "rigid_rotation_gradient_covariance_numerically_assessed": True,
        "rigid_rotation_vorticity_covariance_numerically_assessed": True,
        "offgrid_rigid_velocity_covariance_assessed": True,
        "provider_certified_axis_zero_support_checked": True,
        "provider_certified_outer_spatial_support_checked": True,
        "analytic_divergence_identity_proved_by_this_increment": False,
        "source_exact_frame_vectors_recovered": False,
        "source_exact_orientation_recovered": False,
        "orientation_classification": "repository_autonomous_candidate_design",
        "self_contained_velocity_xyzt_provider": False,
        "complete_ns_residual_assessed": False,
        "oscillation_before_after_ns_residual_compared": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def public_contract() -> dict[str, Any]:
    params = set(inspect.signature(audit_oriented_multiharmonic_field).parameters)
    forbidden = {
        "step", "steps", "fd_step", "residual", "forcing", "pressure", "viscosity",
        "gain", "threshold", "target", "tolerance", "rtol", "atol", "optimizer",
        "frame_azimuth",
    }
    return {
        "audit_inputs": list(inspect.signature(audit_oriented_multiharmonic_field).parameters),
        "forbidden_audit_inputs_present": sorted(params & forbidden),
        "chart_fd_steps_frozen": list(CHART_FD_STEPS),
        "paper_exact": False,
        "pde_validated": False,
    }


__all__ = [
    "CHART_FD_STEPS",
    "audit_oriented_multiharmonic_field",
    "oriented_multiharmonic_diagnostics_contract",
    "public_contract",
]
