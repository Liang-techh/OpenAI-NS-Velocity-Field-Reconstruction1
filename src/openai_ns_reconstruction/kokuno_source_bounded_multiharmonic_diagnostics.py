"""Three-resolution diagnostics for the bounded provider-driven Kokuno curl family.

This Kokuno Agent-2 increment changes no velocity coefficients or source inputs.
It consumes the public Cartesian ``velocity(x,y,z,t)`` surface from K2-OSC-104
and evaluates an implementation-distinct centered-difference diagnostic at a
small frozen set of off-axis source-chart probes.

The diagnostic reports

* Cartesian divergence at three fixed chart-equivalent spatial resolutions;
* Cartesian vorticity magnitude and axial-fraction observations;
* medium/fine vorticity drift as an off-grid/resolution stability observation;
* exact-zero behavior inside the provider-certified axis core; and
* a deterministic off-grid velocity sample.

The finite-difference ladder is frozen in source-chart units and mapped to
physical x/y/z steps with the corrected Q scaling.  No caller supplies steps,
thresholds, residuals, pressure, forcing, optimization targets, or gains.

This is a numerical diagnostic, not an analytic proof that every caller
provider is a curl.  It does not assess an outer support boundary because the
current provider contract does not expose one.  It also does not evaluate a
complete Navier--Stokes residual: the corrected background/forcing/support
provider, matched pressure, and preregistered restricted forcing remain
missing.  Nothing here is paper-exact source recovery.
"""
from __future__ import annotations

import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from . import kokuno_source_bounded_multiharmonic_velocity as _parent
from .kokuno_source_bounded_multiharmonic_velocity import (
    KokunoBoundedMultiHarmonicPhysicalVelocity,
    bounded_multiharmonic_contract,
)
from .kokuno_source_physicalized_localized_harmonic import source_chart_to_physical_rzt

TASK = "K2-OSC-105"
SCHEMA = "kokuno-a2-bounded-multiharmonic-three-resolution-diagnostics-v1"
PARENT_A2_PR = 1141
PARENT_A2_HEAD = "4241b0b4fcdc36e750798133101baaf803ba8c44"
PARENT_A2_SOURCE_BLOB = "06bc24722b4219419e91cd2d5b08cb1294665f24"

# Fixed source-chart-equivalent derivative ladder.  x/y use Q^(1/2)*h and z
# uses Q^D*h, so the public derivative is always with respect to physical
# Cartesian coordinates even though the numerical resolution is specified in
# the normalized source chart.
CHART_FD_STEPS = (2.0e-2, 1.0e-2, 5.0e-3)

# Frozen dimensionless probe layout.  Radii are offsets from the largest
# provider-certified axis-zero core so no centered stencil approaches the
# singular cylindrical axis.  The layout is diagnostic-only and is not a
# source numerical target.
_R_OFFSETS = np.asarray((0.35, 0.48, 0.63, 0.81, 1.02, 1.24, 0.57, 0.93), dtype=float)
_THETA = np.asarray((0.31, 0.79, 1.27, 1.74, 2.19, 2.63, 1.03, 2.91), dtype=float)
_Z = np.asarray((-0.27, -0.16, -0.05, 0.06, 0.17, 0.28, 0.11, -0.22), dtype=float)
_T = np.asarray((0.72, 0.81, 0.90, 0.99, 1.08, 1.17, 0.86, 1.12), dtype=float)

_OFFGRID_SHIFT = {
    "R": 0.037,
    "theta": 0.071,
    "Z": -0.029,
    "T": 0.013,
}


def _git_blob_sha1(path: str | Path) -> str:
    data = Path(path).read_bytes()
    prefix = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(prefix + data).hexdigest()


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _parent_identity() -> None:
    path = getattr(_parent, "__file__", None)
    if not path or _git_blob_sha1(path) != PARENT_A2_SOURCE_BLOB:
        raise RuntimeError("exact K2-OSC-104 parent source blob drifted")
    contract = bounded_multiharmonic_contract()
    if contract.get("finite_complete_curl_harmonic_family_materialized") is not True:
        raise RuntimeError("K2-OSC-104 complete-curl family prerequisite is unavailable")
    if contract.get("self_contained_velocity_xyzt_provider") is not False:
        raise RuntimeError("K2-OSC-104 source-provider truth boundary drifted")


def _validate_field(field: Any) -> KokunoBoundedMultiHarmonicPhysicalVelocity:
    _parent_identity()
    if not isinstance(field, KokunoBoundedMultiHarmonicPhysicalVelocity):
        raise ValueError("field must be KokunoBoundedMultiHarmonicPhysicalVelocity")
    semantic = str(getattr(field, "semantic_sha256", ""))
    if len(semantic) != 64 or any(c not in "0123456789abcdef" for c in semantic):
        raise ValueError("field semantic identity must be a lowercase 64-hex digest")
    if len(field.terms) < 1:
        raise ValueError("bounded multi-harmonic field must contain at least one term")
    return field


def _probe_coordinates(field: KokunoBoundedMultiHarmonicPhysicalVelocity, *, offgrid: bool):
    max_core = max(float(term.provider.axis_zero_radius_chart) for term in field.terms)
    R = max_core + _R_OFFSETS.copy()
    theta = _THETA.copy()
    Z = _Z.copy()
    T = _T.copy()
    if offgrid:
        R += _OFFGRID_SHIFT["R"]
        theta += _OFFGRID_SHIFT["theta"]
        Z += _OFFGRID_SHIFT["Z"]
        T += _OFFGRID_SHIFT["T"]
    r, z, t = source_chart_to_physical_rzt(
        R,
        Z,
        T,
        ell=field.scaling.ell,
        h=field.scaling.h,
    )
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, z, t, R, theta, Z, T


def _velocity(field: KokunoBoundedMultiHarmonicPhysicalVelocity, x, y, z, t) -> np.ndarray:
    value = np.asarray(field.velocity(x, y, z, t), dtype=float)
    expected = np.broadcast(np.asarray(x), np.asarray(y), np.asarray(z), np.asarray(t)).shape + (3,)
    if value.shape != expected:
        raise RuntimeError(f"public velocity returned shape {value.shape}, expected {expected}")
    if not np.all(np.isfinite(value)):
        raise RuntimeError("public velocity returned nonfinite values")
    return value


def _fd_gradient(
    field: KokunoBoundedMultiHarmonicPhysicalVelocity,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    chart_step: float,
) -> np.ndarray:
    scaling = field.scaling
    h_xy = float(chart_step) * float(scaling.radial_scale)
    h_z = float(chart_step) * float(scaling.axial_scale)
    if not all(math.isfinite(v) and v > 0.0 for v in (h_xy, h_z)):
        raise RuntimeError("physical finite-difference step is not finite and positive")

    for coordinate, step in ((x, h_xy), (y, h_xy), (z, h_z)):
        if np.any(coordinate + step == coordinate) or np.any(coordinate - step == coordinate):
            raise RuntimeError("finite-difference perturbation collapsed in binary64")

    du_dx = (_velocity(field, x + h_xy, y, z, t) - _velocity(field, x - h_xy, y, z, t)) / (2.0 * h_xy)
    du_dy = (_velocity(field, x, y + h_xy, z, t) - _velocity(field, x, y - h_xy, z, t)) / (2.0 * h_xy)
    du_dz = (_velocity(field, x, y, z + h_z, t) - _velocity(field, x, y, z - h_z, t)) / (2.0 * h_z)
    # gradient[..., component, derivative_axis]
    return np.stack((du_dx, du_dy, du_dz), axis=-1)


def _derivative_metrics(gradient: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
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
    grad_rms = float(np.sqrt(np.mean(np.sum(gradient * gradient, axis=(1, 2)))))
    div_rms = float(np.sqrt(np.mean(div * div)))
    vort_sq = np.sum(vort * vort, axis=1)
    vort_rms = float(np.sqrt(np.mean(vort_sq)))
    vort_max = float(np.max(np.sqrt(vort_sq)))
    axial_rms = float(np.sqrt(np.mean(vort[:, 2] * vort[:, 2])))
    return div, vort, {
        "divergence_max_abs": float(np.max(np.abs(div))),
        "divergence_rms": div_rms,
        "gradient_rms": grad_rms,
        "normalized_divergence_rms": div_rms / max(grad_rms, np.finfo(float).tiny),
        "vorticity_rms": vort_rms,
        "vorticity_max": vort_max,
        "axial_vorticity_rms_fraction": axial_rms / max(vort_rms, np.finfo(float).tiny),
    }


def _relative_array_drift(a: np.ndarray, b: np.ndarray) -> float:
    num = float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=-1))))
    den = float(np.sqrt(np.mean(np.sum(b * b, axis=-1))))
    return num / max(den, np.finfo(float).tiny)


def audit_bounded_multiharmonic_field(
    field: KokunoBoundedMultiHarmonicPhysicalVelocity,
) -> dict[str, Any]:
    """Run the frozen three-resolution public-velocity diagnostic.

    This function intentionally accepts no derivative step, gate, residual,
    pressure, forcing, viscosity, optimizer, or tolerance argument.
    """

    field = _validate_field(field)
    x, y, z, t, R, theta, Z, T = _probe_coordinates(field, offgrid=False)
    base_velocity = _velocity(field, x, y, z, t)
    speed = np.linalg.norm(base_velocity, axis=-1)

    resolutions: list[dict[str, Any]] = []
    vorticities: list[np.ndarray] = []
    for chart_step in CHART_FD_STEPS:
        gradient = _fd_gradient(field, x, y, z, t, chart_step=chart_step)
        _, vort, metrics = _derivative_metrics(gradient)
        vorticities.append(vort)
        resolutions.append(
            {
                "chart_step": float(chart_step),
                "physical_xy_step": float(chart_step * field.scaling.radial_scale),
                "physical_z_step": float(chart_step * field.scaling.axial_scale),
                **metrics,
            }
        )

    min_core = min(float(term.provider.axis_zero_radius_chart) for term in field.terms)
    r_core, z_core, t_core = source_chart_to_physical_rzt(
        np.asarray((0.0, 0.5 * min_core)),
        np.asarray((0.0, 0.13)),
        np.asarray((0.91, 1.03)),
        ell=field.scaling.ell,
        h=field.scaling.h,
    )
    core_velocity = _velocity(field, r_core, np.zeros_like(r_core), z_core, t_core)
    axis_core_exact_zero = bool(np.array_equal(core_velocity, np.zeros_like(core_velocity)))
    if not axis_core_exact_zero:
        raise RuntimeError("provider-certified common axis core is not exact zero")

    xo, yo, zo, to, Ro, thetao, Zo, To = _probe_coordinates(field, offgrid=True)
    offgrid_velocity = _velocity(field, xo, yo, zo, to)
    offgrid_speed = np.linalg.norm(offgrid_velocity, axis=-1)

    medium_fine_drift = _relative_array_drift(vorticities[1], vorticities[2])
    coarse_fine_drift = _relative_array_drift(vorticities[0], vorticities[2])

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "parent": {
            "pr": PARENT_A2_PR,
            "head": PARENT_A2_HEAD,
            "source_blob_sha1": PARENT_A2_SOURCE_BLOB,
        },
        "field_semantic_sha256": field.semantic_sha256,
        "protocol": {
            "chart_fd_steps": list(CHART_FD_STEPS),
            "probe_count": int(R.size),
            "radius_policy": "largest certified axis-zero core plus frozen offsets",
            "offgrid_shift": dict(_OFFGRID_SHIFT),
            "derivatives_consume_public_velocity_only": True,
            "scientific_acceptance_gate": None,
        },
        "probe_chart": {
            "R": R.tolist(),
            "theta": theta.tolist(),
            "Z": Z.tolist(),
            "T": T.tolist(),
        },
        "velocity_observation": {
            "speed_rms": float(np.sqrt(np.mean(speed * speed))),
            "speed_max": float(np.max(speed)),
            "nonzero_probe_count": int(np.count_nonzero(speed > 0.0)),
        },
        "resolutions": resolutions,
        "resolution_stability_observation": {
            "medium_to_fine_vorticity_relative_rms_drift": medium_fine_drift,
            "coarse_to_fine_vorticity_relative_rms_drift": coarse_fine_drift,
        },
        "axis_support_observation": {
            "common_axis_core_exact_zero": axis_core_exact_zero,
            "common_core_probe_R": [0.0, 0.5 * min_core],
            "outer_support_assessed": False,
            "reason_outer_support_not_assessed": "provider contract exposes axis-zero core but no authenticated outer-support boundary",
        },
        "offgrid_observation": {
            "finite": bool(np.all(np.isfinite(offgrid_velocity))),
            "speed_rms": float(np.sqrt(np.mean(offgrid_speed * offgrid_speed))),
            "speed_max": float(np.max(offgrid_speed)),
            "chart_R": Ro.tolist(),
            "chart_theta": thetao.tolist(),
            "chart_Z": Zo.tolist(),
            "chart_T": To.tolist(),
        },
        "truth_boundary": bounded_multiharmonic_diagnostics_contract(),
    }
    report["diagnostic_sha256"] = _sha256(report)
    return report


def bounded_multiharmonic_diagnostics_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_a2_pr": PARENT_A2_PR,
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_a2_source_blob": PARENT_A2_SOURCE_BLOB,
        "public_velocity_only_cartesian_derivatives_materialized": True,
        "three_resolution_divergence_observation_materialized": True,
        "three_resolution_vorticity_morphology_observation_materialized": True,
        "offgrid_velocity_observation_materialized": True,
        "certified_axis_zero_support_checked": True,
        "outer_support_boundary_assessed": False,
        "analytic_divergence_identity_proved_by_this_increment": False,
        "source_provider_self_contained": False,
        "complete_ns_residual_assessed": False,
        "oscillation_before_after_ns_residual_compared": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def public_contract() -> dict[str, Any]:
    parameters = set(inspect.signature(audit_bounded_multiharmonic_field).parameters)
    forbidden = {
        "step", "steps", "threshold", "tolerance", "rtol", "atol", "residual",
        "forcing", "pressure", "viscosity", "gain", "target", "optimizer", "panels",
    }
    return {
        "audit_inputs": list(inspect.signature(audit_bounded_multiharmonic_field).parameters),
        "forbidden_audit_inputs_present": sorted(parameters & forbidden),
        "chart_fd_steps": list(CHART_FD_STEPS),
        "paper_exact": False,
        "pde_validated": False,
    }


__all__ = [
    "CHART_FD_STEPS",
    "TASK",
    "audit_bounded_multiharmonic_field",
    "bounded_multiharmonic_diagnostics_contract",
    "public_contract",
]
