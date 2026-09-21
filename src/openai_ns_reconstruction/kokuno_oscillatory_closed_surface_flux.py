"""Closed-surface incompressibility flux diagnostic for the frozen A2 oscillation.

This module changes no velocity byte and adds no oscillatory degree of freedom.  It
samples only the already-frozen public ``velocity_osc(x,y,z,t)`` on several
strictly interior cylindrical control surfaces and evaluates

    integral_{partial V} u_osc . n dS.

For a divergence-free field this closed-surface flux is zero.  The control
surfaces are deliberately kept away from the exact compact-support boundary, so
individual face fluxes are nontrivial and cancellation cannot pass merely because
the public support mask returns zero.

Kokuno's corrected 2026-09-09 reconstruction remains structural provenance for
the localized complete-curl organization.  The public-z pullback, concrete C7
support envelope, autonomous mode/background choices, quadrature protocol, and
this diagnostic are repository realizations and are not paper-exact hidden data.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_public_z_pullback_velocity import default_field, velocity_osc

TASK = "KOKUNO-A2-OSCILLATORY-CLOSED-SURFACE-FLUX-076"
SCHEMA = "kokuno-a2-oscillatory-closed-surface-flux-v1"
PARENT_AGENT2_PR = 934
PARENT_AGENT2_HEAD = "c542ae397f3aadf52d7f63fe32d0d19fac866273"
PARENT_SUPPORT_FLUX_BLOB = "ec4e6f30da60fa9a5645d4e938f5ece1830dc4d4"
PUBLIC_Z_VELOCITY_BLOB = "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

# Control volumes are expressed as fractions of the existing public support.
# Every face is strictly interior, unlike the exact-zero support boundary audited
# by parent #934.  The second volume is axially asymmetric to avoid a symmetry-
# only cancellation certificate.
CONTROL_VOLUME_FRACTIONS = (
    {"r_inner": 0.15, "r_outer": 0.85, "z_lower": 0.18, "z_upper": 0.82},
    {"r_inner": 0.25, "r_outer": 0.75, "z_lower": 0.12, "z_upper": 0.70},
)
EVALUATION_TIMES = (0.41, 0.59)
# (theta nodes, radial nodes on caps, axial nodes on cylindrical faces)
RESOLUTION_LEVELS = ((16, 9, 11), (32, 17, 21), (64, 33, 41))
FINE_THETA_SHIFT_CELLS = 0.371

# Frozen before exact-head Actions output.  These are scoped quadrature/flux
# consistency guards, not CR001 PDE gates.
COARSE_RELATIVE_CLOSURE_MAX = 5.0e-2
MEDIUM_RELATIVE_CLOSURE_MAX = 2.0e-2
FINE_RELATIVE_CLOSURE_MAX = 1.0e-2
SHIFTED_FINE_RELATIVE_CLOSURE_MAX = 1.5e-2
FACE_FLUX_DENOMINATOR_FLOOR = 1.0e-8


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _as_velocity(value: Any) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim < 1 or out.shape[-1] != 3 or not np.all(np.isfinite(out)):
        raise RuntimeError("public oscillatory velocity returned an invalid Cartesian array")
    return out


def _trapz_1d(values: np.ndarray, coordinate: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    coordinate = np.asarray(coordinate, dtype=float)
    if values.ndim != 1 or coordinate.ndim != 1 or values.shape != coordinate.shape:
        raise ValueError("trapezoid input must be matching one-dimensional arrays")
    if coordinate.size < 2 or not np.all(np.diff(coordinate) > 0.0):
        raise ValueError("trapezoid coordinate must be strictly increasing")
    delta = np.diff(coordinate)
    return float(np.sum(0.5 * delta * (values[:-1] + values[1:])))


def _support() -> tuple[float, float, float, float]:
    field = default_field()
    return (
        float(field.radial_inner),
        float(field.radial_outer),
        float(field.axial_lower),
        float(field.axial_upper),
    )


def _physical_control_volumes() -> tuple[dict[str, float], ...]:
    r0, r1, z0, z1 = _support()
    dr = r1 - r0
    dz = z1 - z0
    volumes: list[dict[str, float]] = []
    for frac in CONTROL_VOLUME_FRACTIONS:
        volume = {
            "r_inner": r0 + float(frac["r_inner"]) * dr,
            "r_outer": r0 + float(frac["r_outer"]) * dr,
            "z_lower": z0 + float(frac["z_lower"]) * dz,
            "z_upper": z0 + float(frac["z_upper"]) * dz,
        }
        if not (r0 < volume["r_inner"] < volume["r_outer"] < r1):
            raise RuntimeError("registered radial control volume left public support")
        if not (z0 < volume["z_lower"] < volume["z_upper"] < z1):
            raise RuntimeError("registered axial control volume left public support")
        volumes.append(volume)
    return tuple(volumes)


def _closed_cylinder_flux(
    volume: Mapping[str, float],
    time: float,
    resolution: tuple[int, int, int],
    *,
    theta_shift_cells: float = 0.0,
) -> dict[str, float]:
    ntheta, nr, nz = (int(v) for v in resolution)
    if ntheta < 4 or nr < 3 or nz < 3:
        raise ValueError("closed-surface quadrature resolution is too small")
    r_inner = float(volume["r_inner"])
    r_outer = float(volume["r_outer"])
    z_lower = float(volume["z_lower"])
    z_upper = float(volume["z_upper"])

    dtheta = 2.0 * math.pi / ntheta
    theta = (np.arange(ntheta, dtype=float) + float(theta_shift_cells)) * dtheta
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    # Radial cylindrical faces: dS = r dtheta dz.  Periodic theta quadrature
    # uses the uniform trapezoidal rule, which is simply dtheta * sum.
    z_line = np.linspace(z_lower, z_upper, nz, dtype=float)
    zz, th = np.meshgrid(z_line, theta, indexing="ij")
    radial_face_flux: dict[str, float] = {}
    for label, radius, sign in (
        ("radial_inner", r_inner, -1.0),
        ("radial_outer", r_outer, +1.0),
    ):
        xx = radius * np.cos(th)
        yy = radius * np.sin(th)
        tt = np.full_like(xx, float(time))
        uu = _as_velocity(velocity_osc(xx, yy, zz, tt))
        u_r = uu[..., 0] * np.cos(th) + uu[..., 1] * np.sin(th)
        theta_integral = dtheta * np.sum(sign * radius * u_r, axis=1)
        radial_face_flux[label] = _trapz_1d(theta_integral, z_line)

    # Axial caps: dS = r dr dtheta.
    r_line = np.linspace(r_inner, r_outer, nr, dtype=float)
    rr, th_cap = np.meshgrid(r_line, theta, indexing="ij")
    cap_flux: dict[str, float] = {}
    for label, axial, sign in (
        ("axial_lower", z_lower, -1.0),
        ("axial_upper", z_upper, +1.0),
    ):
        xx = rr * np.cos(th_cap)
        yy = rr * np.sin(th_cap)
        zz_cap = np.full_like(xx, axial)
        tt = np.full_like(xx, float(time))
        uu = _as_velocity(velocity_osc(xx, yy, zz_cap, tt))
        theta_integral = dtheta * np.sum(sign * uu[..., 2] * rr, axis=1)
        cap_flux[label] = _trapz_1d(theta_integral, r_line)

    faces = {**radial_face_flux, **cap_flux}
    values = np.asarray(tuple(faces.values()), dtype=float)
    if not np.all(np.isfinite(values)):
        raise RuntimeError("closed-surface face flux became non-finite")
    net = float(np.sum(values))
    denominator = float(np.sum(np.abs(values)))
    relative = float(abs(net) / max(denominator, np.finfo(float).tiny))
    return {
        **faces,
        "net_flux": net,
        "absolute_face_flux_sum": denominator,
        "relative_closure": relative,
        "max_absolute_face_flux": float(np.max(np.abs(values))),
    }


@dataclass(frozen=True)
class OscillatoryClosedSurfaceFluxResult:
    control_volumes: tuple[dict[str, float], ...]
    times: tuple[float, ...]
    levels: tuple[dict[str, Any], ...]
    shifted_fine: dict[str, Any]

    def summary_payload(self) -> dict[str, Any]:
        return {
            "control_volumes": [dict(volume) for volume in self.control_volumes],
            "times": list(self.times),
            "levels": [dict(level) for level in self.levels],
            "shifted_fine": dict(self.shifted_fine),
        }


class KokunoOscillatoryClosedSurfaceFlux:
    """Materialize the frozen interior closed-surface flux receipt."""

    def semantic_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "parent_support_flux_blob": PARENT_SUPPORT_FLUX_BLOB,
            "public_z_velocity_blob": PUBLIC_Z_VELOCITY_BLOB,
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
                "scope": "localized oscillatory complete-curl organization",
            },
            "frozen_control_volume_fractions": [dict(v) for v in CONTROL_VOLUME_FRACTIONS],
            "frozen_times": list(EVALUATION_TIMES),
            "frozen_resolution_levels": [list(level) for level in RESOLUTION_LEVELS],
            "fine_theta_shift_cells": FINE_THETA_SHIFT_CELLS,
            "quadrature_contract": {
                "geometry": "strictly interior closed cylindrical control volumes",
                "radial_faces": "periodic theta trapezoid x ordinary z trapezoid",
                "axial_caps": "periodic theta trapezoid x ordinary radial trapezoid with r Jacobian",
                "quantity": "closed surface integral of u_osc dot n",
                "individual_face_fluxes_required_nontrivial": True,
                "uses_only_public_velocity_values": True,
                "production_jacobian_used": False,
                "production_divergence_used": False,
            },
            "frozen_engineering_gates": {
                "coarse_relative_closure_max": COARSE_RELATIVE_CLOSURE_MAX,
                "medium_relative_closure_max": MEDIUM_RELATIVE_CLOSURE_MAX,
                "fine_relative_closure_max": FINE_RELATIVE_CLOSURE_MAX,
                "shifted_fine_relative_closure_max": SHIFTED_FINE_RELATIVE_CLOSURE_MAX,
                "face_flux_denominator_floor": FACE_FLUX_DENOMINATOR_FLOOR,
            },
            "scientific_boundary": {
                "oscillatory_velocity_changed": False,
                "oscillatory_support_changed": False,
                "new_oscillatory_parameters_added": False,
                "generic_interior_closed_surface_flux_assessed": True,
                "local_derivative_reconstruction_performed": False,
                "volume_divergence_integral_assessed": False,
                "whole_domain_divergence_volume_l2_assessed": False,
                "mean_projection_performed": False,
                "radial_inverse_performed": False,
                "correction_velocity_constructed": False,
                "pressure_included": False,
                "forcing_included": False,
                "complete_ns_residual": False,
                "same_protocol_st006_comparison_valid": False,
                "residual_reduction_claimed": False,
                "paper_exact": False,
                "pde_validated": False,
            },
        }

    @property
    def closed_surface_flux_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    def materialize(self) -> OscillatoryClosedSurfaceFluxResult:
        volumes = _physical_control_volumes()
        levels: list[dict[str, Any]] = []
        for resolution in RESOLUTION_LEVELS:
            cases: list[dict[str, Any]] = []
            for volume_index, volume in enumerate(volumes):
                for time in EVALUATION_TIMES:
                    flux = _closed_cylinder_flux(volume, time, resolution)
                    cases.append({
                        "volume_index": volume_index,
                        "time": float(time),
                        **flux,
                    })
            levels.append({
                "resolution": list(resolution),
                "cases": cases,
                "max_relative_closure": float(max(case["relative_closure"] for case in cases)),
                "min_absolute_face_flux_sum": float(min(case["absolute_face_flux_sum"] for case in cases)),
            })

        fine = RESOLUTION_LEVELS[-1]
        shifted_cases: list[dict[str, Any]] = []
        for volume_index, volume in enumerate(volumes):
            for time in EVALUATION_TIMES:
                flux = _closed_cylinder_flux(
                    volume,
                    time,
                    fine,
                    theta_shift_cells=FINE_THETA_SHIFT_CELLS,
                )
                shifted_cases.append({
                    "volume_index": volume_index,
                    "time": float(time),
                    **flux,
                })
        shifted = {
            "resolution": list(fine),
            "theta_shift_cells": FINE_THETA_SHIFT_CELLS,
            "cases": shifted_cases,
            "max_relative_closure": float(max(case["relative_closure"] for case in shifted_cases)),
            "min_absolute_face_flux_sum": float(min(case["absolute_face_flux_sum"] for case in shifted_cases)),
        }
        return OscillatoryClosedSurfaceFluxResult(
            control_volumes=volumes,
            times=tuple(float(t) for t in EVALUATION_TIMES),
            levels=tuple(levels),
            shifted_fine=shifted,
        )

    def manifest(self) -> dict[str, Any]:
        payload = self.semantic_payload()
        return {"payload": payload, "closed_surface_flux_sha256": _sha256(payload)}

    def save_manifest(self, path: str | Path) -> dict[str, Any]:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        manifest = self.manifest()
        target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return manifest

    @classmethod
    def load_manifest(cls, path: str | Path) -> "KokunoOscillatoryClosedSurfaceFlux":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"payload", "closed_surface_flux_sha256"}:
            raise ValueError("invalid oscillatory closed-surface-flux manifest envelope")
        payload = raw["payload"]
        if not isinstance(payload, dict) or _sha256(payload) != raw["closed_surface_flux_sha256"]:
            raise ValueError("oscillatory closed-surface-flux manifest checksum mismatch")
        obj = cls()
        if obj.semantic_payload() != payload:
            raise ValueError("oscillatory closed-surface-flux manifest semantics mismatch")
        return obj


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(KokunoOscillatoryClosedSurfaceFlux.materialize)
    forbidden = {
        "residual", "defect", "mean", "m0", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "damping", "nu", "viscosity", "spatial_step", "time_step",
        "amplitude", "phase", "scale", "orientation", "support", "scientific_threshold",
        "threshold", "resolution", "quadrature", "control_volume", "theta_shift",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "forbidden_materialize_inputs_present": sorted(forbidden.intersection(signature.parameters)),
        "generic_interior_closed_surface_flux_assessed": True,
        "individual_face_fluxes_required_nontrivial": True,
        "uses_only_public_velocity_values": True,
        "volume_divergence_integral_assessed": False,
        "whole_domain_divergence_volume_l2_assessed": False,
        "mean_projection_performed": False,
        "correction_velocity_constructed": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }
