"""Compact-support boundary flux contract for the frozen Kokuno oscillation.

This module adds no velocity degree of freedom.  It consumes the already-frozen
public complete-curl oscillatory field and makes one previously implicit support
property first class: the public field is exactly zero on and outside the
registered radial/axial support boundary.  Hence its kinetic-energy flux

    F_K = (|u_osc|^2 / 2) u_osc

has zero normal trace there.

The zero trace is a property of the repository realization: a C7 compact
``cos^8`` envelope is applied to the vector-potential coefficient before the
complete curl, and the public evaluator fail-closes to zero outside the strict
support mask.  Kokuno's corrected 2026-09-09 reconstruction remains structural
provenance for the localized complete-curl organization; this support envelope,
public-z pullback, autonomous mode/background, and the present diagnostic are
not paper-exact hidden data.
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

TASK = "KOKUNO-A2-OSCILLATORY-SUPPORT-BOUNDARY-FLUX-075"
SCHEMA = "kokuno-a2-oscillatory-support-boundary-flux-v1"
PARENT_AGENT2_PR = 926
PARENT_AGENT2_HEAD = "22d89643356d99048b1039b6094040e9dbc589c3"
PARENT_ENERGY_FLUX_BLOB = "cce26ef1593856a029d79a3b06901532dcb76c5d"
PUBLIC_Z_VELOCITY_BLOB = "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

BOUNDARY_THETA_COUNT = 24
BOUNDARY_Z_COUNT = 7
BOUNDARY_RADIUS_COUNT = 7
BOUNDARY_TIMES = (0.31, 0.50, 0.69)
OUTSIDE_RELATIVE_OFFSET = 1.0e-6
INWARD_COLLAR_FRACTIONS = (1.0 / 16.0, 1.0 / 32.0, 1.0 / 64.0)


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _as_velocity(value: Any) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim < 1 or out.shape[-1] != 3 or not np.all(np.isfinite(out)):
        raise RuntimeError("public oscillatory velocity returned an invalid Cartesian array")
    return out


def _surface_flux(velocity: np.ndarray, normal: np.ndarray) -> np.ndarray:
    speed_sq = np.sum(velocity * velocity, axis=-1)
    return 0.5 * speed_sq * np.sum(velocity * normal, axis=-1)


@dataclass(frozen=True)
class OscillatorySupportBoundaryFluxResult:
    radial_support: tuple[float, float]
    axial_support: tuple[float, float]
    boundary_sample_count: int
    outside_sample_count: int
    boundary_velocity_max_abs: float
    outside_velocity_max_abs: float
    boundary_normal_energy_flux_max_abs: float
    pointwise_boundary_trace_zero: bool
    pointwise_outside_zero: bool
    support_boundary_energy_flux_integral: float

    def summary_payload(self) -> dict[str, Any]:
        return {
            "radial_support": list(self.radial_support),
            "axial_support": list(self.axial_support),
            "boundary_sample_count": self.boundary_sample_count,
            "outside_sample_count": self.outside_sample_count,
            "boundary_velocity_max_abs": self.boundary_velocity_max_abs,
            "outside_velocity_max_abs": self.outside_velocity_max_abs,
            "boundary_normal_energy_flux_max_abs": self.boundary_normal_energy_flux_max_abs,
            "pointwise_boundary_trace_zero": self.pointwise_boundary_trace_zero,
            "pointwise_outside_zero": self.pointwise_outside_zero,
            "support_boundary_energy_flux_integral": self.support_boundary_energy_flux_integral,
        }


class KokunoOscillatorySupportBoundaryFlux:
    """Materialize the frozen support-boundary zero-flux contract."""

    def _support(self) -> tuple[float, float, float, float]:
        field = default_field()
        return (
            float(field.radial_inner),
            float(field.radial_outer),
            float(field.axial_lower),
            float(field.axial_upper),
        )

    def semantic_payload(self) -> dict[str, Any]:
        r0, r1, z0, z1 = self._support()
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "parent_energy_flux_blob": PARENT_ENERGY_FLUX_BLOB,
            "public_z_velocity_blob": PUBLIC_Z_VELOCITY_BLOB,
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
                "scope": "localized oscillatory complete-curl organization",
            },
            "repository_support_realization": {
                "radial_support": [r0, r1],
                "axial_support": [z0, z1],
                "vector_potential_envelope": "C7 compact cos^8 radial and axial envelopes",
                "strict_public_support_mask": True,
                "support_kept_away_from_axis": r0 > 0.0,
            },
            "frozen_diagnostic_protocol": {
                "boundary_theta_count": BOUNDARY_THETA_COUNT,
                "boundary_z_count": BOUNDARY_Z_COUNT,
                "boundary_radius_count": BOUNDARY_RADIUS_COUNT,
                "boundary_times": list(BOUNDARY_TIMES),
                "outside_relative_offset": OUTSIDE_RELATIVE_OFFSET,
                "independent_inward_collar_fractions": list(INWARD_COLLAR_FRACTIONS),
            },
            "quantity_contract": {
                "kinetic_energy_flux": "(|u_osc|^2/2) u_osc",
                "boundary_normal_flux": "(|u_osc|^2/2) (u_osc dot n)",
                "support_boundary_trace": "u_osc=0 on radial/axial support boundary",
                "outside_trace": "u_osc=0 outside radial/axial support",
                "boundary_flux_consequence": "pointwise zero normal kinetic-energy flux implies zero support-boundary flux integral",
            },
            "scientific_boundary": {
                "oscillatory_velocity_changed": False,
                "oscillatory_support_changed": False,
                "new_oscillatory_parameters_added": False,
                "support_boundary_trace_assessed": True,
                "support_boundary_energy_flux_assessed": True,
                "whole_domain_energy_balance_assessed": False,
                "volume_integral_of_self_advection_power_assessed": False,
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
    def support_flux_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    def materialize(self) -> OscillatorySupportBoundaryFluxResult:
        r0, r1, z0, z1 = self._support()
        theta = np.linspace(-math.pi, math.pi, BOUNDARY_THETA_COUNT, endpoint=False)
        z_line = np.linspace(z0 + 0.2 * (z1 - z0), z1 - 0.2 * (z1 - z0), BOUNDARY_Z_COUNT)
        r_line = np.linspace(r0 + 0.2 * (r1 - r0), r1 - 0.2 * (r1 - r0), BOUNDARY_RADIUS_COUNT)
        times = np.asarray(BOUNDARY_TIMES, dtype=float)

        boundary_velocity: list[np.ndarray] = []
        boundary_normal: list[np.ndarray] = []
        outside_velocity: list[np.ndarray] = []

        # Radial inner/outer cylinders.
        tt, zz, th = np.meshgrid(times, z_line, theta, indexing="ij")
        for radius, sign in ((r0, -1.0), (r1, 1.0)):
            xx = radius * np.cos(th)
            yy = radius * np.sin(th)
            uu = _as_velocity(velocity_osc(xx, yy, zz, tt))
            nn = np.stack((sign * np.cos(th), sign * np.sin(th), np.zeros_like(th)), axis=-1)
            boundary_velocity.append(uu.reshape((-1, 3)))
            boundary_normal.append(nn.reshape((-1, 3)))

            offset = OUTSIDE_RELATIVE_OFFSET * (r1 - r0)
            outside_radius = radius + sign * offset
            uo = _as_velocity(velocity_osc(outside_radius * np.cos(th), outside_radius * np.sin(th), zz, tt))
            outside_velocity.append(uo.reshape((-1, 3)))

        # Axial lower/upper caps, sampled away from the radial edges.
        tt2, rr, th2 = np.meshgrid(times, r_line, theta, indexing="ij")
        for axial, sign in ((z0, -1.0), (z1, 1.0)):
            xx = rr * np.cos(th2)
            yy = rr * np.sin(th2)
            zz_cap = np.full_like(xx, axial)
            uu = _as_velocity(velocity_osc(xx, yy, zz_cap, tt2))
            nn = np.zeros_like(uu)
            nn[..., 2] = sign
            boundary_velocity.append(uu.reshape((-1, 3)))
            boundary_normal.append(nn.reshape((-1, 3)))

            offset = OUTSIDE_RELATIVE_OFFSET * (z1 - z0)
            uo = _as_velocity(velocity_osc(xx, yy, np.full_like(xx, axial + sign * offset), tt2))
            outside_velocity.append(uo.reshape((-1, 3)))

        ub = np.concatenate(boundary_velocity, axis=0)
        nb = np.concatenate(boundary_normal, axis=0)
        uo = np.concatenate(outside_velocity, axis=0)
        flux = _surface_flux(ub, nb)
        boundary_zero = bool(np.array_equal(ub, np.zeros_like(ub)))
        outside_zero = bool(np.array_equal(uo, np.zeros_like(uo)))
        flux_zero = bool(np.array_equal(flux, np.zeros_like(flux)))
        if not (boundary_zero and outside_zero and flux_zero):
            raise RuntimeError("oscillatory public support boundary lost exact zero trace")

        return OscillatorySupportBoundaryFluxResult(
            radial_support=(r0, r1),
            axial_support=(z0, z1),
            boundary_sample_count=int(ub.shape[0]),
            outside_sample_count=int(uo.shape[0]),
            boundary_velocity_max_abs=float(np.max(np.abs(ub))),
            outside_velocity_max_abs=float(np.max(np.abs(uo))),
            boundary_normal_energy_flux_max_abs=float(np.max(np.abs(flux))),
            pointwise_boundary_trace_zero=True,
            pointwise_outside_zero=True,
            support_boundary_energy_flux_integral=0.0,
        )

    def manifest(self) -> dict[str, Any]:
        payload = self.semantic_payload()
        return {"payload": payload, "support_flux_sha256": _sha256(payload)}

    def save_manifest(self, path: str | Path) -> dict[str, Any]:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        manifest = self.manifest()
        target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return manifest

    @classmethod
    def load_manifest(cls, path: str | Path) -> "KokunoOscillatorySupportBoundaryFlux":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"payload", "support_flux_sha256"}:
            raise ValueError("invalid oscillatory support-flux manifest envelope")
        payload = raw["payload"]
        if not isinstance(payload, dict) or _sha256(payload) != raw["support_flux_sha256"]:
            raise ValueError("oscillatory support-flux manifest checksum mismatch")
        obj = cls()
        if obj.semantic_payload() != payload:
            raise ValueError("oscillatory support-flux manifest semantics mismatch")
        return obj


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(KokunoOscillatorySupportBoundaryFlux.materialize)
    forbidden = {
        "residual", "defect", "mean", "m0", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "damping", "nu", "viscosity", "spatial_step", "time_step",
        "amplitude", "phase", "scale", "orientation", "support", "scientific_threshold",
        "threshold", "collar_fraction", "resolution",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "forbidden_materialize_inputs_present": sorted(forbidden.intersection(signature.parameters)),
        "independent_public_velocity_collar_reaudit_required": True,
        "support_boundary_trace_assessed": True,
        "support_boundary_energy_flux_assessed": True,
        "whole_domain_energy_balance_assessed": False,
        "mean_projection_performed": False,
        "correction_velocity_constructed": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }
