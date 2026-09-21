"""Local kinetic-energy transport contract for the frozen Kokuno oscillation.

This module adds no velocity degrees of freedom.  It consumes the already-frozen
public complete-curl oscillatory field and its existing FD6 quadratic
self-advection diagnostic,

    Q = (u_osc . grad) u_osc,

and exposes the local kinetic-energy identity

    u_osc . Q
      = div((|u_osc|^2 / 2) u_osc)
        - (|u_osc|^2 / 2) div(u_osc).

For an exactly divergence-free complete curl the second term vanishes, so the
quadratic self-interaction is locally conservative kinetic-energy transport.
The production object below does not independently differentiate the flux; it
records the algebraic power/remainder from the frozen FD6 seam.  Independent
public-velocity-only FD4 reconstruction belongs to the exact-head workflow.

Kokuno's corrected 2026-09-09 reconstruction remains structural provenance for
the upstream localized complete-curl organization.  The public-z pullback,
autonomous mode/background choices, finite differences, and this diagnostic are
repository realizations and are not paper-exact hidden data.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_public_oscillatory_self_advection import (
    GRADIENT_FD6_STEP,
    evaluate_oscillatory_self_advection_fd6,
)
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-OSCILLATORY-SELF-ADVECTION-ENERGY-FLUX-074"
SCHEMA = "kokuno-a2-oscillatory-self-advection-energy-flux-v1"
PARENT_AGENT2_PR = 919
PARENT_AGENT2_HEAD = "2162b1bfa9d00d6c142848e52f36f5b699722a0c"
PARENT_CANCELLATION_BLOB = "10e3f29ee1a159c8647dda2ebecefb5e9c6952f0"
PUBLIC_Z_VELOCITY_BLOB = "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _as_vector(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim < 1 or out.shape[-1] != 3:
        raise ValueError(f"{name} must have trailing Cartesian axis of length 3")
    if out.size == 0 or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain at least one finite sample")
    return out


@dataclass(frozen=True)
class OscillatorySelfAdvectionEnergyFluxResult:
    velocity: np.ndarray
    self_advection: np.ndarray
    kinetic_energy_density: np.ndarray
    self_advection_power_density: np.ndarray
    divergence: np.ndarray
    compressibility_remainder: np.ndarray
    conservative_flux_divergence_equivalent: np.ndarray
    inherited_spatial_step: float

    def summary_payload(self) -> dict[str, float]:
        power = np.asarray(self.self_advection_power_density, dtype=float)
        remainder = np.asarray(self.compressibility_remainder, dtype=float)
        flux = np.asarray(self.conservative_flux_divergence_equivalent, dtype=float)
        energy = np.asarray(self.kinetic_energy_density, dtype=float)

        def rms(value: np.ndarray) -> float:
            return float(np.sqrt(np.mean(np.asarray(value, dtype=float) ** 2)))

        return {
            "kinetic_energy_density_rms": rms(energy),
            "self_advection_power_rms": rms(power),
            "compressibility_remainder_rms": rms(remainder),
            "conservative_flux_divergence_equivalent_rms": rms(flux),
            "compressibility_to_power_rms_ratio": rms(remainder)
            / max(rms(power), np.finfo(float).tiny),
        }


class KokunoOscillatorySelfAdvectionEnergyFlux:
    """Expose the conservative-form energy contract of the frozen A2 oscillation."""

    def semantic_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "parent_cancellation_blob": PARENT_CANCELLATION_BLOB,
            "public_z_velocity_blob": PUBLIC_Z_VELOCITY_BLOB,
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
                "scope": "localized oscillatory complete-curl organization",
            },
            "quantity_contract": {
                "velocity": "u_osc",
                "quadratic_self_advection": "(u_osc.grad)u_osc",
                "kinetic_energy_density": "|u_osc|^2/2",
                "self_advection_power_density": "u_osc dot ((u_osc.grad)u_osc)",
                "compressibility_remainder": "(|u_osc|^2/2) div(u_osc)",
                "conservative_flux_divergence_equivalent": (
                    "u_osc dot ((u_osc.grad)u_osc) + "
                    "(|u_osc|^2/2) div(u_osc)"
                ),
                "identity": (
                    "div((|u_osc|^2/2)u_osc) = "
                    "u_osc dot ((u_osc.grad)u_osc) + "
                    "(|u_osc|^2/2) div(u_osc)"
                ),
                "production_derivative": "inherits centered Cartesian FD6",
                "production_spatial_step": GRADIENT_FD6_STEP,
                "production_does_not_independently_differentiate_energy_flux": True,
            },
            "scientific_boundary": {
                "oscillatory_velocity_changed": False,
                "new_oscillatory_parameters_added": False,
                "strictly_local_identity_only": True,
                "whole_domain_energy_balance_assessed": False,
                "boundary_flux_integral_assessed": False,
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
    def energy_flux_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> OscillatorySelfAdvectionEnergyFluxResult:
        inherited = evaluate_oscillatory_self_advection_fd6(
            x, y, z, t, spatial_step=GRADIENT_FD6_STEP
        )
        u = _as_vector(inherited["velocity"], "oscillatory velocity")
        q = _as_vector(inherited["self_advection"], "oscillatory self-advection")
        if u.shape != q.shape:
            raise RuntimeError("oscillatory energy-flux velocity/self-advection shape mismatch")

        direct_u = _as_vector(velocity_osc(x, y, z, t), "public oscillatory velocity")
        if direct_u.shape != u.shape or not np.array_equal(direct_u, u):
            raise RuntimeError("oscillatory velocity did not replay exactly across inherited seams")

        divergence = np.asarray(inherited["divergence"], dtype=float)
        expected_scalar_shape = u.shape[:-1]
        if divergence.shape != expected_scalar_shape or not np.all(np.isfinite(divergence)):
            raise RuntimeError("oscillatory divergence has invalid shape or non-finite values")

        energy = 0.5 * np.sum(u * u, axis=-1)
        power = np.sum(u * q, axis=-1)
        remainder = energy * divergence
        conservative_equivalent = power + remainder
        scalars = (energy, power, remainder, conservative_equivalent)
        if not all(value.shape == expected_scalar_shape for value in scalars):
            raise RuntimeError("oscillatory energy-flux scalar shape mismatch")
        if not all(np.all(np.isfinite(value)) for value in scalars):
            raise RuntimeError("oscillatory energy-flux diagnostic became non-finite")

        return OscillatorySelfAdvectionEnergyFluxResult(
            velocity=u,
            self_advection=q,
            kinetic_energy_density=energy,
            self_advection_power_density=power,
            divergence=divergence,
            compressibility_remainder=remainder,
            conservative_flux_divergence_equivalent=conservative_equivalent,
            inherited_spatial_step=float(inherited["spatial_step"]),
        )

    def manifest(self) -> dict[str, Any]:
        payload = self.semantic_payload()
        return {"payload": payload, "energy_flux_sha256": _sha256(payload)}

    def save_manifest(self, path: str | Path) -> dict[str, Any]:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        manifest = self.manifest()
        target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return manifest

    @classmethod
    def load_manifest(cls, path: str | Path) -> "KokunoOscillatorySelfAdvectionEnergyFlux":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"payload", "energy_flux_sha256"}:
            raise ValueError("invalid oscillatory energy-flux manifest envelope")
        payload = raw["payload"]
        if not isinstance(payload, dict) or _sha256(payload) != raw["energy_flux_sha256"]:
            raise ValueError("oscillatory energy-flux manifest checksum mismatch")
        obj = cls()
        if obj.semantic_payload() != payload:
            raise ValueError("oscillatory energy-flux manifest semantics mismatch")
        return obj


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(KokunoOscillatorySelfAdvectionEnergyFlux.evaluate)
    forbidden = {
        "residual", "defect", "mean", "m0", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "damping", "nu", "viscosity", "spatial_step", "time_step",
        "amplitude", "phase", "scale", "orientation", "scientific_threshold", "threshold",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "forbidden_evaluate_inputs_present": sorted(forbidden.intersection(signature.parameters)),
        "public_velocity_only_independent_reaudit_required": True,
        "production_energy_flux_divergence_independently_differentiated": False,
        "mean_projection_performed": False,
        "correction_velocity_constructed": False,
        "whole_domain_energy_balance_assessed": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }
