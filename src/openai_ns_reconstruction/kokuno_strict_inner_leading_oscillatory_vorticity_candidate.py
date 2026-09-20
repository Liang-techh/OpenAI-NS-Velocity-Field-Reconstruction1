"""First-class strict-inner leading + oscillatory candidate with vorticity.

This Agent-2 increment promotes the vorticity already implied by the checksum-bound
strict-inner spatial candidate to a public callable

    vorticity(x,y,z,t) = curl(u_inner_strict + u_osc_complete_curl).

No new derivative realization is introduced: production vorticity is computed from
the existing #874 spatial Jacobian.  The corrected 2026-09-09 Kokuno reader remains
structural provenance for upstream PA.10 / localized complete-curl organization.
The public-z pullback, autonomous oscillatory parameters, finite-difference
oscillatory Jacobian, artifact identity scheme, and this vorticity surface are
repository realizations and are not claimed paper-exact.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_strict_inner_leading_oscillatory_spatial_candidate import (
    KokunoStrictInnerLeadingOscillatorySpatialCandidate,
    StrictInnerSpatialBackend,
)

TASK = "KOKUNO-A2-STRICT-INNER-VORTICITY-CANDIDATE-067"
SCHEMA = "kokuno-a2-strict-inner-vorticity-candidate-v1"
PARENT_AGENT2_PR = 874
PARENT_AGENT2_HEAD = "9cb3869b9cfd2d5b8dcb1da222df74d80e12d0c2"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

_TRUTH_BOUNDARY = {
    "strict_inner_composite_velocity_executable": True,
    "strict_inner_composite_velocity_dt_executable": True,
    "strict_inner_composite_velocity_jacobian_executable": True,
    "strict_inner_composite_vorticity_executable": True,
    "vorticity_derived_from_existing_spatial_candidate": True,
    "new_derivative_realization_introduced": False,
    "oscillatory_complete_curl_path_inherited": True,
    "inner_leading_is_final_corrected_fixed_point": False,
    "global_leading_velocity_materialized": False,
    "outer_join_materialized": False,
    "matched_pressure_included": False,
    "restricted_forcing_included": False,
    "agent3_correction_velocity_included": False,
    "complete_kokuno_candidate_assembled": False,
    "whole_domain_vorticity_morphology_verified": False,
    "three_resolution_morphology_verified": False,
    "complete_ns_momentum_residual_formed": False,
    "heldout_ns_residual_assessed": False,
    "st006_same_protocol_comparison_valid": False,
    "residual_reduction_claimed": False,
    "paper_exact": False,
    "pde_validated": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256_payload(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class StrictInnerVorticityEvaluation:
    velocity: np.ndarray
    velocity_dt: np.ndarray
    velocity_jacobian: np.ndarray
    divergence: np.ndarray
    vorticity: np.ndarray


@dataclass(frozen=True)
class KokunoStrictInnerLeadingOscillatoryVorticityCandidate:
    """Checksum-bound strict-inner candidate exposing curl of the public velocity."""

    inner_leading_backend: StrictInnerSpatialBackend

    @property
    def spatial_candidate(self) -> KokunoStrictInnerLeadingOscillatorySpatialCandidate:
        return KokunoStrictInnerLeadingOscillatorySpatialCandidate(self.inner_leading_backend)

    @property
    def velocity_candidate_sha256(self) -> str:
        return self.spatial_candidate.velocity_candidate_sha256

    @property
    def differentiable_sha256(self) -> str:
        return self.spatial_candidate.differentiable_sha256

    @property
    def spatial_sha256(self) -> str:
        return self.spatial_candidate.spatial_sha256

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.spatial_candidate.velocity(x, y, z, t)

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.spatial_candidate.velocity_dt(x, y, z, t)

    def velocity_jacobian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.spatial_candidate.velocity_jacobian(x, y, z, t)

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> StrictInnerVorticityEvaluation:
        spatial = self.spatial_candidate.evaluate(x, y, z, t)
        vorticity = np.asarray(spatial.vorticity, dtype=float)
        expected = np.asarray(spatial.velocity, dtype=float).shape
        if vorticity.shape != expected:
            raise RuntimeError(f"composite vorticity must have shape {expected}, got {vorticity.shape}")
        if not np.all(np.isfinite(vorticity)):
            raise RuntimeError("composite vorticity became non-finite")
        return StrictInnerVorticityEvaluation(
            velocity=np.asarray(spatial.velocity, dtype=float),
            velocity_dt=np.asarray(spatial.velocity_dt, dtype=float),
            velocity_jacobian=np.asarray(spatial.velocity_jacobian, dtype=float),
            divergence=np.asarray(spatial.divergence, dtype=float),
            vorticity=vorticity,
        )

    def vorticity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return curl(u) from the existing checksum-bound spatial Jacobian."""
        return self.evaluate(x, y, z, t).vorticity

    def vorticity_configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "velocity_candidate_sha256": self.velocity_candidate_sha256,
            "differentiable_sha256": self.differentiable_sha256,
            "spatial_sha256": self.spatial_sha256,
            "composition": "vorticity=curl_from_existing_velocity_jacobian",
            "jacobian_convention": "J[component,axis]=partial_axis velocity_component",
            "curl_components": ["J[2,1]-J[1,2]", "J[0,2]-J[2,0]", "J[1,0]-J[0,1]"],
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
            },
        }

    @property
    def vorticity_sha256(self) -> str:
        return _sha256_payload(self.vorticity_configuration())

    def manifest(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "velocity_candidate_sha256": self.velocity_candidate_sha256,
            "differentiable_sha256": self.differentiable_sha256,
            "spatial_sha256": self.spatial_sha256,
            "vorticity_configuration": self.vorticity_configuration(),
            "vorticity_sha256": self.vorticity_sha256,
            "truth_boundary": dict(_TRUTH_BOUNDARY),
            "scientific_gates": {
                "momentum_max_l2": 1.0e-3,
                "divergence_max_l2": 1.0e-5,
                "free_residual_defined_forcing_forbidden": True,
            },
        }
        payload["manifest_sha256"] = _sha256_payload(payload)
        return payload

    def save_manifest(self, path: str | Path) -> dict[str, Any]:
        payload = self.manifest()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def from_manifest(
        cls,
        inner_leading_backend: StrictInnerSpatialBackend,
        payload: Mapping[str, Any],
    ) -> "KokunoStrictInnerLeadingOscillatoryVorticityCandidate":
        if not isinstance(payload, Mapping):
            raise TypeError("vorticity candidate manifest must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("vorticity candidate manifest schema mismatch")
        manifest = dict(payload)
        supplied_manifest_sha = manifest.pop("manifest_sha256", None)
        if supplied_manifest_sha != _sha256_payload(manifest):
            raise ValueError("vorticity candidate manifest sha256 mismatch")

        obj = cls(inner_leading_backend=inner_leading_backend)
        expected = obj.manifest()
        for key in (
            "velocity_candidate_sha256",
            "differentiable_sha256",
            "spatial_sha256",
            "vorticity_configuration",
            "vorticity_sha256",
            "truth_boundary",
            "scientific_gates",
        ):
            if payload.get(key) != expected[key]:
                raise ValueError(f"vorticity candidate {key} changed")
        return obj

    @classmethod
    def load_manifest(
        cls,
        inner_leading_backend: StrictInnerSpatialBackend,
        path: str | Path,
    ) -> "KokunoStrictInnerLeadingOscillatoryVorticityCandidate":
        return cls.from_manifest(inner_leading_backend, json.loads(Path(path).read_text()))


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(KokunoStrictInnerLeadingOscillatoryVorticityCandidate.vorticity)
    forbidden = {
        "residual", "defect", "mean", "pressure", "forcing", "target", "gain",
        "damping", "delta_a", "delta_y", "nu", "viscosity", "spatial_step", "threshold",
    }
    return {
        "task": TASK,
        "schema": SCHEMA,
        "vorticity_signature": str(signature),
        "forbidden_vorticity_inputs_present": bool(forbidden & set(signature.parameters)),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
