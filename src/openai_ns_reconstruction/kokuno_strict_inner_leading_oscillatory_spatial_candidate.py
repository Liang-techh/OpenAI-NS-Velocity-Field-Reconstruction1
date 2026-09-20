"""First-class strict-inner leading + oscillatory candidate with spatial Jacobian.

Kokuno Agent 2 owns the oscillatory/complete-curl lane.  Agent-2 #866 already
materializes the currently legal strict-inner candidate surface

    velocity    = u_inner_strict + u_osc_complete_curl,
    velocity_dt = d_t(u_inner_strict + u_osc_complete_curl).

This module adds exactly one candidate-facing capability:

    velocity_jacobian[..., component, axis] = partial_axis velocity_component.

The Agent-1 contribution is consumed through its analytic PA.10 strict-inner
``velocity_jacobian`` seam.  The oscillatory contribution reuses A2's frozen
complete-curl velocity and its centered Cartesian FD6 derivative diagnostic at
a fixed repository step.  No caller may tune that derivative step here.

Scientific boundary
-------------------
The corrected 2026-09-09 Kokuno reconstruction is provenance for the upstream
PA.10 formulas and localized complete-curl structure.  The public-z pullback,
repository-autonomous oscillatory parameters, FD6 oscillatory Jacobian, and
this additive Cartesian composition are repository realizations, not
paper-exact hidden data.  This remains strict-inner only: no corrected/global
leading join, Agent-3 correction velocity, matched pressure, restricted
forcing, or complete Navier--Stokes residual is introduced.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping, Protocol, runtime_checkable

import numpy as np

from .kokuno_public_oscillatory_vorticity_diagnostic import evaluate_vorticity_osc_fd6
from .kokuno_strict_inner_leading_oscillatory_differentiable_candidate import (
    KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate,
)

TASK = "KOKUNO-A2-STRICT-INNER-SPATIAL-CANDIDATE-066"
SCHEMA = "kokuno-a2-strict-inner-spatial-candidate-v1"
PARENT_AGENT2_PR = 866
PARENT_AGENT2_HEAD = "68f84128f07b6743d368a9ab7a051e441f5a59b2"
AGENT1_RUNTIME_PR = 848
AGENT1_RUNTIME_HEAD = "300919075c07396e4e466f553be78eaae0707df2"
AGENT1_RUNTIME_BLOB = "5076ffa69748c2d7e68289a155b67288eab6dab4"
AGENT1_SPATIAL_PR = 819
AGENT1_SPATIAL_HEAD = "cb54e6e1a9061cacf78e446dc66bc992f64a0f8d"
AGENT1_SPATIAL_BLOB = "3dbde84d431ba72e847c631543137a5fdfb69459"
OSCILLATORY_JACOBIAN_BLOB = "4ae525bad1e4b83c9dcabc1a97d3931562d099a5"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
OSCILLATORY_SPATIAL_STEP = 1.0e-3

_TRUTH_BOUNDARY = {
    "strict_inner_composite_velocity_executable": True,
    "strict_inner_composite_velocity_dt_executable": True,
    "strict_inner_composite_velocity_jacobian_executable": True,
    "base_velocity_candidate_identity_preserved": True,
    "base_differentiable_candidate_identity_preserved": True,
    "inner_leading_velocity_jacobian_analytic": True,
    "oscillatory_velocity_jacobian_uses_fixed_fd6": True,
    "oscillatory_complete_curl_path_inherited": True,
    "agent1_inner_domain_fail_closed": True,
    "inner_leading_is_final_corrected_fixed_point": False,
    "global_leading_velocity_materialized": False,
    "outer_join_materialized": False,
    "matched_pressure_included": False,
    "restricted_forcing_included": False,
    "agent3_correction_velocity_included": False,
    "complete_kokuno_candidate_assembled": False,
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


def _validate_sha256(value: Any, name: str) -> str:
    text = str(value)
    if len(text) != 64 or any(c not in "0123456789abcdef" for c in text):
        raise ValueError(f"{name} must be one lowercase hexadecimal sha256")
    return text


@runtime_checkable
class StrictInnerSpatialBackend(Protocol):
    """Minimal Agent-1 contract required by the spatial artifact extension."""

    @property
    def field_sha256(self) -> str: ...

    @property
    def temporal_derivative_sha256(self) -> str: ...

    @property
    def spatial_derivative_sha256(self) -> str: ...

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...

    def velocity_jacobian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...


@dataclass(frozen=True)
class StrictInnerSpatialEvaluation:
    velocity: np.ndarray
    velocity_dt: np.ndarray
    inner_leading_jacobian: np.ndarray
    oscillatory_jacobian: np.ndarray
    velocity_jacobian: np.ndarray
    jacobian_additive_closure_abs_max: float
    divergence: np.ndarray
    vorticity: np.ndarray


@dataclass(frozen=True)
class KokunoStrictInnerLeadingOscillatorySpatialCandidate:
    """Checksum-bound strict-inner candidate exposing a spatial Jacobian."""

    inner_leading_backend: StrictInnerSpatialBackend

    def __post_init__(self) -> None:
        for name in ("velocity", "velocity_dt", "velocity_jacobian"):
            if not callable(getattr(self.inner_leading_backend, name, None)):
                raise TypeError(f"inner_leading_backend must expose {name}(x,y,z,t)")
        for name in (
            "field_sha256",
            "temporal_derivative_sha256",
            "spatial_derivative_sha256",
        ):
            _validate_sha256(getattr(self.inner_leading_backend, name, ""), f"inner_leading_backend.{name}")

    @property
    def differentiable_candidate(self) -> KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate:
        return KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate(self.inner_leading_backend)

    @property
    def velocity_candidate_sha256(self) -> str:
        return self.differentiable_candidate.velocity_candidate_sha256

    @property
    def differentiable_sha256(self) -> str:
        return self.differentiable_candidate.differentiable_sha256

    @property
    def inner_spatial_derivative_sha256(self) -> str:
        return _validate_sha256(
            self.inner_leading_backend.spatial_derivative_sha256,
            "inner_leading_backend.spatial_derivative_sha256",
        )

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.differentiable_candidate.velocity(x, y, z, t)

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.differentiable_candidate.velocity_dt(x, y, z, t)

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> StrictInnerSpatialEvaluation:
        timed = self.differentiable_candidate.evaluate(x, y, z, t)
        velocity = np.asarray(timed.velocity, dtype=float)
        expected_jacobian_shape = velocity.shape[:-1] + (3, 3)

        inner = np.asarray(
            self.inner_leading_backend.velocity_jacobian(x, y, z, t), dtype=float
        )
        if inner.shape != expected_jacobian_shape:
            raise ValueError(
                "inner leading velocity_jacobian must have shape "
                f"{expected_jacobian_shape}, got {inner.shape}"
            )
        if not np.all(np.isfinite(inner)):
            raise ValueError("inner leading velocity_jacobian must be finite")

        osc = evaluate_vorticity_osc_fd6(
            x, y, z, t, spatial_step=OSCILLATORY_SPATIAL_STEP
        )
        oscillatory = np.asarray(osc["velocity_gradient_fd6"], dtype=float)
        if oscillatory.shape != expected_jacobian_shape:
            raise RuntimeError("oscillatory FD6 Jacobian returned an unexpected shape")
        if not np.all(np.isfinite(oscillatory)):
            raise RuntimeError("oscillatory FD6 Jacobian became non-finite")

        total = inner + oscillatory
        closure = (
            float(np.max(np.abs(total - inner - oscillatory))) if total.size else 0.0
        )
        if closure > 1.0e-13:
            raise RuntimeError("strict-inner spatial Jacobian additive closure failed")

        divergence = np.trace(total, axis1=-2, axis2=-1)
        vorticity = np.stack(
            (
                total[..., 2, 1] - total[..., 1, 2],
                total[..., 0, 2] - total[..., 2, 0],
                total[..., 1, 0] - total[..., 0, 1],
            ),
            axis=-1,
        )
        return StrictInnerSpatialEvaluation(
            velocity=velocity,
            velocity_dt=np.asarray(timed.velocity_dt, dtype=float),
            inner_leading_jacobian=inner,
            oscillatory_jacobian=oscillatory,
            velocity_jacobian=total,
            jacobian_additive_closure_abs_max=closure,
            divergence=divergence,
            vorticity=vorticity,
        )

    def velocity_jacobian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return J[...,component,axis] with a fixed A2 oscillatory FD6 seam."""
        return self.evaluate(x, y, z, t).velocity_jacobian

    def spatial_configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "base_differentiable_candidate": {
                "velocity_candidate_sha256": self.velocity_candidate_sha256,
                "differentiable_sha256": self.differentiable_sha256,
                "parent_agent2_pr": PARENT_AGENT2_PR,
                "parent_agent2_head": PARENT_AGENT2_HEAD,
            },
            "inner_leading_spatial_derivative": {
                "spatial_derivative_sha256": self.inner_spatial_derivative_sha256,
                "runtime_agent1_pr": AGENT1_RUNTIME_PR,
                "runtime_agent1_head": AGENT1_RUNTIME_HEAD,
                "runtime_agent1_source_blob": AGENT1_RUNTIME_BLOB,
                "source_agent1_spatial_pr": AGENT1_SPATIAL_PR,
                "source_agent1_spatial_head": AGENT1_SPATIAL_HEAD,
                "source_agent1_spatial_blob": AGENT1_SPATIAL_BLOB,
                "semantics": "analytic PA.10 strict-inner Cartesian velocity Jacobian",
            },
            "oscillatory_spatial_derivative": {
                "source_blob": OSCILLATORY_JACOBIAN_BLOB,
                "operator": "centered_cartesian_fd6",
                "fixed_spatial_step": OSCILLATORY_SPATIAL_STEP,
                "semantics": "numerical derivative of frozen complete-curl oscillatory velocity",
            },
            "composition": "velocity_jacobian=jacobian_inner_strict+jacobian_osc_complete_curl",
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
            },
        }

    @property
    def spatial_sha256(self) -> str:
        return _sha256_payload(self.spatial_configuration())

    def manifest(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "velocity_candidate_sha256": self.velocity_candidate_sha256,
            "differentiable_sha256": self.differentiable_sha256,
            "spatial_configuration": self.spatial_configuration(),
            "spatial_sha256": self.spatial_sha256,
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
    ) -> "KokunoStrictInnerLeadingOscillatorySpatialCandidate":
        if not isinstance(payload, Mapping):
            raise TypeError("spatial candidate manifest must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("spatial candidate manifest schema mismatch")
        manifest = dict(payload)
        supplied_manifest_sha = manifest.pop("manifest_sha256", None)
        if supplied_manifest_sha != _sha256_payload(manifest):
            raise ValueError("spatial candidate manifest sha256 mismatch")

        obj = cls(inner_leading_backend=inner_leading_backend)
        if payload.get("velocity_candidate_sha256") != obj.velocity_candidate_sha256:
            raise ValueError("base velocity candidate semantic sha256 mismatch")
        if payload.get("differentiable_sha256") != obj.differentiable_sha256:
            raise ValueError("base differentiable candidate semantic sha256 mismatch")
        if payload.get("spatial_configuration") != obj.spatial_configuration():
            raise ValueError("spatial candidate configuration/provenance changed")
        if payload.get("spatial_sha256") != obj.spatial_sha256:
            raise ValueError("spatial candidate semantic sha256 mismatch")
        if payload.get("truth_boundary") != dict(_TRUTH_BOUNDARY):
            raise ValueError("spatial candidate truth boundary changed")
        if payload.get("scientific_gates") != obj.manifest()["scientific_gates"]:
            raise ValueError("spatial candidate scientific gates changed")
        return obj

    @classmethod
    def load_manifest(
        cls,
        inner_leading_backend: StrictInnerSpatialBackend,
        path: str | Path,
    ) -> "KokunoStrictInnerLeadingOscillatorySpatialCandidate":
        return cls.from_manifest(
            inner_leading_backend, json.loads(Path(path).read_text())
        )


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(
        KokunoStrictInnerLeadingOscillatorySpatialCandidate.velocity_jacobian
    )
    forbidden = {
        "residual",
        "defect",
        "mean",
        "pressure",
        "forcing",
        "target",
        "gain",
        "damping",
        "delta_a",
        "delta_y",
        "nu",
        "viscosity",
        "spatial_step",
        "threshold",
    }
    return {
        "task": TASK,
        "schema": SCHEMA,
        "velocity_jacobian_signature": str(signature),
        "forbidden_velocity_jacobian_inputs_present": bool(forbidden & set(signature.parameters)),
        "oscillatory_spatial_step": OSCILLATORY_SPATIAL_STEP,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
