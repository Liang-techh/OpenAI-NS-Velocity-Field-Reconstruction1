"""First-class strict-inner leading + oscillatory candidate with ``velocity_dt``.

Kokuno Agent 2 owns the oscillatory/complete-curl lane.  PR #857 already
materializes the currently legal first-class velocity artifact

    velocity = u_inner_strict + u_osc_complete_curl.

This module adds exactly one candidate-facing capability: a checksum-bound
fixed-Cartesian ``velocity_dt`` surface.  It deliberately reuses the existing
Agent-2 #821 additive time-derivative seam and the Agent-1 backend's inherited
analytic ``velocity_dt``.  The #857 velocity object and its semantic SHA are
preserved rather than rebuilt or retuned.

The scope is still strict-inner PA.10 only.  Any Agent-1 domain failure
propagates.  There is no outer/global join, matched pressure, restricted
forcing, Agent-3 correction velocity, complete momentum residual, or paper-
exact claim here.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping, Protocol, runtime_checkable

import numpy as np

from .kokuno_public_inner_leading_oscillatory_time_derivative import (
    evaluate_inner_leading_oscillatory_time_derivative,
)
from .kokuno_strict_inner_leading_oscillatory_candidate import (
    KokunoStrictInnerLeadingOscillatoryCandidate,
)

TASK = "KOKUNO-A2-STRICT-INNER-DIFFERENTIABLE-CANDIDATE-065"
SCHEMA = "kokuno-a2-strict-inner-differentiable-candidate-v1"
PARENT_AGENT2_PR = 857
PARENT_AGENT2_HEAD = "a586b7afe4bb47dbfd4b25177586c79d362afff9"
AGENT1_REFERENCE_PR = 848
AGENT1_REFERENCE_HEAD = "300919075c07396e4e466f553be78eaae0707df2"
AGENT1_REFERENCE_BLOB = "5076ffa69748c2d7e68289a155b67288eab6dab4"
A2_TIME_COMPOSITION_BLOB = "c780aea6d55daa62734b7906fe4e9e4f9ef752f0"
OSCILLATORY_TIME_DERIVATIVE_BLOB = "e80fe68ad031197de95a8cb8692c02c6a1005936"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

_TRUTH_BOUNDARY = {
    "strict_inner_composite_velocity_executable": True,
    "strict_inner_composite_velocity_dt_executable": True,
    "base_velocity_candidate_identity_preserved": True,
    "production_velocity_dt_uses_finite_difference": False,
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
class StrictInnerDifferentiableBackend(Protocol):
    """Minimal Agent-1 contract needed by this artifact-layer extension."""

    @property
    def field_sha256(self) -> str: ...

    @property
    def temporal_derivative_sha256(self) -> str: ...

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...


@dataclass(frozen=True)
class StrictInnerDifferentiableEvaluation:
    inner_leading_velocity: np.ndarray
    oscillatory_velocity: np.ndarray
    velocity: np.ndarray
    inner_leading_velocity_dt: np.ndarray
    oscillatory_velocity_dt: np.ndarray
    velocity_dt: np.ndarray
    velocity_replay_abs_max: float
    velocity_dt_additive_closure_abs_max: float


@dataclass(frozen=True)
class KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate:
    """Checksum-bound strict-inner candidate exposing both velocity and velocity_dt."""

    inner_leading_backend: StrictInnerDifferentiableBackend

    def __post_init__(self) -> None:
        if not callable(getattr(self.inner_leading_backend, "velocity", None)):
            raise TypeError("inner_leading_backend must expose velocity(x,y,z,t)")
        if not callable(getattr(self.inner_leading_backend, "velocity_dt", None)):
            raise TypeError("inner_leading_backend must expose velocity_dt(x,y,z,t)")
        _validate_sha256(
            getattr(self.inner_leading_backend, "field_sha256", ""),
            "inner_leading_backend.field_sha256",
        )
        _validate_sha256(
            getattr(self.inner_leading_backend, "temporal_derivative_sha256", ""),
            "inner_leading_backend.temporal_derivative_sha256",
        )

    @property
    def base_candidate(self) -> KokunoStrictInnerLeadingOscillatoryCandidate:
        return KokunoStrictInnerLeadingOscillatoryCandidate(self.inner_leading_backend)

    @property
    def inner_field_sha256(self) -> str:
        return _validate_sha256(
            self.inner_leading_backend.field_sha256,
            "inner_leading_backend.field_sha256",
        )

    @property
    def inner_temporal_derivative_sha256(self) -> str:
        return _validate_sha256(
            self.inner_leading_backend.temporal_derivative_sha256,
            "inner_leading_backend.temporal_derivative_sha256",
        )

    @property
    def velocity_candidate_sha256(self) -> str:
        """Preserve #857's semantic identity for the unchanged velocity callable."""
        return self.base_candidate.candidate_sha256

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Replay the exact #857 velocity surface without changing its semantics."""
        return self.base_candidate.velocity(x, y, z, t)

    def evaluate(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> StrictInnerDifferentiableEvaluation:
        base = self.base_candidate.evaluate(x, y, z, t)
        timed = evaluate_inner_leading_oscillatory_time_derivative(
            self.inner_leading_backend, x, y, z, t
        )

        expected = np.asarray(base.velocity, dtype=float)
        replay = np.asarray(timed.inner_plus_oscillatory_velocity, dtype=float)
        if replay.shape != expected.shape:
            raise RuntimeError("time seam and base candidate returned different shapes")
        velocity_replay_abs_max = (
            float(np.max(np.abs(expected - replay))) if expected.size else 0.0
        )
        if velocity_replay_abs_max > 1.0e-14:
            raise RuntimeError("time seam changed the frozen #857 candidate velocity")

        dt_inner = np.asarray(timed.inner_leading_velocity_dt, dtype=float)
        dt_osc = np.asarray(timed.oscillatory_velocity_dt, dtype=float)
        dt_total = np.asarray(timed.inner_plus_oscillatory_velocity_dt, dtype=float)
        if dt_total.shape != expected.shape:
            raise RuntimeError("composite velocity_dt has the wrong shape")
        if not (
            np.all(np.isfinite(dt_inner))
            and np.all(np.isfinite(dt_osc))
            and np.all(np.isfinite(dt_total))
        ):
            raise RuntimeError("composite velocity_dt became non-finite")
        dt_closure = (
            float(np.max(np.abs(dt_total - dt_inner - dt_osc)))
            if dt_total.size
            else 0.0
        )
        if dt_closure > 1.0e-13:
            raise RuntimeError("strict-inner velocity_dt additive closure failed")

        return StrictInnerDifferentiableEvaluation(
            inner_leading_velocity=np.asarray(base.inner_leading_velocity, dtype=float),
            oscillatory_velocity=np.asarray(base.oscillatory_velocity, dtype=float),
            velocity=expected,
            inner_leading_velocity_dt=dt_inner,
            oscillatory_velocity_dt=dt_osc,
            velocity_dt=dt_total,
            velocity_replay_abs_max=velocity_replay_abs_max,
            velocity_dt_additive_closure_abs_max=dt_closure,
        )

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Fixed-Cartesian derivative of the current strict-inner velocity artifact."""
        return self.evaluate(x, y, z, t).velocity_dt

    def differentiable_configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "base_velocity_candidate": {
                "candidate_sha256": self.velocity_candidate_sha256,
                "field_configuration": self.base_candidate.field_configuration(),
                "parent_agent2_pr": PARENT_AGENT2_PR,
                "parent_agent2_head": PARENT_AGENT2_HEAD,
            },
            "inner_leading_time_derivative": {
                "field_sha256": self.inner_field_sha256,
                "temporal_derivative_sha256": self.inner_temporal_derivative_sha256,
                "reference_agent1_pr": AGENT1_REFERENCE_PR,
                "reference_agent1_head": AGENT1_REFERENCE_HEAD,
                "reference_agent1_source_blob": AGENT1_REFERENCE_BLOB,
                "scope": "PA.10 inner contraction-center only",
            },
            "oscillatory_time_derivative": {
                "source_blob": OSCILLATORY_TIME_DERIVATIVE_BLOB,
                "composition_source_blob": A2_TIME_COMPOSITION_BLOB,
                "semantics": (
                    "exact derivative of the frozen repository-autonomous "
                    "oscillatory time modulation"
                ),
            },
            "composition": (
                "velocity_dt=velocity_dt_inner_strict+"
                "velocity_dt_osc_complete_curl"
            ),
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
            },
        }

    @property
    def differentiable_sha256(self) -> str:
        return _sha256_payload(self.differentiable_configuration())

    def manifest(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "velocity_candidate_sha256": self.velocity_candidate_sha256,
            "differentiable_configuration": self.differentiable_configuration(),
            "differentiable_sha256": self.differentiable_sha256,
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
        inner_leading_backend: StrictInnerDifferentiableBackend,
        payload: Mapping[str, Any],
    ) -> "KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate":
        if not isinstance(payload, Mapping):
            raise TypeError("differentiable candidate manifest must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("differentiable candidate manifest schema mismatch")
        manifest = dict(payload)
        supplied_manifest_sha = manifest.pop("manifest_sha256", None)
        if supplied_manifest_sha != _sha256_payload(manifest):
            raise ValueError("differentiable candidate manifest sha256 mismatch")

        obj = cls(inner_leading_backend=inner_leading_backend)
        if payload.get("velocity_candidate_sha256") != obj.velocity_candidate_sha256:
            raise ValueError("base velocity candidate semantic sha256 mismatch")
        if payload.get("differentiable_configuration") != obj.differentiable_configuration():
            raise ValueError("differentiable candidate configuration/provenance changed")
        if payload.get("differentiable_sha256") != obj.differentiable_sha256:
            raise ValueError("differentiable candidate semantic sha256 mismatch")
        if payload.get("truth_boundary") != dict(_TRUTH_BOUNDARY):
            raise ValueError("differentiable candidate truth boundary changed")
        if payload.get("scientific_gates") != obj.manifest()["scientific_gates"]:
            raise ValueError("differentiable candidate scientific gates changed")
        return obj

    @classmethod
    def load_manifest(
        cls,
        inner_leading_backend: StrictInnerDifferentiableBackend,
        path: str | Path,
    ) -> "KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate":
        return cls.from_manifest(
            inner_leading_backend, json.loads(Path(path).read_text())
        )


def public_contract() -> dict[str, Any]:
    velocity_signature = inspect.signature(
        KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate.velocity
    )
    velocity_dt_signature = inspect.signature(
        KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate.velocity_dt
    )
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "alpha", "damping", "delta_y", "delta_a", "nu", "viscosity",
        "scientific_threshold", "spatial_step", "time_step",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "agent1_reference_pr": AGENT1_REFERENCE_PR,
        "agent1_reference_head": AGENT1_REFERENCE_HEAD,
        "agent1_reference_blob": AGENT1_REFERENCE_BLOB,
        "a2_time_composition_blob": A2_TIME_COMPOSITION_BLOB,
        "oscillatory_time_derivative_blob": OSCILLATORY_TIME_DERIVATIVE_BLOB,
        "forbidden_velocity_inputs_present": sorted(
            forbidden.intersection(velocity_signature.parameters)
        ),
        "forbidden_velocity_dt_inputs_present": sorted(
            forbidden.intersection(velocity_dt_signature.parameters)
        ),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
