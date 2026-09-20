"""Executable strict-inner PA.10 leading + complete-curl oscillatory candidate.

Kokuno Agent 2 owns only the oscillatory/complete-curl lane.  This module is a
small composition layer over an Agent-1 *inner contraction-center* velocity and
the already-frozen Agent-2 localized complete-curl oscillatory field.  It makes
priority-5's currently legal object explicit as one callable

    velocity(x,y,z,t) = u_inner(x,y,z,t) + u_osc(x,y,z,t).

The word ``strict-inner`` is part of the contract.  Agent 1 still owns the
missing corrected/global leading construction and outer join; any Agent-1
domain failure propagates rather than being replaced by a zero extension or an
A2 taper.

Source / realization boundary
-----------------------------
The corrected Kokuno reader dated 2026-09-09 is structural provenance for the
PA.10 leading formulas and the localized complete-curl organization.  The
public-z pullback, frozen autonomous oscillatory parameters, this additive
composition, its manifest/identity scheme, and finite probe receipts are
repository realizations.  This object is not paper-exact, does not include a
matched pressure/restricted forcing/Agent-3 correction, and is not a complete
Navier--Stokes candidate.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping, Protocol, runtime_checkable

import numpy as np

from .kokuno_public_z_pullback_velocity import default_field, velocity_osc

TASK = "KOKUNO-A2-STRICT-INNER-LEADING-OSCILLATORY-CANDIDATE-064"
SCHEMA = "kokuno-a2-strict-inner-leading-oscillatory-candidate-v1"
PARENT_AGENT2_PR = 849
PARENT_AGENT2_HEAD = "616e61abd2420a7c7aff3590a08393360039cd1a"
AGENT1_REFERENCE_PR = 848
AGENT1_REFERENCE_HEAD = "300919075c07396e4e466f553be78eaae0707df2"
AGENT1_REFERENCE_BLOB = "5076ffa69748c2d7e68289a155b67288eab6dab4"
OSCILLATORY_VELOCITY_BLOB = "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

_TRUTH_BOUNDARY = {
    "strict_inner_composite_velocity_executable": True,
    "strict_inner_composite_manifest_executable": True,
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
class StrictInnerVelocityBackend(Protocol):
    """Minimal Agent-1 contract consumed by this A2 candidate composition."""

    @property
    def field_sha256(self) -> str: ...

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...


@dataclass(frozen=True)
class StrictInnerLeadingOscillatoryEvaluation:
    inner_leading_velocity: np.ndarray
    oscillatory_velocity: np.ndarray
    velocity: np.ndarray
    additive_closure_abs_max: float


@dataclass(frozen=True)
class KokunoStrictInnerLeadingOscillatoryCandidate:
    """Callable, checksum-bound strict-inner ``u_inner + u_osc`` candidate view."""

    inner_leading_backend: StrictInnerVelocityBackend

    def __post_init__(self) -> None:
        method = getattr(self.inner_leading_backend, "velocity", None)
        if not callable(method):
            raise TypeError("inner_leading_backend must expose velocity(x,y,z,t)")
        _validate_sha256(
            getattr(self.inner_leading_backend, "field_sha256", ""),
            "inner_leading_backend.field_sha256",
        )

    @property
    def inner_field_sha256(self) -> str:
        return _validate_sha256(
            self.inner_leading_backend.field_sha256,
            "inner_leading_backend.field_sha256",
        )

    @property
    def oscillatory_field_sha256(self) -> str:
        # Inherited property hashes the subclass' full public-z-pullback payload.
        return _validate_sha256(default_field().sha256, "oscillatory field sha256")

    @staticmethod
    def _broadcast_xyzt(x: Any, y: Any, z: Any, t: Any) -> tuple[np.ndarray, ...]:
        arrays = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        if not all(np.all(np.isfinite(a)) for a in arrays):
            raise ValueError("x, y, z and t must be finite and broadcastable")
        return tuple(arrays)

    def _inner_velocity(
        self, x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray
    ) -> np.ndarray:
        value = np.asarray(
            self.inner_leading_backend.velocity(x, y, z, t), dtype=float
        )
        expected = x.shape + (3,)
        if value.shape != expected:
            raise ValueError(
                f"inner leading velocity must have shape {expected}, got {value.shape}"
            )
        if not np.all(np.isfinite(value)):
            raise ValueError("inner leading velocity must be finite")
        return value

    def evaluate(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> StrictInnerLeadingOscillatoryEvaluation:
        """Evaluate the exact additive composition on the current strict-inner domain."""
        xb, yb, zb, tb = self._broadcast_xyzt(x, y, z, t)
        inner = self._inner_velocity(xb, yb, zb, tb)
        osc = np.asarray(velocity_osc(xb, yb, zb, tb), dtype=float)
        expected = xb.shape + (3,)
        if osc.shape != expected or not np.all(np.isfinite(osc)):
            raise RuntimeError("oscillatory velocity is non-finite or has the wrong shape")
        total = inner + osc
        closure = float(np.max(np.abs(total - inner - osc))) if total.size else 0.0
        if closure > 1.0e-14:
            raise RuntimeError("strict-inner leading/oscillatory additive closure failed")
        return StrictInnerLeadingOscillatoryEvaluation(
            inner_leading_velocity=inner,
            oscillatory_velocity=osc,
            velocity=total,
            additive_closure_abs_max=closure,
        )

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Unified currently-legal candidate callable: ``u_inner + u_osc``."""
        return self.evaluate(x, y, z, t).velocity

    def field_configuration(self) -> dict[str, Any]:
        """Semantic identity for this additive candidate composition."""
        return {
            "schema": SCHEMA,
            "composition": "velocity=u_inner_strict+u_osc_complete_curl",
            "inner_leading": {
                "field_sha256": self.inner_field_sha256,
                "reference_agent1_pr": AGENT1_REFERENCE_PR,
                "reference_agent1_head": AGENT1_REFERENCE_HEAD,
                "reference_agent1_source_blob": AGENT1_REFERENCE_BLOB,
                "scope": "PA.10 inner contraction-center only",
            },
            "oscillatory": {
                "field_sha256": self.oscillatory_field_sha256,
                "source_blob": OSCILLATORY_VELOCITY_BLOB,
                "callable": "kokuno_public_z_pullback_velocity.velocity_osc",
                "construction": "localized vector potential -> complete curl",
            },
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
            },
        }

    @property
    def candidate_sha256(self) -> str:
        return _sha256_payload(self.field_configuration())

    def manifest(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "field_configuration": self.field_configuration(),
            "candidate_sha256": self.candidate_sha256,
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
        inner_leading_backend: StrictInnerVelocityBackend,
        payload: Mapping[str, Any],
    ) -> "KokunoStrictInnerLeadingOscillatoryCandidate":
        """Rebind a manifest to a supplied backend and verify semantic identity.

        Agent 2 does not reconstruct or vendor Agent-1 mathematics here.  The
        caller supplies the Agent-1 backend; the manifest must match both that
        backend's field identity and the frozen A2 oscillatory realization.
        """
        if not isinstance(payload, Mapping):
            raise TypeError("candidate manifest must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("candidate manifest schema mismatch")
        manifest = dict(payload)
        supplied_manifest_sha = manifest.pop("manifest_sha256", None)
        if supplied_manifest_sha != _sha256_payload(manifest):
            raise ValueError("candidate manifest sha256 mismatch")

        obj = cls(inner_leading_backend=inner_leading_backend)
        if payload.get("field_configuration") != obj.field_configuration():
            raise ValueError("candidate field configuration/provenance changed")
        if payload.get("candidate_sha256") != obj.candidate_sha256:
            raise ValueError("candidate semantic sha256 mismatch")
        expected_truth = dict(_TRUTH_BOUNDARY)
        if payload.get("truth_boundary") != expected_truth:
            raise ValueError("candidate truth boundary changed")
        if payload.get("scientific_gates") != obj.manifest()["scientific_gates"]:
            raise ValueError("candidate scientific gates changed")
        return obj

    @classmethod
    def load_manifest(
        cls,
        inner_leading_backend: StrictInnerVelocityBackend,
        path: str | Path,
    ) -> "KokunoStrictInnerLeadingOscillatoryCandidate":
        return cls.from_manifest(
            inner_leading_backend, json.loads(Path(path).read_text())
        )


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(KokunoStrictInnerLeadingOscillatoryCandidate.velocity)
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
        "oscillatory_velocity_blob": OSCILLATORY_VELOCITY_BLOB,
        "forbidden_velocity_inputs_present": sorted(
            forbidden.intersection(signature.parameters)
        ),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
