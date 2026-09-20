"""Typed Agent-2 -> Agent-3 handoff for the oscillatory transport increment.

This module adds no new velocity model and performs no cylindrical mean/radial
correction work.  It wraps the already executable strict-inner Agent-2
transport increment from ``kokuno_public_inner_leading_oscillatory_transport_delta``
and gives downstream code a checksum-bound, explicitly scoped interface.

The handed quantity is

    delta_T_osc = d_t u_osc
                  + (u_inner . grad) u_osc
                  + (u_osc . grad) u_inner
                  + (u_osc . grad) u_osc
                  - nu Delta u_osc,

with repository viscosity fixed upstream at ``nu=0.01``.  It is a Cartesian,
pressure/forcing-free strict-inner transport increment.  It is *not* a complete
Navier--Stokes residual, and this module deliberately does not project an m=0
mean, invert a radial operator, or construct a correction velocity.  Those
operations remain Agent-3 ownership.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping, Protocol

import numpy as np

from .kokuno_oscillatory_source_contract import (
    default_kokuno_oscillatory_source_contract,
)
from .kokuno_public_inner_leading_oscillatory_transport_delta import (
    InnerLeadingOscillatoryTransportDeltaResult,
    evaluate_inner_leading_oscillatory_transport_delta,
)


TASK = "KOKUNO-A2-OSCILLATORY-DELTA-HANDOFF-070"
SCHEMA = "kokuno-a2-oscillatory-transport-handoff-v1"
PARENT_AGENT2_PR = 894
PARENT_AGENT2_HEAD = "a467bc175a801ed69c197fc9e38ca68ac8c0b86b"
TRANSPORT_DELTA_PR = 849
TRANSPORT_DELTA_HEAD = "616e61abd2420a7c7aff3590a08393360039cd1a"
AGENT3_CONSUMER_LANE = "cylindrical m=0 mean/radial correction; no A2 reimplementation"


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _is_sha256(value: Any) -> bool:
    if not isinstance(value, str) or len(value) != 64:
        return False
    return all(c in "0123456789abcdef" for c in value)


class StrictInnerTransportBackend(Protocol):
    field_sha256: str

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...
    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...
    def self_advection(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...
    def velocity_laplacian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...


@dataclass(frozen=True)
class OscillatoryTransportHandoffResult:
    """One evaluated raw Cartesian handoff; no mean projection is performed."""

    oscillatory_transport_increment: np.ndarray
    oscillatory_time_increment: np.ndarray
    oscillatory_nonlinear_increment: np.ndarray
    oscillatory_viscous_increment: np.ndarray
    inner_leading_transport: np.ndarray
    inner_plus_oscillatory_transport: np.ndarray
    inner_leading_velocity: np.ndarray
    oscillatory_velocity: np.ndarray
    inner_plus_oscillatory_velocity: np.ndarray
    viscosity: float
    decomposition_closure_abs_max: float


class KokunoOscillatoryTransportHandoff:
    """Checksum-bound raw transport-delta handoff from Agent 2 to Agent 3."""

    def __init__(self, inner_leading_backend: StrictInnerTransportBackend):
        self._backend = inner_leading_backend
        self._validate_backend_contract()
        self._source_contract = default_kokuno_oscillatory_source_contract()

    def _validate_backend_contract(self) -> None:
        field_sha = getattr(self._backend, "field_sha256", None)
        if not _is_sha256(field_sha):
            raise TypeError("inner_leading_backend must expose a lowercase 64-hex field_sha256")
        for name in ("velocity", "velocity_dt", "self_advection", "velocity_laplacian"):
            if not callable(getattr(self._backend, name, None)):
                raise TypeError(
                    "inner_leading_backend must expose callable "
                    f"{name}(x,y,z,t)"
                )

    @property
    def backend_field_sha256(self) -> str:
        return str(self._backend.field_sha256)

    @property
    def source_contract_sha256(self) -> str:
        return self._source_contract.contract_sha256

    def semantic_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "transport_delta_pr": TRANSPORT_DELTA_PR,
            "transport_delta_head": TRANSPORT_DELTA_HEAD,
            "backend_field_sha256": self.backend_field_sha256,
            "source_contract_sha256": self.source_contract_sha256,
            "quantity_contract": {
                "name": "oscillatory_transport_increment",
                "frame": "Cartesian",
                "formula": (
                    "d_t u_osc + (u_inner.grad)u_osc + (u_osc.grad)u_inner + "
                    "(u_osc.grad)u_osc - nu Delta u_osc"
                ),
                "pressure_included": False,
                "forcing_included": False,
                "complete_ns_residual": False,
                "mean_projection_performed": False,
                "radial_inverse_performed": False,
                "correction_velocity_constructed": False,
                "agent3_consumer_lane": AGENT3_CONSUMER_LANE,
            },
            "scientific_boundary": {
                "strict_inner_only": True,
                "global_leading_join_included": False,
                "matched_pressure_included": False,
                "restricted_forcing_included": False,
                "same_protocol_st006_comparison_valid": False,
                "residual_reduction_claimed": False,
                "paper_exact": False,
                "pde_validated": False,
            },
        }

    @property
    def handoff_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> OscillatoryTransportHandoffResult:
        parent = evaluate_inner_leading_oscillatory_transport_delta(
            self._backend, x, y, z, t
        )
        return self._from_parent_result(parent)

    @staticmethod
    def _from_parent_result(
        parent: InnerLeadingOscillatoryTransportDeltaResult,
    ) -> OscillatoryTransportHandoffResult:
        reconstructed = (
            np.asarray(parent.oscillatory_time_increment, dtype=float)
            + np.asarray(parent.oscillatory_nonlinear_increment, dtype=float)
            + np.asarray(parent.oscillatory_viscous_increment, dtype=float)
        )
        delta = np.asarray(parent.oscillatory_transport_increment, dtype=float)
        if reconstructed.shape != delta.shape:
            raise RuntimeError("oscillatory transport decomposition shape mismatch")
        if not np.all(np.isfinite(reconstructed)) or not np.all(np.isfinite(delta)):
            raise RuntimeError("oscillatory transport handoff became non-finite")
        closure = float(np.max(np.abs(reconstructed - delta))) if delta.size else 0.0
        if closure > 5.0e-13:
            raise RuntimeError("oscillatory transport handoff decomposition failed closure")

        return OscillatoryTransportHandoffResult(
            oscillatory_transport_increment=delta,
            oscillatory_time_increment=np.asarray(parent.oscillatory_time_increment, dtype=float),
            oscillatory_nonlinear_increment=np.asarray(parent.oscillatory_nonlinear_increment, dtype=float),
            oscillatory_viscous_increment=np.asarray(parent.oscillatory_viscous_increment, dtype=float),
            inner_leading_transport=np.asarray(parent.inner_leading_transport, dtype=float),
            inner_plus_oscillatory_transport=np.asarray(
                parent.inner_plus_oscillatory_transport, dtype=float
            ),
            inner_leading_velocity=np.asarray(parent.inner_leading_velocity, dtype=float),
            oscillatory_velocity=np.asarray(parent.oscillatory_velocity, dtype=float),
            inner_plus_oscillatory_velocity=np.asarray(
                parent.inner_plus_oscillatory_velocity, dtype=float
            ),
            viscosity=float(parent.viscosity),
            decomposition_closure_abs_max=closure,
        )

    def transport_delta(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.evaluate(x, y, z, t).oscillatory_transport_increment

    def nonlinear_increment(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return the raw nonlinear oscillatory increment for Agent-3 projection.

        This is intentionally still Cartesian and unprojected.  A2 does not infer
        or solve the mean/radial correction from this quantity.
        """
        return self.evaluate(x, y, z, t).oscillatory_nonlinear_increment

    def manifest(self) -> dict[str, Any]:
        self._validate_backend_contract()
        self._source_contract.validate()
        payload = self.semantic_payload()
        return {
            "payload": payload,
            "handoff_sha256": _sha256(payload),
        }

    def save_manifest(self, path: str | Path) -> dict[str, Any]:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        manifest = self.manifest()
        target.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return manifest

    @classmethod
    def load_manifest(
        cls,
        inner_leading_backend: StrictInnerTransportBackend,
        path: str | Path,
    ) -> "KokunoOscillatoryTransportHandoff":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"payload", "handoff_sha256"}:
            raise ValueError("invalid oscillatory transport handoff manifest envelope")
        payload = raw["payload"]
        if not isinstance(payload, dict) or _sha256(payload) != raw["handoff_sha256"]:
            raise ValueError("oscillatory transport handoff manifest checksum mismatch")
        obj = cls(inner_leading_backend)
        if obj.semantic_payload() != payload:
            raise ValueError("oscillatory transport handoff manifest semantics mismatch")
        return obj


def public_contract() -> dict[str, Any]:
    """Static routing contract; useful to Agent 3 without constructing a backend."""
    evaluate_sig = inspect.signature(KokunoOscillatoryTransportHandoff.evaluate)
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "alpha", "damping", "delta_y", "delta_a", "nu",
        "viscosity", "scientific_threshold", "spatial_step", "time_step",
    }
    source = default_kokuno_oscillatory_source_contract()
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "transport_delta_pr": TRANSPORT_DELTA_PR,
        "transport_delta_head": TRANSPORT_DELTA_HEAD,
        "source_contract_sha256": source.contract_sha256,
        "forbidden_evaluate_inputs_present": sorted(
            forbidden.intersection(evaluate_sig.parameters)
        ),
        "handoff_quantity": "raw Cartesian oscillatory transport increment",
        "consumer": "Kokuno Agent 3",
        "mean_projection_performed": False,
        "correction_velocity_constructed": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }
