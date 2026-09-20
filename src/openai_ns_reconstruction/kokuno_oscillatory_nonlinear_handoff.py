"""Typed A2 -> A3 handoff for the three nonlinear oscillatory transport pieces.

The existing Agent-2 transport handoff exposes the aggregate strict-inner
nonlinear increment

    N_osc = (u_inner . grad)u_osc
            + (u_osc . grad)u_inner
            + (u_osc . grad)u_osc.

For mean-defect attribution Agent 3 needs to know whether a projected signal
comes from the two leading/oscillation cross terms or from the genuinely
quadratic oscillation self-interaction.  This module exposes those already
implemented pieces separately without performing any cylindrical projection,
radial inverse, mean correction, pressure/forcing completion, or fitting.

No oscillatory parameter is added or changed.  The corrected 2026-09-09 Kokuno
reconstruction remains structural provenance only; the public-z realization,
finite-difference derivatives and this typed decomposition are repository
realizations and are not paper-exact.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_oscillatory_transport_handoff import (
    KokunoOscillatoryTransportHandoff,
    StrictInnerTransportBackend,
)
from .kokuno_public_inner_leading_oscillatory_full_advection import (
    evaluate_inner_leading_oscillatory_full_advection,
)
from .kokuno_public_inner_leading_oscillatory_transport import (
    INNER_ADVECTION_SPATIAL_STEP,
    OSCILLATORY_ADVECTION_SPATIAL_STEP,
)

TASK = "KOKUNO-A2-OSCILLATORY-NONLINEAR-HANDOFF-071"
SCHEMA = "kokuno-a2-oscillatory-nonlinear-handoff-v1"
PARENT_AGENT2_PR = 901
PARENT_AGENT2_HEAD = "4a78b06afbab95c108c184a32cb5951e66f42991"
PARENT_HANDOFF_BLOB = "57932afda3b65e4cb63e9ae7c8ceeb0aba141cc5"
FULL_ADVECTION_PR = 830
FULL_ADVECTION_BLOB = "e44ccbe993d08951f48cfb5b093ee02ea5642dcf"


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class OscillatoryNonlinearHandoffResult:
    oscillatory_nonlinear_increment: np.ndarray
    inner_advects_oscillation: np.ndarray
    oscillation_advects_inner: np.ndarray
    oscillatory_self_advection: np.ndarray
    mixed_cross_advection: np.ndarray
    decomposition_closure_abs_max: float
    parent_replay_abs_max: float


class KokunoOscillatoryNonlinearHandoff:
    """Expose the frozen nonlinear A2 pieces for downstream A3 attribution."""

    def __init__(self, inner_leading_backend: StrictInnerTransportBackend):
        self._backend = inner_leading_backend
        self._parent = KokunoOscillatoryTransportHandoff(inner_leading_backend)

    @property
    def backend_field_sha256(self) -> str:
        return self._parent.backend_field_sha256

    @property
    def source_contract_sha256(self) -> str:
        return self._parent.source_contract_sha256

    @property
    def parent_handoff_sha256(self) -> str:
        return self._parent.handoff_sha256

    def semantic_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "parent_handoff_blob": PARENT_HANDOFF_BLOB,
            "parent_handoff_sha256": self.parent_handoff_sha256,
            "full_advection_pr": FULL_ADVECTION_PR,
            "full_advection_blob": FULL_ADVECTION_BLOB,
            "backend_field_sha256": self.backend_field_sha256,
            "source_contract_sha256": self.source_contract_sha256,
            "production_steps": {
                "inner_advection_spatial_step": INNER_ADVECTION_SPATIAL_STEP,
                "oscillatory_advection_spatial_step": OSCILLATORY_ADVECTION_SPATIAL_STEP,
            },
            "quantity_contract": {
                "aggregate": "(u_inner.grad)u_osc + (u_osc.grad)u_inner + (u_osc.grad)u_osc",
                "inner_advects_oscillation": "(u_inner.grad)u_osc",
                "oscillation_advects_inner": "(u_osc.grad)u_inner",
                "oscillatory_self_advection": "(u_osc.grad)u_osc",
                "frame": "Cartesian",
                "mean_projection_performed": False,
                "radial_inverse_performed": False,
                "correction_velocity_constructed": False,
                "pressure_included": False,
                "forcing_included": False,
                "complete_ns_residual": False,
            },
            "scientific_boundary": {
                "strict_inner_only": True,
                "no_new_oscillatory_parameters": True,
                "global_leading_join_included": False,
                "same_protocol_st006_comparison_valid": False,
                "residual_reduction_claimed": False,
                "paper_exact": False,
                "pde_validated": False,
            },
        }

    @property
    def handoff_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> OscillatoryNonlinearHandoffResult:
        parent = self._parent.evaluate(x, y, z, t)
        full = evaluate_inner_leading_oscillatory_full_advection(
            self._backend,
            x,
            y,
            z,
            t,
            inner_leading_spatial_step=INNER_ADVECTION_SPATIAL_STEP,
            oscillatory_spatial_step=OSCILLATORY_ADVECTION_SPATIAL_STEP,
        )
        a = np.asarray(full.inner_advects_oscillation, dtype=float)
        b = np.asarray(full.oscillation_advects_inner, dtype=float)
        q = np.asarray(full.oscillatory_self_advection, dtype=float)
        mixed = np.asarray(full.mixed_cross_advection, dtype=float)
        reconstructed = a + b + q
        aggregate = np.asarray(parent.oscillatory_nonlinear_increment, dtype=float)
        if reconstructed.shape != aggregate.shape:
            raise RuntimeError("nonlinear oscillatory handoff shape mismatch")
        arrays = (a, b, q, mixed, reconstructed, aggregate)
        if not all(np.all(np.isfinite(value)) for value in arrays):
            raise RuntimeError("nonlinear oscillatory handoff became non-finite")
        mixed_closure = float(np.max(np.abs(mixed - (a + b)))) if mixed.size else 0.0
        decomposition = float(np.max(np.abs(aggregate - reconstructed))) if aggregate.size else 0.0
        if mixed_closure > 5.0e-13:
            raise RuntimeError("mixed nonlinear oscillatory terms failed closure")
        if decomposition > 5.0e-13:
            raise RuntimeError("nonlinear oscillatory handoff failed decomposition closure")
        return OscillatoryNonlinearHandoffResult(
            oscillatory_nonlinear_increment=aggregate,
            inner_advects_oscillation=a,
            oscillation_advects_inner=b,
            oscillatory_self_advection=q,
            mixed_cross_advection=mixed,
            decomposition_closure_abs_max=decomposition,
            parent_replay_abs_max=decomposition,
        )

    def manifest(self) -> dict[str, Any]:
        payload = self.semantic_payload()
        return {"payload": payload, "handoff_sha256": _sha256(payload)}

    def save_manifest(self, path: str | Path) -> dict[str, Any]:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        manifest = self.manifest()
        target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return manifest

    @classmethod
    def load_manifest(
        cls,
        inner_leading_backend: StrictInnerTransportBackend,
        path: str | Path,
    ) -> "KokunoOscillatoryNonlinearHandoff":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"payload", "handoff_sha256"}:
            raise ValueError("invalid nonlinear handoff manifest envelope")
        payload = raw["payload"]
        if not isinstance(payload, dict) or _sha256(payload) != raw["handoff_sha256"]:
            raise ValueError("nonlinear handoff manifest checksum mismatch")
        obj = cls(inner_leading_backend)
        if obj.semantic_payload() != payload:
            raise ValueError("nonlinear handoff manifest semantics mismatch")
        return obj


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(KokunoOscillatoryNonlinearHandoff.evaluate)
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "damping", "nu", "viscosity", "spatial_step", "time_step",
        "amplitude", "phase", "scale", "orientation", "scientific_threshold",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "forbidden_evaluate_inputs_present": sorted(forbidden.intersection(signature.parameters)),
        "consumer": "Kokuno Agent 3",
        "decomposition": [
            "inner_advects_oscillation",
            "oscillation_advects_inner",
            "oscillatory_self_advection",
        ],
        "mean_projection_performed": False,
        "correction_velocity_constructed": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }
