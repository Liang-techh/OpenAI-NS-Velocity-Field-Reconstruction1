"""Diagnostic attribution for the frozen strict-inner oscillatory nonlinear handoff.

This module does not construct a new oscillatory field and does not perform the
Agent-3 cylindrical mean projection.  It summarizes the three already-frozen
Cartesian nonlinear pieces

    N_i->o = (u_inner . grad) u_osc,
    N_o->i = (u_osc . grad) u_inner,
    N_o->o = (u_osc . grad) u_osc,

with non-cancelling vector-RMS magnitudes and normalized magnitude shares.  The
shares are descriptive routing evidence only: they are not NS residual shares,
not m=0 mean-defect shares, and not acceptance gates.

Kokuno's corrected 2026-09-09 reconstruction remains structural provenance for
the upstream localized complete-curl construction.  The public-z realization,
finite-difference derivative machinery, this attribution protocol, and all
reported numerical shares are repository realizations, not paper-exact data.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_oscillatory_nonlinear_handoff import (
    KokunoOscillatoryNonlinearHandoff,
    StrictInnerTransportBackend,
)

TASK = "KOKUNO-A2-OSCILLATORY-NONLINEAR-ATTRIBUTION-072"
SCHEMA = "kokuno-a2-oscillatory-nonlinear-attribution-v1"
PARENT_AGENT2_PR = 908
PARENT_AGENT2_HEAD = "b88f97cb14ea9849b7bbd7cf9ff9f410161e8c2a"
PARENT_NONLINEAR_HANDOFF_BLOB = "e8d5e2d590ae839296a76268a8c787e609bc7fb1"


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _vector_rms(value: np.ndarray) -> float:
    value = np.asarray(value, dtype=float)
    if value.ndim < 1 or value.shape[-1] != 3:
        raise ValueError("nonlinear attribution expects trailing Cartesian component axis of length 3")
    if value.size == 0:
        raise ValueError("nonlinear attribution requires at least one sample")
    if not np.all(np.isfinite(value)):
        raise RuntimeError("nonlinear attribution received non-finite values")
    return float(np.sqrt(np.mean(np.sum(value * value, axis=-1))))


@dataclass(frozen=True)
class OscillatoryNonlinearAttributionResult:
    inner_advects_oscillation_rms: float
    oscillation_advects_inner_rms: float
    oscillatory_self_advection_rms: float
    mixed_cross_advection_rms: float
    aggregate_nonlinear_increment_rms: float
    inner_advects_oscillation_share: float
    oscillation_advects_inner_share: float
    oscillatory_self_advection_share: float
    mixed_magnitude_share: float
    magnitude_denominator: float
    decomposition_closure_abs_max: float
    mixed_closure_abs_max: float

    def to_payload(self) -> dict[str, float]:
        return {
            "inner_advects_oscillation_rms": self.inner_advects_oscillation_rms,
            "oscillation_advects_inner_rms": self.oscillation_advects_inner_rms,
            "oscillatory_self_advection_rms": self.oscillatory_self_advection_rms,
            "mixed_cross_advection_rms": self.mixed_cross_advection_rms,
            "aggregate_nonlinear_increment_rms": self.aggregate_nonlinear_increment_rms,
            "inner_advects_oscillation_share": self.inner_advects_oscillation_share,
            "oscillation_advects_inner_share": self.oscillation_advects_inner_share,
            "oscillatory_self_advection_share": self.oscillatory_self_advection_share,
            "mixed_magnitude_share": self.mixed_magnitude_share,
            "magnitude_denominator": self.magnitude_denominator,
            "decomposition_closure_abs_max": self.decomposition_closure_abs_max,
            "mixed_closure_abs_max": self.mixed_closure_abs_max,
        }


class KokunoOscillatoryNonlinearAttribution:
    """Summarize raw Cartesian A2 nonlinear pieces without entering A3's lane."""

    def __init__(self, inner_leading_backend: StrictInnerTransportBackend):
        self._parent = KokunoOscillatoryNonlinearHandoff(inner_leading_backend)

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
            "parent_nonlinear_handoff_blob": PARENT_NONLINEAR_HANDOFF_BLOB,
            "parent_handoff_sha256": self.parent_handoff_sha256,
            "backend_field_sha256": self.backend_field_sha256,
            "source_contract_sha256": self.source_contract_sha256,
            "metric_contract": {
                "vector_rms": "sqrt(mean(sum(v_i^2 over Cartesian components)))",
                "magnitude_denominator": "rms(i_to_o)+rms(o_to_i)+rms(o_to_o)",
                "piece_shares": "piece_rms/magnitude_denominator",
                "mixed_magnitude_share": "(rms(i_to_o)+rms(o_to_i))/magnitude_denominator",
                "shares_are_non_cancelling_magnitude_diagnostics": True,
                "shares_are_ns_residual_fractions": False,
                "shares_are_m0_mean_defect_fractions": False,
                "shares_are_acceptance_gates": False,
            },
            "scientific_boundary": {
                "strict_inner_only": True,
                "frame": "Cartesian",
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
    def attribution_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> OscillatoryNonlinearAttributionResult:
        parent = self._parent.evaluate(x, y, z, t)
        a = np.asarray(parent.inner_advects_oscillation, dtype=float)
        b = np.asarray(parent.oscillation_advects_inner, dtype=float)
        q = np.asarray(parent.oscillatory_self_advection, dtype=float)
        mixed = np.asarray(parent.mixed_cross_advection, dtype=float)
        aggregate = np.asarray(parent.oscillatory_nonlinear_increment, dtype=float)
        arrays = (a, b, q, mixed, aggregate)
        if len({value.shape for value in arrays}) != 1:
            raise RuntimeError("nonlinear attribution shape mismatch")

        rms_a = _vector_rms(a)
        rms_b = _vector_rms(b)
        rms_q = _vector_rms(q)
        rms_mixed = _vector_rms(mixed)
        rms_aggregate = _vector_rms(aggregate)
        denominator = rms_a + rms_b + rms_q
        if not np.isfinite(denominator) or denominator <= 0.0:
            raise RuntimeError("nonlinear attribution magnitude denominator is degenerate")

        mixed_closure = float(np.max(np.abs(mixed - (a + b))))
        decomposition = float(np.max(np.abs(aggregate - (a + b + q))))
        if mixed_closure > 5.0e-13:
            raise RuntimeError("mixed nonlinear oscillatory terms failed closure")
        if decomposition > 5.0e-13:
            raise RuntimeError("nonlinear oscillatory attribution failed decomposition closure")

        return OscillatoryNonlinearAttributionResult(
            inner_advects_oscillation_rms=rms_a,
            oscillation_advects_inner_rms=rms_b,
            oscillatory_self_advection_rms=rms_q,
            mixed_cross_advection_rms=rms_mixed,
            aggregate_nonlinear_increment_rms=rms_aggregate,
            inner_advects_oscillation_share=rms_a / denominator,
            oscillation_advects_inner_share=rms_b / denominator,
            oscillatory_self_advection_share=rms_q / denominator,
            mixed_magnitude_share=(rms_a + rms_b) / denominator,
            magnitude_denominator=denominator,
            decomposition_closure_abs_max=decomposition,
            mixed_closure_abs_max=mixed_closure,
        )

    def manifest(self) -> dict[str, Any]:
        payload = self.semantic_payload()
        return {"payload": payload, "attribution_sha256": _sha256(payload)}

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
    ) -> "KokunoOscillatoryNonlinearAttribution":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"payload", "attribution_sha256"}:
            raise ValueError("invalid nonlinear attribution manifest envelope")
        payload = raw["payload"]
        if not isinstance(payload, dict) or _sha256(payload) != raw["attribution_sha256"]:
            raise ValueError("nonlinear attribution manifest checksum mismatch")
        obj = cls(inner_leading_backend)
        if obj.semantic_payload() != payload:
            raise ValueError("nonlinear attribution manifest semantics mismatch")
        return obj


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(KokunoOscillatoryNonlinearAttribution.evaluate)
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
        "consumer": "Kokuno Agent 3",
        "mean_projection_performed": False,
        "correction_velocity_constructed": False,
        "shares_are_acceptance_gates": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }
