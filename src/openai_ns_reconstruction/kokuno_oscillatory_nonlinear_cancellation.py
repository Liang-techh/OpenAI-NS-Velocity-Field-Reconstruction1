"""Cancellation/alignment diagnostics for frozen strict-inner oscillatory nonlinear pieces.

This module adds no oscillatory degrees of freedom.  It consumes the already-
frozen Agent-2 decomposition

    A = (u_inner . grad) u_osc,
    B = (u_osc . grad) u_inner,
    Q = (u_osc . grad) u_osc,

and reports raw-Cartesian L2 alignment/cancellation diagnostics.  These quantities
help downstream Agent 3 distinguish weak pieces from strong pieces that nearly
cancel before/after cylindrical projection.  They are descriptive routing
evidence only: not m=0 mean-defect fractions, not NS residual fractions, and not
scientific acceptance gates.

Kokuno's corrected 2026-09-09 reconstruction remains structural provenance for
the upstream localized complete-curl organization.  The public-z realization,
numerical derivative seams, this diagnostic protocol, and all reported numbers
are repository realizations and are not paper-exact hidden data.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_oscillatory_nonlinear_attribution import KokunoOscillatoryNonlinearAttribution
from .kokuno_oscillatory_nonlinear_handoff import (
    KokunoOscillatoryNonlinearHandoff,
    StrictInnerTransportBackend,
)

TASK = "KOKUNO-A2-OSCILLATORY-NONLINEAR-CANCELLATION-073"
SCHEMA = "kokuno-a2-oscillatory-nonlinear-cancellation-v1"
PARENT_AGENT2_PR = 914
PARENT_AGENT2_HEAD = "b8159807d53c362bbdd691b56bdcd4c03456cff0"
PARENT_ATTRIBUTION_BLOB = "ff74ef1d12da5b5e5bb4291e949d0150c1b59148"


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _as_vectors(value: Any) -> np.ndarray:
    value = np.asarray(value, dtype=float)
    if value.ndim < 1 or value.shape[-1] != 3:
        raise ValueError("nonlinear cancellation diagnostic expects trailing Cartesian axis of length 3")
    if value.size == 0:
        raise ValueError("nonlinear cancellation diagnostic requires at least one sample")
    if not np.all(np.isfinite(value)):
        raise RuntimeError("nonlinear cancellation diagnostic received non-finite values")
    return value


def _inner(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.mean(np.sum(a * b, axis=-1)))


def _norm(value: np.ndarray) -> float:
    return float(np.sqrt(max(_inner(value, value), 0.0)))


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    na = _norm(a)
    nb = _norm(b)
    if na <= 0.0 or nb <= 0.0:
        raise RuntimeError("directional alignment is undefined for a zero nonlinear piece")
    value = _inner(a, b) / (na * nb)
    if not np.isfinite(value):
        raise RuntimeError("nonlinear alignment is non-finite")
    # Roundoff can exceed the mathematical interval by a few ulps.
    return float(np.clip(value, -1.0, 1.0))


@dataclass(frozen=True)
class OscillatoryNonlinearCancellationResult:
    inner_to_osc_vs_osc_to_inner_cosine: float
    inner_to_osc_vs_self_cosine: float
    osc_to_inner_vs_self_cosine: float
    mixed_vs_self_cosine: float
    mixed_cancellation_ratio: float
    aggregate_cancellation_ratio: float
    gram_reconstruction_relative_error: float
    parent_attribution_replay_relative_error: float

    def to_payload(self) -> dict[str, float]:
        return {
            "inner_to_osc_vs_osc_to_inner_cosine": self.inner_to_osc_vs_osc_to_inner_cosine,
            "inner_to_osc_vs_self_cosine": self.inner_to_osc_vs_self_cosine,
            "osc_to_inner_vs_self_cosine": self.osc_to_inner_vs_self_cosine,
            "mixed_vs_self_cosine": self.mixed_vs_self_cosine,
            "mixed_cancellation_ratio": self.mixed_cancellation_ratio,
            "aggregate_cancellation_ratio": self.aggregate_cancellation_ratio,
            "gram_reconstruction_relative_error": self.gram_reconstruction_relative_error,
            "parent_attribution_replay_relative_error": self.parent_attribution_replay_relative_error,
        }


class KokunoOscillatoryNonlinearCancellation:
    """Measure raw-Cartesian directional cancellation without entering Agent 3's lane."""

    def __init__(self, inner_leading_backend: StrictInnerTransportBackend):
        self._handoff = KokunoOscillatoryNonlinearHandoff(inner_leading_backend)
        self._attribution = KokunoOscillatoryNonlinearAttribution(inner_leading_backend)
        if self._handoff.handoff_sha256 != self._attribution.parent_handoff_sha256:
            raise RuntimeError("nonlinear handoff/attribution identity mismatch")

    @property
    def backend_field_sha256(self) -> str:
        return self._handoff.backend_field_sha256

    @property
    def source_contract_sha256(self) -> str:
        return self._handoff.source_contract_sha256

    @property
    def parent_handoff_sha256(self) -> str:
        return self._handoff.handoff_sha256

    @property
    def parent_attribution_sha256(self) -> str:
        return self._attribution.attribution_sha256

    def semantic_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "parent_attribution_blob": PARENT_ATTRIBUTION_BLOB,
            "parent_attribution_sha256": self.parent_attribution_sha256,
            "parent_handoff_sha256": self.parent_handoff_sha256,
            "backend_field_sha256": self.backend_field_sha256,
            "source_contract_sha256": self.source_contract_sha256,
            "metric_contract": {
                "inner_product": "mean(sum(a_i*b_i over Cartesian components))",
                "norm": "sqrt(inner_product(v,v))",
                "pairwise_cosines": "inner(a,b)/(norm(a)*norm(b))",
                "mixed_cancellation_ratio": "norm(A+B)/(norm(A)+norm(B))",
                "aggregate_cancellation_ratio": "norm(A+B+Q)/(norm(A)+norm(B)+norm(Q))",
                "gram_identity": "norm(A+B+Q)^2 from three norms plus pairwise inner products",
                "ratios_are_raw_cartesian_cancellation_diagnostics": True,
                "ratios_are_m0_mean_defect_fractions": False,
                "ratios_are_ns_residual_fractions": False,
                "ratios_are_acceptance_gates": False,
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
    def cancellation_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> OscillatoryNonlinearCancellationResult:
        raw = self._handoff.evaluate(x, y, z, t)
        parent = self._attribution.evaluate(x, y, z, t)

        a = _as_vectors(raw.inner_advects_oscillation)
        b = _as_vectors(raw.oscillation_advects_inner)
        q = _as_vectors(raw.oscillatory_self_advection)
        mixed = _as_vectors(raw.mixed_cross_advection)
        aggregate = _as_vectors(raw.oscillatory_nonlinear_increment)
        if len({value.shape for value in (a, b, q, mixed, aggregate)}) != 1:
            raise RuntimeError("nonlinear cancellation shape mismatch")

        na, nb, nq = _norm(a), _norm(b), _norm(q)
        nmixed, nagg = _norm(mixed), _norm(aggregate)
        if min(na, nb, nq) <= 0.0:
            raise RuntimeError("nonlinear cancellation diagnostic requires nonzero constituent pieces")

        mixed_denominator = na + nb
        aggregate_denominator = na + nb + nq
        mixed_ratio = nmixed / mixed_denominator
        aggregate_ratio = nagg / aggregate_denominator
        if not (np.isfinite(mixed_ratio) and np.isfinite(aggregate_ratio)):
            raise RuntimeError("nonlinear cancellation ratio is non-finite")
        if mixed_ratio < -1.0e-12 or mixed_ratio > 1.0 + 1.0e-12:
            raise RuntimeError("mixed cancellation ratio violates the triangle inequality")
        if aggregate_ratio < -1.0e-12 or aggregate_ratio > 1.0 + 1.0e-12:
            raise RuntimeError("aggregate cancellation ratio violates the triangle inequality")

        ip_ab = _inner(a, b)
        ip_aq = _inner(a, q)
        ip_bq = _inner(b, q)
        gram_sq = na * na + nb * nb + nq * nq + 2.0 * (ip_ab + ip_aq + ip_bq)
        gram_sq = max(gram_sq, 0.0)
        gram_norm = float(np.sqrt(gram_sq))
        gram_error = abs(gram_norm - nagg) / max(nagg, aggregate_denominator, 1.0e-30)

        parent_values = np.asarray(
            [
                parent.inner_advects_oscillation_rms,
                parent.oscillation_advects_inner_rms,
                parent.oscillatory_self_advection_rms,
                parent.mixed_cross_advection_rms,
                parent.aggregate_nonlinear_increment_rms,
            ],
            dtype=float,
        )
        replay_values = np.asarray([na, nb, nq, nmixed, nagg], dtype=float)
        parent_replay = float(
            np.max(np.abs(parent_values - replay_values) / np.maximum(np.abs(parent_values), 1.0e-30))
        )

        return OscillatoryNonlinearCancellationResult(
            inner_to_osc_vs_osc_to_inner_cosine=_cosine(a, b),
            inner_to_osc_vs_self_cosine=_cosine(a, q),
            osc_to_inner_vs_self_cosine=_cosine(b, q),
            mixed_vs_self_cosine=_cosine(mixed, q),
            mixed_cancellation_ratio=float(np.clip(mixed_ratio, 0.0, 1.0)),
            aggregate_cancellation_ratio=float(np.clip(aggregate_ratio, 0.0, 1.0)),
            gram_reconstruction_relative_error=gram_error,
            parent_attribution_replay_relative_error=parent_replay,
        )

    def manifest(self) -> dict[str, Any]:
        payload = self.semantic_payload()
        return {"payload": payload, "cancellation_sha256": _sha256(payload)}

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
    ) -> "KokunoOscillatoryNonlinearCancellation":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"payload", "cancellation_sha256"}:
            raise ValueError("invalid nonlinear cancellation manifest envelope")
        payload = raw["payload"]
        if not isinstance(payload, dict) or _sha256(payload) != raw["cancellation_sha256"]:
            raise ValueError("nonlinear cancellation manifest checksum mismatch")
        obj = cls(inner_leading_backend)
        if obj.semantic_payload() != payload:
            raise ValueError("nonlinear cancellation manifest semantics mismatch")
        return obj


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(KokunoOscillatoryNonlinearCancellation.evaluate)
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
        "cancellation_metrics_are_acceptance_gates": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }
