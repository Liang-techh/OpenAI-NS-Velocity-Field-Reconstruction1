"""Agent-3 cylindrical m=0 projection of the typed A2 oscillatory transport handoff.

The exact Agent-2 #901 handoff exposes the strict-inner Cartesian increment

    delta T_osc = d_t u_osc
                  + (u_inner . grad) u_osc
                  + (u_osc . grad) u_inner
                  + (u_osc . grad) u_osc
                  - nu Delta u_osc,              nu = 0.01.

This module performs the Agent-3-owned rotating cylindrical mean projection of
that increment and its time/nonlinear/viscous decomposition.  It additionally
checks the independent upstream identity

    delta T_osc = T_(inner+osc) - T_inner.

It does not add pressure, forcing, a global leading join, a radial inverse, or a
correction velocity.  Consequently this witness is not a complete NS defect and
is not authorized as a finite-cycle correction target.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_oscillatory_correction_cross_mean import (
    _project_cartesian_ring_to_cylindrical_mean,
)
from .kokuno_oscillatory_quadratic_mean import (
    ANGULAR_ORDERS,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    _git_blob_sha,
    _relative_rms_difference,
    _validate_rings,
)

TASK = "KOKUNO-A3-OSCILLATORY-TRANSPORT-MEAN-084"
SCHEMA = "kokuno-a3-oscillatory-transport-mean-v1"
PARENT_AGENT3_PR = 895
PARENT_AGENT3_HEAD = "d305b64edcdde4f16dd8a48908f53fa23c2b9baf"

AGENT2_PR = 901
AGENT2_HEAD = "4a78b06afbab95c108c184a32cb5951e66f42991"
AGENT2_MODULE = "openai_ns_reconstruction.kokuno_oscillatory_transport_handoff"
AGENT2_CLASS = "KokunoOscillatoryTransportHandoff"
AGENT2_SOURCE_BLOB_SHA = "57932afda3b65e4cb63e9ae7c8ceeb0aba141cc5"
AGENT2_PARENT_TRANSPORT_PR = 849
AGENT2_PARENT_TRANSPORT_HEAD = "616e61abd2420a7c7aff3590a08393360039cd1a"

REPOSITORY_VISCOSITY = 0.01
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

OscillatoryTransportProvider = Callable[
    [np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float],
]


@dataclass(frozen=True)
class ExactAgent2OscillatoryTransportHandoff:
    """Bind the exact checksum-aware Agent-2 #901 handoff implementation."""

    handoff: object
    source_file: str
    source_blob_sha: str
    handoff_sha256: str
    source_contract_sha256: str

    @classmethod
    def bind(cls, handoff: object) -> "ExactAgent2OscillatoryTransportHandoff":
        handoff_type = type(handoff)
        if handoff_type.__module__ != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 oscillatory handoff module identity")
        if handoff_type.__name__ != AGENT2_CLASS:
            raise ValueError("unexpected Agent-2 oscillatory handoff class identity")
        if not callable(getattr(handoff, "evaluate", None)):
            raise TypeError("Agent-2 oscillatory handoff must expose evaluate(x,y,z,t)")
        if not callable(getattr(handoff, "semantic_payload", None)):
            raise TypeError("Agent-2 oscillatory handoff must expose semantic_payload()")
        source = inspect.getsourcefile(handoff_type)
        if source is None:
            raise ValueError("Agent-2 oscillatory handoff source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 oscillatory handoff source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError("Agent-2 handoff source blob does not match pinned PR #901")

        payload = handoff.semantic_payload()
        quantity = payload.get("quantity_contract", {})
        boundary = payload.get("scientific_boundary", {})
        if payload.get("parent_agent2_head") != "a467bc175a801ed69c197fc9e38ca68ac8c0b86b":
            raise ValueError("Agent-2 parent identity drifted")
        if payload.get("transport_delta_head") != AGENT2_PARENT_TRANSPORT_HEAD:
            raise ValueError("Agent-2 transport-delta identity drifted")
        if quantity.get("frame") != "Cartesian":
            raise ValueError("Agent-2 handoff frame drifted")
        if quantity.get("complete_ns_residual") is not False:
            raise ValueError("Agent-2 handoff improperly promotes a complete NS residual")
        if quantity.get("mean_projection_performed") is not False:
            raise ValueError("Agent-2 handoff crossed the Agent-3 mean ownership boundary")
        if boundary.get("strict_inner_only") is not True:
            raise ValueError("Agent-2 handoff strict-inner scope drifted")
        if boundary.get("matched_pressure_included") is not False:
            raise ValueError("Agent-2 handoff unexpectedly includes pressure")
        if boundary.get("restricted_forcing_included") is not False:
            raise ValueError("Agent-2 handoff unexpectedly includes forcing")

        handoff_sha = getattr(handoff, "handoff_sha256", None)
        source_contract_sha = getattr(handoff, "source_contract_sha256", None)
        for label, value in (
            ("handoff_sha256", handoff_sha),
            ("source_contract_sha256", source_contract_sha),
        ):
            if not isinstance(value, str) or len(value) != 64 or any(
                c not in "0123456789abcdef" for c in value
            ):
                raise ValueError(f"Agent-2 {label} must be lowercase 64-hex")
        return cls(
            handoff=handoff,
            source_file=str(path),
            source_blob_sha=blob_sha,
            handoff_sha256=str(handoff_sha),
            source_contract_sha256=str(source_contract_sha),
        )

    def components(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
        result = self.handoff.evaluate(x, y, z, t)
        return (
            np.asarray(getattr(result, "oscillatory_time_increment"), dtype=float),
            np.asarray(getattr(result, "oscillatory_nonlinear_increment"), dtype=float),
            np.asarray(getattr(result, "oscillatory_viscous_increment"), dtype=float),
            np.asarray(getattr(result, "oscillatory_transport_increment"), dtype=float),
            np.asarray(getattr(result, "inner_leading_transport"), dtype=float),
            np.asarray(getattr(result, "inner_plus_oscillatory_transport"), dtype=float),
            float(getattr(result, "viscosity")),
        )

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent2_pr": AGENT2_PR,
            "agent2_expected_head": AGENT2_HEAD,
            "module": AGENT2_MODULE,
            "class": AGENT2_CLASS,
            "source_blob_sha": self.source_blob_sha,
            "source_blob_matches_pinned": self.source_blob_sha == AGENT2_SOURCE_BLOB_SHA,
            "handoff_sha256": self.handoff_sha256,
            "source_contract_sha256": self.source_contract_sha256,
        }


@dataclass(frozen=True)
class OscillatoryTransportMeanWitness:
    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_time_increment_cylindrical: np.ndarray
    mean_nonlinear_increment_cylindrical: np.ndarray
    mean_viscous_increment_cylindrical: np.ndarray
    mean_transport_increment_cylindrical: np.ndarray
    mean_inner_transport_cylindrical: np.ndarray
    mean_inner_plus_oscillatory_transport_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_transport_mean_relative_differences: tuple[float, ...]
    pointwise_decomposition_closure_absolute_max: float
    pointwise_before_after_closure_absolute_max: float
    projected_decomposition_closure_absolute_max: float
    projected_before_after_closure_absolute_max: float
    transport_increment_mean_rms: float
    full_ring_transport_increment_rms: float
    transport_increment_mean_to_full_rms_ratio: float
    viscosity: float
    backend: ExactAgent2OscillatoryTransportHandoff | None
    backend_kind: str

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "radius": np.asarray(self.radius, dtype=float).tolist(),
            "z": np.asarray(self.z, dtype=float).tolist(),
            "t": np.asarray(self.t, dtype=float).tolist(),
            "cylindrical_components": ["radial", "tangential", "axial"],
            "mean_time_increment_cylindrical": np.asarray(self.mean_time_increment_cylindrical).tolist(),
            "mean_nonlinear_increment_cylindrical": np.asarray(self.mean_nonlinear_increment_cylindrical).tolist(),
            "mean_viscous_increment_cylindrical": np.asarray(self.mean_viscous_increment_cylindrical).tolist(),
            "mean_transport_increment_cylindrical": np.asarray(self.mean_transport_increment_cylindrical).tolist(),
            "mean_inner_transport_cylindrical": np.asarray(self.mean_inner_transport_cylindrical).tolist(),
            "mean_inner_plus_oscillatory_transport_cylindrical": np.asarray(
                self.mean_inner_plus_oscillatory_transport_cylindrical
            ).tolist(),
            "angular_orders": list(self.angular_orders),
            "successive_transport_mean_relative_differences": list(
                self.successive_transport_mean_relative_differences
            ),
            "pointwise_decomposition_closure_absolute_max": self.pointwise_decomposition_closure_absolute_max,
            "pointwise_before_after_closure_absolute_max": self.pointwise_before_after_closure_absolute_max,
            "projected_decomposition_closure_absolute_max": self.projected_decomposition_closure_absolute_max,
            "projected_before_after_closure_absolute_max": self.projected_before_after_closure_absolute_max,
            "transport_increment_mean_rms": self.transport_increment_mean_rms,
            "full_ring_transport_increment_rms": self.full_ring_transport_increment_rms,
            "transport_increment_mean_to_full_rms_ratio": self.transport_increment_mean_to_full_rms_ratio,
            "viscosity": self.viscosity,
            "backend_kind": self.backend_kind,
            "agent2_backend": None if self.backend is None else self.backend.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _validate_vector(value: Any, expected_shape: tuple[int, ...], label: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != expected_shape:
        raise ValueError(f"{label} returned shape {array.shape}, expected {expected_shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{label} returned non-finite values")
    return array


def _mean_for_order(
    provider: OscillatoryTransportProvider,
    radius: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, float, float]:
    theta = 2.0 * np.pi * np.arange(order, dtype=float) / float(order)
    c = np.cos(theta)
    s = np.sin(theta)
    xq = radius[..., None] * c
    yq = radius[..., None] * s
    zq = z[..., None] + np.zeros_like(theta)
    tq = t[..., None] + np.zeros_like(theta)
    expected = radius.shape + (order, 3)

    raw = provider(xq, yq, zq, tq)
    if not isinstance(raw, tuple) or len(raw) != 7:
        raise TypeError("oscillatory transport provider must return seven entries")
    time_part = _validate_vector(raw[0], expected, "time increment")
    nonlinear = _validate_vector(raw[1], expected, "nonlinear increment")
    viscous = _validate_vector(raw[2], expected, "viscous increment")
    delta = _validate_vector(raw[3], expected, "transport increment")
    inner = _validate_vector(raw[4], expected, "inner transport")
    combined = _validate_vector(raw[5], expected, "inner-plus-oscillatory transport")
    viscosity = float(raw[6])
    if not np.isfinite(viscosity) or viscosity != REPOSITORY_VISCOSITY:
        raise ValueError("oscillatory handoff viscosity must remain exactly 0.01")

    decomposition_closure = float(np.max(np.abs(delta - time_part - nonlinear - viscous)))
    before_after_closure = float(np.max(np.abs(delta - (combined - inner))))

    mean_time = _project_cartesian_ring_to_cylindrical_mean(time_part, c, s)
    mean_nonlinear = _project_cartesian_ring_to_cylindrical_mean(nonlinear, c, s)
    mean_viscous = _project_cartesian_ring_to_cylindrical_mean(viscous, c, s)
    mean_delta = _project_cartesian_ring_to_cylindrical_mean(delta, c, s)
    mean_inner = _project_cartesian_ring_to_cylindrical_mean(inner, c, s)
    mean_combined = _project_cartesian_ring_to_cylindrical_mean(combined, c, s)
    projected_decomposition = float(
        np.max(np.abs(mean_delta - mean_time - mean_nonlinear - mean_viscous))
    )
    projected_before_after = float(
        np.max(np.abs(mean_delta - (mean_combined - mean_inner)))
    )
    full_ring_rms = float(np.sqrt(np.mean(np.sum(delta * delta, axis=-1))))
    return (
        mean_time,
        mean_nonlinear,
        mean_viscous,
        mean_delta,
        mean_inner,
        mean_combined,
        decomposition_closure,
        before_after_closure,
        projected_decomposition,
        projected_before_after,
        full_ring_rms,
        viscosity,
    )


def _materialize_from_provider(
    provider: OscillatoryTransportProvider,
    radius: Any,
    z: Any,
    t: Any,
    *,
    backend: ExactAgent2OscillatoryTransportHandoff | None,
    backend_kind: str,
) -> OscillatoryTransportMeanWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    per_order = [_mean_for_order(provider, rb, zb, tb, order) for order in ANGULAR_ORDERS]
    totals = [entry[3] for entry in per_order]
    successive = tuple(
        _relative_rms_difference(totals[i + 1], totals[i])
        for i in range(len(totals) - 1)
    )
    final = per_order[-1]
    mean_delta = final[3]
    mean_rms = float(np.sqrt(np.mean(np.sum(mean_delta * mean_delta, axis=-1))))
    full_rms = float(final[10])
    ratio = mean_rms / max(full_rms, np.finfo(float).tiny)
    return OscillatoryTransportMeanWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_time_increment_cylindrical=np.array(final[0], copy=True),
        mean_nonlinear_increment_cylindrical=np.array(final[1], copy=True),
        mean_viscous_increment_cylindrical=np.array(final[2], copy=True),
        mean_transport_increment_cylindrical=np.array(final[3], copy=True),
        mean_inner_transport_cylindrical=np.array(final[4], copy=True),
        mean_inner_plus_oscillatory_transport_cylindrical=np.array(final[5], copy=True),
        angular_orders=ANGULAR_ORDERS,
        successive_transport_mean_relative_differences=successive,
        pointwise_decomposition_closure_absolute_max=float(max(v[6] for v in per_order)),
        pointwise_before_after_closure_absolute_max=float(max(v[7] for v in per_order)),
        projected_decomposition_closure_absolute_max=float(max(v[8] for v in per_order)),
        projected_before_after_closure_absolute_max=float(max(v[9] for v in per_order)),
        transport_increment_mean_rms=mean_rms,
        full_ring_transport_increment_rms=full_rms,
        transport_increment_mean_to_full_rms_ratio=float(ratio),
        viscosity=float(final[11]),
        backend=backend,
        backend_kind=backend_kind,
    )


def materialize_oscillatory_transport_mean(
    handoff_backend: ExactAgent2OscillatoryTransportHandoff,
    radius: Any,
    z: Any,
    t: Any,
) -> OscillatoryTransportMeanWitness:
    """Project exact A2 #901 raw Cartesian transport increment onto cylindrical m=0."""
    if not isinstance(handoff_backend, ExactAgent2OscillatoryTransportHandoff):
        raise TypeError("handoff_backend must be ExactAgent2OscillatoryTransportHandoff")

    def provider(x: np.ndarray, y: np.ndarray, zz: np.ndarray, tt: np.ndarray):
        return handoff_backend.components(x, y, zz, tt)

    return _materialize_from_provider(
        provider,
        radius,
        z,
        t,
        backend=handoff_backend,
        backend_kind="exact-agent2-901-oscillatory-transport-handoff",
    )


def truth_boundary() -> dict[str, object]:
    return {
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "repository_viscosity": REPOSITORY_VISCOSITY,
        "agent2_handoff_head": AGENT2_HEAD,
        "agent2_parent_transport_head": AGENT2_PARENT_TRANSPORT_HEAD,
        "strict_inner_only": True,
        "oscillatory_transport_increment_materialized": True,
        "cylindrical_m0_projection_materialized": True,
        "time_increment_included": True,
        "nonlinear_increment_included": True,
        "viscous_increment_included": True,
        "before_after_transport_difference_checked": True,
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "global_corrected_leading_join_materialized": False,
        "complete_ns_defect": False,
        "radial_inverse_performed": False,
        "scoped_transport_mean_authorized_as_correction_target": False,
        "mean_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_mean_allowed": False,
        "caller_supplied_pressure_allowed": False,
        "caller_supplied_forcing_allowed": False,
        "caller_supplied_viscosity_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
