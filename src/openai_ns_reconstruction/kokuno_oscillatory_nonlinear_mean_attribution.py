"""Attribute the current typed oscillatory nonlinear transport at cylindrical m=0.

Agent 2 PR #908 exposes the #901-bound strict-inner nonlinear increment as

    N_osc =
        (u_inner . grad) u_osc
      + (u_osc . grad) u_inner
      + (u_osc . grad) u_osc.

This module performs only the Agent-3-owned rotating cylindrical m=0 projection.
It keeps the two leading/oscillation mixed pieces separate from the genuinely
quadratic oscillatory self-interaction and checks that their projected sum
replays the aggregate nonlinear mean.

The result is a strict-inner nonlinear attribution witness.  It contains no
pressure, forcing, global leading join, radial inverse, correction velocity or
complete Navier--Stokes residual.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import inspect
import json
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

TASK = "KOKUNO-A3-OSCILLATORY-NONLINEAR-MEAN-ATTRIBUTION-086"
SCHEMA = "kokuno-a3-oscillatory-nonlinear-mean-attribution-v1"
PARENT_AGENT3_PR = 909
PARENT_AGENT3_HEAD = "cd53700a4644b1c10cb008fdca068d52e99cdcb2"

AGENT2_PR = 908
AGENT2_HEAD = "b88f97cb14ea9849b7bbd7cf9ff9f410161e8c2a"
AGENT2_MODULE = "openai_ns_reconstruction.kokuno_oscillatory_nonlinear_handoff"
AGENT2_CLASS = "KokunoOscillatoryNonlinearHandoff"
AGENT2_SOURCE_BLOB_SHA = "e8d5e2d590ae839296a76268a8c787e609bc7fb1"
AGENT2_PARENT_PR = 901
AGENT2_PARENT_HEAD = "4a78b06afbab95c108c184a32cb5951e66f42991"
AGENT2_PARENT_HANDOFF_BLOB = "57932afda3b65e4cb63e9ae7c8ceeb0aba141cc5"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

POINTWISE_CLOSURE_GATE = 5.0e-13
PROJECTED_CLOSURE_GATE = 1.0e-12
ANGULAR_CONVERGENCE_GATE = 1.0e-5
NONZERO_MEAN_FLOOR = 1.0e-12

NonlinearProvider = Callable[
    [np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
]


@dataclass(frozen=True)
class ExactAgent2OscillatoryNonlinearHandoff:
    """Bind the exact checksum-aware Agent-2 #908 nonlinear handoff."""

    handoff: object
    source_file: str
    source_blob_sha: str
    handoff_sha256: str
    parent_handoff_sha256: str
    source_contract_sha256: str

    @classmethod
    def bind(cls, handoff: object) -> "ExactAgent2OscillatoryNonlinearHandoff":
        handoff_type = type(handoff)
        if handoff_type.__module__ != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 nonlinear handoff module identity")
        if handoff_type.__name__ != AGENT2_CLASS:
            raise ValueError("unexpected Agent-2 nonlinear handoff class identity")
        if not callable(getattr(handoff, "evaluate", None)):
            raise TypeError("Agent-2 nonlinear handoff must expose evaluate(x,y,z,t)")
        if not callable(getattr(handoff, "semantic_payload", None)):
            raise TypeError("Agent-2 nonlinear handoff must expose semantic_payload()")

        source = inspect.getsourcefile(handoff_type)
        if source is None:
            raise ValueError("Agent-2 nonlinear handoff source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 nonlinear handoff source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError("Agent-2 nonlinear handoff source blob does not match pinned PR #908")

        payload = handoff.semantic_payload()
        if payload.get("parent_agent2_pr") != AGENT2_PARENT_PR:
            raise ValueError("Agent-2 nonlinear handoff parent PR drifted")
        if payload.get("parent_agent2_head") != AGENT2_PARENT_HEAD:
            raise ValueError("Agent-2 nonlinear handoff parent head drifted")
        if payload.get("parent_handoff_blob") != AGENT2_PARENT_HANDOFF_BLOB:
            raise ValueError("Agent-2 nonlinear handoff parent source blob drifted")

        quantity = payload.get("quantity_contract", {})
        scientific = payload.get("scientific_boundary", {})
        if quantity.get("frame") != "Cartesian":
            raise ValueError("Agent-2 nonlinear handoff frame drifted")
        if quantity.get("mean_projection_performed") is not False:
            raise ValueError("Agent-2 nonlinear handoff crossed Agent-3 mean ownership")
        if quantity.get("radial_inverse_performed") is not False:
            raise ValueError("Agent-2 nonlinear handoff crossed Agent-3 radial ownership")
        if quantity.get("correction_velocity_constructed") is not False:
            raise ValueError("Agent-2 nonlinear handoff crossed Agent-3 correction ownership")
        if quantity.get("complete_ns_residual") is not False:
            raise ValueError("Agent-2 nonlinear handoff improperly promotes a full residual")
        if scientific.get("strict_inner_only") is not True:
            raise ValueError("Agent-2 nonlinear handoff strict-inner scope drifted")

        values = {
            "handoff_sha256": getattr(handoff, "handoff_sha256", None),
            "parent_handoff_sha256": getattr(handoff, "parent_handoff_sha256", None),
            "source_contract_sha256": getattr(handoff, "source_contract_sha256", None),
        }
        for label, value in values.items():
            if (
                not isinstance(value, str)
                or len(value) != 64
                or any(c not in "0123456789abcdef" for c in value)
            ):
                raise ValueError(f"Agent-2 {label} must be lowercase 64-hex")

        return cls(
            handoff=handoff,
            source_file=str(path),
            source_blob_sha=blob_sha,
            handoff_sha256=str(values["handoff_sha256"]),
            parent_handoff_sha256=str(values["parent_handoff_sha256"]),
            source_contract_sha256=str(values["source_contract_sha256"]),
        )

    def components(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        result = self.handoff.evaluate(x, y, z, t)
        return (
            np.asarray(getattr(result, "inner_advects_oscillation"), dtype=float),
            np.asarray(getattr(result, "oscillation_advects_inner"), dtype=float),
            np.asarray(getattr(result, "oscillatory_self_advection"), dtype=float),
            np.asarray(getattr(result, "mixed_cross_advection"), dtype=float),
            np.asarray(getattr(result, "oscillatory_nonlinear_increment"), dtype=float),
        )

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent2_pr": AGENT2_PR,
            "agent2_expected_head": AGENT2_HEAD,
            "agent2_parent_pr": AGENT2_PARENT_PR,
            "agent2_parent_head": AGENT2_PARENT_HEAD,
            "module": AGENT2_MODULE,
            "class": AGENT2_CLASS,
            "source_blob_sha": self.source_blob_sha,
            "source_blob_matches_pinned": self.source_blob_sha == AGENT2_SOURCE_BLOB_SHA,
            "handoff_sha256": self.handoff_sha256,
            "parent_handoff_sha256": self.parent_handoff_sha256,
            "source_contract_sha256": self.source_contract_sha256,
        }


@dataclass(frozen=True)
class OscillatoryNonlinearMeanAttributionWitness:
    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_inner_advects_oscillation_cylindrical: np.ndarray
    mean_oscillation_advects_inner_cylindrical: np.ndarray
    mean_mixed_cross_cylindrical: np.ndarray
    mean_oscillatory_self_advection_cylindrical: np.ndarray
    mean_aggregate_nonlinear_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_aggregate_mean_relative_differences: tuple[float, ...]
    successive_mixed_mean_relative_differences: tuple[float, ...]
    successive_quadratic_mean_relative_differences: tuple[float, ...]
    pointwise_three_piece_closure_absolute_max: float
    pointwise_mixed_closure_absolute_max: float
    projected_three_piece_closure_absolute_max: float
    projected_mixed_closure_absolute_max: float
    aggregate_mean_rms: float
    mixed_mean_rms: float
    quadratic_mean_rms: float
    full_ring_aggregate_rms: float
    full_ring_mixed_rms: float
    full_ring_quadratic_rms: float
    aggregate_mean_to_full_rms_ratio: float
    mixed_mean_to_full_rms_ratio: float
    quadratic_mean_to_full_rms_ratio: float
    mixed_mean_to_aggregate_mean_rms_ratio: float
    quadratic_mean_to_aggregate_mean_rms_ratio: float
    backend: ExactAgent2OscillatoryNonlinearHandoff | None
    backend_kind: str

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "radius": np.asarray(self.radius, dtype=float).tolist(),
            "z": np.asarray(self.z, dtype=float).tolist(),
            "t": np.asarray(self.t, dtype=float).tolist(),
            "cylindrical_components": ["radial", "tangential", "axial"],
            "mean_inner_advects_oscillation_cylindrical": np.asarray(
                self.mean_inner_advects_oscillation_cylindrical, dtype=float
            ).tolist(),
            "mean_oscillation_advects_inner_cylindrical": np.asarray(
                self.mean_oscillation_advects_inner_cylindrical, dtype=float
            ).tolist(),
            "mean_mixed_cross_cylindrical": np.asarray(
                self.mean_mixed_cross_cylindrical, dtype=float
            ).tolist(),
            "mean_oscillatory_self_advection_cylindrical": np.asarray(
                self.mean_oscillatory_self_advection_cylindrical, dtype=float
            ).tolist(),
            "mean_aggregate_nonlinear_cylindrical": np.asarray(
                self.mean_aggregate_nonlinear_cylindrical, dtype=float
            ).tolist(),
            "angular_orders": list(self.angular_orders),
            "successive_aggregate_mean_relative_differences": list(
                self.successive_aggregate_mean_relative_differences
            ),
            "successive_mixed_mean_relative_differences": list(
                self.successive_mixed_mean_relative_differences
            ),
            "successive_quadratic_mean_relative_differences": list(
                self.successive_quadratic_mean_relative_differences
            ),
            "pointwise_three_piece_closure_absolute_max": self.pointwise_three_piece_closure_absolute_max,
            "pointwise_mixed_closure_absolute_max": self.pointwise_mixed_closure_absolute_max,
            "projected_three_piece_closure_absolute_max": self.projected_three_piece_closure_absolute_max,
            "projected_mixed_closure_absolute_max": self.projected_mixed_closure_absolute_max,
            "aggregate_mean_rms": self.aggregate_mean_rms,
            "mixed_mean_rms": self.mixed_mean_rms,
            "quadratic_mean_rms": self.quadratic_mean_rms,
            "full_ring_aggregate_rms": self.full_ring_aggregate_rms,
            "full_ring_mixed_rms": self.full_ring_mixed_rms,
            "full_ring_quadratic_rms": self.full_ring_quadratic_rms,
            "aggregate_mean_to_full_rms_ratio": self.aggregate_mean_to_full_rms_ratio,
            "mixed_mean_to_full_rms_ratio": self.mixed_mean_to_full_rms_ratio,
            "quadratic_mean_to_full_rms_ratio": self.quadratic_mean_to_full_rms_ratio,
            "mixed_mean_to_aggregate_mean_rms_ratio": self.mixed_mean_to_aggregate_mean_rms_ratio,
            "quadratic_mean_to_aggregate_mean_rms_ratio": self.quadratic_mean_to_aggregate_mean_rms_ratio,
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


def _vector_rms(value: np.ndarray) -> float:
    array = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def _mean_for_order(
    provider: NonlinearProvider,
    radius: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    order: int,
) -> tuple[np.ndarray, ...]:
    theta = 2.0 * np.pi * np.arange(order, dtype=float) / float(order)
    c = np.cos(theta)
    s = np.sin(theta)
    xq = radius[..., None] * c
    yq = radius[..., None] * s
    zq = z[..., None] + np.zeros_like(theta)
    tq = t[..., None] + np.zeros_like(theta)
    expected = radius.shape + (order, 3)

    raw = provider(xq, yq, zq, tq)
    if not isinstance(raw, tuple) or len(raw) != 5:
        raise TypeError("nonlinear provider must return five arrays")
    a = _validate_vector(raw[0], expected, "inner_advects_oscillation")
    b = _validate_vector(raw[1], expected, "oscillation_advects_inner")
    q = _validate_vector(raw[2], expected, "oscillatory_self_advection")
    mixed = _validate_vector(raw[3], expected, "mixed_cross_advection")
    aggregate = _validate_vector(raw[4], expected, "oscillatory_nonlinear_increment")

    pointwise_mixed = float(np.max(np.abs(mixed - (a + b))))
    pointwise_three = float(np.max(np.abs(aggregate - (a + b + q))))
    if pointwise_mixed > POINTWISE_CLOSURE_GATE:
        raise ValueError("pointwise mixed nonlinear closure failed")
    if pointwise_three > POINTWISE_CLOSURE_GATE:
        raise ValueError("pointwise three-piece nonlinear closure failed")

    mean_a = _project_cartesian_ring_to_cylindrical_mean(a, c, s)
    mean_b = _project_cartesian_ring_to_cylindrical_mean(b, c, s)
    mean_q = _project_cartesian_ring_to_cylindrical_mean(q, c, s)
    mean_mixed = _project_cartesian_ring_to_cylindrical_mean(mixed, c, s)
    mean_aggregate = _project_cartesian_ring_to_cylindrical_mean(aggregate, c, s)

    projected_mixed = float(np.max(np.abs(mean_mixed - (mean_a + mean_b))))
    projected_three = float(
        np.max(np.abs(mean_aggregate - (mean_a + mean_b + mean_q)))
    )
    if projected_mixed > PROJECTED_CLOSURE_GATE:
        raise ValueError("projected mixed nonlinear closure failed")
    if projected_three > PROJECTED_CLOSURE_GATE:
        raise ValueError("projected three-piece nonlinear closure failed")

    return (
        mean_a,
        mean_b,
        mean_q,
        mean_mixed,
        mean_aggregate,
        np.sqrt(np.mean(np.sum(q * q, axis=-1), axis=-1)),
        np.sqrt(np.mean(np.sum(mixed * mixed, axis=-1), axis=-1)),
        np.sqrt(np.mean(np.sum(aggregate * aggregate, axis=-1), axis=-1)),
        pointwise_mixed,
        pointwise_three,
        projected_mixed,
        projected_three,
    )


def _materialize_from_provider(
    provider: NonlinearProvider,
    radius: Any,
    z: Any,
    t: Any,
    *,
    backend: ExactAgent2OscillatoryNonlinearHandoff | None,
    backend_kind: str,
) -> OscillatoryNonlinearMeanAttributionWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    per_order = [_mean_for_order(provider, rb, zb, tb, order) for order in ANGULAR_ORDERS]

    aggregates = [entry[4] for entry in per_order]
    mixeds = [entry[3] for entry in per_order]
    quadratics = [entry[2] for entry in per_order]
    aggregate_successive = tuple(
        _relative_rms_difference(aggregates[i + 1], aggregates[i])
        for i in range(len(aggregates) - 1)
    )
    mixed_successive = tuple(
        _relative_rms_difference(mixeds[i + 1], mixeds[i])
        for i in range(len(mixeds) - 1)
    )
    quadratic_successive = tuple(
        _relative_rms_difference(quadratics[i + 1], quadratics[i])
        for i in range(len(quadratics) - 1)
    )

    final = per_order[-1]
    aggregate_mean_rms = _vector_rms(final[4])
    mixed_mean_rms = _vector_rms(final[3])
    quadratic_mean_rms = _vector_rms(final[2])
    full_q = float(np.sqrt(np.mean(np.asarray(final[5], dtype=float) ** 2)))
    full_mixed = float(np.sqrt(np.mean(np.asarray(final[6], dtype=float) ** 2)))
    full_aggregate = float(np.sqrt(np.mean(np.asarray(final[7], dtype=float) ** 2)))
    tiny = np.finfo(float).tiny

    return OscillatoryNonlinearMeanAttributionWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_inner_advects_oscillation_cylindrical=np.array(final[0], copy=True),
        mean_oscillation_advects_inner_cylindrical=np.array(final[1], copy=True),
        mean_mixed_cross_cylindrical=np.array(final[3], copy=True),
        mean_oscillatory_self_advection_cylindrical=np.array(final[2], copy=True),
        mean_aggregate_nonlinear_cylindrical=np.array(final[4], copy=True),
        angular_orders=ANGULAR_ORDERS,
        successive_aggregate_mean_relative_differences=aggregate_successive,
        successive_mixed_mean_relative_differences=mixed_successive,
        successive_quadratic_mean_relative_differences=quadratic_successive,
        pointwise_three_piece_closure_absolute_max=float(max(entry[9] for entry in per_order)),
        pointwise_mixed_closure_absolute_max=float(max(entry[8] for entry in per_order)),
        projected_three_piece_closure_absolute_max=float(max(entry[11] for entry in per_order)),
        projected_mixed_closure_absolute_max=float(max(entry[10] for entry in per_order)),
        aggregate_mean_rms=aggregate_mean_rms,
        mixed_mean_rms=mixed_mean_rms,
        quadratic_mean_rms=quadratic_mean_rms,
        full_ring_aggregate_rms=full_aggregate,
        full_ring_mixed_rms=full_mixed,
        full_ring_quadratic_rms=full_q,
        aggregate_mean_to_full_rms_ratio=aggregate_mean_rms / max(full_aggregate, tiny),
        mixed_mean_to_full_rms_ratio=mixed_mean_rms / max(full_mixed, tiny),
        quadratic_mean_to_full_rms_ratio=quadratic_mean_rms / max(full_q, tiny),
        mixed_mean_to_aggregate_mean_rms_ratio=mixed_mean_rms / max(aggregate_mean_rms, tiny),
        quadratic_mean_to_aggregate_mean_rms_ratio=quadratic_mean_rms / max(aggregate_mean_rms, tiny),
        backend=backend,
        backend_kind=backend_kind,
    )


def materialize_oscillatory_nonlinear_mean_attribution(
    handoff_backend: ExactAgent2OscillatoryNonlinearHandoff,
    radius: Any,
    z: Any,
    t: Any,
) -> OscillatoryNonlinearMeanAttributionWitness:
    """Project exact A2 #908 nonlinear pieces without accepting surrogate means."""

    if not isinstance(handoff_backend, ExactAgent2OscillatoryNonlinearHandoff):
        raise TypeError("handoff_backend must be ExactAgent2OscillatoryNonlinearHandoff")
    if handoff_backend.source_blob_sha != AGENT2_SOURCE_BLOB_SHA:
        raise ValueError("Agent-2 nonlinear handoff source provenance drifted")
    return _materialize_from_provider(
        handoff_backend.components,
        radius,
        z,
        t,
        backend=handoff_backend,
        backend_kind="exact-agent2-pr908-current-typed-nonlinear-attribution",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_oscillatory_nonlinear_mean_attribution)
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "inverse", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "angular_order",
        "spatial_step", "time_step", "viscosity", "nu", "delta_y", "delta_a",
        "normalized_score", "scientific_threshold",
    }
    return {
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_pr": AGENT2_PR,
        "agent2_head": AGENT2_HEAD,
        "agent2_source_blob": AGENT2_SOURCE_BLOB_SHA,
        "agent2_transport_parent_pr": AGENT2_PARENT_PR,
        "agent2_transport_parent_head": AGENT2_PARENT_HEAD,
        "angular_orders": list(ANGULAR_ORDERS),
        "strict_inner_only": True,
        "current_typed_nonlinear_decomposition_consumed": True,
        "inner_advects_oscillation_mean_materialized": True,
        "oscillation_advects_inner_mean_materialized": True,
        "mixed_cross_mean_materialized": True,
        "oscillatory_quadratic_self_mean_materialized": True,
        "aggregate_nonlinear_mean_materialized": True,
        "three_piece_mean_attribution_checked": True,
        "agent2_curl_or_advection_reimplemented_by_agent3": False,
        "radial_inverse_performed_in_this_increment": False,
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "global_corrected_leading_join_materialized": False,
        "complete_ns_defect": False,
        "scoped_nonlinear_mean_authorized_as_correction_target": False,
        "mean_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "caller_supplied_mean_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


def _cartesian_from_cylindrical(
    cylindrical: np.ndarray,
    theta: np.ndarray,
) -> np.ndarray:
    radial = cylindrical[..., 0]
    tangential = cylindrical[..., 1]
    axial = cylindrical[..., 2]
    c = np.cos(theta)
    s = np.sin(theta)
    return np.stack(
        (radial * c - tangential * s, radial * s + tangential * c, axial),
        axis=-1,
    )


def _analytic_provider(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    del z, t
    theta = np.arctan2(np.asarray(y, dtype=float), np.asarray(x, dtype=float))
    a0 = np.asarray([0.21, -0.06, 0.14], dtype=float)
    b0 = np.asarray([-0.08, 0.02, 0.11], dtype=float)
    q0 = np.asarray([0.17, -0.03, 0.09], dtype=float)

    m2c = np.cos(2.0 * theta)
    m2s = np.sin(2.0 * theta)
    a_cyl = np.stack(
        (
            a0[0] + 0.03 * m2c,
            a0[1] + 0.02 * m2s,
            a0[2] + 0.01 * m2c,
        ),
        axis=-1,
    )
    b_cyl = np.stack(
        (
            b0[0] - 0.02 * m2s,
            b0[1] + 0.01 * m2c,
            b0[2] - 0.015 * m2s,
        ),
        axis=-1,
    )
    q_cyl = np.stack(
        (
            q0[0] + 0.025 * m2c,
            q0[1] - 0.015 * m2s,
            q0[2] + 0.02 * m2c,
        ),
        axis=-1,
    )
    a = _cartesian_from_cylindrical(a_cyl, theta)
    b = _cartesian_from_cylindrical(b_cyl, theta)
    q = _cartesian_from_cylindrical(q_cyl, theta)
    mixed = a + b
    aggregate = mixed + q
    return a, b, q, mixed, aggregate


def analytic_regression_receipt() -> dict[str, object]:
    witness = _materialize_from_provider(
        _analytic_provider,
        radius=np.asarray([0.35, 0.65, 0.95]),
        z=np.asarray([-0.2, 0.0, 0.25]),
        t=np.asarray([0.44, 0.50, 0.56]),
        backend=None,
        backend_kind="analytic-projection-regression-only",
    )
    expected_a = np.broadcast_to(np.asarray([0.21, -0.06, 0.14]), (3, 3))
    expected_b = np.broadcast_to(np.asarray([-0.08, 0.02, 0.11]), (3, 3))
    expected_q = np.broadcast_to(np.asarray([0.17, -0.03, 0.09]), (3, 3))
    expected_total = expected_a + expected_b + expected_q
    errors = {
        "inner_advects_oscillation_abs_max": float(
            np.max(np.abs(witness.mean_inner_advects_oscillation_cylindrical - expected_a))
        ),
        "oscillation_advects_inner_abs_max": float(
            np.max(np.abs(witness.mean_oscillation_advects_inner_cylindrical - expected_b))
        ),
        "quadratic_self_abs_max": float(
            np.max(np.abs(witness.mean_oscillatory_self_advection_cylindrical - expected_q))
        ),
        "aggregate_abs_max": float(
            np.max(np.abs(witness.mean_aggregate_nonlinear_cylindrical - expected_total))
        ),
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_pr": AGENT2_PR,
        "agent2_expected_head": AGENT2_HEAD,
        "agent2_source_blob": AGENT2_SOURCE_BLOB_SHA,
        "angular_orders": list(ANGULAR_ORDERS),
        "registered_means": {
            "inner_advects_oscillation": [0.21, -0.06, 0.14],
            "oscillation_advects_inner": [-0.08, 0.02, 0.11],
            "mixed_cross": [0.13, -0.04, 0.25],
            "oscillatory_self": [0.17, -0.03, 0.09],
            "aggregate": [0.30, -0.07, 0.34],
        },
        "projection_errors": errors,
        "witness": witness.to_receipt(),
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "truth_boundary": truth_boundary(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = analytic_regression_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
