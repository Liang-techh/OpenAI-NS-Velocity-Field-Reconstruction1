"""Agent-3 cylindrical mean projection of Agent-2 correction transport increment.

This module owns only the m=0 mean-projection seam for the exact Agent-2
correction-induced local transport increment

    delta T_A2 = d_t(delta u)
               + (u_osc . grad) delta u
               + (delta u . grad) u_osc
               + (delta u . grad) delta u
               - nu Delta(delta u).

Agent 2 owns every raw correction/oscillatory derivative in that expression.
Agent 3 only projects the already-materialized Cartesian vector terms into the
rotating cylindrical frame and checks that the projected components close back
to the projected total.  This object is still not a complete Navier--Stokes
defect: leading-field cross terms, pressure gradient and restricted forcing are
not included.
"""
from __future__ import annotations

import argparse
import inspect
import json
from dataclasses import dataclass
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

TASK = "KOKUNO-A3-CORRECTION-TRANSPORT-MEAN-067"
SCHEMA = "kokuno-a3-correction-transport-mean-v1"
PARENT_AGENT3_PR = 797
PARENT_AGENT3_HEAD = "5e329af22540ca791a61482ab777769eb337b1d9"

AGENT2_PR = 804
AGENT2_HEAD = "9bd1623f238a7827019a336aee435216fc38d393"
AGENT2_MODULE = "openai_ns_reconstruction.kokuno_public_correction_transport_increment"
AGENT2_FUNCTION = "evaluate_correction_transport_increment"
AGENT2_SOURCE_BLOB_SHA = "d8becee43f0954fe0b476c187388591109cf0228"
AGENT2_VISCOSITY = 0.01

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

COMPONENT_KEYS = (
    "correction_velocity_dt",
    "cross_advection",
    "correction_self_advection",
    "correction_viscous_term",
)
TOTAL_KEY = "correction_transport_increment"

TransportEvaluator = Callable[..., dict[str, Any]]


@dataclass(frozen=True)
class ExactAgent2CorrectionTransportBackend:
    """Executable binding to the exact Agent-2 #804 transport-composition bytes."""

    evaluator: TransportEvaluator
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(
        cls, evaluator: TransportEvaluator
    ) -> "ExactAgent2CorrectionTransportBackend":
        if not callable(evaluator):
            raise TypeError("Agent-2 correction transport evaluator must be callable")
        if getattr(evaluator, "__module__", None) != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 correction transport module identity")
        if getattr(evaluator, "__name__", None) != AGENT2_FUNCTION:
            raise ValueError("unexpected Agent-2 correction transport function identity")
        source = inspect.getsourcefile(evaluator)
        if source is None:
            raise ValueError("Agent-2 correction transport source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 correction transport source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError(
                "Agent-2 correction transport source blob does not match pinned PR #804"
            )
        return cls(evaluator=evaluator, source_file=str(path), source_blob_sha=blob_sha)

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent2_pr": AGENT2_PR,
            "agent2_expected_head": AGENT2_HEAD,
            "module": AGENT2_MODULE,
            "function": AGENT2_FUNCTION,
            "source_blob_sha": self.source_blob_sha,
            "source_blob_matches_pinned": self.source_blob_sha == AGENT2_SOURCE_BLOB_SHA,
        }


@dataclass(frozen=True)
class CorrectionTransportMeanWitness:
    """Typed m=0 projection of the exact A2 correction transport increment."""

    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_correction_transport_cylindrical: np.ndarray
    mean_correction_velocity_dt_cylindrical: np.ndarray
    mean_cross_advection_cylindrical: np.ndarray
    mean_correction_self_advection_cylindrical: np.ndarray
    mean_correction_viscous_term_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_mean_relative_differences: tuple[float, ...]
    pointwise_component_closure_absolute_max: float
    projected_component_closure_absolute_max: float
    projected_component_closure_relative_rms: float
    correction_transport_mean_rms: float
    full_ring_correction_transport_rms: float
    correction_transport_mean_to_full_rms_ratio: float
    viscosity: float
    backend: ExactAgent2CorrectionTransportBackend | None
    backend_kind: str
    correction_kind: str

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "radius": np.asarray(self.radius, dtype=float).tolist(),
            "z": np.asarray(self.z, dtype=float).tolist(),
            "t": np.asarray(self.t, dtype=float).tolist(),
            "cylindrical_components": ["radial", "tangential", "axial"],
            "mean_correction_transport_cylindrical": np.asarray(
                self.mean_correction_transport_cylindrical, dtype=float
            ).tolist(),
            "mean_correction_velocity_dt_cylindrical": np.asarray(
                self.mean_correction_velocity_dt_cylindrical, dtype=float
            ).tolist(),
            "mean_cross_advection_cylindrical": np.asarray(
                self.mean_cross_advection_cylindrical, dtype=float
            ).tolist(),
            "mean_correction_self_advection_cylindrical": np.asarray(
                self.mean_correction_self_advection_cylindrical, dtype=float
            ).tolist(),
            "mean_correction_viscous_term_cylindrical": np.asarray(
                self.mean_correction_viscous_term_cylindrical, dtype=float
            ).tolist(),
            "angular_orders": list(self.angular_orders),
            "successive_mean_relative_differences": list(
                self.successive_mean_relative_differences
            ),
            "pointwise_component_closure_absolute_max": (
                self.pointwise_component_closure_absolute_max
            ),
            "projected_component_closure_absolute_max": (
                self.projected_component_closure_absolute_max
            ),
            "projected_component_closure_relative_rms": (
                self.projected_component_closure_relative_rms
            ),
            "correction_transport_mean_rms": self.correction_transport_mean_rms,
            "full_ring_correction_transport_rms": self.full_ring_correction_transport_rms,
            "correction_transport_mean_to_full_rms_ratio": (
                self.correction_transport_mean_to_full_rms_ratio
            ),
            "viscosity": self.viscosity,
            "backend_kind": self.backend_kind,
            "correction_kind": self.correction_kind,
            "backend": None if self.backend is None else self.backend.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _validate_handoff_contract(diagnostic: dict[str, Any]) -> None:
    if float(diagnostic.get("viscosity", np.nan)) != AGENT2_VISCOSITY:
        raise ValueError("Agent-2 correction transport viscosity is not frozen at 0.01")
    handoff = diagnostic.get("handoff_contract")
    if not isinstance(handoff, dict):
        raise TypeError("Agent-2 correction transport handoff_contract is missing")
    required = {
        "consumer_lane": "Kokuno Agent 3 mean/radial bookkeeping",
        "quantity_kind": "A2 correction-induced local transport increment",
        "complete_ns_defect": False,
        "includes_correction_time_derivative": True,
        "includes_oscillation_correction_cross_advection": True,
        "includes_correction_self_advection": True,
        "includes_correction_viscosity": True,
        "includes_leading_cross_terms": False,
        "includes_pressure_gradient": False,
        "includes_restricted_forcing": False,
    }
    for key, value in required.items():
        if handoff.get(key) != value:
            raise ValueError(f"Agent-2 handoff contract drifted at {key}")


def _mean_for_order(
    evaluator: TransportEvaluator,
    correction: object,
    radius: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    order: int,
) -> tuple[dict[str, np.ndarray], np.ndarray, float, float]:
    theta = 2.0 * np.pi * np.arange(order, dtype=float) / float(order)
    c = np.cos(theta)
    s = np.sin(theta)
    xq = radius[..., None] * c
    yq = radius[..., None] * s
    zq = z[..., None] + np.zeros_like(theta)
    tq = t[..., None] + np.zeros_like(theta)

    diagnostic = evaluator(correction, xq, yq, zq, tq)
    if not isinstance(diagnostic, dict):
        raise TypeError("Agent-2 evaluator must return a dict")
    _validate_handoff_contract(diagnostic)

    expected_shape = radius.shape + (order, 3)
    vectors: dict[str, np.ndarray] = {}
    for key in (*COMPONENT_KEYS, TOTAL_KEY):
        if key not in diagnostic:
            raise TypeError(f"Agent-2 evaluator must return {key}")
        value = np.asarray(diagnostic[key], dtype=float)
        if value.shape != expected_shape:
            raise ValueError(f"Agent-2 {key} payload has an unexpected ring shape")
        if not np.all(np.isfinite(value)):
            raise ValueError(f"Agent-2 {key} payload contains non-finite values")
        vectors[key] = value

    component_sum = sum((vectors[key] for key in COMPONENT_KEYS), np.zeros_like(vectors[TOTAL_KEY]))
    pointwise_closure = vectors[TOTAL_KEY] - component_sum
    pointwise_closure_abs_max = float(np.max(np.abs(pointwise_closure)))

    means = {
        key: _project_cartesian_ring_to_cylindrical_mean(value, c, s)
        for key, value in vectors.items()
    }
    full_ring_rms = np.sqrt(
        np.mean(np.sum(vectors[TOTAL_KEY] * vectors[TOTAL_KEY], axis=-1), axis=-1)
    )
    return means, full_ring_rms, pointwise_closure_abs_max, AGENT2_VISCOSITY


def _materialize_with_evaluator(
    evaluator: TransportEvaluator,
    correction: object,
    radius: Any,
    z: Any,
    t: Any,
    *,
    backend: ExactAgent2CorrectionTransportBackend | None,
    backend_kind: str,
) -> CorrectionTransportMeanWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    means_by_order: list[dict[str, np.ndarray]] = []
    full_rms: np.ndarray | None = None
    pointwise_closure_abs_max = 0.0
    viscosity = AGENT2_VISCOSITY

    for order in ANGULAR_ORDERS:
        means, ring_rms, closure_abs_max, viscosity = _mean_for_order(
            evaluator, correction, rb, zb, tb, order
        )
        means_by_order.append(means)
        full_rms = ring_rms
        pointwise_closure_abs_max = max(pointwise_closure_abs_max, closure_abs_max)

    total_means = [means[TOTAL_KEY] for means in means_by_order]
    successive = tuple(
        _relative_rms_difference(total_means[index + 1], total_means[index])
        for index in range(len(total_means) - 1)
    )
    final = means_by_order[-1]
    final_total = final[TOTAL_KEY]
    projected_component_sum = sum(
        (final[key] for key in COMPONENT_KEYS), np.zeros_like(final_total)
    )
    projected_closure = final_total - projected_component_sum
    projected_closure_abs_max = float(np.max(np.abs(projected_closure)))
    projected_closure_rms = float(
        np.sqrt(np.mean(np.sum(projected_closure * projected_closure, axis=-1)))
    )
    total_mean_rms = float(
        np.sqrt(np.mean(np.sum(final_total * final_total, axis=-1)))
    )
    projected_closure_relative_rms = projected_closure_rms / max(
        total_mean_rms, np.finfo(float).tiny
    )

    assert full_rms is not None
    full_ring_rms_scalar = float(np.sqrt(np.mean(full_rms * full_rms)))
    ratio = total_mean_rms / max(full_ring_rms_scalar, np.finfo(float).tiny)
    correction_kind = f"{type(correction).__module__}.{type(correction).__qualname__}"

    return CorrectionTransportMeanWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_correction_transport_cylindrical=np.array(final[TOTAL_KEY], copy=True),
        mean_correction_velocity_dt_cylindrical=np.array(
            final["correction_velocity_dt"], copy=True
        ),
        mean_cross_advection_cylindrical=np.array(final["cross_advection"], copy=True),
        mean_correction_self_advection_cylindrical=np.array(
            final["correction_self_advection"], copy=True
        ),
        mean_correction_viscous_term_cylindrical=np.array(
            final["correction_viscous_term"], copy=True
        ),
        angular_orders=ANGULAR_ORDERS,
        successive_mean_relative_differences=successive,
        pointwise_component_closure_absolute_max=pointwise_closure_abs_max,
        projected_component_closure_absolute_max=projected_closure_abs_max,
        projected_component_closure_relative_rms=float(projected_closure_relative_rms),
        correction_transport_mean_rms=total_mean_rms,
        full_ring_correction_transport_rms=full_ring_rms_scalar,
        correction_transport_mean_to_full_rms_ratio=float(ratio),
        viscosity=float(viscosity),
        backend=backend,
        backend_kind=backend_kind,
        correction_kind=correction_kind,
    )


def materialize_agent2_correction_transport_mean(
    backend: ExactAgent2CorrectionTransportBackend,
    correction: object,
    radius: Any,
    z: Any,
    t: Any,
) -> CorrectionTransportMeanWitness:
    """Project exact pinned A2 #804 correction transport onto cylindrical m=0.

    The caller cannot supply a residual, defect, mean, stress, inverse, pressure,
    forcing, target, gain, viscosity, derivative resolution, damping, scientific
    threshold, raw ``delta_y``, raw ``delta_a`` or alpha.
    """
    if not isinstance(backend, ExactAgent2CorrectionTransportBackend):
        raise TypeError("backend must be ExactAgent2CorrectionTransportBackend")
    if backend.source_blob_sha != AGENT2_SOURCE_BLOB_SHA:
        raise ValueError("Agent-2 correction transport backend provenance drifted")
    return _materialize_with_evaluator(
        backend.evaluator,
        correction,
        radius,
        z,
        t,
        backend=backend,
        backend_kind="exact_agent2_pr804_correction_transport_increment",
    )


def truth_boundary() -> dict[str, object]:
    return {
        "agent2_correction_transport_reimplemented_by_agent3": False,
        "cylindrical_m0_correction_transport_projection_executable": True,
        "componentwise_projected_closure_checked": True,
        "projection_uses_naive_cartesian_vector_mean": False,
        "correction_transport_mean_is_complete_ns_defect": False,
        "includes_correction_time_derivative": True,
        "includes_oscillation_correction_cross_advection": True,
        "includes_correction_self_advection": True,
        "includes_correction_viscosity": True,
        "leading_cross_terms_included": False,
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "real_agent3_delta_a_bound": False,
        "mean_correction_velocity_materialized_from_this_witness": False,
        "real_full_candidate_bound": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_mean_values_allowed": False,
        "caller_supplied_viscosity_allowed": False,
        "caller_supplied_derivative_resolution_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "pde_validated": False,
        "paper_exact": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


def _cylindrical_vector(
    theta: np.ndarray,
    mean: tuple[float, float, float],
    phase: float,
) -> np.ndarray:
    c = np.cos(theta)
    s = np.sin(theta)
    radial = mean[0] + 0.031 * np.cos(2.0 * theta + phase)
    tangential = mean[1] + 0.023 * np.sin(2.0 * theta - 0.7 * phase)
    axial = mean[2] + 0.019 * np.cos(2.0 * theta + 0.4 + phase)
    return np.stack(
        [radial * c - tangential * s, radial * s + tangential * c, axial],
        axis=-1,
    )


def _analytic_regression_evaluator(
    correction: object,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
) -> dict[str, Any]:
    del correction, z, t
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    theta = np.arctan2(y, x)

    dt = _cylindrical_vector(theta, (0.08, -0.02, 0.03), 0.1)
    cross = _cylindrical_vector(theta, (0.25, -0.05, 0.30), 0.5)
    self_term = _cylindrical_vector(theta, (0.12, -0.04, 0.18), 0.9)
    viscous = _cylindrical_vector(theta, (-0.03, 0.01, -0.05), 1.3)
    total = dt + cross + self_term + viscous
    return {
        "correction_velocity_dt": dt,
        "cross_advection": cross,
        "correction_self_advection": self_term,
        "correction_viscous_term": viscous,
        "correction_transport_increment": total,
        "viscosity": AGENT2_VISCOSITY,
        "handoff_contract": {
            "consumer_lane": "Kokuno Agent 3 mean/radial bookkeeping",
            "quantity_kind": "A2 correction-induced local transport increment",
            "complete_ns_defect": False,
            "includes_correction_time_derivative": True,
            "includes_oscillation_correction_cross_advection": True,
            "includes_correction_self_advection": True,
            "includes_correction_viscosity": True,
            "includes_leading_cross_terms": False,
            "includes_pressure_gradient": False,
            "includes_restricted_forcing": False,
        },
    }


def analytic_regression_receipt() -> dict[str, object]:
    witness = _materialize_with_evaluator(
        _analytic_regression_evaluator,
        correction=object(),
        radius=np.asarray([0.45, 0.80, 1.15]),
        z=np.asarray([-0.6, 0.0, 0.7]),
        t=np.asarray([0.44, 0.50, 0.56]),
        backend=None,
        backend_kind="analytic_projection_regression_only",
    )
    expected_components = {
        "mean_correction_velocity_dt_cylindrical": np.asarray([0.08, -0.02, 0.03]),
        "mean_cross_advection_cylindrical": np.asarray([0.25, -0.05, 0.30]),
        "mean_correction_self_advection_cylindrical": np.asarray([0.12, -0.04, 0.18]),
        "mean_correction_viscous_term_cylindrical": np.asarray([-0.03, 0.01, -0.05]),
        "mean_correction_transport_cylindrical": np.asarray([0.42, -0.10, 0.46]),
    }
    errors = {}
    for field, expected in expected_components.items():
        observed = np.asarray(getattr(witness, field), dtype=float)
        errors[field] = float(
            np.max(np.abs(observed - np.broadcast_to(expected, observed.shape)))
        )
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_pr": AGENT2_PR,
        "agent2_expected_head": AGENT2_HEAD,
        "agent2_source_blob_sha": AGENT2_SOURCE_BLOB_SHA,
        "source_reader": {
            "repository": SOURCE_READER_REPO,
            "head": SOURCE_READER_HEAD,
            "date": SOURCE_READER_DATE,
        },
        "projection_ladder": list(ANGULAR_ORDERS),
        "analytic_projection_absolute_max_errors": errors,
        "analytic_witness": witness.to_receipt(),
        "real_agent2_numeric_correction_transport_mean_observed_in_this_receipt": False,
        "truth_boundary": truth_boundary(),
    }


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    receipt = analytic_regression_receipt()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(output)


if __name__ == "__main__":
    _main()
