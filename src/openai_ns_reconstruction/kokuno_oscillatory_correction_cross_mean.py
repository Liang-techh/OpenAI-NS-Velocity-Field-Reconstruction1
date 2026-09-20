"""Agent-3 cylindrical mean projection of Agent-2 oscillation/correction cross advection.

This module owns only the mean-projection seam for the mixed nonlinear term

    N_cross = (u_osc . grad) delta_u + (delta_u . grad) u_osc.

Agent 2 owns the oscillatory field, correction field, spatial Jacobians and raw
mixed-advection evaluator. Agent 3 projects that already-materialized vector
field onto the axisymmetric m=0 cylindrical coefficient needed by later
mean/radial correction bookkeeping. The result is a nonlinear subterm, not a
complete Navier--Stokes defect and not evidence of residual reduction.
"""
from __future__ import annotations

import argparse
import inspect
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_oscillatory_quadratic_mean import (
    ANGULAR_ORDERS,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    _git_blob_sha,
    _relative_rms_difference,
    _validate_rings,
)

TASK = "KOKUNO-A3-OSCILLATORY-CORRECTION-CROSS-MEAN-065"
SCHEMA = "kokuno-a3-oscillatory-correction-cross-mean-v1"
PARENT_AGENT3_PR = 780
PARENT_AGENT3_HEAD = "4bbfb3ac60ebfc235a6ef691589199f8e80cfdae"

AGENT2_PR = 788
AGENT2_HEAD = "c1e45c12d1842d8e401a4b54ee926576aacfb7f5"
AGENT2_MODULE = (
    "openai_ns_reconstruction.kokuno_public_oscillatory_correction_cross_advection"
)
AGENT2_FUNCTION = "evaluate_oscillatory_correction_cross_advection"
AGENT2_SOURCE_BLOB_SHA = "8c096a8a1bf91e33f08d1794a61d46d8badd9bcb"
AGENT2_OSCILLATORY_SPATIAL_STEP = 1.0e-3
AGENT2_CORRECTION_SPATIAL_STEP = 1.0e-3

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

CrossAdvectionEvaluator = Callable[..., dict[str, Any]]


@dataclass(frozen=True)
class ExactAgent2CrossAdvectionBackend:
    """Executable binding to the exact Agent-2 mixed-advection source bytes."""

    evaluator: CrossAdvectionEvaluator
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(cls, evaluator: CrossAdvectionEvaluator) -> "ExactAgent2CrossAdvectionBackend":
        if not callable(evaluator):
            raise TypeError("Agent-2 cross-advection evaluator must be callable")
        if getattr(evaluator, "__module__", None) != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 cross-advection module identity")
        if getattr(evaluator, "__name__", None) != AGENT2_FUNCTION:
            raise ValueError("unexpected Agent-2 cross-advection function identity")
        source = inspect.getsourcefile(evaluator)
        if source is None:
            raise ValueError("Agent-2 cross-advection source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 cross-advection source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError("Agent-2 cross-advection source blob does not match pinned PR #788")
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
class OscillatoryCorrectionCrossMeanWitness:
    """Typed m=0 projection of the A2 oscillation/correction mixed nonlinearity."""

    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_oscillation_advects_correction_cylindrical: np.ndarray
    mean_correction_advects_oscillation_cylindrical: np.ndarray
    mean_cross_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_cross_mean_relative_differences: tuple[float, ...]
    termwise_mean_closure_absolute_max: float
    cross_mean_rms: float
    full_ring_cross_rms: float
    cross_mean_to_full_rms_ratio: float
    backend: ExactAgent2CrossAdvectionBackend | None
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
            "mean_oscillation_advects_correction_cylindrical": np.asarray(
                self.mean_oscillation_advects_correction_cylindrical, dtype=float
            ).tolist(),
            "mean_correction_advects_oscillation_cylindrical": np.asarray(
                self.mean_correction_advects_oscillation_cylindrical, dtype=float
            ).tolist(),
            "mean_cross_cylindrical": np.asarray(
                self.mean_cross_cylindrical, dtype=float
            ).tolist(),
            "angular_orders": list(self.angular_orders),
            "successive_cross_mean_relative_differences": list(
                self.successive_cross_mean_relative_differences
            ),
            "termwise_mean_closure_absolute_max": self.termwise_mean_closure_absolute_max,
            "cross_mean_rms": self.cross_mean_rms,
            "full_ring_cross_rms": self.full_ring_cross_rms,
            "cross_mean_to_full_rms_ratio": self.cross_mean_to_full_rms_ratio,
            "backend_kind": self.backend_kind,
            "correction_kind": self.correction_kind,
            "backend": None if self.backend is None else self.backend.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _project_cartesian_ring_to_cylindrical_mean(
    vector: np.ndarray, c: np.ndarray, s: np.ndarray
) -> np.ndarray:
    value = np.asarray(vector, dtype=float)
    if value.shape[-2:] != (c.size, 3):
        raise ValueError("ring vector must have shape (..., angular_order, 3)")
    radial = value[..., 0] * c + value[..., 1] * s
    tangential = -value[..., 0] * s + value[..., 1] * c
    axial = value[..., 2]
    return np.stack(
        [
            np.mean(radial, axis=-1),
            np.mean(tangential, axis=-1),
            np.mean(axial, axis=-1),
        ],
        axis=-1,
    )


def _mean_for_order(
    evaluator: CrossAdvectionEvaluator,
    correction: object,
    radius: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    theta = 2.0 * np.pi * np.arange(order, dtype=float) / float(order)
    c = np.cos(theta)
    s = np.sin(theta)
    xq = radius[..., None] * c
    yq = radius[..., None] * s
    zq = z[..., None] + np.zeros_like(theta)
    tq = t[..., None] + np.zeros_like(theta)

    diagnostic = evaluator(
        correction,
        xq,
        yq,
        zq,
        tq,
        oscillatory_spatial_step=AGENT2_OSCILLATORY_SPATIAL_STEP,
        correction_spatial_step=AGENT2_CORRECTION_SPATIAL_STEP,
    )
    required = (
        "oscillation_advects_correction",
        "correction_advects_oscillation",
        "cross_advection",
    )
    if not isinstance(diagnostic, dict) or not all(key in diagnostic for key in required):
        raise TypeError("Agent-2 evaluator must return both mixed terms and cross_advection")

    expected_shape = radius.shape + (order, 3)
    term_oc = np.asarray(diagnostic[required[0]], dtype=float)
    term_co = np.asarray(diagnostic[required[1]], dtype=float)
    cross = np.asarray(diagnostic[required[2]], dtype=float)
    if term_oc.shape != expected_shape or term_co.shape != expected_shape or cross.shape != expected_shape:
        raise ValueError("Agent-2 mixed-advection payload has an unexpected ring shape")
    if not all(np.all(np.isfinite(a)) for a in (term_oc, term_co, cross)):
        raise ValueError("Agent-2 mixed-advection payload contains non-finite values")

    scale = max(float(np.max(np.abs(cross))), 1.0)
    closure = float(np.max(np.abs(cross - (term_oc + term_co))))
    if closure > 1.0e-12 * scale:
        raise ValueError("Agent-2 mixed-advection termwise closure failed")

    mean_oc = _project_cartesian_ring_to_cylindrical_mean(term_oc, c, s)
    mean_co = _project_cartesian_ring_to_cylindrical_mean(term_co, c, s)
    mean_cross = _project_cartesian_ring_to_cylindrical_mean(cross, c, s)
    mean_closure = float(np.max(np.abs(mean_cross - (mean_oc + mean_co))))
    full_ring_rms = np.sqrt(np.mean(np.sum(cross * cross, axis=-1), axis=-1))
    return mean_oc, mean_co, mean_cross, full_ring_rms, mean_closure


def _materialize_with_evaluator(
    evaluator: CrossAdvectionEvaluator,
    correction: object,
    radius: Any,
    z: Any,
    t: Any,
    *,
    backend: ExactAgent2CrossAdvectionBackend | None,
    backend_kind: str,
) -> OscillatoryCorrectionCrossMeanWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    means_oc: list[np.ndarray] = []
    means_co: list[np.ndarray] = []
    means_cross: list[np.ndarray] = []
    full_rms: np.ndarray | None = None
    mean_closures: list[float] = []
    for order in ANGULAR_ORDERS:
        mean_oc, mean_co, mean_cross, ring_rms, closure = _mean_for_order(
            evaluator, correction, rb, zb, tb, order
        )
        means_oc.append(mean_oc)
        means_co.append(mean_co)
        means_cross.append(mean_cross)
        full_rms = ring_rms
        mean_closures.append(closure)

    successive = tuple(
        _relative_rms_difference(means_cross[index + 1], means_cross[index])
        for index in range(len(means_cross) - 1)
    )
    final_oc = means_oc[-1]
    final_co = means_co[-1]
    final_cross = means_cross[-1]
    assert full_rms is not None
    cross_mean_rms = float(np.sqrt(np.mean(np.sum(final_cross * final_cross, axis=-1))))
    full_ring_cross_rms = float(np.sqrt(np.mean(full_rms * full_rms)))
    ratio = cross_mean_rms / max(full_ring_cross_rms, np.finfo(float).tiny)
    correction_kind = f"{type(correction).__module__}.{type(correction).__qualname__}"

    return OscillatoryCorrectionCrossMeanWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_oscillation_advects_correction_cylindrical=np.array(final_oc, copy=True),
        mean_correction_advects_oscillation_cylindrical=np.array(final_co, copy=True),
        mean_cross_cylindrical=np.array(final_cross, copy=True),
        angular_orders=ANGULAR_ORDERS,
        successive_cross_mean_relative_differences=successive,
        termwise_mean_closure_absolute_max=float(max(mean_closures)),
        cross_mean_rms=cross_mean_rms,
        full_ring_cross_rms=full_ring_cross_rms,
        cross_mean_to_full_rms_ratio=float(ratio),
        backend=backend,
        backend_kind=backend_kind,
        correction_kind=correction_kind,
    )


def materialize_agent2_oscillatory_correction_cross_mean(
    backend: ExactAgent2CrossAdvectionBackend,
    correction: object,
    radius: Any,
    z: Any,
    t: Any,
) -> OscillatoryCorrectionCrossMeanWitness:
    """Project the exact pinned A2 mixed term onto its cylindrical m=0 coefficient.

    The correction object is passed only to the exact A2 evaluator. The caller
    cannot provide a residual, defect, mean value, stress, inverse target, gain,
    damping choice, scientific threshold, raw ``delta_y`` or raw ``delta_a``.
    """
    if not isinstance(backend, ExactAgent2CrossAdvectionBackend):
        raise TypeError("backend must be ExactAgent2CrossAdvectionBackend")
    if backend.source_blob_sha != AGENT2_SOURCE_BLOB_SHA:
        raise ValueError("Agent-2 cross-advection backend source provenance drifted")
    return _materialize_with_evaluator(
        backend.evaluator,
        correction,
        radius,
        z,
        t,
        backend=backend,
        backend_kind="exact_agent2_pr788_oscillation_correction_cross_advection",
    )


def truth_boundary() -> dict[str, object]:
    return {
        "agent2_cross_advection_reimplemented_by_agent3": False,
        "cylindrical_m0_cross_projection_executable": True,
        "projection_uses_naive_cartesian_vector_mean": False,
        "cross_mean_is_complete_ns_defect": False,
        "oscillatory_self_advection_mean_included_here": False,
        "correction_self_advection_included": False,
        "leading_cross_terms_included": False,
        "pressure_gradient_included": False,
        "viscous_term_included": False,
        "restricted_forcing_included": False,
        "real_agent3_delta_a_bound": False,
        "mean_correction_velocity_materialized_from_this_witness": False,
        "real_full_candidate_bound": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_mean_values_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "pde_validated": False,
        "paper_exact": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


def _analytic_regression_evaluator(
    correction: object,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    oscillatory_spatial_step: float,
    correction_spatial_step: float,
) -> dict[str, np.ndarray]:
    del correction, z, t
    if oscillatory_spatial_step != AGENT2_OSCILLATORY_SPATIAL_STEP:
        raise ValueError("unexpected frozen oscillatory spatial step")
    if correction_spatial_step != AGENT2_CORRECTION_SPATIAL_STEP:
        raise ValueError("unexpected frozen correction spatial step")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    theta = np.arctan2(y, x)
    c = np.cos(theta)
    s = np.sin(theta)

    def cylindrical_vector(mean: tuple[float, float, float], phase: float) -> np.ndarray:
        radial = mean[0] + 0.07 * np.cos(2.0 * theta + phase)
        tangential = mean[1] + 0.05 * np.sin(2.0 * theta - phase)
        axial = mean[2] + 0.03 * np.cos(2.0 * theta + 0.5 * phase)
        return np.stack(
            [radial * c - tangential * s, radial * s + tangential * c, axial],
            axis=-1,
        )

    term_oc = cylindrical_vector((0.40, -0.10, 0.20), 0.2)
    term_co = cylindrical_vector((-0.15, 0.05, 0.10), -0.4)
    return {
        "oscillation_advects_correction": term_oc,
        "correction_advects_oscillation": term_co,
        "cross_advection": term_oc + term_co,
    }


def analytic_regression_receipt() -> dict[str, object]:
    witness = _materialize_with_evaluator(
        _analytic_regression_evaluator,
        correction=object(),
        radius=np.asarray([0.45, 0.8, 1.15]),
        z=np.asarray([-0.6, 0.0, 0.7]),
        t=np.asarray([0.43, 0.50, 0.57]),
        backend=None,
        backend_kind="analytic_projection_regression_only",
    )
    expected_oc = np.broadcast_to(np.asarray([0.40, -0.10, 0.20]), (3, 3))
    expected_co = np.broadcast_to(np.asarray([-0.15, 0.05, 0.10]), (3, 3))
    expected_cross = expected_oc + expected_co
    error = float(
        max(
            np.max(np.abs(witness.mean_oscillation_advects_correction_cylindrical - expected_oc)),
            np.max(np.abs(witness.mean_correction_advects_oscillation_cylindrical - expected_co)),
            np.max(np.abs(witness.mean_cross_cylindrical - expected_cross)),
        )
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
        "analytic_projection_absolute_max_error": error,
        "analytic_witness": witness.to_receipt(),
        "real_agent2_numeric_cross_mean_observed_in_this_receipt": False,
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
