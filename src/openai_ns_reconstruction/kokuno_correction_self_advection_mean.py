"""Agent-3 cylindrical mean projection of Agent-2 correction self-advection.

This module owns only the mean-projection seam for

    N_corr = (delta_u . grad) delta_u.

Agent 2 owns the signed complete-curl correction field, its spatial Jacobian,
and the raw correction self-advection evaluator. Agent 3 projects that already
materialized vector field onto the axisymmetric m=0 cylindrical coefficient
needed by later mean/radial correction bookkeeping. The result is a nonlinear
subterm, not a complete Navier--Stokes defect and not evidence of residual
reduction.
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

TASK = "KOKUNO-A3-CORRECTION-SELF-ADVECTION-MEAN-066"
SCHEMA = "kokuno-a3-correction-self-advection-mean-v1"
PARENT_AGENT3_PR = 789
PARENT_AGENT3_HEAD = "053845282c8f01e53adb417a0a0cf66f7fc2a402"

AGENT2_PR = 795
AGENT2_HEAD = "af99abe67feb45e29614972ab791f227ea004aa6"
AGENT2_MODULE = "openai_ns_reconstruction.kokuno_public_correction_self_advection"
AGENT2_FUNCTION = "evaluate_correction_self_advection_fd4"
AGENT2_SOURCE_BLOB_SHA = "44135ef640a1f076628c9e63f5997fe6400f18f1"
AGENT2_CORRECTION_SPATIAL_STEP = 1.0e-3

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

SelfAdvectionEvaluator = Callable[..., dict[str, Any]]


@dataclass(frozen=True)
class ExactAgent2CorrectionSelfAdvectionBackend:
    """Executable binding to the exact Agent-2 correction self-advection bytes."""

    evaluator: SelfAdvectionEvaluator
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(
        cls, evaluator: SelfAdvectionEvaluator
    ) -> "ExactAgent2CorrectionSelfAdvectionBackend":
        if not callable(evaluator):
            raise TypeError("Agent-2 correction self-advection evaluator must be callable")
        if getattr(evaluator, "__module__", None) != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 correction self-advection module identity")
        if getattr(evaluator, "__name__", None) != AGENT2_FUNCTION:
            raise ValueError("unexpected Agent-2 correction self-advection function identity")
        source = inspect.getsourcefile(evaluator)
        if source is None:
            raise ValueError("Agent-2 correction self-advection source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 correction self-advection source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError(
                "Agent-2 correction self-advection source blob does not match pinned PR #795"
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
class CorrectionSelfAdvectionMeanWitness:
    """Typed m=0 projection of the A2 pure-correction nonlinearity."""

    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_correction_self_advection_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_mean_relative_differences: tuple[float, ...]
    correction_self_mean_rms: float
    full_ring_correction_self_rms: float
    correction_self_mean_to_full_rms_ratio: float
    backend: ExactAgent2CorrectionSelfAdvectionBackend | None
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
            "mean_correction_self_advection_cylindrical": np.asarray(
                self.mean_correction_self_advection_cylindrical, dtype=float
            ).tolist(),
            "angular_orders": list(self.angular_orders),
            "successive_mean_relative_differences": list(
                self.successive_mean_relative_differences
            ),
            "correction_self_mean_rms": self.correction_self_mean_rms,
            "full_ring_correction_self_rms": self.full_ring_correction_self_rms,
            "correction_self_mean_to_full_rms_ratio": (
                self.correction_self_mean_to_full_rms_ratio
            ),
            "backend_kind": self.backend_kind,
            "correction_kind": self.correction_kind,
            "backend": None if self.backend is None else self.backend.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _mean_for_order(
    evaluator: SelfAdvectionEvaluator,
    correction: object,
    radius: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    order: int,
) -> tuple[np.ndarray, np.ndarray]:
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
        spatial_step=AGENT2_CORRECTION_SPATIAL_STEP,
    )
    if not isinstance(diagnostic, dict) or "correction_self_advection" not in diagnostic:
        raise TypeError("Agent-2 evaluator must return correction_self_advection")

    expected_shape = radius.shape + (order, 3)
    self_advection = np.asarray(diagnostic["correction_self_advection"], dtype=float)
    if self_advection.shape != expected_shape:
        raise ValueError("Agent-2 correction self-advection payload has an unexpected ring shape")
    if not np.all(np.isfinite(self_advection)):
        raise ValueError("Agent-2 correction self-advection payload contains non-finite values")

    mean_value = _project_cartesian_ring_to_cylindrical_mean(self_advection, c, s)
    full_ring_rms = np.sqrt(
        np.mean(np.sum(self_advection * self_advection, axis=-1), axis=-1)
    )
    return mean_value, full_ring_rms


def _materialize_with_evaluator(
    evaluator: SelfAdvectionEvaluator,
    correction: object,
    radius: Any,
    z: Any,
    t: Any,
    *,
    backend: ExactAgent2CorrectionSelfAdvectionBackend | None,
    backend_kind: str,
) -> CorrectionSelfAdvectionMeanWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    means: list[np.ndarray] = []
    full_rms: np.ndarray | None = None
    for order in ANGULAR_ORDERS:
        mean_value, ring_rms = _mean_for_order(
            evaluator, correction, rb, zb, tb, order
        )
        means.append(mean_value)
        full_rms = ring_rms

    successive = tuple(
        _relative_rms_difference(means[index + 1], means[index])
        for index in range(len(means) - 1)
    )
    final_mean = means[-1]
    assert full_rms is not None
    correction_self_mean_rms = float(
        np.sqrt(np.mean(np.sum(final_mean * final_mean, axis=-1)))
    )
    full_ring_correction_self_rms = float(np.sqrt(np.mean(full_rms * full_rms)))
    ratio = correction_self_mean_rms / max(
        full_ring_correction_self_rms, np.finfo(float).tiny
    )
    correction_kind = f"{type(correction).__module__}.{type(correction).__qualname__}"

    return CorrectionSelfAdvectionMeanWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_correction_self_advection_cylindrical=np.array(final_mean, copy=True),
        angular_orders=ANGULAR_ORDERS,
        successive_mean_relative_differences=successive,
        correction_self_mean_rms=correction_self_mean_rms,
        full_ring_correction_self_rms=full_ring_correction_self_rms,
        correction_self_mean_to_full_rms_ratio=float(ratio),
        backend=backend,
        backend_kind=backend_kind,
        correction_kind=correction_kind,
    )


def materialize_agent2_correction_self_advection_mean(
    backend: ExactAgent2CorrectionSelfAdvectionBackend,
    correction: object,
    radius: Any,
    z: Any,
    t: Any,
) -> CorrectionSelfAdvectionMeanWitness:
    """Project the exact pinned A2 pure-correction nonlinear term onto m=0.

    The correction object is passed only to the exact A2 evaluator. The caller
    cannot provide a residual, defect, mean value, stress, inverse target, gain,
    damping choice, scientific threshold, raw ``delta_y`` or raw ``delta_a``.
    """
    if not isinstance(backend, ExactAgent2CorrectionSelfAdvectionBackend):
        raise TypeError("backend must be ExactAgent2CorrectionSelfAdvectionBackend")
    if backend.source_blob_sha != AGENT2_SOURCE_BLOB_SHA:
        raise ValueError("Agent-2 correction self-advection backend provenance drifted")
    return _materialize_with_evaluator(
        backend.evaluator,
        correction,
        radius,
        z,
        t,
        backend=backend,
        backend_kind="exact_agent2_pr795_correction_self_advection",
    )


def truth_boundary() -> dict[str, object]:
    return {
        "agent2_correction_self_advection_reimplemented_by_agent3": False,
        "cylindrical_m0_correction_self_projection_executable": True,
        "projection_uses_naive_cartesian_vector_mean": False,
        "correction_self_mean_is_complete_ns_defect": False,
        "oscillatory_self_advection_mean_included_here": False,
        "oscillation_correction_cross_mean_included_here": False,
        "correction_nonlinear_mean_aggregate_complete": False,
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
    spatial_step: float,
) -> dict[str, np.ndarray]:
    del correction, z, t
    if spatial_step != AGENT2_CORRECTION_SPATIAL_STEP:
        raise ValueError("unexpected frozen correction spatial step")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    theta = np.arctan2(y, x)
    c = np.cos(theta)
    s = np.sin(theta)

    radial = 0.12 + 0.07 * np.cos(2.0 * theta + 0.3)
    tangential = -0.04 + 0.05 * np.sin(2.0 * theta - 0.2)
    axial = 0.18 + 0.03 * np.cos(2.0 * theta + 0.6)
    vector = np.stack(
        [radial * c - tangential * s, radial * s + tangential * c, axial],
        axis=-1,
    )
    return {"correction_self_advection": vector}


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
    expected = np.broadcast_to(np.asarray([0.12, -0.04, 0.18]), (3, 3))
    error = float(
        np.max(np.abs(witness.mean_correction_self_advection_cylindrical - expected))
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
        "real_agent2_numeric_correction_self_mean_observed_in_this_receipt": False,
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
