"""Agent-3 cylindrical m=0 projection of A2 inner-leading/oscillatory cross advection.

Agent 2 PR #812 exposes the raw Cartesian mixed convective term

    N_inner,osc = (u_inner . grad) u_osc + (u_osc . grad) u_inner

for the currently executable Agent-1 PA.10 *inner contraction-center* Cartesian
velocity.  Agent 3 owns only the rotating cylindrical mean needed by later
mean/radial bookkeeping.  The inner field is not the corrected fixed point and
has no outer/global join, matched pressure, restricted forcing or full NS
admission, so this witness is intentionally a scoped subterm rather than a
complete defect.
"""
from __future__ import annotations

import argparse
import inspect
import json
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
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

TASK = "KOKUNO-A3-INNER-LEADING-OSCILLATORY-CROSS-MEAN-068"
SCHEMA = "kokuno-a3-inner-leading-oscillatory-cross-mean-v1"
PARENT_AGENT3_PR = 805
PARENT_AGENT3_HEAD = "dc4eb0791eb6c31413d5d5ef1cd7d71c0411e0ae"

AGENT2_PR = 812
AGENT2_HEAD = "39a55b6f80a4bf9dfa6e0d6bee7825ef0cf7e600"
AGENT2_MODULE = (
    "openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_cross_advection"
)
AGENT2_FUNCTION = "evaluate_inner_leading_oscillatory_cross_advection"
AGENT2_SOURCE_BLOB_SHA = "6d37b8a2db30e29a4ee199a0a96443e77b3d49ec"
AGENT2_INNER_LEADING_FD4_STEP = 1.0e-3
AGENT2_OSCILLATORY_FD6_STEP = 1.0e-3

AGENT1_PR = 803
AGENT1_HEAD = "0666fff748d7c2659774b0aacba70042d045aab5"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

CrossEvaluator = Callable[..., Any]


@dataclass(frozen=True)
class ExactAgent2InnerLeadingOscillatoryCrossBackend:
    """Bind the exact Agent-2 #812 raw cross-advection implementation bytes."""

    evaluator: CrossEvaluator
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(
        cls, evaluator: CrossEvaluator
    ) -> "ExactAgent2InnerLeadingOscillatoryCrossBackend":
        if not callable(evaluator):
            raise TypeError("Agent-2 inner-leading/oscillatory evaluator must be callable")
        if getattr(evaluator, "__module__", None) != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 inner-leading/oscillatory module identity")
        if getattr(evaluator, "__name__", None) != AGENT2_FUNCTION:
            raise ValueError("unexpected Agent-2 inner-leading/oscillatory function identity")
        source = inspect.getsourcefile(evaluator)
        if source is None:
            raise ValueError("Agent-2 inner-leading/oscillatory source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 inner-leading/oscillatory source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError("Agent-2 source blob does not match pinned PR #812")
        return cls(evaluator=evaluator, source_file=str(path), source_blob_sha=blob_sha)

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent2_pr": AGENT2_PR,
            "agent2_expected_head": AGENT2_HEAD,
            "agent1_inner_velocity_pr": AGENT1_PR,
            "agent1_inner_velocity_expected_head": AGENT1_HEAD,
            "module": AGENT2_MODULE,
            "function": AGENT2_FUNCTION,
            "source_blob_sha": self.source_blob_sha,
            "source_blob_matches_pinned": self.source_blob_sha == AGENT2_SOURCE_BLOB_SHA,
        }


@dataclass(frozen=True)
class InnerLeadingOscillatoryCrossMeanWitness:
    """Typed cylindrical mean of the scoped A1-inner/A2-oscillatory cross term."""

    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_inner_advects_oscillation_cylindrical: np.ndarray
    mean_oscillation_advects_inner_cylindrical: np.ndarray
    mean_cross_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_cross_mean_relative_differences: tuple[float, ...]
    termwise_mean_closure_absolute_max: float
    cross_mean_rms: float
    full_ring_cross_rms: float
    cross_mean_to_full_rms_ratio: float
    backend: ExactAgent2InnerLeadingOscillatoryCrossBackend | None
    backend_kind: str
    inner_leading_backend_kind: str

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
            "inner_leading_backend_kind": self.inner_leading_backend_kind,
            "backend": None if self.backend is None else self.backend.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _validate_result_steps(result: Any) -> None:
    lead_step = float(getattr(result, "inner_leading_spatial_step", np.nan))
    osc_step = float(getattr(result, "oscillatory_spatial_step", np.nan))
    if lead_step != AGENT2_INNER_LEADING_FD4_STEP:
        raise ValueError("Agent-2 inner-leading FD4 step drifted from frozen 1e-3")
    if osc_step != AGENT2_OSCILLATORY_FD6_STEP:
        raise ValueError("Agent-2 oscillatory FD6 step drifted from frozen 1e-3")


def _mean_for_order(
    evaluator: CrossEvaluator,
    inner_leading_backend: object,
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

    result = evaluator(
        inner_leading_backend,
        xq,
        yq,
        zq,
        tq,
        inner_leading_spatial_step=AGENT2_INNER_LEADING_FD4_STEP,
        oscillatory_spatial_step=AGENT2_OSCILLATORY_FD6_STEP,
    )
    _validate_result_steps(result)

    required = (
        "inner_advects_oscillation",
        "oscillation_advects_inner",
        "cross_advection",
    )
    if not all(hasattr(result, key) for key in required):
        raise TypeError("Agent-2 result must expose both mixed terms and cross_advection")

    expected_shape = radius.shape + (order, 3)
    term_io = np.asarray(getattr(result, required[0]), dtype=float)
    term_oi = np.asarray(getattr(result, required[1]), dtype=float)
    cross = np.asarray(getattr(result, required[2]), dtype=float)
    if (
        term_io.shape != expected_shape
        or term_oi.shape != expected_shape
        or cross.shape != expected_shape
    ):
        raise ValueError("Agent-2 inner-leading/oscillatory payload has an unexpected ring shape")
    if not all(np.all(np.isfinite(a)) for a in (term_io, term_oi, cross)):
        raise ValueError("Agent-2 inner-leading/oscillatory payload contains non-finite values")

    scale = max(float(np.max(np.abs(cross))), 1.0)
    pointwise_closure = float(np.max(np.abs(cross - term_io - term_oi)))
    if pointwise_closure > 1.0e-12 * scale:
        raise ValueError("Agent-2 inner-leading/oscillatory termwise closure failed")

    mean_io = _project_cartesian_ring_to_cylindrical_mean(term_io, c, s)
    mean_oi = _project_cartesian_ring_to_cylindrical_mean(term_oi, c, s)
    mean_cross = _project_cartesian_ring_to_cylindrical_mean(cross, c, s)
    mean_closure = float(np.max(np.abs(mean_cross - mean_io - mean_oi)))
    full_ring_rms = np.sqrt(np.mean(np.sum(cross * cross, axis=-1), axis=-1))
    return mean_io, mean_oi, mean_cross, full_ring_rms, mean_closure


def _materialize_with_evaluator(
    evaluator: CrossEvaluator,
    inner_leading_backend: object,
    radius: Any,
    z: Any,
    t: Any,
    *,
    backend: ExactAgent2InnerLeadingOscillatoryCrossBackend | None,
    backend_kind: str,
) -> InnerLeadingOscillatoryCrossMeanWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    means_io: list[np.ndarray] = []
    means_oi: list[np.ndarray] = []
    means_cross: list[np.ndarray] = []
    full_rms: np.ndarray | None = None
    mean_closures: list[float] = []

    for order in ANGULAR_ORDERS:
        mean_io, mean_oi, mean_cross, ring_rms, mean_closure = _mean_for_order(
            evaluator, inner_leading_backend, rb, zb, tb, order
        )
        means_io.append(mean_io)
        means_oi.append(mean_oi)
        means_cross.append(mean_cross)
        full_rms = ring_rms
        mean_closures.append(mean_closure)

    successive = tuple(
        _relative_rms_difference(means_cross[index + 1], means_cross[index])
        for index in range(len(means_cross) - 1)
    )
    final_io = means_io[-1]
    final_oi = means_oi[-1]
    final_cross = means_cross[-1]
    assert full_rms is not None
    cross_mean_rms = float(np.sqrt(np.mean(np.sum(final_cross * final_cross, axis=-1))))
    full_ring_cross_rms = float(np.sqrt(np.mean(full_rms * full_rms)))
    ratio = cross_mean_rms / max(full_ring_cross_rms, np.finfo(float).tiny)
    inner_kind = (
        f"{type(inner_leading_backend).__module__}."
        f"{type(inner_leading_backend).__qualname__}"
    )

    return InnerLeadingOscillatoryCrossMeanWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_inner_advects_oscillation_cylindrical=np.array(final_io, copy=True),
        mean_oscillation_advects_inner_cylindrical=np.array(final_oi, copy=True),
        mean_cross_cylindrical=np.array(final_cross, copy=True),
        angular_orders=ANGULAR_ORDERS,
        successive_cross_mean_relative_differences=successive,
        termwise_mean_closure_absolute_max=float(max(mean_closures)),
        cross_mean_rms=cross_mean_rms,
        full_ring_cross_rms=full_ring_cross_rms,
        cross_mean_to_full_rms_ratio=float(ratio),
        backend=backend,
        backend_kind=backend_kind,
        inner_leading_backend_kind=inner_kind,
    )


def materialize_agent2_inner_leading_oscillatory_cross_mean(
    backend: ExactAgent2InnerLeadingOscillatoryCrossBackend,
    inner_leading_backend: object,
    radius: Any,
    z: Any,
    t: Any,
) -> InnerLeadingOscillatoryCrossMeanWitness:
    """Project pinned A2 #812 inner-leading/oscillatory cross advection to m=0.

    The caller cannot supply residual/defect values, mean values, radial stress,
    inverse data, pressure, forcing, target/gain/damping, derivative resolution,
    angular order, viscosity, raw ``delta_y``/``delta_a`` or scientific gates.
    """
    if not isinstance(backend, ExactAgent2InnerLeadingOscillatoryCrossBackend):
        raise TypeError("backend must be ExactAgent2InnerLeadingOscillatoryCrossBackend")
    if backend.source_blob_sha != AGENT2_SOURCE_BLOB_SHA:
        raise ValueError("Agent-2 inner-leading/oscillatory backend provenance drifted")
    if not hasattr(inner_leading_backend, "velocity") or not callable(
        getattr(inner_leading_backend, "velocity")
    ):
        raise TypeError("inner_leading_backend must expose velocity(x,y,z,t)")
    return _materialize_with_evaluator(
        backend.evaluator,
        inner_leading_backend,
        radius,
        z,
        t,
        backend=backend,
        backend_kind="exact_agent2_pr812_inner_leading_oscillatory_cross_advection",
    )


def truth_boundary() -> dict[str, object]:
    return {
        "agent2_inner_leading_oscillatory_cross_reimplemented_by_agent3": False,
        "cylindrical_m0_inner_leading_oscillatory_cross_projection_executable": True,
        "projection_uses_naive_cartesian_vector_mean": False,
        "uses_agent1_pa10_inner_contraction_center_only": True,
        "inner_leading_is_final_corrected_fixed_point": False,
        "outer_global_leading_join_materialized": False,
        "global_leading_velocity_materialized": False,
        "inner_leading_oscillatory_cross_mean_is_complete_ns_defect": False,
        "inner_leading_self_advection_included_here": False,
        "inner_leading_velocity_dt_included_here": False,
        "oscillatory_self_advection_mean_included_here": False,
        "correction_transport_mean_included_here": False,
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
        "caller_supplied_derivative_resolution_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "pde_validated": False,
        "paper_exact": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


def _analytic_regression_evaluator(
    inner_leading_backend: object,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    inner_leading_spatial_step: float,
    oscillatory_spatial_step: float,
) -> Any:
    del inner_leading_backend, z, t
    if inner_leading_spatial_step != AGENT2_INNER_LEADING_FD4_STEP:
        raise ValueError("unexpected frozen inner-leading FD4 step")
    if oscillatory_spatial_step != AGENT2_OSCILLATORY_FD6_STEP:
        raise ValueError("unexpected frozen oscillatory FD6 step")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    theta = np.arctan2(y, x)
    c = np.cos(theta)
    s = np.sin(theta)

    def cylindrical_vector(mean: tuple[float, float, float], phase: float) -> np.ndarray:
        radial = mean[0] + 0.06 * np.cos(2.0 * theta + phase)
        tangential = mean[1] + 0.04 * np.sin(2.0 * theta - phase)
        axial = mean[2] + 0.03 * np.cos(2.0 * theta + 0.5 * phase)
        return np.stack(
            [radial * c - tangential * s, radial * s + tangential * c, axial],
            axis=-1,
        )

    term_io = cylindrical_vector((0.21, -0.06, 0.14), 0.3)
    term_oi = cylindrical_vector((-0.08, 0.02, 0.11), -0.5)
    return SimpleNamespace(
        inner_advects_oscillation=term_io,
        oscillation_advects_inner=term_oi,
        cross_advection=term_io + term_oi,
        inner_leading_spatial_step=AGENT2_INNER_LEADING_FD4_STEP,
        oscillatory_spatial_step=AGENT2_OSCILLATORY_FD6_STEP,
    )


def analytic_regression_receipt() -> dict[str, object]:
    class DummyInner:
        def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
            xb, yb, zb, tb = np.broadcast_arrays(
                np.asarray(x, dtype=float),
                np.asarray(y, dtype=float),
                np.asarray(z, dtype=float),
                np.asarray(t, dtype=float),
            )
            return np.stack((xb * 0.0 + 1.0, yb * 0.0, zb * 0.0 + tb * 0.0), axis=-1)

    witness = _materialize_with_evaluator(
        _analytic_regression_evaluator,
        DummyInner(),
        radius=np.asarray([0.45, 0.80, 1.15]),
        z=np.asarray([-0.5, 0.0, 0.6]),
        t=np.asarray([0.46, 0.50, 0.54]),
        backend=None,
        backend_kind="analytic_projection_regression_only",
    )
    expected_io = np.broadcast_to(np.asarray([0.21, -0.06, 0.14]), (3, 3))
    expected_oi = np.broadcast_to(np.asarray([-0.08, 0.02, 0.11]), (3, 3))
    expected_cross = expected_io + expected_oi
    error = float(
        max(
            np.max(np.abs(witness.mean_inner_advects_oscillation_cylindrical - expected_io)),
            np.max(np.abs(witness.mean_oscillation_advects_inner_cylindrical - expected_oi)),
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
        "agent1_inner_velocity_pr": AGENT1_PR,
        "agent1_inner_velocity_expected_head": AGENT1_HEAD,
        "source_reader": {
            "repository": SOURCE_READER_REPO,
            "head": SOURCE_READER_HEAD,
            "date": SOURCE_READER_DATE,
        },
        "projection_ladder": list(ANGULAR_ORDERS),
        "analytic_projection_absolute_max_error": error,
        "analytic_witness": witness.to_receipt(),
        "real_agent2_numeric_inner_leading_cross_mean_observed_in_this_receipt": False,
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
