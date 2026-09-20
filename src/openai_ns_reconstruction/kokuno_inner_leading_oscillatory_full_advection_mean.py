"""Agent-3 cylindrical m=0 projection of strict-inner leading+oscillatory advection.

Agent 2 PR #830 exposes the complete convective composition for the executable
Agent-1 PA.10 inner contraction-center velocity plus the frozen Agent-2
oscillatory field:

    ((u_inner + u_osc) . grad)(u_inner + u_osc)
      = (u_inner . grad)u_inner
      + (u_inner . grad)u_osc
      + (u_osc . grad)u_inner
      + (u_osc . grad)u_osc.

Agent 3 owns only the rotating cylindrical mean used by downstream mean/radial
bookkeeping. This is still an inner-center scoped convective block, not a
complete Navier--Stokes defect: the corrected/global leading join, time and
viscous completion, matched pressure, restricted forcing, and correction-cycle
full candidate remain separate.
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

TASK = "KOKUNO-A3-INNER-LEADING-OSCILLATORY-FULL-ADVECTION-MEAN-070"
SCHEMA = "kokuno-a3-inner-leading-oscillatory-full-advection-mean-v1"
PARENT_AGENT3_PR = 822
PARENT_AGENT3_HEAD = "16e89b139c14c45d783f5b843d9f6748945a2b3e"

AGENT2_PR = 830
AGENT2_HEAD = "3479bed94ac6ef5c67547227f81a224ad43fb409"
AGENT2_MODULE = (
    "openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_full_advection"
)
AGENT2_FUNCTION = "evaluate_inner_leading_oscillatory_full_advection"
AGENT2_SOURCE_BLOB_SHA = "e44ccbe993d08951f48cfb5b093ee02ea5642dcf"
AGENT2_INNER_LEADING_FD4_STEP = 1.0e-3
AGENT2_OSCILLATORY_FD6_STEP = 1.0e-3

AGENT1_SELF_ADVECTION_PR = 829
AGENT1_SELF_ADVECTION_HEAD = "177190d770e0306e116674c5a2e1cbf977628760"
AGENT1_SELF_ADVECTION_BLOB_SHA = "b75f3b08ef3f1853550d8dfdef912980c15287c3"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

AdvectionEvaluator = Callable[..., Any]


@dataclass(frozen=True)
class ExactAgent2InnerLeadingOscillatoryFullAdvectionBackend:
    """Bind exact Agent-2 #830 full-convective implementation bytes."""

    evaluator: AdvectionEvaluator
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(
        cls, evaluator: AdvectionEvaluator
    ) -> "ExactAgent2InnerLeadingOscillatoryFullAdvectionBackend":
        if not callable(evaluator):
            raise TypeError("Agent-2 full-advection evaluator must be callable")
        if getattr(evaluator, "__module__", None) != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 full-advection module identity")
        if getattr(evaluator, "__name__", None) != AGENT2_FUNCTION:
            raise ValueError("unexpected Agent-2 full-advection function identity")
        source = inspect.getsourcefile(evaluator)
        if source is None:
            raise ValueError("Agent-2 full-advection source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 full-advection source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError("Agent-2 source blob does not match pinned PR #830")
        return cls(evaluator=evaluator, source_file=str(path), source_blob_sha=blob_sha)

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent2_pr": AGENT2_PR,
            "agent2_expected_head": AGENT2_HEAD,
            "agent1_self_advection_pr": AGENT1_SELF_ADVECTION_PR,
            "agent1_self_advection_expected_head": AGENT1_SELF_ADVECTION_HEAD,
            "agent1_self_advection_blob_sha": AGENT1_SELF_ADVECTION_BLOB_SHA,
            "module": AGENT2_MODULE,
            "function": AGENT2_FUNCTION,
            "source_blob_sha": self.source_blob_sha,
            "source_blob_matches_pinned": self.source_blob_sha == AGENT2_SOURCE_BLOB_SHA,
        }


@dataclass(frozen=True)
class InnerLeadingOscillatoryFullAdvectionMeanWitness:
    """Typed cylindrical m=0 witness for the strict-inner convective block."""

    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_inner_leading_self_advection_cylindrical: np.ndarray
    mean_inner_advects_oscillation_cylindrical: np.ndarray
    mean_oscillation_advects_inner_cylindrical: np.ndarray
    mean_mixed_cross_advection_cylindrical: np.ndarray
    mean_oscillatory_self_advection_cylindrical: np.ndarray
    mean_inner_plus_oscillatory_full_advection_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_total_advection_mean_relative_differences: tuple[float, ...]
    mixed_term_mean_closure_absolute_max: float
    full_decomposition_mean_closure_absolute_max: float
    total_advection_mean_rms: float
    full_ring_total_advection_rms: float
    total_advection_mean_to_full_rms_ratio: float
    backend: ExactAgent2InnerLeadingOscillatoryFullAdvectionBackend | None
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
            "mean_inner_leading_self_advection_cylindrical": np.asarray(
                self.mean_inner_leading_self_advection_cylindrical, dtype=float
            ).tolist(),
            "mean_inner_advects_oscillation_cylindrical": np.asarray(
                self.mean_inner_advects_oscillation_cylindrical, dtype=float
            ).tolist(),
            "mean_oscillation_advects_inner_cylindrical": np.asarray(
                self.mean_oscillation_advects_inner_cylindrical, dtype=float
            ).tolist(),
            "mean_mixed_cross_advection_cylindrical": np.asarray(
                self.mean_mixed_cross_advection_cylindrical, dtype=float
            ).tolist(),
            "mean_oscillatory_self_advection_cylindrical": np.asarray(
                self.mean_oscillatory_self_advection_cylindrical, dtype=float
            ).tolist(),
            "mean_inner_plus_oscillatory_full_advection_cylindrical": np.asarray(
                self.mean_inner_plus_oscillatory_full_advection_cylindrical, dtype=float
            ).tolist(),
            "angular_orders": list(self.angular_orders),
            "successive_total_advection_mean_relative_differences": list(
                self.successive_total_advection_mean_relative_differences
            ),
            "mixed_term_mean_closure_absolute_max": self.mixed_term_mean_closure_absolute_max,
            "full_decomposition_mean_closure_absolute_max": (
                self.full_decomposition_mean_closure_absolute_max
            ),
            "total_advection_mean_rms": self.total_advection_mean_rms,
            "full_ring_total_advection_rms": self.full_ring_total_advection_rms,
            "total_advection_mean_to_full_rms_ratio": (
                self.total_advection_mean_to_full_rms_ratio
            ),
            "backend_kind": self.backend_kind,
            "inner_leading_backend_kind": self.inner_leading_backend_kind,
            "backend": None if self.backend is None else self.backend.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _validate_result_steps(result: Any) -> None:
    inner_step = float(getattr(result, "inner_leading_spatial_step", np.nan))
    osc_step = float(getattr(result, "oscillatory_spatial_step", np.nan))
    if inner_step != AGENT2_INNER_LEADING_FD4_STEP:
        raise ValueError("Agent-2 inner-leading FD4 step drifted from frozen 1e-3")
    if osc_step != AGENT2_OSCILLATORY_FD6_STEP:
        raise ValueError("Agent-2 oscillatory FD6 step drifted from frozen 1e-3")


def _mean_for_order(
    evaluator: AdvectionEvaluator,
    inner_leading_backend: object,
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

    names = (
        "inner_leading_self_advection",
        "inner_advects_oscillation",
        "oscillation_advects_inner",
        "mixed_cross_advection",
        "oscillatory_self_advection",
        "inner_plus_oscillatory_self_advection",
    )
    if not all(hasattr(result, name) for name in names):
        raise TypeError("Agent-2 result does not expose the frozen full-advection decomposition")

    expected_shape = radius.shape + (order, 3)
    values = [np.asarray(getattr(result, name), dtype=float) for name in names]
    if any(value.shape != expected_shape for value in values):
        raise ValueError("Agent-2 full-advection payload has an unexpected ring shape")
    if not all(np.all(np.isfinite(value)) for value in values):
        raise ValueError("Agent-2 full-advection payload contains non-finite values")

    inner_self, inner_osc, osc_inner, mixed, osc_self, total = values
    scale = max(float(np.max(np.abs(total))), 1.0)
    mixed_pointwise = float(np.max(np.abs(mixed - inner_osc - osc_inner)))
    full_pointwise = float(np.max(np.abs(total - inner_self - mixed - osc_self)))
    if mixed_pointwise > 1.0e-12 * scale:
        raise ValueError("Agent-2 mixed-advection closure failed")
    if full_pointwise > 1.0e-12 * scale:
        raise ValueError("Agent-2 full-advection decomposition closure failed")

    mean_inner_self = _project_cartesian_ring_to_cylindrical_mean(inner_self, c, s)
    mean_inner_osc = _project_cartesian_ring_to_cylindrical_mean(inner_osc, c, s)
    mean_osc_inner = _project_cartesian_ring_to_cylindrical_mean(osc_inner, c, s)
    mean_mixed = _project_cartesian_ring_to_cylindrical_mean(mixed, c, s)
    mean_osc_self = _project_cartesian_ring_to_cylindrical_mean(osc_self, c, s)
    mean_total = _project_cartesian_ring_to_cylindrical_mean(total, c, s)

    mixed_mean_closure = float(
        np.max(np.abs(mean_mixed - mean_inner_osc - mean_osc_inner))
    )
    full_mean_closure = float(
        np.max(np.abs(mean_total - mean_inner_self - mean_mixed - mean_osc_self))
    )
    full_ring_rms = np.sqrt(np.mean(np.sum(total * total, axis=-1), axis=-1))
    return (
        mean_inner_self,
        mean_inner_osc,
        mean_osc_inner,
        mean_mixed,
        mean_osc_self,
        mean_total,
        full_ring_rms,
        np.asarray(mixed_mean_closure),
        np.asarray(full_mean_closure),
    )


def _materialize_with_evaluator(
    evaluator: AdvectionEvaluator,
    inner_leading_backend: object,
    radius: Any,
    z: Any,
    t: Any,
    *,
    backend: ExactAgent2InnerLeadingOscillatoryFullAdvectionBackend | None,
    backend_kind: str,
) -> InnerLeadingOscillatoryFullAdvectionMeanWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    per_order: list[tuple[np.ndarray, ...]] = []
    for order in ANGULAR_ORDERS:
        per_order.append(
            _mean_for_order(evaluator, inner_leading_backend, rb, zb, tb, order)
        )

    totals = [entry[5] for entry in per_order]
    successive = tuple(
        _relative_rms_difference(totals[index + 1], totals[index])
        for index in range(len(totals) - 1)
    )
    final = per_order[-1]
    (
        mean_inner_self,
        mean_inner_osc,
        mean_osc_inner,
        mean_mixed,
        mean_osc_self,
        mean_total,
        full_rms,
        _,
        _,
    ) = final

    total_mean_rms = float(np.sqrt(np.mean(np.sum(mean_total * mean_total, axis=-1))))
    full_ring_rms = float(np.sqrt(np.mean(full_rms * full_rms)))
    ratio = total_mean_rms / max(full_ring_rms, np.finfo(float).tiny)
    mixed_closure = float(max(float(entry[7]) for entry in per_order))
    full_closure = float(max(float(entry[8]) for entry in per_order))
    inner_kind = (
        f"{type(inner_leading_backend).__module__}."
        f"{type(inner_leading_backend).__qualname__}"
    )

    return InnerLeadingOscillatoryFullAdvectionMeanWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_inner_leading_self_advection_cylindrical=np.array(mean_inner_self, copy=True),
        mean_inner_advects_oscillation_cylindrical=np.array(mean_inner_osc, copy=True),
        mean_oscillation_advects_inner_cylindrical=np.array(mean_osc_inner, copy=True),
        mean_mixed_cross_advection_cylindrical=np.array(mean_mixed, copy=True),
        mean_oscillatory_self_advection_cylindrical=np.array(mean_osc_self, copy=True),
        mean_inner_plus_oscillatory_full_advection_cylindrical=np.array(
            mean_total, copy=True
        ),
        angular_orders=ANGULAR_ORDERS,
        successive_total_advection_mean_relative_differences=successive,
        mixed_term_mean_closure_absolute_max=mixed_closure,
        full_decomposition_mean_closure_absolute_max=full_closure,
        total_advection_mean_rms=total_mean_rms,
        full_ring_total_advection_rms=full_ring_rms,
        total_advection_mean_to_full_rms_ratio=float(ratio),
        backend=backend,
        backend_kind=backend_kind,
        inner_leading_backend_kind=inner_kind,
    )


def materialize_agent2_inner_leading_oscillatory_full_advection_mean(
    backend: ExactAgent2InnerLeadingOscillatoryFullAdvectionBackend,
    inner_leading_backend: object,
    radius: Any,
    z: Any,
    t: Any,
) -> InnerLeadingOscillatoryFullAdvectionMeanWitness:
    """Project exact A2 #830 strict-inner full advection onto cylindrical m=0."""
    if not isinstance(backend, ExactAgent2InnerLeadingOscillatoryFullAdvectionBackend):
        raise TypeError(
            "backend must be ExactAgent2InnerLeadingOscillatoryFullAdvectionBackend"
        )
    if backend.source_blob_sha != AGENT2_SOURCE_BLOB_SHA:
        raise ValueError("Agent-2 full-advection backend provenance drifted")
    for method_name in ("velocity", "self_advection"):
        if not hasattr(inner_leading_backend, method_name) or not callable(
            getattr(inner_leading_backend, method_name)
        ):
            raise TypeError(
                f"inner_leading_backend must expose {method_name}(x,y,z,t)"
            )
    return _materialize_with_evaluator(
        backend.evaluator,
        inner_leading_backend,
        radius,
        z,
        t,
        backend=backend,
        backend_kind="exact_agent2_pr830_inner_leading_oscillatory_full_advection",
    )


def truth_boundary() -> dict[str, object]:
    return {
        "agent2_full_advection_reimplemented_by_agent3": False,
        "cylindrical_m0_strict_inner_full_advection_projection_executable": True,
        "projection_uses_naive_cartesian_vector_mean": False,
        "uses_agent1_pa10_inner_contraction_center_only": True,
        "inner_leading_is_final_corrected_fixed_point": False,
        "outer_global_leading_join_materialized": False,
        "global_leading_velocity_materialized": False,
        "strict_inner_full_advection_mean_is_complete_ns_defect": False,
        "inner_leading_self_advection_included_here": True,
        "inner_leading_oscillatory_cross_advection_included_here": True,
        "oscillatory_self_advection_mean_included_here": True,
        "inner_plus_oscillatory_velocity_dt_mean_included_here": False,
        "base_viscous_term_included_here": False,
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
        "openai_field_identified": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


def public_contract() -> dict[str, object]:
    signature = inspect.signature(
        materialize_agent2_inner_leading_oscillatory_full_advection_mean
    )
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "normalized_score", "alpha", "damping", "angular_order",
        "derivative_step", "spatial_step", "nu", "viscosity", "delta_y", "delta_a",
        "scientific_threshold",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_pr": AGENT2_PR,
        "agent2_expected_head": AGENT2_HEAD,
        "agent2_source_blob_sha": AGENT2_SOURCE_BLOB_SHA,
        "agent1_self_advection_pr": AGENT1_SELF_ADVECTION_PR,
        "agent1_self_advection_expected_head": AGENT1_SELF_ADVECTION_HEAD,
        "agent1_self_advection_blob_sha": AGENT1_SELF_ADVECTION_BLOB_SHA,
        "source_reader": {
            "repository": SOURCE_READER_REPO,
            "head": SOURCE_READER_HEAD,
            "date": SOURCE_READER_DATE,
        },
        "angular_orders": list(ANGULAR_ORDERS),
        "frozen_agent2_derivative_steps": {
            "inner_leading_fd4": AGENT2_INNER_LEADING_FD4_STEP,
            "oscillatory_fd6": AGENT2_OSCILLATORY_FD6_STEP,
        },
        "forbidden_public_inputs_present": sorted(
            forbidden.intersection(signature.parameters)
        ),
        "final_gates_unchanged": {
            "normalized_momentum_max_l2": FINAL_NORMALIZED_MOMENTUM_GATE,
            "normalized_divergence_max_l2": FINAL_NORMALIZED_DIVERGENCE_GATE,
        },
        "truth_boundary": truth_boundary(),
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

    def cyl(mean: tuple[float, float, float], phase: float) -> np.ndarray:
        radial = mean[0] + 0.03 * np.cos(2.0 * theta + phase)
        tangential = mean[1] + 0.02 * np.sin(2.0 * theta - phase)
        axial = mean[2] + 0.015 * np.cos(2.0 * theta + 0.5 * phase)
        return np.stack(
            [radial * c - tangential * s, radial * s + tangential * c, axial],
            axis=-1,
        )

    inner_self = cyl((0.16, -0.03, 0.12), 0.2)
    inner_osc = cyl((0.21, -0.06, 0.14), 0.4)
    osc_inner = cyl((-0.08, 0.02, 0.11), -0.5)
    osc_self = cyl((0.07, 0.04, -0.02), 0.8)
    mixed = inner_osc + osc_inner
    total = inner_self + mixed + osc_self
    return SimpleNamespace(
        inner_leading_self_advection=inner_self,
        inner_advects_oscillation=inner_osc,
        oscillation_advects_inner=osc_inner,
        mixed_cross_advection=mixed,
        oscillatory_self_advection=osc_self,
        inner_plus_oscillatory_self_advection=total,
        inner_leading_spatial_step=AGENT2_INNER_LEADING_FD4_STEP,
        oscillatory_spatial_step=AGENT2_OSCILLATORY_FD6_STEP,
    )


def analytic_regression_receipt() -> dict[str, object]:
    class DummyInner:
        def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
            xb, yb, zb, tb = np.broadcast_arrays(
                np.asarray(x, dtype=float), np.asarray(y, dtype=float),
                np.asarray(z, dtype=float), np.asarray(t, dtype=float),
            )
            return np.stack(
                (xb * 0.0 + 1.0, yb * 0.0, zb * 0.0 + tb * 0.0), axis=-1
            )

        def self_advection(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
            xb, yb, zb, tb = np.broadcast_arrays(
                np.asarray(x, dtype=float), np.asarray(y, dtype=float),
                np.asarray(z, dtype=float), np.asarray(t, dtype=float),
            )
            return np.stack(
                (0.0 * xb, 0.0 * yb, 0.0 * zb + 0.0 * tb), axis=-1
            )

    witness = _materialize_with_evaluator(
        _analytic_regression_evaluator,
        DummyInner(),
        radius=np.asarray([0.45, 0.80, 1.15]),
        z=np.asarray([-0.5, 0.0, 0.6]),
        t=np.asarray([0.46, 0.50, 0.54]),
        backend=None,
        backend_kind="analytic_projection_regression_only",
    )
    means = {
        "inner_self": np.asarray([0.16, -0.03, 0.12]),
        "inner_osc": np.asarray([0.21, -0.06, 0.14]),
        "osc_inner": np.asarray([-0.08, 0.02, 0.11]),
        "osc_self": np.asarray([0.07, 0.04, -0.02]),
    }
    means["mixed"] = means["inner_osc"] + means["osc_inner"]
    means["total"] = means["inner_self"] + means["mixed"] + means["osc_self"]
    expected = {key: np.broadcast_to(value, (3, 3)) for key, value in means.items()}
    observed = {
        "inner_self": witness.mean_inner_leading_self_advection_cylindrical,
        "inner_osc": witness.mean_inner_advects_oscillation_cylindrical,
        "osc_inner": witness.mean_oscillation_advects_inner_cylindrical,
        "mixed": witness.mean_mixed_cross_advection_cylindrical,
        "osc_self": witness.mean_oscillatory_self_advection_cylindrical,
        "total": witness.mean_inner_plus_oscillatory_full_advection_cylindrical,
    }
    error = float(max(np.max(np.abs(observed[key] - expected[key])) for key in expected))
    return {
        "schema": SCHEMA,
        "task": TASK,
        "analytic_projection_absolute_max_error": error,
        "registered_total_mean": means["total"].tolist(),
        "analytic_witness": witness.to_receipt(),
        "real_agent2_numeric_full_advection_mean_observed_in_this_receipt": False,
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
