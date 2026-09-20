"""Agent-3 cylindrical m=0 projection of A2 inner-leading/oscillatory time derivative.

Agent 2 PR #821 exposes the additive strict-inner time derivative

    d_t u_inner+osc = d_t u_inner + d_t u_osc

for the executable Agent-1 PA.10 inner contraction-center Cartesian velocity and
Agent-2's frozen oscillatory field. Agent 3 owns only the rotating cylindrical
mean used by later mean/radial bookkeeping. This witness is intentionally a
scoped time-derivative subterm: the inner field is not the corrected/global
leading field and pressure, forcing, base viscous terms and the full NS defect
are not assembled here.
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

TASK = "KOKUNO-A3-INNER-LEADING-OSCILLATORY-TIME-MEAN-069"
SCHEMA = "kokuno-a3-inner-leading-oscillatory-time-mean-v1"
PARENT_AGENT3_PR = 813
PARENT_AGENT3_HEAD = "617ea2e83b2ab68f02e052307a9067d91ae505cc"

AGENT2_PR = 821
AGENT2_HEAD = "110c7150df73114c482e9f1673ef3a928b95053e"
AGENT2_MODULE = (
    "openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_time_derivative"
)
AGENT2_FUNCTION = "evaluate_inner_leading_oscillatory_time_derivative"
AGENT2_SOURCE_BLOB_SHA = "c780aea6d55daa62734b7906fe4e9e4f9ef752f0"

AGENT1_PR = 811
AGENT1_HEAD = "6f16e837c3a95a6e54af5c2cb7725d094f667026"
AGENT1_DT_SOURCE_BLOB_SHA = "c1d576820416945443bbc24cc931cb64e9e032d9"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

TimeEvaluator = Callable[..., Any]


@dataclass(frozen=True)
class ExactAgent2InnerLeadingOscillatoryTimeBackend:
    """Bind the exact Agent-2 #821 additive time-derivative implementation bytes."""

    evaluator: TimeEvaluator
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(
        cls, evaluator: TimeEvaluator
    ) -> "ExactAgent2InnerLeadingOscillatoryTimeBackend":
        if not callable(evaluator):
            raise TypeError("Agent-2 inner-leading/oscillatory time evaluator must be callable")
        if getattr(evaluator, "__module__", None) != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 inner-leading/oscillatory time module identity")
        if getattr(evaluator, "__name__", None) != AGENT2_FUNCTION:
            raise ValueError("unexpected Agent-2 inner-leading/oscillatory time function identity")
        source = inspect.getsourcefile(evaluator)
        if source is None:
            raise ValueError("Agent-2 inner-leading/oscillatory time source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 inner-leading/oscillatory time source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError("Agent-2 source blob does not match pinned PR #821")
        return cls(evaluator=evaluator, source_file=str(path), source_blob_sha=blob_sha)

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent2_pr": AGENT2_PR,
            "agent2_expected_head": AGENT2_HEAD,
            "agent1_inner_velocity_dt_pr": AGENT1_PR,
            "agent1_inner_velocity_dt_expected_head": AGENT1_HEAD,
            "module": AGENT2_MODULE,
            "function": AGENT2_FUNCTION,
            "source_blob_sha": self.source_blob_sha,
            "source_blob_matches_pinned": self.source_blob_sha == AGENT2_SOURCE_BLOB_SHA,
        }


@dataclass(frozen=True)
class InnerLeadingOscillatoryTimeMeanWitness:
    """Typed cylindrical mean of strict-inner A1-leading + A2-osc velocity_dt."""

    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_inner_leading_velocity_dt_cylindrical: np.ndarray
    mean_oscillatory_velocity_dt_cylindrical: np.ndarray
    mean_inner_plus_oscillatory_velocity_dt_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_total_dt_mean_relative_differences: tuple[float, ...]
    additive_mean_closure_absolute_max: float
    total_dt_mean_rms: float
    full_ring_total_dt_rms: float
    total_dt_mean_to_full_rms_ratio: float
    backend: ExactAgent2InnerLeadingOscillatoryTimeBackend | None
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
            "mean_inner_leading_velocity_dt_cylindrical": np.asarray(
                self.mean_inner_leading_velocity_dt_cylindrical, dtype=float
            ).tolist(),
            "mean_oscillatory_velocity_dt_cylindrical": np.asarray(
                self.mean_oscillatory_velocity_dt_cylindrical, dtype=float
            ).tolist(),
            "mean_inner_plus_oscillatory_velocity_dt_cylindrical": np.asarray(
                self.mean_inner_plus_oscillatory_velocity_dt_cylindrical, dtype=float
            ).tolist(),
            "angular_orders": list(self.angular_orders),
            "successive_total_dt_mean_relative_differences": list(
                self.successive_total_dt_mean_relative_differences
            ),
            "additive_mean_closure_absolute_max": self.additive_mean_closure_absolute_max,
            "total_dt_mean_rms": self.total_dt_mean_rms,
            "full_ring_total_dt_rms": self.full_ring_total_dt_rms,
            "total_dt_mean_to_full_rms_ratio": self.total_dt_mean_to_full_rms_ratio,
            "backend_kind": self.backend_kind,
            "inner_leading_backend_kind": self.inner_leading_backend_kind,
            "backend": None if self.backend is None else self.backend.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _mean_for_order(
    evaluator: TimeEvaluator,
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

    result = evaluator(inner_leading_backend, xq, yq, zq, tq)
    required = (
        "inner_leading_velocity_dt",
        "oscillatory_velocity_dt",
        "inner_plus_oscillatory_velocity_dt",
    )
    if not all(hasattr(result, key) for key in required):
        raise TypeError("Agent-2 result must expose inner, oscillatory and summed velocity_dt")

    expected_shape = radius.shape + (order, 3)
    dt_inner = np.asarray(getattr(result, required[0]), dtype=float)
    dt_osc = np.asarray(getattr(result, required[1]), dtype=float)
    dt_total = np.asarray(getattr(result, required[2]), dtype=float)
    if (
        dt_inner.shape != expected_shape
        or dt_osc.shape != expected_shape
        or dt_total.shape != expected_shape
    ):
        raise ValueError("Agent-2 inner-leading/oscillatory time payload has an unexpected ring shape")
    if not all(np.all(np.isfinite(a)) for a in (dt_inner, dt_osc, dt_total)):
        raise ValueError("Agent-2 inner-leading/oscillatory time payload contains non-finite values")

    scale = max(float(np.max(np.abs(dt_total))), 1.0)
    pointwise_closure = float(np.max(np.abs(dt_total - dt_inner - dt_osc)))
    if pointwise_closure > 1.0e-12 * scale:
        raise ValueError("Agent-2 inner-leading/oscillatory velocity_dt additive closure failed")

    mean_inner = _project_cartesian_ring_to_cylindrical_mean(dt_inner, c, s)
    mean_osc = _project_cartesian_ring_to_cylindrical_mean(dt_osc, c, s)
    mean_total = _project_cartesian_ring_to_cylindrical_mean(dt_total, c, s)
    mean_closure = float(np.max(np.abs(mean_total - mean_inner - mean_osc)))
    full_ring_rms = np.sqrt(np.mean(np.sum(dt_total * dt_total, axis=-1), axis=-1))
    return mean_inner, mean_osc, mean_total, full_ring_rms, mean_closure


def _materialize_with_evaluator(
    evaluator: TimeEvaluator,
    inner_leading_backend: object,
    radius: Any,
    z: Any,
    t: Any,
    *,
    backend: ExactAgent2InnerLeadingOscillatoryTimeBackend | None,
    backend_kind: str,
) -> InnerLeadingOscillatoryTimeMeanWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    means_inner: list[np.ndarray] = []
    means_osc: list[np.ndarray] = []
    means_total: list[np.ndarray] = []
    full_rms: np.ndarray | None = None
    mean_closures: list[float] = []

    for order in ANGULAR_ORDERS:
        mean_inner, mean_osc, mean_total, ring_rms, mean_closure = _mean_for_order(
            evaluator, inner_leading_backend, rb, zb, tb, order
        )
        means_inner.append(mean_inner)
        means_osc.append(mean_osc)
        means_total.append(mean_total)
        full_rms = ring_rms
        mean_closures.append(mean_closure)

    successive = tuple(
        _relative_rms_difference(means_total[index + 1], means_total[index])
        for index in range(len(means_total) - 1)
    )
    final_inner = means_inner[-1]
    final_osc = means_osc[-1]
    final_total = means_total[-1]
    assert full_rms is not None
    total_dt_mean_rms = float(np.sqrt(np.mean(np.sum(final_total * final_total, axis=-1))))
    full_ring_total_dt_rms = float(np.sqrt(np.mean(full_rms * full_rms)))
    ratio = total_dt_mean_rms / max(full_ring_total_dt_rms, np.finfo(float).tiny)
    inner_kind = (
        f"{type(inner_leading_backend).__module__}."
        f"{type(inner_leading_backend).__qualname__}"
    )

    return InnerLeadingOscillatoryTimeMeanWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_inner_leading_velocity_dt_cylindrical=np.array(final_inner, copy=True),
        mean_oscillatory_velocity_dt_cylindrical=np.array(final_osc, copy=True),
        mean_inner_plus_oscillatory_velocity_dt_cylindrical=np.array(final_total, copy=True),
        angular_orders=ANGULAR_ORDERS,
        successive_total_dt_mean_relative_differences=successive,
        additive_mean_closure_absolute_max=float(max(mean_closures)),
        total_dt_mean_rms=total_dt_mean_rms,
        full_ring_total_dt_rms=full_ring_total_dt_rms,
        total_dt_mean_to_full_rms_ratio=float(ratio),
        backend=backend,
        backend_kind=backend_kind,
        inner_leading_backend_kind=inner_kind,
    )


def materialize_agent2_inner_leading_oscillatory_time_mean(
    backend: ExactAgent2InnerLeadingOscillatoryTimeBackend,
    inner_leading_backend: object,
    radius: Any,
    z: Any,
    t: Any,
) -> InnerLeadingOscillatoryTimeMeanWitness:
    """Project pinned A2 #821 strict-inner leading+oscillatory velocity_dt to m=0.

    The caller cannot supply residual/defect values, mean values, radial stress,
    pressure, forcing, target/gain/damping, derivative resolution, angular order,
    viscosity, raw ``delta_y``/``delta_a`` or scientific gates.
    """
    if not isinstance(backend, ExactAgent2InnerLeadingOscillatoryTimeBackend):
        raise TypeError("backend must be ExactAgent2InnerLeadingOscillatoryTimeBackend")
    if backend.source_blob_sha != AGENT2_SOURCE_BLOB_SHA:
        raise ValueError("Agent-2 inner-leading/oscillatory time backend provenance drifted")
    for method_name in ("velocity", "velocity_dt"):
        if not hasattr(inner_leading_backend, method_name) or not callable(
            getattr(inner_leading_backend, method_name)
        ):
            raise TypeError(f"inner_leading_backend must expose {method_name}(x,y,z,t)")
    return _materialize_with_evaluator(
        backend.evaluator,
        inner_leading_backend,
        radius,
        z,
        t,
        backend=backend,
        backend_kind="exact_agent2_pr821_inner_leading_oscillatory_time_derivative",
    )


def truth_boundary() -> dict[str, object]:
    return {
        "agent2_inner_leading_oscillatory_time_reimplemented_by_agent3": False,
        "cylindrical_m0_inner_leading_oscillatory_velocity_dt_projection_executable": True,
        "projection_uses_naive_cartesian_vector_mean": False,
        "uses_agent1_pa10_inner_contraction_center_only": True,
        "inner_leading_is_final_corrected_fixed_point": False,
        "outer_global_leading_join_materialized": False,
        "global_leading_velocity_materialized": False,
        "inner_leading_oscillatory_velocity_dt_mean_is_complete_ns_defect": False,
        "inner_leading_self_advection_included_here": False,
        "inner_leading_oscillatory_cross_advection_included_here": False,
        "oscillatory_self_advection_mean_included_here": False,
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
    }


def public_contract() -> dict[str, object]:
    signature = inspect.signature(materialize_agent2_inner_leading_oscillatory_time_mean)
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure", "forcing",
        "target", "gain", "normalized_score", "alpha", "damping", "angular_order",
        "derivative_step", "time_step", "nu", "viscosity", "delta_y", "delta_a",
        "scientific_threshold",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_pr": AGENT2_PR,
        "agent2_head": AGENT2_HEAD,
        "agent2_source_blob_sha": AGENT2_SOURCE_BLOB_SHA,
        "agent1_pr": AGENT1_PR,
        "agent1_head": AGENT1_HEAD,
        "agent1_velocity_dt_source_blob_sha": AGENT1_DT_SOURCE_BLOB_SHA,
        "source_provenance": {
            "repository": SOURCE_READER_REPO,
            "head": SOURCE_READER_HEAD,
            "date": SOURCE_READER_DATE,
        },
        "angular_orders": list(ANGULAR_ORDERS),
        "forbidden_public_inputs_present": sorted(forbidden.intersection(signature.parameters)),
        "final_gates_unchanged": {
            "normalized_momentum_max_l2": FINAL_NORMALIZED_MOMENTUM_GATE,
            "normalized_divergence_max_l2": FINAL_NORMALIZED_DIVERGENCE_GATE,
        },
        "truth_boundary": truth_boundary(),
    }


def _mechanics_fixture() -> dict[str, object]:
    class DummyInner:
        def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
            x, y, z, t = np.broadcast_arrays(
                np.asarray(x, dtype=float), np.asarray(y, dtype=float),
                np.asarray(z, dtype=float), np.asarray(t, dtype=float)
            )
            return np.stack([0.2 + 0.0*x, -0.1 + 0.0*y, 0.3 + 0.0*z], axis=-1)

        def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
            x, y, z, t = np.broadcast_arrays(
                np.asarray(x, dtype=float), np.asarray(y, dtype=float),
                np.asarray(z, dtype=float), np.asarray(t, dtype=float)
            )
            theta = np.arctan2(y, x)
            c = np.cos(theta)
            s = np.sin(theta)
            radial = 0.09 + 0.03*np.cos(2.0*theta)
            tangential = -0.03 + 0.02*np.sin(2.0*theta)
            axial = 0.04 + 0.01*np.cos(2.0*theta)
            return np.stack([radial*c - tangential*s, radial*s + tangential*c, axial], axis=-1)

    def evaluator(inner: object, x: Any, y: Any, z: Any, t: Any) -> Any:
        dt_inner = np.asarray(inner.velocity_dt(x, y, z, t), dtype=float)
        theta = np.arctan2(np.asarray(y, dtype=float), np.asarray(x, dtype=float))
        c = np.cos(theta)
        s = np.sin(theta)
        radial = -0.02 + 0.015*np.sin(2.0*theta)
        tangential = 0.05 + 0.01*np.cos(2.0*theta)
        axial = 0.07 - 0.02*np.sin(2.0*theta)
        dt_osc = np.stack([radial*c - tangential*s, radial*s + tangential*c, axial], axis=-1)
        return type("TimeResult", (), {
            "inner_leading_velocity_dt": dt_inner,
            "oscillatory_velocity_dt": dt_osc,
            "inner_plus_oscillatory_velocity_dt": dt_inner + dt_osc,
        })()

    witness = _materialize_with_evaluator(
        evaluator,
        DummyInner(),
        np.asarray([0.7, 1.1], dtype=float),
        np.asarray([-0.2, 0.25], dtype=float),
        np.asarray([0.47, 0.53], dtype=float),
        backend=None,
        backend_kind="mechanics_fixture",
    )
    expected_inner = np.asarray([0.09, -0.03, 0.04], dtype=float)
    expected_osc = np.asarray([-0.02, 0.05, 0.07], dtype=float)
    expected_total = expected_inner + expected_osc
    error = float(np.max(np.abs(
        witness.mean_inner_plus_oscillatory_velocity_dt_cylindrical - expected_total
    )))
    return {
        "expected_inner_mean": expected_inner.tolist(),
        "expected_oscillatory_mean": expected_osc.tolist(),
        "expected_total_mean": expected_total.tolist(),
        "recovered_total_mean": witness.mean_inner_plus_oscillatory_velocity_dt_cylindrical.tolist(),
        "max_abs_error": error,
        "additive_mean_closure_absolute_max": witness.additive_mean_closure_absolute_max,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    payload = {
        "contract": public_contract(),
        "mechanics_fixture": _mechanics_fixture(),
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
