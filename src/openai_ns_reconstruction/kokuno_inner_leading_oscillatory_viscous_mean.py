"""Agent-3 cylindrical m=0 projection of the strict-inner base viscous block.

This module consumes, without reimplementing their derivative physics,

* Agent 1 PR #839: ``KokunoPA10CartesianCenterLaplacian.velocity_laplacian``
  for the current PA.10 inner contraction-center velocity; and
* Agent 2 PR #708: ``evaluate_velocity_laplacian_fd6`` for the frozen
  oscillatory velocity.

It forms only

    -nu (Delta u_inner + Delta u_osc),   nu = 0.01,

and projects the two summands and their sum onto the rotating cylindrical
m=0 frame.  This is the missing viscous ingredient for the *strict-inner*
transport bookkeeping assembled by Agent 3 #822/#831.  It is not a complete
Navier--Stokes defect: the final corrected/global leading join, matched
pressure, restricted forcing and real correction cycle remain separate.
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

TASK = "KOKUNO-A3-INNER-LEADING-OSCILLATORY-VISCOUS-MEAN-071"
SCHEMA = "kokuno-a3-inner-leading-oscillatory-viscous-mean-v1"
PARENT_AGENT3_PR = 831
PARENT_AGENT3_HEAD = "5a04f14a13aef019cf83fa03bc597c8e68fba11f"

AGENT1_PR = 839
AGENT1_HEAD = "e96ee90144976a992b62d76d4361a94eb16bd91e"
AGENT1_MODULE = "openai_ns_reconstruction.kokuno_pa10_cartesian_center_laplacian"
AGENT1_CLASS = "KokunoPA10CartesianCenterLaplacian"
AGENT1_SOURCE_BLOB_SHA = "74b0e185e8e0bbb08695d493e0a091c305a03006"
AGENT1_LAPLACIAN_STEP = 1.0e-3

AGENT2_PR = 708
AGENT2_HEAD = "9f8bae37c4d3b2559bee6db050655702e4c058e9"
AGENT2_MODULE = "openai_ns_reconstruction.kokuno_public_oscillatory_spatial_laplacian"
AGENT2_FUNCTION = "evaluate_velocity_laplacian_fd6"
AGENT2_SOURCE_BLOB_SHA = "7dfd3dabd28147dffac6ff77ec5457c17689bfe1"
# Inherit the finest preregistered #708 Laplacian step.  It is frozen here and
# is not a caller knob.
AGENT2_LAPLACIAN_STEP = 4.5e-3

REPOSITORY_VISCOSITY = 0.01
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

OscillatoryLaplacianEvaluator = Callable[..., Any]
LaplacianProvider = Callable[[Any, Any, Any, Any], np.ndarray]


@dataclass(frozen=True)
class ExactAgent1InnerLaplacianBackend:
    """Bind the exact Agent-1 #839 inner-center Laplacian implementation."""

    candidate: object
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(cls, candidate: object) -> "ExactAgent1InnerLaplacianBackend":
        candidate_type = type(candidate)
        if candidate_type.__module__ != AGENT1_MODULE:
            raise ValueError("unexpected Agent-1 Laplacian module identity")
        if candidate_type.__name__ != AGENT1_CLASS:
            raise ValueError("unexpected Agent-1 Laplacian class identity")
        velocity_laplacian = getattr(candidate, "velocity_laplacian", None)
        if not callable(velocity_laplacian):
            raise TypeError("Agent-1 backend must expose velocity_laplacian")
        source = inspect.getsourcefile(candidate_type)
        if source is None:
            raise ValueError("Agent-1 Laplacian source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-1 Laplacian source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT1_SOURCE_BLOB_SHA:
            raise ValueError("Agent-1 source blob does not match pinned PR #839")
        return cls(candidate=candidate, source_file=str(path), source_blob_sha=blob_sha)

    def laplacian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        value = np.asarray(self.candidate.velocity_laplacian(x, y, z, t), dtype=float)
        return value

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent1_pr": AGENT1_PR,
            "agent1_expected_head": AGENT1_HEAD,
            "module": AGENT1_MODULE,
            "class": AGENT1_CLASS,
            "source_blob_sha": self.source_blob_sha,
            "source_blob_matches_pinned": self.source_blob_sha == AGENT1_SOURCE_BLOB_SHA,
            "production_spatial_step": AGENT1_LAPLACIAN_STEP,
        }


@dataclass(frozen=True)
class ExactAgent2OscillatoryLaplacianBackend:
    """Bind the exact Agent-2 #708 oscillatory Laplacian evaluator."""

    evaluator: OscillatoryLaplacianEvaluator
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(
        cls, evaluator: OscillatoryLaplacianEvaluator
    ) -> "ExactAgent2OscillatoryLaplacianBackend":
        if not callable(evaluator):
            raise TypeError("Agent-2 oscillatory Laplacian evaluator must be callable")
        if getattr(evaluator, "__module__", None) != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 oscillatory Laplacian module identity")
        if getattr(evaluator, "__name__", None) != AGENT2_FUNCTION:
            raise ValueError("unexpected Agent-2 oscillatory Laplacian function identity")
        source = inspect.getsourcefile(evaluator)
        if source is None:
            raise ValueError("Agent-2 oscillatory Laplacian source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 oscillatory Laplacian source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError("Agent-2 source blob does not match pinned PR #708")
        return cls(evaluator=evaluator, source_file=str(path), source_blob_sha=blob_sha)

    def laplacian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        result = self.evaluator(
            x,
            y,
            z,
            t,
            spatial_step=AGENT2_LAPLACIAN_STEP,
        )
        if not isinstance(result, dict):
            raise TypeError("Agent-2 Laplacian evaluator must return a dictionary")
        step = float(result.get("spatial_step", np.nan))
        if step != AGENT2_LAPLACIAN_STEP:
            raise ValueError("Agent-2 oscillatory Laplacian step drifted")
        if result.get("spatial_operator") != "centered_cartesian_fd6_second_derivative":
            raise ValueError("Agent-2 oscillatory spatial operator identity drifted")
        return np.asarray(result.get("laplacian"), dtype=float)

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent2_pr": AGENT2_PR,
            "agent2_expected_head": AGENT2_HEAD,
            "module": AGENT2_MODULE,
            "function": AGENT2_FUNCTION,
            "source_blob_sha": self.source_blob_sha,
            "source_blob_matches_pinned": self.source_blob_sha == AGENT2_SOURCE_BLOB_SHA,
            "spatial_step": AGENT2_LAPLACIAN_STEP,
        }


@dataclass(frozen=True)
class InnerLeadingOscillatoryViscousMeanWitness:
    """Typed cylindrical mean witness for ``-nu Delta(u_inner+u_osc)``."""

    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_inner_viscous_cylindrical: np.ndarray
    mean_oscillatory_viscous_cylindrical: np.ndarray
    mean_total_viscous_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_total_viscous_mean_relative_differences: tuple[float, ...]
    projected_additive_closure_absolute_max: float
    total_viscous_mean_rms: float
    full_ring_total_viscous_rms: float
    total_viscous_mean_to_full_rms_ratio: float
    viscosity: float
    inner_backend: ExactAgent1InnerLaplacianBackend | None
    oscillatory_backend: ExactAgent2OscillatoryLaplacianBackend | None
    backend_kind: str

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "radius": np.asarray(self.radius, dtype=float).tolist(),
            "z": np.asarray(self.z, dtype=float).tolist(),
            "t": np.asarray(self.t, dtype=float).tolist(),
            "cylindrical_components": ["radial", "tangential", "axial"],
            "mean_inner_viscous_cylindrical": np.asarray(
                self.mean_inner_viscous_cylindrical, dtype=float
            ).tolist(),
            "mean_oscillatory_viscous_cylindrical": np.asarray(
                self.mean_oscillatory_viscous_cylindrical, dtype=float
            ).tolist(),
            "mean_total_viscous_cylindrical": np.asarray(
                self.mean_total_viscous_cylindrical, dtype=float
            ).tolist(),
            "angular_orders": list(self.angular_orders),
            "successive_total_viscous_mean_relative_differences": list(
                self.successive_total_viscous_mean_relative_differences
            ),
            "projected_additive_closure_absolute_max": (
                self.projected_additive_closure_absolute_max
            ),
            "total_viscous_mean_rms": self.total_viscous_mean_rms,
            "full_ring_total_viscous_rms": self.full_ring_total_viscous_rms,
            "total_viscous_mean_to_full_rms_ratio": (
                self.total_viscous_mean_to_full_rms_ratio
            ),
            "viscosity": self.viscosity,
            "backend_kind": self.backend_kind,
            "inner_backend": None if self.inner_backend is None else self.inner_backend.to_receipt(),
            "oscillatory_backend": (
                None if self.oscillatory_backend is None else self.oscillatory_backend.to_receipt()
            ),
            "truth_boundary": truth_boundary(),
        }


def _evaluate_provider(
    provider: LaplacianProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    expected_shape: tuple[int, ...],
    *,
    label: str,
) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    if value.shape != expected_shape:
        raise ValueError(f"{label} returned shape {value.shape}, expected {expected_shape}")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{label} returned non-finite values")
    return value


def _mean_for_order(
    inner_laplacian_provider: LaplacianProvider,
    oscillatory_laplacian_provider: LaplacianProvider,
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
    expected = radius.shape + (order, 3)

    inner_laplacian = _evaluate_provider(
        inner_laplacian_provider, xq, yq, zq, tq, expected, label="inner Laplacian"
    )
    oscillatory_laplacian = _evaluate_provider(
        oscillatory_laplacian_provider,
        xq,
        yq,
        zq,
        tq,
        expected,
        label="oscillatory Laplacian",
    )

    inner_viscous = -REPOSITORY_VISCOSITY * inner_laplacian
    oscillatory_viscous = -REPOSITORY_VISCOSITY * oscillatory_laplacian
    total_viscous = inner_viscous + oscillatory_viscous

    mean_inner = _project_cartesian_ring_to_cylindrical_mean(inner_viscous, c, s)
    mean_osc = _project_cartesian_ring_to_cylindrical_mean(oscillatory_viscous, c, s)
    mean_total = _project_cartesian_ring_to_cylindrical_mean(total_viscous, c, s)
    closure = float(np.max(np.abs(mean_total - mean_inner - mean_osc)))
    full_ring_rms = float(
        np.sqrt(np.mean(np.sum(total_viscous * total_viscous, axis=-1)))
    )
    return mean_inner, mean_osc, mean_total, np.asarray(full_ring_rms), closure


def _materialize_from_laplacian_providers(
    inner_laplacian_provider: LaplacianProvider,
    oscillatory_laplacian_provider: LaplacianProvider,
    radius: Any,
    z: Any,
    t: Any,
    *,
    inner_backend: ExactAgent1InnerLaplacianBackend | None,
    oscillatory_backend: ExactAgent2OscillatoryLaplacianBackend | None,
    backend_kind: str,
) -> InnerLeadingOscillatoryViscousMeanWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    per_order = [
        _mean_for_order(
            inner_laplacian_provider,
            oscillatory_laplacian_provider,
            rb,
            zb,
            tb,
            order,
        )
        for order in ANGULAR_ORDERS
    ]
    totals = [entry[2] for entry in per_order]
    successive = tuple(
        _relative_rms_difference(totals[index + 1], totals[index])
        for index in range(len(totals) - 1)
    )
    mean_inner, mean_osc, mean_total, full_ring_rms_array, _ = per_order[-1]
    total_mean_rms = float(np.sqrt(np.mean(np.sum(mean_total * mean_total, axis=-1))))
    full_ring_total_rms = float(full_ring_rms_array)
    ratio = total_mean_rms / max(full_ring_total_rms, np.finfo(float).tiny)
    closure = float(max(entry[4] for entry in per_order))

    return InnerLeadingOscillatoryViscousMeanWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_inner_viscous_cylindrical=np.array(mean_inner, copy=True),
        mean_oscillatory_viscous_cylindrical=np.array(mean_osc, copy=True),
        mean_total_viscous_cylindrical=np.array(mean_total, copy=True),
        angular_orders=ANGULAR_ORDERS,
        successive_total_viscous_mean_relative_differences=successive,
        projected_additive_closure_absolute_max=closure,
        total_viscous_mean_rms=total_mean_rms,
        full_ring_total_viscous_rms=full_ring_total_rms,
        total_viscous_mean_to_full_rms_ratio=float(ratio),
        viscosity=REPOSITORY_VISCOSITY,
        inner_backend=inner_backend,
        oscillatory_backend=oscillatory_backend,
        backend_kind=backend_kind,
    )


def materialize_inner_leading_oscillatory_viscous_mean(
    inner_backend: ExactAgent1InnerLaplacianBackend,
    oscillatory_backend: ExactAgent2OscillatoryLaplacianBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> InnerLeadingOscillatoryViscousMeanWitness:
    """Project exact A1/A2 base viscous terms onto the frozen A3 mean ladder."""
    if not isinstance(inner_backend, ExactAgent1InnerLaplacianBackend):
        raise TypeError("inner_backend must be ExactAgent1InnerLaplacianBackend")
    if not isinstance(oscillatory_backend, ExactAgent2OscillatoryLaplacianBackend):
        raise TypeError("oscillatory_backend must be ExactAgent2OscillatoryLaplacianBackend")
    return _materialize_from_laplacian_providers(
        inner_backend.laplacian,
        oscillatory_backend.laplacian,
        radius,
        z,
        t,
        inner_backend=inner_backend,
        oscillatory_backend=oscillatory_backend,
        backend_kind="exact-agent1-839-plus-agent2-708",
    )


def truth_boundary() -> dict[str, object]:
    return {
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "repository_viscosity": REPOSITORY_VISCOSITY,
        "viscosity_caller_tunable": False,
        "agent1_inner_laplacian_head": AGENT1_HEAD,
        "agent2_oscillatory_laplacian_head": AGENT2_HEAD,
        "inner_center_scoped_only": True,
        "inner_plus_oscillatory_base_viscous_term_included": True,
        "time_derivative_included_here": False,
        "convective_term_included_here": False,
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "correction_transport_included": False,
        "complete_ns_defect": False,
        "source_center_is_final_corrected_fixed_point": False,
        "global_corrected_leading_join_materialized": False,
        "real_agent3_delta_a_bound": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


def _cylindrical_field(
    radial: float,
    tangential: float,
    axial: float,
    mode_scale: float,
) -> LaplacianProvider:
    """Manufactured Cartesian ring field with registered cylindrical m=0 mean."""
    def provider(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        xb, yb, zb, tb = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        del zb, tb
        theta = np.arctan2(yb, xb)
        c = np.cos(theta)
        s = np.sin(theta)
        m2 = mode_scale * np.cos(2.0 * theta)
        mr = radial + m2
        mt = tangential - 0.5 * m2
        mz = axial + 0.75 * mode_scale * np.sin(2.0 * theta)
        return np.stack((mr * c - mt * s, mr * s + mt * c, mz), axis=-1)
    return provider


def build_mechanics_report() -> dict[str, object]:
    inner = _cylindrical_field(1.0, -2.0, 3.0, 0.7)
    osc = _cylindrical_field(0.5, 1.0, -1.0, -0.4)
    witness = _materialize_from_laplacian_providers(
        inner,
        osc,
        np.array([0.7, 1.1]),
        np.array([-0.2, 0.3]),
        np.array([0.48, 0.52]),
        inner_backend=None,
        oscillatory_backend=None,
        backend_kind="manufactured-mechanics-only",
    )
    receipt = witness.to_receipt()
    receipt["mechanics_only"] = True
    receipt["expected_inner_viscous_mean"] = [-0.01, 0.02, -0.03]
    receipt["expected_oscillatory_viscous_mean"] = [-0.005, -0.01, 0.01]
    receipt["expected_total_viscous_mean"] = [-0.015, 0.01, -0.02]
    return receipt


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_mechanics_report()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
