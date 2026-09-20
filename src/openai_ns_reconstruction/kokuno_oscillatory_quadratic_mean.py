"""Agent-3 cylindrical mean projection of Agent-2 oscillatory self-advection.

This module owns only the mean-projection seam needed by the Kokuno mean/radial
correction lane.  It deliberately does *not* reimplement Agent-2's frozen
oscillatory velocity, spatial gradient, or quadratic contraction.

The exact Agent-2 public evaluator produces

    N_osc = (u_osc . grad) u_osc.

For each fixed cylindrical ring ``(r,z,t)``, Agent 3 needs the m=0 vector
coefficient in the rotating cylindrical frame, not the naive Cartesian vector
average.  With ``e_r=(cos(theta),sin(theta),0)`` and
``e_theta=(-sin(theta),cos(theta),0)`` this module computes

    <N_r>     = mean_theta N_osc . e_r,
    <N_theta> = mean_theta N_osc . e_theta,
    <N_z>     = mean_theta N_osc . e_z.

These three profiles are a real nonlinear oscillatory *subterm* that may feed a
future mean-defect/correction assembly.  They are not the complete Navier--Stokes
residual, do not include leading/cross/pressure/viscous/forcing terms, and must
not be used as a surrogate success metric for the final PDE gate.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np

TASK = "KOKUNO-A3-OSCILLATORY-QUADRATIC-MEAN-064"
SCHEMA = "kokuno-a3-oscillatory-quadratic-mean-v1"
PARENT_AGENT3_PR = 770
PARENT_AGENT3_HEAD = "076dd9e7ca1b53d439761d5b726d4028ad51d587"

AGENT2_PR = 779
AGENT2_HEAD = "a20878fb781fe75698ab22cfbe5c36a78f33ddc5"
AGENT2_MODULE = (
    "openai_ns_reconstruction.kokuno_public_oscillatory_self_advection"
)
AGENT2_FUNCTION = "evaluate_oscillatory_self_advection_fd6"
AGENT2_SOURCE_BLOB_SHA = "5a21ed184eaf6de0f506f709dede36a046e4f2e5"
AGENT2_SPATIAL_STEP = 1.0e-3

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

# Frozen numerical projection ladder.  The caller cannot retune it.
ANGULAR_ORDERS = (32, 64, 128)
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5

SelfAdvectionEvaluator = Callable[..., dict[str, Any]]


def _git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


@dataclass(frozen=True)
class ExactAgent2SelfAdvectionBackend:
    """Executable binding to the exact Agent-2 self-advection source bytes.

    The caller supplies only the imported callable.  Module/function identity and
    the Git blob SHA of its source file are checked here; head/provenance values
    are frozen constants rather than caller-supplied labels.
    """

    evaluator: SelfAdvectionEvaluator
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(cls, evaluator: SelfAdvectionEvaluator) -> "ExactAgent2SelfAdvectionBackend":
        if not callable(evaluator):
            raise TypeError("Agent-2 self-advection evaluator must be callable")
        if getattr(evaluator, "__module__", None) != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 self-advection module identity")
        if getattr(evaluator, "__name__", None) != AGENT2_FUNCTION:
            raise ValueError("unexpected Agent-2 self-advection function identity")
        source = inspect.getsourcefile(evaluator)
        if source is None:
            raise ValueError("Agent-2 evaluator source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 evaluator source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError("Agent-2 self-advection source blob does not match pinned PR #779")
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
class OscillatoryQuadraticMeanWitness:
    """Typed m=0 cylindrical projection of the frozen A2 nonlinear subterm."""

    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_mean_relative_differences: tuple[float, ...]
    mean_rms: float
    full_ring_self_advection_rms: float
    mean_to_full_rms_ratio: float
    backend: ExactAgent2SelfAdvectionBackend | None
    backend_kind: str

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "radius": np.asarray(self.radius, dtype=float).tolist(),
            "z": np.asarray(self.z, dtype=float).tolist(),
            "t": np.asarray(self.t, dtype=float).tolist(),
            "cylindrical_components": ["radial", "tangential", "axial"],
            "mean_cylindrical": np.asarray(self.mean_cylindrical, dtype=float).tolist(),
            "angular_orders": list(self.angular_orders),
            "successive_mean_relative_differences": list(
                self.successive_mean_relative_differences
            ),
            "mean_rms": self.mean_rms,
            "full_ring_self_advection_rms": self.full_ring_self_advection_rms,
            "mean_to_full_rms_ratio": self.mean_to_full_rms_ratio,
            "backend_kind": self.backend_kind,
            "backend": None if self.backend is None else self.backend.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _validate_rings(radius: Any, z: Any, t: Any) -> tuple[np.ndarray, ...]:
    rb, zb, tb = np.broadcast_arrays(
        np.asarray(radius, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    if not (np.all(np.isfinite(rb)) and np.all(np.isfinite(zb)) and np.all(np.isfinite(tb))):
        raise ValueError("radius, z and t must be finite and broadcastable")
    if np.any(rb <= 0.0) or np.any(rb >= 2.0):
        raise ValueError("ring radius must satisfy 0 < r < 2")
    if np.any(np.abs(zb) >= 2.0):
        raise ValueError("ring z must satisfy |z| < 2")
    if np.any(tb < 0.25) or np.any(tb > 0.75):
        raise ValueError("ring time must lie in [0.25, 0.75]")
    return rb, zb, tb


def _mean_for_order(
    evaluator: SelfAdvectionEvaluator,
    radius: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    order: int,
) -> tuple[np.ndarray, np.ndarray]:
    theta = 2.0 * math.pi * np.arange(order, dtype=float) / float(order)
    c = np.cos(theta)
    s = np.sin(theta)

    xq = radius[..., None] * c
    yq = radius[..., None] * s
    zq = z[..., None] + np.zeros_like(theta)
    tq = t[..., None] + np.zeros_like(theta)

    diagnostic = evaluator(
        xq,
        yq,
        zq,
        tq,
        spatial_step=AGENT2_SPATIAL_STEP,
    )
    if not isinstance(diagnostic, dict) or "self_advection" not in diagnostic:
        raise TypeError("Agent-2 evaluator must return a dict containing self_advection")
    adv = np.asarray(diagnostic["self_advection"], dtype=float)
    expected_shape = radius.shape + (order, 3)
    if adv.shape != expected_shape:
        raise ValueError(
            f"Agent-2 self_advection shape {adv.shape} != expected {expected_shape}"
        )
    if not np.all(np.isfinite(adv)):
        raise ValueError("Agent-2 self_advection contains non-finite values")

    radial = adv[..., 0] * c + adv[..., 1] * s
    tangential = -adv[..., 0] * s + adv[..., 1] * c
    axial = adv[..., 2]
    mean = np.stack(
        [
            np.mean(radial, axis=-1),
            np.mean(tangential, axis=-1),
            np.mean(axial, axis=-1),
        ],
        axis=-1,
    )
    full_ring_rms = np.sqrt(np.mean(np.sum(adv * adv, axis=-1), axis=-1))
    return mean, full_ring_rms


def _relative_rms_difference(reference: np.ndarray, value: np.ndarray) -> float:
    ref = np.asarray(reference, dtype=float)
    val = np.asarray(value, dtype=float)
    numerator = float(np.sqrt(np.mean((val - ref) ** 2)))
    denominator = max(float(np.sqrt(np.mean(ref**2))), np.finfo(float).tiny)
    return numerator / denominator


def _materialize_with_evaluator(
    evaluator: SelfAdvectionEvaluator,
    radius: Any,
    z: Any,
    t: Any,
    *,
    backend: ExactAgent2SelfAdvectionBackend | None,
    backend_kind: str,
) -> OscillatoryQuadraticMeanWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    means: list[np.ndarray] = []
    full_rms: np.ndarray | None = None
    for order in ANGULAR_ORDERS:
        mean, ring_rms = _mean_for_order(evaluator, rb, zb, tb, order)
        means.append(mean)
        full_rms = ring_rms

    successive = tuple(
        _relative_rms_difference(means[index + 1], means[index])
        for index in range(len(means) - 1)
    )
    final_mean = means[-1]
    assert full_rms is not None
    mean_rms = float(np.sqrt(np.mean(np.sum(final_mean * final_mean, axis=-1))))
    full_ring = float(np.sqrt(np.mean(full_rms * full_rms)))
    ratio = mean_rms / max(full_ring, np.finfo(float).tiny)

    return OscillatoryQuadraticMeanWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_cylindrical=np.array(final_mean, copy=True),
        angular_orders=ANGULAR_ORDERS,
        successive_mean_relative_differences=successive,
        mean_rms=mean_rms,
        full_ring_self_advection_rms=full_ring,
        mean_to_full_rms_ratio=float(ratio),
        backend=backend,
        backend_kind=backend_kind,
    )


def materialize_agent2_oscillatory_quadratic_mean(
    backend: ExactAgent2SelfAdvectionBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> OscillatoryQuadraticMeanWitness:
    """Project the exact pinned Agent-2 quadratic term onto its m=0 ring mean.

    No residual, defect, pressure, forcing, stress, correction target, gain,
    damping choice, validation score, or scientific threshold is accepted.
    """

    if not isinstance(backend, ExactAgent2SelfAdvectionBackend):
        raise TypeError("backend must be ExactAgent2SelfAdvectionBackend")
    if backend.source_blob_sha != AGENT2_SOURCE_BLOB_SHA:
        raise ValueError("Agent-2 backend source provenance drifted")
    return _materialize_with_evaluator(
        backend.evaluator,
        radius,
        z,
        t,
        backend=backend,
        backend_kind="exact_agent2_pr779_self_advection",
    )


def truth_boundary() -> dict[str, object]:
    return {
        "agent2_self_advection_reimplemented_by_agent3": False,
        "cylindrical_m0_projection_executable": True,
        "projection_uses_naive_cartesian_vector_mean": False,
        "quadratic_mean_is_complete_ns_defect": False,
        "leading_cross_terms_included": False,
        "pressure_gradient_included": False,
        "viscous_term_included": False,
        "restricted_forcing_included": False,
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
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> dict[str, np.ndarray]:
    # Pure operator regression: an axisymmetric cylindrical mean plus m=2 modes.
    del z, t
    if spatial_step != AGENT2_SPATIAL_STEP:
        raise ValueError("unexpected frozen spatial step")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    theta = np.arctan2(y, x)
    c = np.cos(theta)
    s = np.sin(theta)
    radial = 2.0 + 0.3 * np.cos(2.0 * theta)
    tangential = -0.5 + 0.2 * np.sin(2.0 * theta)
    axial = 0.25 + 0.1 * np.cos(2.0 * theta)
    adv = np.stack(
        [radial * c - tangential * s, radial * s + tangential * c, axial],
        axis=-1,
    )
    return {"self_advection": adv}


def analytic_regression_receipt() -> dict[str, object]:
    witness = _materialize_with_evaluator(
        _analytic_regression_evaluator,
        radius=np.asarray([0.45, 0.8, 1.15]),
        z=np.asarray([-0.6, 0.0, 0.7]),
        t=np.asarray([0.42, 0.50, 0.58]),
        backend=None,
        backend_kind="analytic_projection_regression_only",
    )
    expected = np.broadcast_to(np.asarray([2.0, -0.5, 0.25]), (3, 3))
    error = float(np.max(np.abs(witness.mean_cylindrical - expected)))
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
        "agent2_spatial_step": AGENT2_SPATIAL_STEP,
        "analytic_expected_cylindrical_mean": [2.0, -0.5, 0.25],
        "analytic_projection_absolute_max_error": error,
        "analytic_witness": witness.to_receipt(),
        "real_agent2_numeric_mean_observed_in_this_receipt": False,
        "truth_boundary": truth_boundary(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = analytic_regression_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True)
    if args.output is None:
        print(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
        print(args.output)


if __name__ == "__main__":
    main()
