"""Agent-3 cylindrical m=0 projection of the strict-inner transport block.

This module consumes the exact Agent-2 PR #840 transport composition for the
current Agent-1 PA.10 inner contraction-center velocity plus the frozen
oscillatory velocity,

    T = d_t u + (u . grad) u - nu Delta u,    nu = 0.01,
    u = u_inner + u_osc,

and projects the time, advection, viscous, and total Cartesian vectors onto the
rotating cylindrical m=0 frame.  It does not reconstruct Agent-2 derivatives
or curl machinery and does not add pressure, forcing, a global leading join,
or an Agent-3 correction velocity.  Therefore the witness is a strict-inner
pressure/forcing-free transport mean, not a complete Navier--Stokes defect.
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

TASK = "KOKUNO-A3-STRICT-INNER-TRANSPORT-MEAN-072"
SCHEMA = "kokuno-a3-strict-inner-transport-mean-v1"
PARENT_AGENT3_PR = 841
PARENT_AGENT3_HEAD = "d8b0e792d2a263bb97754d99ee6be998611c90c3"

AGENT2_PR = 840
AGENT2_HEAD = "016e3c152d7d11e4fedcd79b819df337f8942983"
AGENT2_MODULE = "openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_transport"
AGENT2_FUNCTION = "evaluate_inner_leading_oscillatory_transport"
AGENT2_SOURCE_BLOB_SHA = "d6082ebd346ba8150321bf3d62f104c107a70300"

AGENT1_PR = 839
AGENT1_HEAD = "e96ee90144976a992b62d76d4361a94eb16bd91e"
AGENT1_MODULE = "openai_ns_reconstruction.kokuno_pa10_cartesian_center_laplacian"
AGENT1_CLASS = "KokunoPA10CartesianCenterLaplacian"
AGENT1_SOURCE_BLOB_SHA = "74b0e185e8e0bbb08695d493e0a091c305a03006"

REPOSITORY_VISCOSITY = 0.01
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

TransportEvaluator = Callable[..., Any]
TransportProvider = Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]]


@dataclass(frozen=True)
class ExactAgent1InnerTransportBackend:
    """Bind the exact Agent-1 #839 strict-inner backend object."""

    candidate: object
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(cls, candidate: object) -> "ExactAgent1InnerTransportBackend":
        candidate_type = type(candidate)
        if candidate_type.__module__ != AGENT1_MODULE:
            raise ValueError("unexpected Agent-1 strict-inner module identity")
        if candidate_type.__name__ != AGENT1_CLASS:
            raise ValueError("unexpected Agent-1 strict-inner class identity")
        for name in ("velocity", "velocity_dt", "self_advection", "velocity_laplacian"):
            if not callable(getattr(candidate, name, None)):
                raise TypeError(f"Agent-1 backend must expose callable {name}")
        source = inspect.getsourcefile(candidate_type)
        if source is None:
            raise ValueError("Agent-1 strict-inner source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-1 strict-inner source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT1_SOURCE_BLOB_SHA:
            raise ValueError("Agent-1 source blob does not match pinned PR #839")
        return cls(candidate=candidate, source_file=str(path), source_blob_sha=blob_sha)

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent1_pr": AGENT1_PR,
            "agent1_expected_head": AGENT1_HEAD,
            "module": AGENT1_MODULE,
            "class": AGENT1_CLASS,
            "source_blob_sha": self.source_blob_sha,
            "source_blob_matches_pinned": self.source_blob_sha == AGENT1_SOURCE_BLOB_SHA,
        }


@dataclass(frozen=True)
class ExactAgent2StrictInnerTransportBackend:
    """Bind the exact Agent-2 #840 strict-inner transport evaluator."""

    evaluator: TransportEvaluator
    source_file: str
    source_blob_sha: str

    @classmethod
    def bind(
        cls, evaluator: TransportEvaluator
    ) -> "ExactAgent2StrictInnerTransportBackend":
        if not callable(evaluator):
            raise TypeError("Agent-2 strict-inner transport evaluator must be callable")
        if getattr(evaluator, "__module__", None) != AGENT2_MODULE:
            raise ValueError("unexpected Agent-2 strict-inner transport module identity")
        if getattr(evaluator, "__name__", None) != AGENT2_FUNCTION:
            raise ValueError("unexpected Agent-2 strict-inner transport function identity")
        source = inspect.getsourcefile(evaluator)
        if source is None:
            raise ValueError("Agent-2 strict-inner transport source file is unavailable")
        path = Path(source).resolve()
        if not path.is_file():
            raise ValueError("Agent-2 strict-inner transport source file does not exist")
        blob_sha = _git_blob_sha(path.read_bytes())
        if blob_sha != AGENT2_SOURCE_BLOB_SHA:
            raise ValueError("Agent-2 source blob does not match pinned PR #840")
        return cls(evaluator=evaluator, source_file=str(path), source_blob_sha=blob_sha)

    def components(
        self,
        inner_backend: ExactAgent1InnerTransportBackend,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
        if not isinstance(inner_backend, ExactAgent1InnerTransportBackend):
            raise TypeError("inner_backend must be ExactAgent1InnerTransportBackend")
        result = self.evaluator(inner_backend.candidate, x, y, z, t)
        viscosity = float(getattr(result, "viscosity", np.nan))
        if viscosity != REPOSITORY_VISCOSITY:
            raise ValueError("Agent-2 strict-inner transport viscosity drifted")
        time_part = np.asarray(
            getattr(result, "inner_plus_oscillatory_velocity_dt"), dtype=float
        )
        advection = np.asarray(
            getattr(result, "inner_plus_oscillatory_self_advection"), dtype=float
        )
        viscous = np.asarray(getattr(result, "viscous_term"), dtype=float)
        transport = np.asarray(getattr(result, "transport"), dtype=float)
        return time_part, advection, viscous, transport, viscosity

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
class StrictInnerTransportMeanWitness:
    """Typed cylindrical mean witness for strict-inner ``u_t+u.grad u-nu Delta u``."""

    radius: np.ndarray
    z: np.ndarray
    t: np.ndarray
    mean_time_derivative_cylindrical: np.ndarray
    mean_advection_cylindrical: np.ndarray
    mean_viscous_cylindrical: np.ndarray
    mean_transport_cylindrical: np.ndarray
    angular_orders: tuple[int, ...]
    successive_transport_mean_relative_differences: tuple[float, ...]
    pointwise_component_closure_absolute_max: float
    projected_component_closure_absolute_max: float
    transport_mean_rms: float
    full_ring_transport_rms: float
    transport_mean_to_full_rms_ratio: float
    viscosity: float
    agent1_backend: ExactAgent1InnerTransportBackend | None
    agent2_backend: ExactAgent2StrictInnerTransportBackend | None
    backend_kind: str

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "radius": np.asarray(self.radius, dtype=float).tolist(),
            "z": np.asarray(self.z, dtype=float).tolist(),
            "t": np.asarray(self.t, dtype=float).tolist(),
            "cylindrical_components": ["radial", "tangential", "axial"],
            "mean_time_derivative_cylindrical": np.asarray(
                self.mean_time_derivative_cylindrical, dtype=float
            ).tolist(),
            "mean_advection_cylindrical": np.asarray(
                self.mean_advection_cylindrical, dtype=float
            ).tolist(),
            "mean_viscous_cylindrical": np.asarray(
                self.mean_viscous_cylindrical, dtype=float
            ).tolist(),
            "mean_transport_cylindrical": np.asarray(
                self.mean_transport_cylindrical, dtype=float
            ).tolist(),
            "angular_orders": list(self.angular_orders),
            "successive_transport_mean_relative_differences": list(
                self.successive_transport_mean_relative_differences
            ),
            "pointwise_component_closure_absolute_max": (
                self.pointwise_component_closure_absolute_max
            ),
            "projected_component_closure_absolute_max": (
                self.projected_component_closure_absolute_max
            ),
            "transport_mean_rms": self.transport_mean_rms,
            "full_ring_transport_rms": self.full_ring_transport_rms,
            "transport_mean_to_full_rms_ratio": self.transport_mean_to_full_rms_ratio,
            "viscosity": self.viscosity,
            "backend_kind": self.backend_kind,
            "agent1_backend": (
                None if self.agent1_backend is None else self.agent1_backend.to_receipt()
            ),
            "agent2_backend": (
                None if self.agent2_backend is None else self.agent2_backend.to_receipt()
            ),
            "truth_boundary": truth_boundary(),
        }


def _validate_vector(
    value: Any, expected_shape: tuple[int, ...], *, label: str
) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != expected_shape:
        raise ValueError(f"{label} returned shape {array.shape}, expected {expected_shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{label} returned non-finite values")
    return array


def _mean_for_order(
    provider: TransportProvider,
    radius: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, float, float]:
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
        raise TypeError("transport provider must return five entries")
    time_part = _validate_vector(raw[0], expected, label="time derivative")
    advection = _validate_vector(raw[1], expected, label="advection")
    viscous = _validate_vector(raw[2], expected, label="viscous term")
    transport = _validate_vector(raw[3], expected, label="transport")
    viscosity = float(raw[4])
    if not np.isfinite(viscosity) or viscosity != REPOSITORY_VISCOSITY:
        raise ValueError("transport provider viscosity must remain exactly 0.01")

    assembled = time_part + advection + viscous
    pointwise_closure = float(np.max(np.abs(transport - assembled)))

    mean_time = _project_cartesian_ring_to_cylindrical_mean(time_part, c, s)
    mean_advection = _project_cartesian_ring_to_cylindrical_mean(advection, c, s)
    mean_viscous = _project_cartesian_ring_to_cylindrical_mean(viscous, c, s)
    mean_transport = _project_cartesian_ring_to_cylindrical_mean(transport, c, s)
    projected_closure = float(
        np.max(np.abs(mean_transport - mean_time - mean_advection - mean_viscous))
    )
    full_ring_rms = float(
        np.sqrt(np.mean(np.sum(transport * transport, axis=-1)))
    )
    return (
        mean_time,
        mean_advection,
        mean_viscous,
        mean_transport,
        pointwise_closure,
        projected_closure,
        full_ring_rms,
        viscosity,
    )


def _materialize_from_transport_provider(
    provider: TransportProvider,
    radius: Any,
    z: Any,
    t: Any,
    *,
    agent1_backend: ExactAgent1InnerTransportBackend | None,
    agent2_backend: ExactAgent2StrictInnerTransportBackend | None,
    backend_kind: str,
) -> StrictInnerTransportMeanWitness:
    rb, zb, tb = _validate_rings(radius, z, t)
    per_order = [
        _mean_for_order(provider, rb, zb, tb, order) for order in ANGULAR_ORDERS
    ]
    totals = [entry[3] for entry in per_order]
    successive = tuple(
        _relative_rms_difference(totals[index + 1], totals[index])
        for index in range(len(totals) - 1)
    )
    mean_time, mean_advection, mean_viscous, mean_transport, _, _, full_rms, viscosity = per_order[-1]
    mean_rms = float(
        np.sqrt(np.mean(np.sum(mean_transport * mean_transport, axis=-1)))
    )
    ratio = mean_rms / max(float(full_rms), np.finfo(float).tiny)
    pointwise_closure = float(max(entry[4] for entry in per_order))
    projected_closure = float(max(entry[5] for entry in per_order))

    return StrictInnerTransportMeanWitness(
        radius=np.array(rb, copy=True),
        z=np.array(zb, copy=True),
        t=np.array(tb, copy=True),
        mean_time_derivative_cylindrical=np.array(mean_time, copy=True),
        mean_advection_cylindrical=np.array(mean_advection, copy=True),
        mean_viscous_cylindrical=np.array(mean_viscous, copy=True),
        mean_transport_cylindrical=np.array(mean_transport, copy=True),
        angular_orders=ANGULAR_ORDERS,
        successive_transport_mean_relative_differences=successive,
        pointwise_component_closure_absolute_max=pointwise_closure,
        projected_component_closure_absolute_max=projected_closure,
        transport_mean_rms=mean_rms,
        full_ring_transport_rms=float(full_rms),
        transport_mean_to_full_rms_ratio=float(ratio),
        viscosity=float(viscosity),
        agent1_backend=agent1_backend,
        agent2_backend=agent2_backend,
        backend_kind=backend_kind,
    )


def materialize_strict_inner_transport_mean(
    transport_backend: ExactAgent2StrictInnerTransportBackend,
    inner_backend: ExactAgent1InnerTransportBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> StrictInnerTransportMeanWitness:
    """Project the exact #840 raw strict-inner transport onto the frozen A3 ladder."""
    if not isinstance(transport_backend, ExactAgent2StrictInnerTransportBackend):
        raise TypeError("transport_backend must be ExactAgent2StrictInnerTransportBackend")
    if not isinstance(inner_backend, ExactAgent1InnerTransportBackend):
        raise TypeError("inner_backend must be ExactAgent1InnerTransportBackend")

    def provider(x: np.ndarray, y: np.ndarray, zz: np.ndarray, tt: np.ndarray):
        return transport_backend.components(inner_backend, x, y, zz, tt)

    return _materialize_from_transport_provider(
        provider,
        radius,
        z,
        t,
        agent1_backend=inner_backend,
        agent2_backend=transport_backend,
        backend_kind="exact-agent2-840-plus-agent1-839",
    )


def truth_boundary() -> dict[str, object]:
    return {
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "repository_viscosity": REPOSITORY_VISCOSITY,
        "viscosity_caller_tunable": False,
        "agent1_inner_head": AGENT1_HEAD,
        "agent2_transport_head": AGENT2_HEAD,
        "inner_center_scoped_only": True,
        "time_derivative_included": True,
        "convective_term_included": True,
        "base_viscous_term_included": True,
        "strict_inner_pressure_forcing_free_transport_mean_materialized": True,
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "correction_transport_included": False,
        "complete_ns_defect": False,
        "source_center_is_final_corrected_fixed_point": False,
        "global_corrected_leading_join_materialized": False,
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


def _cylindrical_field(
    radial: float, tangential: float, axial: float, mode_scale: float
) -> Callable[[Any, Any, Any, Any], np.ndarray]:
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
        mt = tangential - 0.4 * m2
        mz = axial + 0.6 * mode_scale * np.sin(2.0 * theta)
        return np.stack((mr * c - mt * s, mr * s + mt * c, mz), axis=-1)
    return provider


def _mechanics_provider() -> TransportProvider:
    time_provider = _cylindrical_field(0.07, 0.02, 0.11, 0.23)
    advection_provider = _cylindrical_field(0.36, -0.03, 0.35, -0.31)
    viscous_provider = _cylindrical_field(-0.015, 0.01, -0.02, 0.17)

    def provider(x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray):
        time_part = time_provider(x, y, z, t)
        advection = advection_provider(x, y, z, t)
        viscous = viscous_provider(x, y, z, t)
        transport = time_part + advection + viscous
        return time_part, advection, viscous, transport, REPOSITORY_VISCOSITY

    return provider


def build_mechanics_report() -> dict[str, object]:
    witness = _materialize_from_transport_provider(
        _mechanics_provider(),
        np.array([0.7, 1.1]),
        np.array([-0.2, 0.3]),
        np.array([0.48, 0.52]),
        agent1_backend=None,
        agent2_backend=None,
        backend_kind="manufactured-mechanics-only",
    )
    receipt = witness.to_receipt()
    receipt["mechanics_only"] = True
    receipt["expected_time_mean"] = [0.07, 0.02, 0.11]
    receipt["expected_advection_mean"] = [0.36, -0.03, 0.35]
    receipt["expected_viscous_mean"] = [-0.015, 0.01, -0.02]
    receipt["expected_transport_mean"] = [0.415, 0.0, 0.44]
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
