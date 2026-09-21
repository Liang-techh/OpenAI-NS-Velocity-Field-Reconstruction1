"""Current-lineage oscillatory-induced nonlinear mean attribution.

This Agent-3 increment consumes the exact Agent-2 #970 partial field

    u_partial = u_lead,current + u_osc

and the exact Agent-2 #960 parameter-free oscillatory Jacobian.  Agent 3 adds only
one fixed numerical realization needed by the current interface: a centered
Cartesian derivative of the already-authenticated current leading velocity.  It
then forms

    A = (u_lead . grad) u_osc,
    B = (u_osc . grad) u_lead,
    Q = (u_osc . grad) u_osc,
    N = A + B + Q,

and reuses the existing Agent-3 cylindrical m=0 projector from PR #915.  No
pressure, forcing, radial inverse, correction velocity or complete NS residual is
introduced here.  The fixed leading derivative is a repository numerical
realization, not a paper-exact derivative or an independent PDE validation.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np

from .kokuno_oscillatory_nonlinear_mean_attribution import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    OscillatoryNonlinearMeanAttributionWitness,
    _materialize_from_provider,
)

TASK = "KOKUNO-A3-CURRENT-PARTIAL-NONLINEAR-MEAN-102"
SCHEMA = "kokuno-a3-current-partial-nonlinear-mean-v1"
PARENT_AGENT3_PR = 967
PARENT_AGENT3_HEAD = "4259def27bc2395a046dd8bcb778fb3fa036d64d"

AGENT2_COMPOSITE_PR = 970
AGENT2_COMPOSITE_HEAD = "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4"
AGENT2_COMPOSITE_MODULE = (
    "openai_ns_reconstruction.kokuno_current_partial_leading_oscillatory_velocity"
)
AGENT2_COMPOSITE_CLASS = "_CurrentPartialLeadingOscillatoryField"
AGENT2_COMPOSITE_SOURCE_BLOB = "532435705bd5db341e82371180da8c04c3dd6c14"

AGENT2_DIFFERENTIAL_PR = 960
AGENT2_DIFFERENTIAL_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
AGENT2_DIFFERENTIAL_MODULE = "openai_ns_reconstruction.kokuno_oscillatory_batch_differentials"
AGENT2_DIFFERENTIAL_FUNCTION = "evaluate_oscillatory_batch_differentials"
AGENT2_DIFFERENTIAL_SOURCE_BLOB = "12df3bf6c949baeebf609b366f2973e7f515a157"

AGENT1_LEADING_PR = 965
AGENT1_LEADING_HEAD = "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"

SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

# Deliberately fixed and not a public input.  Matching the A2 #960 spatial scale
# makes the mixed-term diagnostic easier to compare without creating a tuning knob.
FIXED_LEADING_SPATIAL_STEP = 1.0e-3
COMPOSITION_CLOSURE_GATE = 2.0e-12


def _git_blob_sha1(path: str | Path) -> str:
    raw = Path(path).read_bytes()
    header = f"blob {len(raw)}\0".encode("ascii")
    return hashlib.sha1(header + raw).hexdigest()


def _valid_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(c in "0123456789abcdef" for c in value)
    )


def _broadcast_xyzt(
    x: Any, y: Any, z: Any, t: Any
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    try:
        xx, yy, zz, tt = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
    except ValueError as exc:
        raise ValueError("x, y, z, t must be mutually broadcastable") from exc
    if not all(np.all(np.isfinite(v)) for v in (xx, yy, zz, tt)):
        raise ValueError("x, y, z, t must be finite")
    return xx, yy, zz, tt


def _validate_vector(value: Any, shape: tuple[int, ...], label: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != shape + (3,):
        raise RuntimeError(f"{label} returned shape {out.shape}, expected {shape + (3,)}")
    if not np.all(np.isfinite(out)):
        raise RuntimeError(f"{label} returned non-finite values")
    return out


def _validate_jacobian(value: Any, shape: tuple[int, ...], label: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != shape + (3, 3):
        raise RuntimeError(
            f"{label} returned shape {out.shape}, expected {shape + (3, 3)}"
        )
    if not np.all(np.isfinite(out)):
        raise RuntimeError(f"{label} returned non-finite values")
    return out


def _leading_jacobian_fixed_centered(
    leading_velocity: Callable[[Any, Any, Any, Any], Any],
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    """Centered Cartesian derivative of the current A1 leading velocity.

    The convention is ``J[..., component, axis] = partial_axis u_component``.
    The step is frozen above and is intentionally absent from the public API.
    """
    h = FIXED_LEADING_SPATIAL_STEP
    coords = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]
    shape = x.shape
    jac = np.empty(shape + (3, 3), dtype=float)
    for axis in range(3):
        plus = [np.array(v, copy=True) for v in coords]
        minus = [np.array(v, copy=True) for v in coords]
        plus[axis] += h
        minus[axis] -= h
        up = _validate_vector(
            leading_velocity(plus[0], plus[1], plus[2], t), shape, "leading velocity +h"
        )
        um = _validate_vector(
            leading_velocity(minus[0], minus[1], minus[2], t), shape, "leading velocity -h"
        )
        jac[..., :, axis] = (up - um) / (2.0 * h)
    return jac


def _components_from_interfaces(
    leading_velocity: Callable[[Any, Any, Any, Any], Any],
    composite_velocity: Callable[[Any, Any, Any, Any], Any],
    oscillatory_differentials: Callable[[Any, Any], Any],
    x: Any,
    y: Any,
    z: Any,
    t: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Internal formula adapter; public callers cannot supply these interfaces."""
    xx, yy, zz, tt = _broadcast_xyzt(x, y, z, t)
    shape = xx.shape
    points = np.stack((xx, yy, zz), axis=-1)

    leading = _validate_vector(leading_velocity(xx, yy, zz, tt), shape, "leading velocity")
    total = _validate_vector(composite_velocity(xx, yy, zz, tt), shape, "partial composite velocity")
    osc_eval = oscillatory_differentials(points, tt)
    oscillatory = _validate_vector(getattr(osc_eval, "velocity", None), shape, "oscillatory velocity")
    osc_jac = _validate_jacobian(
        getattr(osc_eval, "velocity_jacobian", None), shape, "oscillatory Jacobian"
    )

    composition_error = float(np.max(np.abs(total - (leading + oscillatory)), initial=0.0))
    if composition_error > COMPOSITION_CLOSURE_GATE:
        raise RuntimeError("A2 #970 composite no longer closes as leading + oscillatory")

    lead_jac = _leading_jacobian_fixed_centered(leading_velocity, xx, yy, zz, tt)
    inner_advects_osc = np.einsum("...ij,...j->...i", osc_jac, leading)
    osc_advects_inner = np.einsum("...ij,...j->...i", lead_jac, oscillatory)
    osc_self = np.einsum("...ij,...j->...i", osc_jac, oscillatory)
    mixed = inner_advects_osc + osc_advects_inner
    aggregate = mixed + osc_self
    return inner_advects_osc, osc_advects_inner, osc_self, mixed, aggregate


@dataclass(frozen=True)
class ExactCurrentPartialNonlinearBackend:
    """Checksum-bound bridge from A2 #970/#960 into the A3 mean projector."""

    composite_field: object
    differential_function: Callable[[Any, Any], Any]
    composite_semantic_sha256: str
    differential_semantic_sha256: str
    composite_source_blob: str
    differential_source_blob: str

    @classmethod
    def bind(
        cls,
        composite_field: object,
        differential_function: Callable[[Any, Any], Any],
    ) -> "ExactCurrentPartialNonlinearBackend":
        field_type = type(composite_field)
        if field_type.__module__ != AGENT2_COMPOSITE_MODULE:
            raise ValueError("unexpected A2 #970 composite module identity")
        if field_type.__name__ != AGENT2_COMPOSITE_CLASS:
            raise ValueError("unexpected A2 #970 composite class identity")
        if not callable(getattr(composite_field, "velocity", None)):
            raise TypeError("A2 #970 composite must expose velocity(x,y,z,t)")
        if not callable(getattr(composite_field, "semantic_payload", None)):
            raise TypeError("A2 #970 composite must expose semantic_payload()")
        leading = getattr(composite_field, "leading_backend", None)
        if leading is None or not callable(getattr(leading, "velocity", None)):
            raise TypeError("A2 #970 composite must retain its authenticated leading backend")

        field_source = inspect.getsourcefile(field_type)
        if field_source is None:
            raise ValueError("A2 #970 composite source is unavailable")
        field_blob = _git_blob_sha1(field_source)
        if field_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
            raise ValueError("A2 #970 composite source blob drifted")

        if not callable(differential_function):
            raise TypeError("A2 #960 differential function must be callable")
        if differential_function.__module__ != AGENT2_DIFFERENTIAL_MODULE:
            raise ValueError("unexpected A2 #960 differential module identity")
        if differential_function.__name__ != AGENT2_DIFFERENTIAL_FUNCTION:
            raise ValueError("unexpected A2 #960 differential function identity")
        diff_source = inspect.getsourcefile(differential_function)
        if diff_source is None:
            raise ValueError("A2 #960 differential source is unavailable")
        diff_blob = _git_blob_sha1(diff_source)
        if diff_blob != AGENT2_DIFFERENTIAL_SOURCE_BLOB:
            raise ValueError("A2 #960 differential source blob drifted")

        payload = composite_field.semantic_payload()
        parent = payload.get("parent_agent2", {})
        a1 = payload.get("agent1_leading", {})
        truth = payload.get("truth_boundary", {})
        if parent.get("pr") != AGENT2_DIFFERENTIAL_PR or parent.get("head") != AGENT2_DIFFERENTIAL_HEAD:
            raise ValueError("A2 #970 parent differential lineage drifted")
        if a1.get("pr") != AGENT1_LEADING_PR or a1.get("head") != AGENT1_LEADING_HEAD:
            raise ValueError("A2 #970 Agent-1 leading lineage drifted")
        required_false = (
            "velocity_beyond_Xh_materialized",
            "outer_global_leading_velocity_materialized",
            "matched_global_pressure_materialized",
            "restricted_forcing_materialized",
            "heldout_ns_residual_assessed",
            "same_protocol_comparable_to_st006",
            "residual_reduction_claimed",
            "pde_validated",
        )
        if truth.get("leading_plus_oscillatory_partial_velocity_materialized") is not True:
            raise ValueError("A2 #970 partial composition truth boundary drifted")
        for key in required_false:
            if truth.get(key) is not False:
                raise ValueError(f"A2 #970 truth boundary drifted at {key}")

        composite_semantic = str(getattr(composite_field, "semantic_sha256", ""))
        if not _valid_sha256(composite_semantic):
            raise ValueError("A2 #970 composite semantic identity is malformed")

        diff_module = inspect.getmodule(differential_function)
        diff_semantic_fn = getattr(diff_module, "batch_differentials_sha256", None)
        if not callable(diff_semantic_fn):
            raise ValueError("A2 #960 differential semantic identity is unavailable")
        diff_semantic = str(diff_semantic_fn())
        if not _valid_sha256(diff_semantic):
            raise ValueError("A2 #960 differential semantic identity is malformed")
        if parent.get("batch_differentials_sha256") != diff_semantic:
            raise ValueError("A2 #970/#960 differential semantic identities disagree")

        return cls(
            composite_field=composite_field,
            differential_function=differential_function,
            composite_semantic_sha256=composite_semantic,
            differential_semantic_sha256=diff_semantic,
            composite_source_blob=field_blob,
            differential_source_blob=diff_blob,
        )

    def components(
        self, x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return _components_from_interfaces(
            self.composite_field.leading_backend.velocity,
            self.composite_field.velocity,
            self.differential_function,
            x,
            y,
            z,
            t,
        )

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent2_composite_pr": AGENT2_COMPOSITE_PR,
            "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
            "agent2_composite_source_blob": self.composite_source_blob,
            "agent2_composite_semantic_sha256": self.composite_semantic_sha256,
            "agent2_differential_pr": AGENT2_DIFFERENTIAL_PR,
            "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
            "agent2_differential_source_blob": self.differential_source_blob,
            "agent2_differential_semantic_sha256": self.differential_semantic_sha256,
            "agent1_leading_pr": AGENT1_LEADING_PR,
            "agent1_leading_head": AGENT1_LEADING_HEAD,
            "fixed_leading_spatial_step": FIXED_LEADING_SPATIAL_STEP,
            "leading_derivative_realization": "centered_cartesian_fd2_repository_fixed",
        }


def materialize_current_partial_nonlinear_mean_attribution(
    backend: ExactCurrentPartialNonlinearBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> OscillatoryNonlinearMeanAttributionWitness:
    """Project current-lineage nonlinear pieces without accepting surrogate defects."""
    if not isinstance(backend, ExactCurrentPartialNonlinearBackend):
        raise TypeError("backend must be ExactCurrentPartialNonlinearBackend")
    if backend.composite_source_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
        raise ValueError("A2 #970 composite provenance drifted")
    if backend.differential_source_blob != AGENT2_DIFFERENTIAL_SOURCE_BLOB:
        raise ValueError("A2 #960 differential provenance drifted")
    return _materialize_from_provider(
        backend.components,
        radius,
        z,
        t,
        backend=backend,  # witness only requires a to_receipt() backend at runtime
        backend_kind="exact-current-a2-pr970-partial-nonlinear-attribution",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_partial_nonlinear_mean_attribution)
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "inverse", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "angular_order",
        "spatial_step", "derivative_step", "time_step", "viscosity", "nu",
        "normalized_score", "scientific_threshold", "correction",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_composite_pr": AGENT2_COMPOSITE_PR,
        "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
        "agent2_differential_pr": AGENT2_DIFFERENTIAL_PR,
        "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
        "current_partial_leading_plus_oscillation_consumed": True,
        "current_partial_mixed_nonlinear_mean_materialized": True,
        "current_partial_quadratic_nonlinear_mean_materialized": True,
        "current_partial_aggregate_nonlinear_mean_materialized": True,
        "current_leading_spatial_derivative_is_fixed_repository_fd": True,
        "caller_tunable_leading_spatial_step": False,
        "agent2_oscillatory_jacobian_reimplemented_by_agent3": False,
        "radial_inverse_performed_in_this_increment": False,
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "velocity_beyond_Xh_materialized": False,
        "outer_global_leading_velocity_materialized": False,
        "complete_ns_defect": False,
        "scoped_nonlinear_mean_authorized_as_correction_target": False,
        "mean_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "caller_supplied_surrogate_defect_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "paper_exact": False,
        "pde_validated": False,
        "blowup_proved": False,
    }
