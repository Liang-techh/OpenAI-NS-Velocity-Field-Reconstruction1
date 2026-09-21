"""Batch spatial differentials for the current partial Kokuno composite velocity.

This Agent-2 increment extends the exact current partial Cartesian composition

    u_partial = u_lead,current-through-X_h + u_osc,frozen-complete-curl

from PR #970 with one parameter-free spatial differential surface.  The oscillatory
Jacobian is consumed directly from Agent-2 PR #960.  The currently available Agent-1
leading velocity has no public Jacobian, so its spatial derivative is realized here by
one frozen centered Cartesian FD6 operator using the same repository step as #960.
The sum therefore remains a repository numerical realization, not a paper-exact
Kokuno derivative and not an independent PDE validation.

No mean projection, radial inverse, correction velocity, pressure, forcing or complete
Navier--Stokes residual is introduced here.  Calls remain fail-closed outside the
currently materialized Agent-1 X_h domain.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from typing import Any, Callable, Mapping

import numpy as np

from .kokuno_current_partial_leading_oscillatory_velocity import (
    default_field as _default_composite_field,
)
from .kokuno_oscillatory_batch_differentials import (
    FIXED_SPATIAL_STEP,
    batch_differentials_sha256,
    evaluate_oscillatory_batch_differentials,
)

TASK = "K2-OSC-081"
SCHEMA = "kokuno-a2-current-partial-composite-batch-differentials-v1"
RECEIPT_SCHEMA = "kokuno-a2-current-partial-composite-batch-differentials-receipt-v1"

PARENT_AGENT2_PR = 970
PARENT_AGENT2_HEAD = "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4"
OSCILLATORY_DIFFERENTIAL_PR = 960
OSCILLATORY_DIFFERENTIAL_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
AGENT1_LEADING_PR = 965
AGENT1_LEADING_HEAD = "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

# Frozen repository numerical realization.  This is deliberately not a public input.
LEADING_FD6_SPATIAL_STEP = FIXED_SPATIAL_STEP
COMPOSITION_CLOSURE_GATE = 2.0e-12
DIRECT_FD6_REPLAY_GATE = 2.0e-8
EXACT_CLOSURE_GATE = 2.0e-13
OSCILLATORY_VORTICITY_RMS_FLOOR = 1.0e-12

_RECEIPT_POINTS = np.asarray(
    [
        (0.32, 0.11, -0.20),
        (0.41, -0.17, -0.08),
        (-0.36, 0.24, 0.05),
        (-0.52, -0.16, 0.16),
        (0.58, 0.21, 0.22),
        (0.47, -0.31, -0.14),
    ],
    dtype=float,
)
_RECEIPT_TIMES = np.asarray((0.31, 0.39, 0.47, 0.55, 0.63, 0.71), dtype=float)
_AXIS_POINTS = np.asarray(
    ((0.0, 0.0, -0.10), (0.0, 0.0, 0.0), (0.0, 0.0, 0.12)), dtype=float
)
_AXIS_TIMES = np.asarray((0.37, 0.51, 0.67), dtype=float)
# These points stay safely inside the current leading domain while lying below the
# frozen oscillatory annulus's positive inner radius.  A +/-3h FD6 stencil remains
# below that radial support boundary as well.
_OSC_EXTERIOR_POINTS = np.asarray(
    ((0.050, 0.000, -0.08), (0.000, -0.050, 0.02), (-0.035, 0.035, 0.09)), dtype=float
)
_OSC_EXTERIOR_TIMES = np.asarray((0.35, 0.49, 0.65), dtype=float)


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _broadcast_time(time: Any, shape: tuple[int, ...]) -> np.ndarray:
    value = np.asarray(time, dtype=float)
    try:
        out = np.broadcast_to(value, shape)
    except ValueError as exc:
        raise ValueError("time must be broadcastable to points.shape[:-1]") from exc
    if not np.all(np.isfinite(out)):
        raise ValueError("time must be finite")
    return np.asarray(out, dtype=float)


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


def _curl_from_jacobian(jacobian: np.ndarray) -> np.ndarray:
    j = np.asarray(jacobian, dtype=float)
    if j.shape[-2:] != (3, 3):
        raise ValueError("jacobian must have trailing shape (3,3)")
    return np.stack(
        (
            j[..., 2, 1] - j[..., 1, 2],
            j[..., 0, 2] - j[..., 2, 0],
            j[..., 1, 0] - j[..., 0, 1],
        ),
        axis=-1,
    )


def _fd6_jacobian(
    velocity_function: Callable[[Any, Any, Any, Any], Any],
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    """Centered Cartesian FD6 Jacobian with the frozen repository step.

    Convention: ``J[..., component, axis] = partial_axis velocity_component``.
    Any upstream domain failure at a shifted stencil point is intentionally allowed to
    propagate; the partial candidate is not extended by fallback or extrapolation.
    """
    coords = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]
    shape = coords[0].shape
    if not (coords[1].shape == shape and coords[2].shape == shape and t.shape == shape):
        raise ValueError("x, y, z, t must have identical broadcast shapes")

    # Standard sixth-order centered first derivative:
    # [-f(-3h)+9f(-2h)-45f(-h)+45f(h)-9f(2h)+f(3h)]/(60h).
    offsets_and_weights = ((-3, -1.0), (-2, 9.0), (-1, -45.0), (1, 45.0), (2, -9.0), (3, 1.0))
    h = float(LEADING_FD6_SPATIAL_STEP)
    jacobian = np.empty(shape + (3, 3), dtype=float)
    for axis in range(3):
        accum = np.zeros(shape + (3,), dtype=float)
        for offset, weight in offsets_and_weights:
            shifted = [np.array(v, copy=True) for v in coords]
            shifted[axis] += offset * h
            sample = _validate_vector(
                velocity_function(shifted[0], shifted[1], shifted[2], t),
                shape,
                f"velocity FD6 axis={axis} offset={offset}",
            )
            accum += weight * sample
        jacobian[..., :, axis] = accum / (60.0 * h)
    if not np.all(np.isfinite(jacobian)):
        raise RuntimeError("FD6 Jacobian became non-finite")
    return jacobian


@dataclass(frozen=True)
class CurrentPartialCompositeBatchDifferentialEvaluation:
    velocity: np.ndarray
    velocity_jacobian: np.ndarray
    divergence: np.ndarray
    vorticity: np.ndarray
    leading_velocity: np.ndarray
    leading_velocity_jacobian: np.ndarray
    oscillatory_velocity: np.ndarray
    oscillatory_velocity_jacobian: np.ndarray
    oscillatory_interior_mask: np.ndarray


def _evaluate_with_interfaces(
    composite_field: Any,
    oscillatory_differentials: Callable[[Any, Any], Any],
    points: Any,
    time: Any,
) -> CurrentPartialCompositeBatchDifferentialEvaluation:
    xyz = np.asarray(points, dtype=float)
    if xyz.ndim < 1 or xyz.shape[-1] != 3:
        raise ValueError("points must have shape [...,3]")
    if not np.all(np.isfinite(xyz)):
        raise ValueError("points must be finite")
    shape = xyz.shape[:-1]
    tt = _broadcast_time(time, shape)

    # Call the A2 oscillatory differential surface first.  It validates the full time
    # array before spatial masking and therefore preserves the frozen A2 time contract.
    osc = oscillatory_differentials(xyz, tt)
    oscillatory_velocity = _validate_vector(
        getattr(osc, "velocity", None), shape, "oscillatory velocity"
    )
    oscillatory_jacobian = _validate_jacobian(
        getattr(osc, "velocity_jacobian", None), shape, "oscillatory Jacobian"
    )
    interior_mask = np.asarray(getattr(osc, "interior_mask", None), dtype=bool)
    if interior_mask.shape != shape:
        raise RuntimeError("oscillatory interior mask has unexpected shape")

    if not callable(getattr(composite_field, "velocity", None)):
        raise TypeError("composite field must expose velocity(x,y,z,t)")
    leading_backend = getattr(composite_field, "leading_backend", None)
    if leading_backend is None or not callable(getattr(leading_backend, "velocity", None)):
        raise TypeError("composite field must retain its authenticated leading backend")

    xx, yy, zz = (xyz[..., 0], xyz[..., 1], xyz[..., 2])
    leading_velocity = _validate_vector(
        leading_backend.velocity(xx, yy, zz, tt), shape, "leading velocity"
    )
    velocity = _validate_vector(
        composite_field.velocity(xx, yy, zz, tt), shape, "partial composite velocity"
    )
    composition_error = float(
        np.max(np.abs(velocity - (leading_velocity + oscillatory_velocity)), initial=0.0)
    )
    if composition_error > COMPOSITION_CLOSURE_GATE:
        raise RuntimeError("partial composite no longer closes as leading + oscillatory")

    leading_jacobian = _fd6_jacobian(
        leading_backend.velocity, xx, yy, zz, tt
    )
    velocity_jacobian = leading_jacobian + oscillatory_jacobian
    divergence = np.trace(velocity_jacobian, axis1=-2, axis2=-1)
    vorticity = _curl_from_jacobian(velocity_jacobian)
    if not all(
        np.all(np.isfinite(v))
        for v in (velocity_jacobian, divergence, vorticity)
    ):
        raise RuntimeError("partial composite differential evaluation became non-finite")

    return CurrentPartialCompositeBatchDifferentialEvaluation(
        velocity=velocity,
        velocity_jacobian=velocity_jacobian,
        divergence=divergence,
        vorticity=vorticity,
        leading_velocity=leading_velocity,
        leading_velocity_jacobian=leading_jacobian,
        oscillatory_velocity=oscillatory_velocity,
        oscillatory_velocity_jacobian=oscillatory_jacobian,
        oscillatory_interior_mask=interior_mask,
    )


def evaluate_current_partial_composite_batch_differentials(
    points: Any, time: Any
) -> CurrentPartialCompositeBatchDifferentialEvaluation:
    """Return velocity/Jacobian/divergence/vorticity for the current partial composite.

    The only public inputs are the spacetime samples.  The leading FD6 step and every
    oscillatory parameter remain frozen.  The exact A1 domain failure beyond current
    ``X_h`` is preserved rather than silently extrapolated.
    """
    return _evaluate_with_interfaces(
        _default_composite_field(),
        evaluate_oscillatory_batch_differentials,
        points,
        time,
    )


def semantic_payload(composite_field: Any | None = None) -> dict[str, Any]:
    field = _default_composite_field() if composite_field is None else composite_field
    composite_semantic = str(getattr(field, "semantic_sha256", ""))
    if len(composite_semantic) != 64 or any(c not in "0123456789abcdef" for c in composite_semantic):
        raise RuntimeError("parent partial-composite semantic identity is malformed")
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2": {
            "pr": PARENT_AGENT2_PR,
            "head": PARENT_AGENT2_HEAD,
            "composite_semantic_sha256": composite_semantic,
        },
        "oscillatory_differential": {
            "pr": OSCILLATORY_DIFFERENTIAL_PR,
            "head": OSCILLATORY_DIFFERENTIAL_HEAD,
            "semantic_sha256": batch_differentials_sha256(),
            "fixed_spatial_step": FIXED_SPATIAL_STEP,
        },
        "agent1_leading": {
            "pr": AGENT1_LEADING_PR,
            "head": AGENT1_LEADING_HEAD,
            "spatial_derivative": "centered_cartesian_fd6_repository_fixed",
            "fixed_spatial_step": LEADING_FD6_SPATIAL_STEP,
        },
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
            "scope": "localized waves / complete curls plus public leading-coordinate structure",
        },
        "truth_boundary": {
            "current_partial_leading_plus_oscillatory_velocity_materialized": True,
            "current_partial_composite_spatial_differentials_materialized": True,
            "velocity_beyond_Xh_materialized": False,
            "outer_global_leading_velocity_materialized": False,
            "global_compact_support_completed": False,
            "full_concrete_oscillatory_runtime_digest_bound": False,
            "identity_preserving_composite_save_load_available": False,
            "mean_projection_performed": False,
            "radial_inverse_performed": False,
            "correction_velocity_constructed": False,
            "matched_global_pressure_materialized": False,
            "restricted_forcing_materialized": False,
            "complete_velocity_pressure_forcing_api": False,
            "complete_ns_residual": False,
            "same_protocol_comparable_to_st006": False,
            "residual_reduction_claimed": False,
            "velocity_export_ready": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "pde_validated": False,
        },
    }


def current_partial_composite_differentials_sha256() -> str:
    return _sha256(semantic_payload())


def public_contract() -> dict[str, Any]:
    params = inspect.signature(evaluate_current_partial_composite_batch_differentials).parameters
    forbidden = {
        "amplitude", "phase", "phase_offset", "scale", "orientation", "support",
        "spatial_step", "derivative_step", "residual", "defect", "target", "forcing",
        "pressure", "viscosity", "nu", "gain", "damping", "threshold", "mean",
        "stress", "inverse", "correction", "angular_order",
    }
    return {
        "public_inputs": list(params),
        "forbidden_inputs_present": sorted(forbidden.intersection(params)),
        "shape_preserving_batch_api": True,
        "fixed_leading_fd6_step": LEADING_FD6_SPATIAL_STEP,
        "fixed_oscillatory_fd6_step": FIXED_SPATIAL_STEP,
        "caller_tunable_derivative_step": False,
        "leading_profile_reimplemented": False,
        "complete_curl_reimplemented": False,
        "mean_projection_performed": False,
        "radial_inverse_performed": False,
        "correction_velocity_constructed": False,
        "pressure_or_forcing_added": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def _diagnostic_for_evaluation(
    evaluation: CurrentPartialCompositeBatchDifferentialEvaluation,
) -> dict[str, float]:
    curl = _curl_from_jacobian(evaluation.velocity_jacobian)
    trace = np.trace(evaluation.velocity_jacobian, axis1=-2, axis2=-1)
    return {
        "composition_closure_max_abs": float(
            np.max(
                np.abs(
                    evaluation.velocity
                    - (evaluation.leading_velocity + evaluation.oscillatory_velocity)
                ),
                initial=0.0,
            )
        ),
        "jacobian_stage_closure_max_abs": float(
            np.max(
                np.abs(
                    evaluation.velocity_jacobian
                    - (
                        evaluation.leading_velocity_jacobian
                        + evaluation.oscillatory_velocity_jacobian
                    )
                ),
                initial=0.0,
            )
        ),
        "divergence_trace_closure_max_abs": float(
            np.max(np.abs(evaluation.divergence - trace), initial=0.0)
        ),
        "vorticity_curl_closure_max_abs": float(
            np.max(np.abs(evaluation.vorticity - curl), initial=0.0)
        ),
    }


def materialize_current_partial_composite_differentials_receipt() -> dict[str, Any]:
    field = _default_composite_field()
    evaluation = _evaluate_with_interfaces(
        field,
        evaluate_oscillatory_batch_differentials,
        _RECEIPT_POINTS,
        _RECEIPT_TIMES,
    )
    direct_jacobian = _fd6_jacobian(
        field.velocity,
        _RECEIPT_POINTS[:, 0],
        _RECEIPT_POINTS[:, 1],
        _RECEIPT_POINTS[:, 2],
        _RECEIPT_TIMES,
    )
    axis = _evaluate_with_interfaces(
        field,
        evaluate_oscillatory_batch_differentials,
        _AXIS_POINTS,
        _AXIS_TIMES,
    )
    exterior = _evaluate_with_interfaces(
        field,
        evaluate_oscillatory_batch_differentials,
        _OSC_EXTERIOR_POINTS,
        _OSC_EXTERIOR_TIMES,
    )

    diagnostic = _diagnostic_for_evaluation(evaluation)
    diagnostic.update(
        {
            "interior_probe_count": int(_RECEIPT_POINTS.shape[0]),
            "axis_probe_count": int(_AXIS_POINTS.shape[0]),
            "oscillatory_exterior_probe_count": int(_OSC_EXTERIOR_POINTS.shape[0]),
            "direct_composite_fd6_replay_max_abs": float(
                np.max(np.abs(evaluation.velocity_jacobian - direct_jacobian), initial=0.0)
            ),
            "oscillatory_vorticity_vector_rms": float(
                np.sqrt(
                    np.mean(
                        np.sum(
                            _curl_from_jacobian(evaluation.oscillatory_velocity_jacobian) ** 2,
                            axis=-1,
                        )
                    )
                )
            ),
            "axis_oscillatory_velocity_abs_max": float(
                np.max(np.abs(axis.oscillatory_velocity), initial=0.0)
            ),
            "axis_oscillatory_jacobian_abs_max": float(
                np.max(np.abs(axis.oscillatory_velocity_jacobian), initial=0.0)
            ),
            "axis_total_minus_leading_jacobian_abs_max": float(
                np.max(
                    np.abs(axis.velocity_jacobian - axis.leading_velocity_jacobian), initial=0.0
                )
            ),
            "exterior_oscillatory_velocity_abs_max": float(
                np.max(np.abs(exterior.oscillatory_velocity), initial=0.0)
            ),
            "exterior_oscillatory_jacobian_abs_max": float(
                np.max(np.abs(exterior.oscillatory_velocity_jacobian), initial=0.0)
            ),
            "exterior_total_minus_leading_jacobian_abs_max": float(
                np.max(
                    np.abs(exterior.velocity_jacobian - exterior.leading_velocity_jacobian),
                    initial=0.0,
                )
            ),
        }
    )

    beyond_fail_closed = False
    try:
        evaluate_current_partial_composite_batch_differentials(
            np.asarray([[1.0e6, 0.0, 0.0]], dtype=float), 0.5
        )
    except (ValueError, RuntimeError):
        beyond_fail_closed = True
    diagnostic["beyond_current_Xh_fail_closed"] = beyond_fail_closed

    return {
        "schema": RECEIPT_SCHEMA,
        "semantic_sha256": current_partial_composite_differentials_sha256(),
        "semantic_payload": semantic_payload(field),
        "public_contract": public_contract(),
        "diagnostic": diagnostic,
        "frozen_mechanical_gates": {
            "composition_closure_max_abs": COMPOSITION_CLOSURE_GATE,
            "jacobian_stage_closure_max_abs": EXACT_CLOSURE_GATE,
            "divergence_trace_closure_max_abs": EXACT_CLOSURE_GATE,
            "vorticity_curl_closure_max_abs": EXACT_CLOSURE_GATE,
            "direct_composite_fd6_replay_max_abs": DIRECT_FD6_REPLAY_GATE,
            "oscillatory_vorticity_vector_rms_min": OSCILLATORY_VORTICITY_RMS_FLOOR,
            "axis_oscillatory_velocity_abs_max": 0.0,
            "axis_oscillatory_jacobian_abs_max": 0.0,
            "axis_total_minus_leading_jacobian_abs_max": 0.0,
            "exterior_oscillatory_velocity_abs_max": 0.0,
            "exterior_oscillatory_jacobian_abs_max": 0.0,
            "exterior_total_minus_leading_jacobian_abs_max": 0.0,
            "beyond_current_Xh_fail_closed": True,
        },
    }


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != RECEIPT_SCHEMA:
        raise AssertionError("unexpected partial-composite differential receipt schema")
    if receipt.get("semantic_sha256") != current_partial_composite_differentials_sha256():
        raise AssertionError("partial-composite differential semantic identity drifted")
    contract = receipt.get("public_contract", {})
    if contract.get("public_inputs") != ["points", "time"]:
        raise AssertionError("public input surface drifted")
    if contract.get("forbidden_inputs_present") != []:
        raise AssertionError("scientific tuning input leaked into the public API")

    diagnostic = receipt.get("diagnostic", {})
    gates = receipt.get("frozen_mechanical_gates", {})
    for key in (
        "composition_closure_max_abs",
        "jacobian_stage_closure_max_abs",
        "divergence_trace_closure_max_abs",
        "vorticity_curl_closure_max_abs",
        "direct_composite_fd6_replay_max_abs",
        "axis_oscillatory_velocity_abs_max",
        "axis_oscillatory_jacobian_abs_max",
        "axis_total_minus_leading_jacobian_abs_max",
        "exterior_oscillatory_velocity_abs_max",
        "exterior_oscillatory_jacobian_abs_max",
        "exterior_total_minus_leading_jacobian_abs_max",
    ):
        value = float(diagnostic[key])
        if not np.isfinite(value) or value > float(gates[key]):
            raise AssertionError(f"frozen mechanical gate failed: {key}={value}")
    signal = float(diagnostic["oscillatory_vorticity_vector_rms"])
    if not np.isfinite(signal) or signal < float(gates["oscillatory_vorticity_vector_rms_min"]):
        raise AssertionError("oscillatory vorticity signal became trivial")
    if diagnostic.get("beyond_current_Xh_fail_closed") is not True:
        raise AssertionError("partial composite differential surface no longer fails closed beyond X_h")

    truth = receipt.get("semantic_payload", {}).get("truth_boundary", {})
    required_false = (
        "velocity_beyond_Xh_materialized",
        "outer_global_leading_velocity_materialized",
        "global_compact_support_completed",
        "full_concrete_oscillatory_runtime_digest_bound",
        "identity_preserving_composite_save_load_available",
        "mean_projection_performed",
        "radial_inverse_performed",
        "correction_velocity_constructed",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_velocity_pressure_forcing_api",
        "complete_ns_residual",
        "same_protocol_comparable_to_st006",
        "residual_reduction_claimed",
        "velocity_export_ready",
        "paper_exact",
        "openai_field_identified",
        "pde_validated",
    )
    for key in required_false:
        if truth.get(key) is not False:
            raise AssertionError(f"truth boundary drifted at {key}")
