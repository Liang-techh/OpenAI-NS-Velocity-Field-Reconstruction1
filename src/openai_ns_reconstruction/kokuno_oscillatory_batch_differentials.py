"""Axis-safe batch spatial differentials for the frozen Kokuno A2 oscillation.

This increment extends the already-frozen #954 point-cloud velocity adapter with a
parameter-free differential surface.  It does not alter the oscillatory field.
Strict-support interior points reuse A2's existing centered Cartesian FD6
vorticity diagnostic at the fixed repository step already consumed by #874;
axis/support-exterior points are returned as exact zero before that differential
routine is entered.

The corrected 2026-09-09 Kokuno reconstruction remains structural provenance for
localized waves / complete curls.  The public-z pullback, autonomous parameters,
FD6 spatial derivative and this batch adapter are repository realizations and
are not claimed paper-exact.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from typing import Any, Mapping

import numpy as np

from .kokuno_oscillatory_batch_axis_safety import (
    _broadcast_time,
    _support,
    batch_axis_safety_sha256,
    velocity_osc_batch,
)
from .kokuno_public_oscillatory_vorticity_diagnostic import evaluate_vorticity_osc_fd6

TASK = "KOKUNO-A2-BATCH-DIFFERENTIALS-079"
SCHEMA = "kokuno-a2-oscillatory-batch-differentials-v1"
PARENT_AGENT2_PR = 954
PARENT_AGENT2_HEAD = "e43c32c3d0258f635e8b1f999d89c38524355856"
PARENT_BATCH_BLOB = "599baa190ec742d0e532e5517078534124e68d6b"
OSCILLATORY_DIAGNOSTIC_PR = 643
FIXED_STEP_CONSUMER_PR = 874
OSCILLATORY_JACOBIAN_SOURCE_BLOB = "4ae525bad1e4b83c9dcabc1a97d3931562d099a5"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

# This is deliberately not a public input.  It is the existing fixed A2 FD6
# spatial step used by the strict-inner spatial candidate (#874).
FIXED_SPATIAL_STEP = 1.0e-3
SCALAR_REPLAY_ATOL = 1.0e-12
INTERIOR_VORTICITY_RMS_FLOOR = 1.0e-12

_RECEIPT_POINTS = np.asarray(
    [
        (0.42, 0.19, -1.10),
        (0.63, -0.27, -0.45),
        (-0.51, 0.38, 0.15),
        (-0.74, -0.22, 0.65),
        (0.91, 0.31, 1.10),
        (-1.02, 0.24, -0.80),
        (0.36, -0.44, 0.92),
        (-0.58, -0.49, -1.32),
    ],
    dtype=float,
)
_RECEIPT_TIMES = np.asarray((0.31, 0.37, 0.43, 0.49, 0.55, 0.61, 0.67, 0.71))


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


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


@dataclass(frozen=True)
class OscillatoryBatchDifferentialEvaluation:
    velocity: np.ndarray
    velocity_jacobian: np.ndarray
    divergence: np.ndarray
    vorticity: np.ndarray
    interior_mask: np.ndarray


def evaluate_oscillatory_batch_differentials(
    points: Any, time: Any
) -> OscillatoryBatchDifferentialEvaluation:
    """Return frozen ``u_osc``, FD6 Jacobian/divergence/curl on ``points[...,3]``.

    No derivative step or oscillatory parameter is caller-tunable.  The parent
    #954 time contract is checked before spatial masking.  The symmetry axis,
    support faces and support exterior return exact zero for every field and do
    not enter the cylindrical public evaluator or FD6 differential routine.
    """
    xyz = np.asarray(points, dtype=float)
    if xyz.ndim < 1 or xyz.shape[-1] != 3:
        raise ValueError("points must have shape [...,3]")
    if not np.all(np.isfinite(xyz)):
        raise ValueError("points must be finite")

    lead_shape = xyz.shape[:-1]
    tt = _broadcast_time(time, lead_shape)
    velocity = np.asarray(velocity_osc_batch(xyz, tt), dtype=float)
    if velocity.shape != lead_shape + (3,) or not np.all(np.isfinite(velocity)):
        raise RuntimeError("parent batch oscillatory velocity returned an invalid shape/value")

    flat_xyz = xyz.reshape((-1, 3))
    flat_t = tt.reshape(-1)
    flat_velocity = velocity.reshape((-1, 3))
    jacobian = np.zeros((flat_xyz.shape[0], 3, 3), dtype=float)
    divergence = np.zeros(flat_xyz.shape[0], dtype=float)
    vorticity = np.zeros((flat_xyz.shape[0], 3), dtype=float)

    if flat_xyz.shape[0] == 0:
        return OscillatoryBatchDifferentialEvaluation(
            velocity=velocity,
            velocity_jacobian=jacobian.reshape(lead_shape + (3, 3)),
            divergence=divergence.reshape(lead_shape),
            vorticity=vorticity.reshape(lead_shape + (3,)),
            interior_mask=np.zeros(lead_shape, dtype=bool),
        )

    r0, r1, z0, z1 = _support()
    radius = np.hypot(flat_xyz[:, 0], flat_xyz[:, 1])
    inside = (
        (radius > r0)
        & (radius < r1)
        & (flat_xyz[:, 2] > z0)
        & (flat_xyz[:, 2] < z1)
    )

    if np.any(inside):
        p = flat_xyz[inside]
        diag = evaluate_vorticity_osc_fd6(
            p[:, 0], p[:, 1], p[:, 2], flat_t[inside],
            spatial_step=FIXED_SPATIAL_STEP,
        )
        diag_velocity = np.asarray(diag["velocity"], dtype=float)
        diag_jacobian = np.asarray(diag["velocity_gradient_fd6"], dtype=float)
        diag_divergence = np.asarray(diag["divergence_fd6"], dtype=float)
        diag_vorticity = np.asarray(diag["vorticity_fd6"], dtype=float)
        count = p.shape[0]
        expected = {
            "velocity": (count, 3),
            "jacobian": (count, 3, 3),
            "divergence": (count,),
            "vorticity": (count, 3),
        }
        actual = {
            "velocity": diag_velocity.shape,
            "jacobian": diag_jacobian.shape,
            "divergence": diag_divergence.shape,
            "vorticity": diag_vorticity.shape,
        }
        if actual != expected:
            raise RuntimeError(f"oscillatory FD6 diagnostic returned unexpected shapes: {actual}")
        if not all(
            np.all(np.isfinite(a))
            for a in (diag_velocity, diag_jacobian, diag_divergence, diag_vorticity)
        ):
            raise RuntimeError("oscillatory FD6 diagnostic became non-finite")
        if np.max(np.abs(diag_velocity - flat_velocity[inside]), initial=0.0) > SCALAR_REPLAY_ATOL:
            raise RuntimeError("FD6 diagnostic velocity disagrees with parent batch velocity")
        jacobian[inside] = diag_jacobian
        divergence[inside] = diag_divergence
        vorticity[inside] = diag_vorticity

    return OscillatoryBatchDifferentialEvaluation(
        velocity=velocity,
        velocity_jacobian=jacobian.reshape(lead_shape + (3, 3)),
        divergence=divergence.reshape(lead_shape),
        vorticity=vorticity.reshape(lead_shape + (3,)),
        interior_mask=inside.reshape(lead_shape),
    )


def semantic_payload() -> dict[str, Any]:
    r0, r1, z0, z1 = _support()
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "parent_batch_blob": PARENT_BATCH_BLOB,
        "parent_batch_semantic_sha256": batch_axis_safety_sha256(),
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
            "scope": "localized oscillatory fields and complete-curl organization",
        },
        "differential_realization": {
            "source_diagnostic_pr": OSCILLATORY_DIAGNOSTIC_PR,
            "fixed_step_consumer_pr": FIXED_STEP_CONSUMER_PR,
            "oscillatory_jacobian_source_blob": OSCILLATORY_JACOBIAN_SOURCE_BLOB,
            "operator": "centered_cartesian_fd6",
            "fixed_spatial_step": FIXED_SPATIAL_STEP,
            "caller_tunable_step": False,
            "jacobian_convention": "J[...,component,axis]=partial_axis velocity_component",
            "vorticity": ["J[2,1]-J[1,2]", "J[0,2]-J[2,0]", "J[1,0]-J[0,1]"],
        },
        "strict_support": {
            "radial_inner": r0,
            "radial_outer": r1,
            "axial_lower": z0,
            "axial_upper": z1,
            "axis_and_exterior_exact_zero": True,
        },
        "truth_boundary": {
            "source_formulas_changed": False,
            "oscillatory_velocity_changed": False,
            "new_derivative_realization_introduced": False,
            "batch_derivative_adapter_introduced": True,
            "global_leading_velocity_materialized": False,
            "mean_projection_performed": False,
            "radial_inverse_performed": False,
            "correction_velocity_constructed": False,
            "pressure_or_forcing_added": False,
            "complete_ns_residual": False,
            "residual_reduction_claimed": False,
            "paper_exact": False,
            "pde_validated": False,
        },
    }


def batch_differentials_sha256() -> str:
    return _sha256(semantic_payload())


def public_contract() -> dict[str, Any]:
    params = inspect.signature(evaluate_oscillatory_batch_differentials).parameters
    forbidden = {
        "amplitude", "phase", "phase_offset", "scale", "orientation", "support",
        "spatial_step", "derivative_step", "residual", "target", "forcing", "pressure",
        "viscosity", "nu", "gain", "threshold", "mean", "stress", "inverse",
    }
    return {
        "public_inputs": list(params),
        "forbidden_inputs_present": sorted(forbidden.intersection(params)),
        "shape_preserving_batch_api": True,
        "fixed_spatial_step": FIXED_SPATIAL_STEP,
        "caller_tunable_step": False,
        "axis_safe_by_parent_strict_support_mask": True,
        "complete_curl_reimplemented": False,
        "new_derivative_realization_introduced": False,
        "mean_projection_performed": False,
        "radial_inverse_performed": False,
        "correction_velocity_constructed": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def materialize_batch_differentials_receipt() -> dict[str, Any]:
    """Materialize deterministic scalar-replay and axis/support-safety checks."""
    batch = evaluate_oscillatory_batch_differentials(_RECEIPT_POINTS, _RECEIPT_TIMES)
    scalar = [
        evaluate_vorticity_osc_fd6(x, y, z, t, spatial_step=FIXED_SPATIAL_STEP)
        for (x, y, z), t in zip(_RECEIPT_POINTS, _RECEIPT_TIMES, strict=True)
    ]
    scalar_jacobian = np.stack(
        [np.asarray(item["velocity_gradient_fd6"], dtype=float).reshape(3, 3) for item in scalar]
    )
    scalar_vorticity = np.stack(
        [np.asarray(item["vorticity_fd6"], dtype=float).reshape(3) for item in scalar]
    )
    scalar_divergence = np.asarray(
        [float(np.asarray(item["divergence_fd6"], dtype=float)) for item in scalar]
    )

    r0, r1, z0, z1 = _support()
    tiny = np.nextafter(0.0, 1.0)
    masked_points = np.asarray(
        [
            (0.0, 0.0, 0.0),
            (tiny, -tiny, 0.25),
            (0.5 * r0, 0.0, -0.4),
            (r0, 0.0, 0.0),
            (r1, 0.0, 0.0),
            (0.5 * (r0 + r1), 0.0, z0),
            (0.5 * (r0 + r1), 0.0, z1),
            (1.01 * r1, 0.0, 0.3),
        ],
        dtype=float,
    )
    masked = evaluate_oscillatory_batch_differentials(masked_points, 0.47)

    curl_replay = _curl_from_jacobian(batch.velocity_jacobian)
    trace_replay = np.trace(batch.velocity_jacobian, axis1=-2, axis2=-1)
    masked_abs_max = max(
        float(np.max(np.abs(masked.velocity), initial=0.0)),
        float(np.max(np.abs(masked.velocity_jacobian), initial=0.0)),
        float(np.max(np.abs(masked.divergence), initial=0.0)),
        float(np.max(np.abs(masked.vorticity), initial=0.0)),
    )
    diagnostic = {
        "interior_probe_count": int(_RECEIPT_POINTS.shape[0]),
        "masked_probe_count": int(masked_points.shape[0]),
        "jacobian_scalar_replay_max_abs": float(np.max(np.abs(batch.velocity_jacobian - scalar_jacobian))),
        "vorticity_scalar_replay_max_abs": float(np.max(np.abs(batch.vorticity - scalar_vorticity))),
        "divergence_scalar_replay_max_abs": float(np.max(np.abs(batch.divergence - scalar_divergence))),
        "vorticity_curl_closure_max_abs": float(np.max(np.abs(batch.vorticity - curl_replay))),
        "divergence_trace_closure_max_abs": float(np.max(np.abs(batch.divergence - trace_replay))),
        "masked_all_fields_abs_max": masked_abs_max,
        "interior_vorticity_vector_rms": float(np.sqrt(np.mean(np.sum(batch.vorticity * batch.vorticity, axis=-1)))),
        "batch_velocity_shape": list(batch.velocity.shape),
        "batch_jacobian_shape": list(batch.velocity_jacobian.shape),
        "batch_vorticity_shape": list(batch.vorticity.shape),
    }
    return {
        "schema": "kokuno-a2-oscillatory-batch-differentials-receipt-v1",
        "semantic_sha256": batch_differentials_sha256(),
        "semantic_payload": semantic_payload(),
        "public_contract": public_contract(),
        "diagnostic": diagnostic,
        "frozen_mechanical_gates": {
            "scalar_replay_max_abs": SCALAR_REPLAY_ATOL,
            "curl_and_trace_closure_max_abs": 0.0,
            "masked_all_fields_abs_max": 0.0,
            "interior_vorticity_vector_rms_min": INTERIOR_VORTICITY_RMS_FLOOR,
        },
    }
