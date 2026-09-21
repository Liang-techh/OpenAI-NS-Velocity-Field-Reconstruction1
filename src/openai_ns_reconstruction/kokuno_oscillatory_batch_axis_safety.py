"""Batch/axis-safe public evaluator for the frozen Kokuno A2 oscillation.

This module does not alter the oscillatory construction.  It wraps the existing
public ``velocity_osc(x,y,z,t)`` in a point-cloud API and makes one support fact
explicit: the registered oscillatory vector-potential coefficient has strict
annular/axial compact support, so points on the axis and outside that support are
exactly zero and need not enter the interior cylindrical realization.

The corrected 2026-09-09 Kokuno reconstruction remains structural provenance for
localized waves and complete curls.  The concrete public support, public-z
pullback, autonomous mode data, and this batch adapter are repository
realizations; they are not recovered paper-exact hidden values.
"""
from __future__ import annotations

import hashlib
import inspect
import json
from typing import Any, Mapping

import numpy as np

from .kokuno_public_z_pullback_velocity import default_field, velocity_osc

TASK = "KOKUNO-A2-BATCH-AXIS-SAFETY-078"
SCHEMA = "kokuno-a2-oscillatory-batch-axis-safety-v1"
PARENT_AGENT2_PR = 948
PARENT_AGENT2_HEAD = "fc89770a9ee2ed1897e53ca83f57b995e7b1250d"
PUBLIC_Z_VELOCITY_BLOB = "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

SCALAR_REPLAY_ATOL = 1.0e-12
INTERIOR_SIGNAL_FLOOR = 1.0e-12

# Fixed public-coordinate probes used only for the deterministic receipt.  They
# do not tune the field and are not an NS validation sample.
_RECEIPT_INTERIOR_POINTS = np.asarray(
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
_RECEIPT_INTERIOR_TIMES = np.asarray((0.31, 0.37, 0.43, 0.49, 0.55, 0.61, 0.67, 0.71))


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _support() -> tuple[float, float, float, float]:
    field = default_field()
    values = (
        float(field.radial_inner),
        float(field.radial_outer),
        float(field.axial_lower),
        float(field.axial_upper),
    )
    if not np.all(np.isfinite(values)):
        raise RuntimeError("oscillatory public support is non-finite")
    r0, r1, z0, z1 = values
    if not (0.0 < r0 < r1 and z0 < z1):
        raise RuntimeError("oscillatory public support is malformed")
    return values


def _broadcast_time(time: Any, shape: tuple[int, ...]) -> np.ndarray:
    raw = np.asarray(time, dtype=float)
    try:
        out = np.broadcast_to(raw, shape)
    except ValueError as exc:
        raise ValueError("time must broadcast to points.shape[:-1]") from exc
    if not np.all(np.isfinite(out)):
        raise ValueError("time must be finite")
    return np.asarray(out, dtype=float)


def velocity_osc_batch(points: Any, time: Any) -> np.ndarray:
    """Evaluate the frozen public oscillation on one ``[...,3]`` point cloud.

    The adapter is deliberately parameter-free.  Points outside the registered
    strict support, including the entire symmetry axis, are assigned exact zero
    before the interior public evaluator is called.  Strictly interior points are
    evaluated in one vectorized call to the pre-existing ``velocity_osc``.
    """
    xyz = np.asarray(points, dtype=float)
    if xyz.ndim < 1 or xyz.shape[-1] != 3:
        raise ValueError("points must have shape [...,3]")
    if not np.all(np.isfinite(xyz)):
        raise ValueError("points must be finite")

    lead_shape = xyz.shape[:-1]
    tt = _broadcast_time(time, lead_shape)
    flat_xyz = xyz.reshape((-1, 3))
    flat_t = tt.reshape(-1)
    out = np.zeros((flat_xyz.shape[0], 3), dtype=float)
    if flat_xyz.shape[0] == 0:
        return out.reshape(lead_shape + (3,))

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
        values = np.asarray(velocity_osc(p[:, 0], p[:, 1], p[:, 2], flat_t[inside]), dtype=float)
        if values.shape != (p.shape[0], 3) or not np.all(np.isfinite(values)):
            raise RuntimeError("public oscillatory velocity returned an invalid batch")
        out[inside] = values
    return out.reshape(lead_shape + (3,))


def semantic_payload() -> dict[str, Any]:
    r0, r1, z0, z1 = _support()
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "public_z_velocity_blob": PUBLIC_Z_VELOCITY_BLOB,
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
            "scope": "localized oscillatory fields and complete-curl organization",
        },
        "batch_contract": {
            "input": "finite points[...,3] plus finite broadcastable time",
            "output": "finite velocity[...,3]",
            "strict_support": {
                "radial_inner": r0,
                "radial_outer": r1,
                "axial_lower": z0,
                "axial_upper": z1,
            },
            "outside_support_behavior": "exact zero without entering interior cylindrical evaluator",
            "inside_support_behavior": "one vectorized call to existing public velocity_osc",
            "axis_safe_by_support_mask": True,
            "new_oscillatory_parameters": False,
        },
        "truth_boundary": {
            "source_formulas_changed": False,
            "oscillatory_velocity_changed": False,
            "support_changed": False,
            "leading_profile_constructed": False,
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


def batch_axis_safety_sha256() -> str:
    return _sha256(semantic_payload())


def public_contract() -> dict[str, Any]:
    params = inspect.signature(velocity_osc_batch).parameters
    forbidden = {
        "amplitude",
        "phase",
        "phase_offset",
        "scale",
        "orientation",
        "support",
        "residual",
        "target",
        "forcing",
        "pressure",
        "viscosity",
        "nu",
        "gain",
        "threshold",
    }
    return {
        "public_inputs": list(params),
        "forbidden_inputs_present": sorted(forbidden.intersection(params)),
        "vectorized_batch_api": True,
        "shape_preserving": True,
        "time_broadcast_supported": True,
        "axis_safe_by_strict_support_mask": True,
        "interior_uses_existing_public_velocity": True,
        "complete_curl_reimplemented": False,
        "mean_projection_performed": False,
        "radial_inverse_performed": False,
        "correction_velocity_constructed": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def materialize_batch_axis_safety_receipt() -> dict[str, Any]:
    """Materialize fixed scalar-parity and axis/support-safety diagnostics."""
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
    masked = velocity_osc_batch(masked_points, 0.47)
    interior = velocity_osc_batch(_RECEIPT_INTERIOR_POINTS, _RECEIPT_INTERIOR_TIMES)

    scalar = np.stack(
        [
            np.asarray(velocity_osc(x, y, z, t), dtype=float).reshape(3)
            for (x, y, z), t in zip(_RECEIPT_INTERIOR_POINTS, _RECEIPT_INTERIOR_TIMES)
        ],
        axis=0,
    )
    if not np.all(np.isfinite(scalar)):
        raise RuntimeError("scalar public oscillatory replay became non-finite")

    scalar_replay_max_abs = float(np.max(np.abs(interior - scalar)))
    masked_speed_max = float(np.max(np.linalg.norm(masked, axis=-1)))
    interior_vector_rms = float(np.sqrt(np.mean(np.sum(interior * interior, axis=-1))))
    payload = semantic_payload()
    return {
        "schema": "kokuno-a2-oscillatory-batch-axis-safety-receipt-v1",
        "semantic_sha256": _sha256(payload),
        "semantic_payload": payload,
        "public_contract": public_contract(),
        "diagnostic": {
            "interior_probe_count": int(_RECEIPT_INTERIOR_POINTS.shape[0]),
            "masked_probe_count": int(masked_points.shape[0]),
            "scalar_replay_max_abs": scalar_replay_max_abs,
            "masked_speed_max": masked_speed_max,
            "interior_vector_rms": interior_vector_rms,
            "batch_shape": list(interior.shape),
        },
        "frozen_mechanical_gates": {
            "scalar_replay_max_abs": SCALAR_REPLAY_ATOL,
            "masked_speed_max": 0.0,
            "interior_vector_rms_min": INTERIOR_SIGNAL_FLOOR,
        },
    }
