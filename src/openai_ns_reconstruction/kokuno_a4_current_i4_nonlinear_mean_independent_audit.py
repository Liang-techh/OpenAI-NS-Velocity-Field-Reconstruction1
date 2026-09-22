"""Implementation-distinct Agent-4 audit of the exact current-I4 nonlinear mean.

The production artifact is Agent-3 PR #1101.  This audit does not use its
leading FD2 / oscillatory-differential construction path.  After an exact
Agent-2 PR #1080 save/load roundtrip (performed by the scientific workflow),
the reference sees only public Cartesian velocity(x,y,z,t) from the total
field and its identity-bound leading backend.  It reconstructs derivatives
with centered five-point FD4 and projects rotating cylindrical means with
nonuniform Gauss-Legendre quadrature.

This is scoped nonlinear-attribution consistency evidence only.  It is not a
complete Navier-Stokes residual, a correction authorization, or PDE admission.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

TASK = "K4-VAL-113"
SCHEMA = "kokuno-a4-current-i4-nonlinear-mean-independent-audit-v1"

PARENT_AGENT3_PR = 1101
PARENT_AGENT3_HEAD = "6c9a4e7b0dec4b3669fb1dc60297b18f10163e0b"
PARENT_AGENT3_SOURCE_BLOB = "10d067570ddba834db3f19c9859c2efdad16b3e8"
PARENT_BACKEND_KIND = "exact-current-a2-pr1080-through-i4-nonlinear-attribution"

AGENT2_PR = 1080
AGENT2_HEAD = "c40d8ddecd2971544a6e07dab093436b423cf326"
AGENT2_MODULE = "openai_ns_reconstruction.kokuno_current_i4_leading_oscillatory_identity"
AGENT2_CLASS = "CurrentI4LeadingOscillatoryField"
AGENT2_SOURCE_BLOB = "2a0a5aa5966b02da856bdcf51940f3c186042802"

AGENT1_PR = 1079
AGENT1_HEAD = "b06742ca6e189499192ede3cce40f62cdc1e35ca"
AGENT1_MODULE = "openai_ns_reconstruction.kokuno_pa16_current_cartesian_i4_leading_preservation"
AGENT1_CLASS = "KokunoPA16CurrentCartesianI4LeadingPreservation"
AGENT1_SOURCE_BLOB = "6f04ce0a856b44430402576dad88438da90d1ebb"

FROZEN_SEED = 9173851
FROZEN_SPATIAL_STEPS = (0.012, 0.006, 0.003)
FROZEN_ANGULAR_ORDERS = (20, 40, 80)
FROZEN_AXIS_RADII = (1.0e-4, 0.005, 0.02)
FROZEN_I4_FRACTIONS = (0.23, 0.61, 0.84)
FROZEN_I4_ETAS = (-0.27, 0.11, 0.34)
FROZEN_I4_TIMES = (0.41, 0.56, 0.68)
FROZEN_INNER_RING_COUNT = 6

FINE_REL_RMS_GATE = 5.0e-2
FINE_REL_WEIGHTED_L2_GATE = 5.0e-2
FINE_REL_MAX_GATE = 1.5e-1
AXIS_NEAR_NORMALIZED_ERROR_GATE = 1.5e-1
FINE_PAIR_REL_RMS_GATE = 5.0e-2
RESOLUTION_DEGRADATION_FACTOR = 1.25
RESOLUTION_NUMERICAL_FLOOR = 2.0e-5
PIECE_CLOSURE_RELATIVE_GATE = 5.0e-10
NONTRIVIAL_MEAN_RMS_FLOOR = 1.0e-12
LOW_SIGNAL_SCALE_FLOOR = 1.0e-10
LOW_SIGNAL_ABSOLUTE_GATE = 1.0e-10

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5

PIECES = ("A", "B", "mixed", "quadratic", "aggregate")
SCIENTIFIC_PIECES = ("mixed", "quadratic", "aggregate")
RECEIPT_KEYS = {
    "A": "mean_inner_advects_oscillation_cylindrical",
    "B": "mean_oscillation_advects_inner_cylindrical",
    "mixed": "mean_mixed_cross_cylindrical",
    "quadratic": "mean_oscillatory_self_advection_cylindrical",
    "aggregate": "mean_aggregate_nonlinear_cylindrical",
}


@dataclass(frozen=True)
class RingSpec:
    label: str
    kind: str
    radius: float
    z: float
    t: float


def _git_blob_sha1(path: str | Path) -> str:
    data = Path(path).read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _valid_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(c in "0123456789abcdef" for c in value)
    )


def frozen_inner_ring_specs() -> tuple[RingSpec, ...]:
    rng = np.random.default_rng(FROZEN_SEED)
    return tuple(
        RingSpec(
            label=f"inner_{i}",
            kind="inner_offgrid",
            radius=float(rng.uniform(0.355, 0.525)),
            z=float(rng.uniform(-0.145, -0.015)),
            t=float(rng.uniform(0.335, 0.445)),
        )
        for i in range(FROZEN_INNER_RING_COUNT)
    )


def frozen_axis_near_specs() -> tuple[RingSpec, ...]:
    return tuple(
        RingSpec(
            label=f"axis_near_{i}",
            kind="positive_radius_axis_near",
            radius=float(r),
            z=(-0.035 + 0.035 * i),
            t=0.39,
        )
        for i, r in enumerate(FROZEN_AXIS_RADII)
    )


def frozen_i4_specs(field: object) -> tuple[RingSpec, ...]:
    x0 = float(getattr(field, "X_I4_start"))
    x1 = float(getattr(field, "X_I4_end"))
    d = float(getattr(field, "D"))
    if not (np.isfinite(x0) and np.isfinite(x1) and np.isfinite(d) and 0.0 < x0 < x1):
        raise ValueError("authenticated current-I4 geometry is invalid")
    log0, log1 = math.log(x0), math.log(x1)
    out: list[RingSpec] = []
    for i, (frac, eta, time) in enumerate(
        zip(FROZEN_I4_FRACTIONS, FROZEN_I4_ETAS, FROZEN_I4_TIMES)
    ):
        X = math.exp(log0 + float(frac) * (log1 - log0))
        if not x0 < X < x1:
            raise ValueError("frozen I4 probe escaped strict I4 interior")
        q = (1.0 - float(time)) / (1.0 - float(eta) ** 2)
        if not np.isfinite(q) or q <= 0.0:
            raise ValueError("frozen I4 probe has invalid similarity q")
        radius = math.sqrt(2.0 * q * X)
        z = (q**d) * float(eta)
        if not (np.isfinite(radius) and np.isfinite(z)):
            raise ValueError("frozen I4 probe is not representable in float64")
        out.append(
            RingSpec(
                label=f"i4_{i}",
                kind="strict_interior_i4",
                radius=float(radius),
                z=float(z),
                t=float(time),
            )
        )
    return tuple(out)


def _as_velocity(value: Any, n: int, label: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.shape != (n, 3) or not np.all(np.isfinite(arr)):
        raise ValueError(f"{label} must return finite shape ({n},3), got {arr.shape}")
    return arr


def _velocity(
    field: object, x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray
) -> np.ndarray:
    fn = getattr(field, "velocity", None)
    if not callable(fn):
        raise TypeError("field must expose velocity(x,y,z,t)")
    return _as_velocity(fn(x, y, z, t), len(x), type(field).__name__)


def _jacobian_fd4(
    field: object,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    h: float,
) -> np.ndarray:
    if not (np.isfinite(h) and h > 0.0):
        raise ValueError("FD4 step must be finite and positive")
    base = (np.asarray(x, float), np.asarray(y, float), np.asarray(z, float))
    n = len(x)
    J = np.empty((n, 3, 3), dtype=float)
    for axis in range(3):
        values: dict[int, np.ndarray] = {}
        for mult in (-2, -1, 1, 2):
            coords = [a.copy() for a in base]
            coords[axis] = coords[axis] + mult * h
            values[mult] = _velocity(field, coords[0], coords[1], coords[2], t)
        J[:, :, axis] = (
            values[-2] - 8.0 * values[-1] + 8.0 * values[1] - values[2]
        ) / (12.0 * h)
    return J


def _cylindrical_components(v: np.ndarray, theta: np.ndarray) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    out = np.empty_like(v)
    out[:, 0] = v[:, 0] * c + v[:, 1] * s
    out[:, 1] = -v[:, 0] * s + v[:, 1] * c
    out[:, 2] = v[:, 2]
    return out


def independent_ring_means(
    total_field: object,
    leading_field: object,
    spec: RingSpec,
    *,
    spatial_step: float,
    angular_order: int,
) -> dict[str, np.ndarray]:
    if angular_order < 8:
        raise ValueError("angular order is too small for the frozen audit")
    nodes, weights = np.polynomial.legendre.leggauss(int(angular_order))
    theta = math.pi * (nodes + 1.0)
    avg_weights = 0.5 * weights
    r = float(spec.radius)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = np.full_like(x, float(spec.z))
    t = np.full_like(x, float(spec.t))

    u_total = _velocity(total_field, x, y, z, t)
    u_lead = _velocity(leading_field, x, y, z, t)
    u_osc = u_total - u_lead
    j_total = _jacobian_fd4(total_field, x, y, z, t, spatial_step)
    j_lead = _jacobian_fd4(leading_field, x, y, z, t, spatial_step)
    j_osc = j_total - j_lead

    A = np.einsum("nij,nj->ni", j_osc, u_lead)
    B = np.einsum("nij,nj->ni", j_lead, u_osc)
    Q = np.einsum("nij,nj->ni", j_osc, u_osc)
    cart = {
        "A": A,
        "B": B,
        "mixed": A + B,
        "quadratic": Q,
        "aggregate": A + B + Q,
    }
    return {
        piece: np.sum(avg_weights[:, None] * _cylindrical_components(values, theta), axis=0)
        for piece, values in cart.items()
    }


def _extract_receipt_vector(receipt: Mapping[str, Any], key: str) -> np.ndarray:
    arr = np.squeeze(np.asarray(receipt.get(key), dtype=float))
    if arr.shape != (3,) or not np.all(np.isfinite(arr)):
        raise ValueError(f"production receipt {key} must be one finite cylindrical vector")
    return arr


def _validate_receipt(
    receipt: Mapping[str, Any], spec: RingSpec
) -> dict[str, np.ndarray]:
    if receipt.get("backend_kind") != PARENT_BACKEND_KIND:
        raise ValueError("production receipt backend kind drifted")
    rr = np.asarray(receipt.get("radius"), dtype=float).reshape(-1)
    zz = np.asarray(receipt.get("z"), dtype=float).reshape(-1)
    tt = np.asarray(receipt.get("t"), dtype=float).reshape(-1)
    if rr.size != 1 or zz.size != 1 or tt.size != 1:
        raise ValueError("A4 expects exactly one ring per public production receipt")
    tol_r = max(2e-13, 4e-15 * abs(spec.radius))
    if abs(float(rr[0]) - spec.radius) > tol_r:
        raise ValueError(f"production radius drifted for {spec.label}")
    if abs(float(zz[0]) - spec.z) > 2e-13 or abs(float(tt[0]) - spec.t) > 2e-13:
        raise ValueError(f"production z/t drifted for {spec.label}")

    backend = receipt.get("agent2_backend")
    if not isinstance(backend, Mapping):
        raise ValueError("production receipt lost Agent-2 backend identity")
    expected = {
        "agent2_composite_pr": AGENT2_PR,
        "agent2_composite_head": AGENT2_HEAD,
        "agent1_leading_pr": AGENT1_PR,
        "agent1_leading_head": AGENT1_HEAD,
    }
    for key, value in expected.items():
        if backend.get(key) != value:
            raise ValueError(f"production receipt identity drifted at {key}")
    if float(receipt.get("pointwise_three_piece_closure_absolute_max")) > 5.0e-13:
        raise ValueError("production pointwise A+B+Q closure failed")
    if float(receipt.get("projected_three_piece_closure_absolute_max")) > 1.0e-12:
        raise ValueError("production projected A+B+Q closure failed")
    return {piece: _extract_receipt_vector(receipt, key) for piece, key in RECEIPT_KEYS.items()}


def _validate_a3_truth_boundary(boundary: Mapping[str, Any]) -> None:
    required_true = (
        "current_I4_leading_plus_oscillatory_identity_consumed",
        "current_I4_mixed_nonlinear_mean_materialized",
        "current_I4_quadratic_oscillatory_mean_materialized",
        "current_I4_aggregate_nonlinear_mean_materialized",
        "source_I4_reserved_mean_correction_interval_recorded",
    )
    required_false = (
        "source_I3_positive_order_correction_materialized",
        "source_I4_five_row_mean_correction_materialized",
        "current_I4_correction_velocity_materialized",
        "post_I4_or_pulse_leading_identity_consumed",
        "outer_global_velocity_consumed",
        "radial_inverse_performed_in_this_increment",
        "compact_radial_stress_materialized_in_this_increment",
        "radial_force_materialized_in_this_increment",
        "scoped_nonlinear_mean_authorized_as_complete_ns_defect",
        "scoped_nonlinear_mean_authorized_as_correction_target",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    )
    for key in required_true:
        if boundary.get(key) is not True:
            raise ValueError(f"A3 #1101 truth boundary drifted at {key}")
    for key in required_false:
        if boundary.get(key) is not False:
            raise ValueError(f"A3 #1101 truth boundary drifted at {key}")
    if boundary.get("forbidden_scientific_controls_exposed") is not False:
        raise ValueError("A3 #1101 unexpectedly exposes a scientific tuning control")
    if float(boundary.get("final_normalized_momentum_gate")) != FINAL_NORMALIZED_MOMENTUM_GATE:
        raise ValueError("final momentum gate drifted")
    if float(boundary.get("final_normalized_divergence_gate")) != FINAL_NORMALIZED_DIVERGENCE_GATE:
        raise ValueError("final divergence gate drifted")


def validate_candidate_identity(field: object) -> object:
    ftype = type(field)
    if ftype.__module__ != AGENT2_MODULE or ftype.__name__ != AGENT2_CLASS:
        raise ValueError("unexpected A2 #1080 candidate runtime identity")
    source = inspect.getsourcefile(ftype)
    if source is None or _git_blob_sha1(source) != AGENT2_SOURCE_BLOB:
        raise ValueError("A2 #1080 source blob drifted")
    semantic = getattr(field, "semantic_sha256", None)
    if not _valid_sha256(semantic):
        raise ValueError("A2 #1080 semantic identity is malformed")
    leading = getattr(field, "leading_backend", None)
    if leading is None:
        raise ValueError("A2 #1080 candidate lost its leading backend")
    ltype = type(leading)
    if ltype.__module__ != AGENT1_MODULE or ltype.__name__ != AGENT1_CLASS:
        raise ValueError("unexpected A1 #1079 leading runtime identity")
    lsource = inspect.getsourcefile(ltype)
    if lsource is None or _git_blob_sha1(lsource) != AGENT1_SOURCE_BLOB:
        raise ValueError("A1 #1079 source blob drifted")
    return leading


def _stack_piece(level: Sequence[Mapping[str, np.ndarray]], piece: str) -> np.ndarray:
    return np.stack([np.asarray(item[piece], float) for item in level], axis=0)


def _rms(a: np.ndarray) -> float:
    x = np.asarray(a, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def _relative_rms(a: np.ndarray, b: np.ndarray) -> float:
    x, y = np.asarray(a, float), np.asarray(b, float)
    return _rms(x - y) / max(_rms(x), _rms(y), np.finfo(float).tiny)


def _comparison_metrics(
    production: np.ndarray, independent: np.ndarray, radii: np.ndarray
) -> dict[str, Any]:
    p, q, r = np.asarray(production, float), np.asarray(independent, float), np.asarray(radii, float)
    if p.shape != q.shape or p.ndim != 2 or p.shape[1] != 3 or r.shape != (p.shape[0],):
        raise ValueError("comparison arrays/radii have incompatible shapes")
    d = p - q
    weights = np.maximum(r, np.finfo(float).tiny)
    weights = weights / np.max(weights)
    weights = weights / np.sum(weights)

    def wl2(v: np.ndarray) -> float:
        return float(np.sqrt(np.sum(weights * np.sum(v * v, axis=1))))

    p_l2, q_l2, d_l2 = wl2(p), wl2(q), wl2(d)
    rms_scale = max(_rms(p), _rms(q), np.finfo(float).tiny)
    max_scale = max(float(np.max(np.abs(p))), float(np.max(np.abs(q))), np.finfo(float).tiny)
    index = np.unravel_index(int(np.argmax(np.abs(d))), d.shape)
    return {
        "relative_rms": _rms(d) / rms_scale,
        "relative_weighted_l2": d_l2 / max(p_l2, q_l2, np.finfo(float).tiny),
        "relative_max": float(np.max(np.abs(d))) / max_scale,
        "absolute_max_error": float(np.max(np.abs(d))),
        "signal_scale": max(float(np.max(np.abs(p))), float(np.max(np.abs(q)))),
        "production_rms": _rms(p),
        "independent_rms": _rms(q),
        "production_weighted_l2": p_l2,
        "independent_weighted_l2": q_l2,
        "worst_ring_index": int(index[0]),
        "worst_component": ("radial", "tangential", "axial")[int(index[1])],
        "worst_signed_error": float(d[index]),
    }


def _piece_closure(level: Sequence[Mapping[str, np.ndarray]]) -> float:
    mixed = _stack_piece(level, "mixed")
    q = _stack_piece(level, "quadratic")
    agg = _stack_piece(level, "aggregate")
    scale = max(_rms(agg), _rms(mixed), _rms(q), np.finfo(float).tiny)
    return _rms(agg - mixed - q) / scale


def _ab_closure(level: Sequence[Mapping[str, np.ndarray]]) -> float:
    A, B, mixed = (
        _stack_piece(level, "A"),
        _stack_piece(level, "B"),
        _stack_piece(level, "mixed"),
    )
    scale = max(_rms(A), _rms(B), _rms(mixed), np.finfo(float).tiny)
    return _rms(mixed - A - B) / scale


def _resolution_report(
    levels: Mapping[str, Sequence[Mapping[str, np.ndarray]]]
) -> dict[str, Any]:
    labels = list(levels)
    if len(labels) != 3:
        raise ValueError("exactly three resolution levels are required")
    out: dict[str, Any] = {"labels": labels, "pieces": {}}
    for piece in SCIENTIFIC_PIECES:
        a = _stack_piece(levels[labels[0]], piece)
        b = _stack_piece(levels[labels[1]], piece)
        c = _stack_piece(levels[labels[2]], piece)
        out["pieces"][piece] = {
            "coarse_to_medium_relative_rms": _relative_rms(a, b),
            "medium_to_fine_relative_rms": _relative_rms(b, c),
            "coarse_to_medium_absolute_max": float(np.max(np.abs(a - b))),
            "medium_to_fine_absolute_max": float(np.max(np.abs(b - c))),
            "fine_signal_scale": float(np.max(np.abs(c))),
        }
    return out


def _axis_metrics(
    production_level: Sequence[Mapping[str, np.ndarray]],
    independent_level: Sequence[Mapping[str, np.ndarray]],
) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for piece in SCIENTIFIC_PIECES:
        p, q = _stack_piece(production_level, piece), _stack_piece(independent_level, piece)
        dmax = float(np.max(np.abs(p - q)))
        scale = max(
            float(np.max(np.abs(p))),
            float(np.max(np.abs(q))),
            LOW_SIGNAL_SCALE_FLOOR,
        )
        out[piece] = {
            "absolute_max_error": dmax,
            "normalized_error": dmax / scale,
            "signal_scale": max(float(np.max(np.abs(p))), float(np.max(np.abs(q)))),
        }
    return out


def audit_current_i4_nonlinear_mean(
    field: object,
    production_receipts: Mapping[str, Mapping[str, Any]],
    a3_truth_boundary: Mapping[str, Any],
) -> dict[str, Any]:
    leading = validate_candidate_identity(field)
    _validate_a3_truth_boundary(a3_truth_boundary)

    inner_specs = frozen_inner_ring_specs()
    axis_specs = frozen_axis_near_specs()
    i4_specs = frozen_i4_specs(field)
    all_specs = inner_specs + axis_specs + i4_specs
    if set(production_receipts) != {s.label for s in all_specs}:
        raise ValueError("production receipt labels do not match frozen A4 probe set")
    production = {
        spec.label: _validate_receipt(production_receipts[spec.label], spec)
        for spec in all_specs
    }

    def levels_for(specs: Sequence[RingSpec], mode: str) -> dict[str, list[dict[str, np.ndarray]]]:
        out: dict[str, list[dict[str, np.ndarray]]] = {}
        if mode == "spatial":
            for h in FROZEN_SPATIAL_STEPS:
                out[f"h={h:.6g}"] = [
                    independent_ring_means(
                        field, leading, spec,
                        spatial_step=h, angular_order=FROZEN_ANGULAR_ORDERS[-1],
                    )
                    for spec in specs
                ]
        elif mode == "angular":
            for ntheta in FROZEN_ANGULAR_ORDERS:
                out[f"n={ntheta}"] = [
                    independent_ring_means(
                        field, leading, spec,
                        spatial_step=FROZEN_SPATIAL_STEPS[-1], angular_order=ntheta,
                    )
                    for spec in specs
                ]
        else:
            raise ValueError("unknown resolution family")
        return out

    inner_spatial = levels_for(inner_specs, "spatial")
    inner_angular = levels_for(inner_specs, "angular")
    i4_spatial = levels_for(i4_specs, "spatial")
    i4_angular = levels_for(i4_specs, "angular")
    fine_inner = inner_angular[f"n={FROZEN_ANGULAR_ORDERS[-1]}"]
    fine_i4 = i4_angular[f"n={FROZEN_ANGULAR_ORDERS[-1]}"]

    prod_inner = [production[s.label] for s in inner_specs]
    prod_i4 = [production[s.label] for s in i4_specs]
    inner_r = np.asarray([s.radius for s in inner_specs], dtype=float)
    i4_r = np.asarray([s.radius for s in i4_specs], dtype=float)

    inner_comparisons = {
        piece: _comparison_metrics(
            _stack_piece(prod_inner, piece), _stack_piece(fine_inner, piece), inner_r
        )
        for piece in PIECES
    }
    i4_comparisons = {
        piece: _comparison_metrics(
            _stack_piece(prod_i4, piece), _stack_piece(fine_i4, piece), i4_r
        )
        for piece in PIECES
    }
    for comparisons, specs in ((inner_comparisons, inner_specs), (i4_comparisons, i4_specs)):
        for piece in PIECES:
            comparisons[piece]["worst_ring_label"] = specs[
                comparisons[piece]["worst_ring_index"]
            ].label

    axis_independent = [
        independent_ring_means(
            field, leading, spec,
            spatial_step=FROZEN_SPATIAL_STEPS[-1],
            angular_order=FROZEN_ANGULAR_ORDERS[-1],
        )
        for spec in axis_specs
    ]
    axis_production = [production[s.label] for s in axis_specs]

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "agent2_pr": AGENT2_PR,
        "agent2_head": AGENT2_HEAD,
        "agent1_pr": AGENT1_PR,
        "agent1_head": AGENT1_HEAD,
        "candidate_semantic_sha256": str(getattr(field, "semantic_sha256")),
        "protocol": {
            "seed": FROZEN_SEED,
            "spatial_operator": "centered_cartesian_fd4_five_point",
            "spatial_steps": list(FROZEN_SPATIAL_STEPS),
            "angular_operator": "gauss_legendre_nonuniform_rotating_cylindrical_mean",
            "angular_orders": list(FROZEN_ANGULAR_ORDERS),
            "inner_ring_count": len(inner_specs),
            "axis_near_radii": list(FROZEN_AXIS_RADII),
            "i4_stage_fractions": list(FROZEN_I4_FRACTIONS),
            "i4_nontriviality_required": False,
            "low_signal_scale_floor": LOW_SIGNAL_SCALE_FLOOR,
            "low_signal_absolute_gate": LOW_SIGNAL_ABSOLUTE_GATE,
        },
        "inner_specs": [spec.__dict__ for spec in inner_specs],
        "axis_near_specs": [spec.__dict__ for spec in axis_specs],
        "i4_specs": [spec.__dict__ for spec in i4_specs],
        "inner_production_vs_independent": inner_comparisons,
        "i4_production_vs_independent": i4_comparisons,
        "inner_spatial_resolution": _resolution_report(inner_spatial),
        "inner_angular_resolution": _resolution_report(inner_angular),
        "i4_spatial_resolution": _resolution_report(i4_spatial),
        "i4_angular_resolution": _resolution_report(i4_angular),
        "inner_independent_piece_closure_relative_rms": _piece_closure(fine_inner),
        "inner_independent_A_plus_B_closure_relative_rms": _ab_closure(fine_inner),
        "i4_independent_piece_closure_relative_rms": _piece_closure(fine_i4),
        "i4_independent_A_plus_B_closure_relative_rms": _ab_closure(fine_i4),
        "axis_near": _axis_metrics(axis_production, axis_independent),
        "inner_production_quadratic_rms": _rms(_stack_piece(prod_inner, "quadratic")),
        "inner_independent_quadratic_rms": _rms(_stack_piece(fine_inner, "quadratic")),
        "inner_production_aggregate_rms": _rms(_stack_piece(prod_inner, "aggregate")),
        "inner_independent_aggregate_rms": _rms(_stack_piece(fine_inner, "aggregate")),
        "scientific_boundary": {
            "complete_ns_defect": False,
            "matched_pressure_materialized": False,
            "restricted_forcing_materialized": False,
            "correction_velocity_materialized": False,
            "leading_only_full_ns_residual_assessed": False,
            "leading_plus_oscillatory_full_ns_residual_assessed": False,
            "after_correction_residual_assessed": False,
            "heldout_normalized_ns_residual_assessed": False,
            "same_protocol_comparable_to_st006": False,
            "pde_validated": False,
            "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
            "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        },
    }
    report["audit_pass"] = scoped_gates_pass(report)
    return report


def _resolution_piece_pass(metrics: Mapping[str, Any]) -> bool:
    scale = float(metrics["fine_signal_scale"])
    cm = float(metrics["coarse_to_medium_relative_rms"])
    mf = float(metrics["medium_to_fine_relative_rms"])
    if scale < LOW_SIGNAL_SCALE_FLOOR:
        return float(metrics["medium_to_fine_absolute_max"]) <= LOW_SIGNAL_ABSOLUTE_GATE
    if mf > FINE_PAIR_REL_RMS_GATE:
        return False
    if cm > RESOLUTION_NUMERICAL_FLOOR and mf > RESOLUTION_NUMERICAL_FLOOR:
        if mf > RESOLUTION_DEGRADATION_FACTOR * cm:
            return False
    return True


def _comparison_pass(metrics: Mapping[str, Any]) -> bool:
    if float(metrics["signal_scale"]) < LOW_SIGNAL_SCALE_FLOOR:
        return float(metrics["absolute_max_error"]) <= LOW_SIGNAL_ABSOLUTE_GATE
    return (
        float(metrics["relative_rms"]) <= FINE_REL_RMS_GATE
        and float(metrics["relative_weighted_l2"]) <= FINE_REL_WEIGHTED_L2_GATE
        and float(metrics["relative_max"]) <= FINE_REL_MAX_GATE
    )


def scoped_gates_pass(report: Mapping[str, Any]) -> bool:
    try:
        for piece in SCIENTIFIC_PIECES:
            if not _comparison_pass(report["inner_production_vs_independent"][piece]):
                return False
            if not _comparison_pass(report["i4_production_vs_independent"][piece]):
                return False
            if float(report["axis_near"][piece]["normalized_error"]) > AXIS_NEAR_NORMALIZED_ERROR_GATE:
                return False
        for family in (
            "inner_spatial_resolution",
            "inner_angular_resolution",
            "i4_spatial_resolution",
            "i4_angular_resolution",
        ):
            for piece in SCIENTIFIC_PIECES:
                if not _resolution_piece_pass(report[family]["pieces"][piece]):
                    return False
        for key in (
            "inner_independent_piece_closure_relative_rms",
            "inner_independent_A_plus_B_closure_relative_rms",
            "i4_independent_piece_closure_relative_rms",
            "i4_independent_A_plus_B_closure_relative_rms",
        ):
            if float(report[key]) > PIECE_CLOSURE_RELATIVE_GATE:
                return False
        if min(
            float(report["inner_production_quadratic_rms"]),
            float(report["inner_independent_quadratic_rms"]),
            float(report["inner_production_aggregate_rms"]),
            float(report["inner_independent_aggregate_rms"]),
        ) < NONTRIVIAL_MEAN_RMS_FLOOR:
            return False
        boundary = report["scientific_boundary"]
        for key in (
            "complete_ns_defect",
            "matched_pressure_materialized",
            "restricted_forcing_materialized",
            "correction_velocity_materialized",
            "leading_only_full_ns_residual_assessed",
            "leading_plus_oscillatory_full_ns_residual_assessed",
            "after_correction_residual_assessed",
            "heldout_normalized_ns_residual_assessed",
            "same_protocol_comparable_to_st006",
            "pde_validated",
        ):
            if boundary.get(key) is not False:
                return False
        return (
            float(boundary["final_normalized_momentum_gate"]) == FINAL_NORMALIZED_MOMENTUM_GATE
            and float(boundary["final_normalized_divergence_gate"]) == FINAL_NORMALIZED_DIVERGENCE_GATE
        )
    except (KeyError, TypeError, ValueError, OverflowError):
        return False


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(audit_current_i4_nonlinear_mean)
    forbidden = {
        "threshold", "tolerance", "residual", "forcing", "pressure", "gain",
        "damping", "viscosity", "nu", "spatial_step", "angular_order",
        "normalized_score", "correction",
    }
    return {
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(forbidden.intersection(signature.parameters)),
        "implementation_distinct_spatial_operator": "centered_cartesian_fd4_five_point",
        "implementation_distinct_angular_operator": "gauss_legendre_nonuniform",
        "pde_validated": False,
    }
