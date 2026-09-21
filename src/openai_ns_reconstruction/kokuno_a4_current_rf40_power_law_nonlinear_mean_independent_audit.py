"""Independent A4 audit of the current RF40 power-law nonlinear m=0 mean.

This module audits the public Agent-3 PR #1028 nonlinear-mean receipt without
using its construction derivative path.  The reference operator sees only the
save/reloaded Agent-2 PR #1010 public Cartesian total velocity and its
identity-bound public leading velocity.  It reconstructs Cartesian Jacobians
with a centered five-point fourth-order stencil and projects rotating
cylindrical means with nonuniform Gauss-Legendre angular quadrature.

Agent-3 production uses the repository-fixed leading FD2 path together with
Agent-2 PR #960 oscillatory differentials and an equispaced 32/64/128 angular
projector.  Therefore the numerical derivative and angular integration paths
here are implementation-distinct.

This remains correction-side attribution consistency evidence only.  It does
not materialize pressure, restricted forcing, a complete Navier--Stokes
defect, an authorized correction target, a Cartesian correction velocity, or
held-out full-PDE residual evidence.
"""
from __future__ import annotations

import hashlib
import inspect
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

TASK = "K4-VAL-105"
SCHEMA = "kokuno-a4-current-rf40-power-law-nonlinear-mean-independent-audit-v1"

PARENT_AGENT3_PR = 1028
PARENT_AGENT3_HEAD = "4b7a4400dcf1ee0be6bb07f23c2f8af6f893de88"
PARENT_AGENT3_SOURCE_BLOB = "4ec98a6adbb7afd16d6595cc24072fe18e62fa87"
PARENT_BACKEND_KIND = (
    "exact-current-a2-pr1010-through-rf40-power-law-nonlinear-attribution"
)

AGENT2_PR = 1010
AGENT2_HEAD = "e36d9da4b7f037e998f5b1658f8c0ea291a76b80"
AGENT2_MODULE = (
    "openai_ns_reconstruction.kokuno_current_rf40_power_law_leading_oscillatory_identity"
)
AGENT2_CLASS = "CurrentRF40PowerLawLeadingOscillatoryField"
AGENT2_SOURCE_BLOB = "4f1dc0566c041b9de20f18b4ae9c49d9fd93c58d"

AGENT1_PR = 1005
AGENT1_HEAD = "2c76ebdc41d6c566f43a1305034ba2ff9dce410b"
AGENT1_MODULE = "openai_ns_reconstruction.kokuno_pa16_current_cartesian_rf40_power_law"
AGENT1_CLASS = "KokunoPA16CurrentCartesianRF40PowerLaw"
AGENT1_SOURCE_BLOB = "36a0183352a2005bd0b5f477594e05a2e4fe3868"

FROZEN_SEED = 9173771
FROZEN_SPATIAL_STEPS = (0.012, 0.006, 0.003)
FROZEN_ANGULAR_ORDERS = (20, 40, 80)
FROZEN_AXIS_RADII = (1.0e-4, 0.005, 0.02)
FROZEN_POWER_STAGE_FRACTIONS = (0.22, 0.53, 0.81)
FROZEN_POWER_ETAS = (-0.31, 0.09, 0.27)
FROZEN_POWER_TIMES = (0.35, 0.53, 0.69)
FROZEN_INNER_RING_COUNT = 6

# Scoped attribution-consistency gates only; these are not PDE gates.
FINE_REL_RMS_GATE = 5.0e-2
FINE_REL_WEIGHTED_L2_GATE = 5.0e-2
FINE_REL_MAX_GATE = 1.5e-1
AXIS_NEAR_NORMALIZED_ERROR_GATE = 1.5e-1
FINE_PAIR_REL_RMS_GATE = 5.0e-2
RESOLUTION_DEGRADATION_FACTOR = 1.25
RESOLUTION_NUMERICAL_FLOOR = 2.0e-5
PIECE_CLOSURE_RELATIVE_GATE = 5.0e-10
NONTRIVIAL_MEAN_RMS_FLOOR = 1.0e-12
POWER_LAW_ZERO_SUPPORT_ABS_GATE = 1.0e-10

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
    """Fresh deterministic off-grid rings around the registered inner support."""
    rng = np.random.default_rng(FROZEN_SEED)
    specs: list[RingSpec] = []
    for i in range(FROZEN_INNER_RING_COUNT):
        specs.append(
            RingSpec(
                label=f"inner_{i}",
                kind="inner",
                radius=float(rng.uniform(0.36, 0.52)),
                z=float(rng.uniform(-0.14, -0.02)),
                t=float(rng.uniform(0.34, 0.44)),
            )
        )
    return tuple(specs)


def frozen_axis_near_specs() -> tuple[RingSpec, ...]:
    return tuple(
        RingSpec(
            label=f"axis_near_{i}",
            kind="axis_near",
            radius=float(r),
            z=-0.04 + 0.04 * i,
            t=0.39,
        )
        for i, r in enumerate(FROZEN_AXIS_RADII)
    )


def frozen_power_law_specs(field: object) -> tuple[RingSpec, ...]:
    """Construct strict X3<X<X4 rings from public candidate geometry."""
    X3 = float(getattr(field, "X_3"))
    X4 = float(getattr(field, "X_4"))
    Tw = float(getattr(field, "T_w"))
    D = float(getattr(field, "D"))
    specs: list[RingSpec] = []
    for i, (frac, eta, time) in enumerate(
        zip(FROZEN_POWER_STAGE_FRACTIONS, FROZEN_POWER_ETAS, FROZEN_POWER_TIMES)
    ):
        X = X3 * math.exp(float(frac) * Tw)
        if not X3 < X < X4:
            raise ValueError("frozen power-law ring escaped strict X3<X<X4")
        q = (1.0 - float(time)) / (1.0 - float(eta) ** 2)
        if q <= 0.0:
            raise ValueError("frozen power-law ring has nonpositive q")
        r = math.sqrt(2.0 * q * X)
        z = (q**D) * float(eta)
        specs.append(
            RingSpec(
                label=f"power_law_{i}", kind="power_law", radius=r, z=z, t=float(time)
            )
        )
    return tuple(specs)


def _as_velocity(value: Any, n: int, label: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.shape != (n, 3) or not np.all(np.isfinite(arr)):
        raise ValueError(f"{label} must return finite shape ({n},3), got {arr.shape}")
    return arr


def _velocity(field: object, x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray) -> np.ndarray:
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
    """Five-point centered Cartesian Jacobian J[i,j]=partial_j u_i."""
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
        deriv = (
            values[-2]
            - 8.0 * values[-1]
            + 8.0 * values[1]
            - values[2]
        ) / (12.0 * h)
        J[:, :, axis] = deriv
    return J


def _cylindrical_components(v: np.ndarray, theta: np.ndarray) -> np.ndarray:
    c = np.cos(theta)
    s = np.sin(theta)
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
    """Independent A/B/Q/N means using public velocities, FD4 and GL quadrature."""
    if angular_order < 8:
        raise ValueError("angular order is too small for the preregistered audit")
    nodes, weights = np.polynomial.legendre.leggauss(int(angular_order))
    theta = math.pi * (nodes + 1.0)
    avg_weights = 0.5 * weights  # (1/2pi) int_0^2pi = 1/2 int_-1^1
    r = float(spec.radius)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = np.full_like(x, float(spec.z))
    t = np.full_like(x, float(spec.t))

    u_total = _velocity(total_field, x, y, z, t)
    u_lead = _velocity(leading_field, x, y, z, t)
    u_osc = u_total - u_lead
    J_total = _jacobian_fd4(total_field, x, y, z, t, spatial_step)
    J_lead = _jacobian_fd4(leading_field, x, y, z, t, spatial_step)
    J_osc = J_total - J_lead

    A = np.einsum("nij,nj->ni", J_osc, u_lead)
    B = np.einsum("nij,nj->ni", J_lead, u_osc)
    Q = np.einsum("nij,nj->ni", J_osc, u_osc)
    mixed = A + B
    aggregate = mixed + Q
    cart = {"A": A, "B": B, "mixed": mixed, "quadratic": Q, "aggregate": aggregate}
    result: dict[str, np.ndarray] = {}
    for piece, values in cart.items():
        cyl = _cylindrical_components(values, theta)
        result[piece] = np.sum(avg_weights[:, None] * cyl, axis=0)
    return result


def _extract_receipt_vector(receipt: Mapping[str, Any], key: str) -> np.ndarray:
    arr = np.asarray(receipt.get(key), dtype=float)
    arr = np.squeeze(arr)
    if arr.shape != (3,) or not np.all(np.isfinite(arr)):
        raise ValueError(f"production receipt {key} must contain one finite cylindrical vector")
    return arr


def _validate_receipt(receipt: Mapping[str, Any], spec: RingSpec) -> dict[str, np.ndarray]:
    if receipt.get("backend_kind") != PARENT_BACKEND_KIND:
        raise ValueError("production receipt backend kind drifted")
    rr = np.asarray(receipt.get("radius"), dtype=float).reshape(-1)
    zz = np.asarray(receipt.get("z"), dtype=float).reshape(-1)
    tt = np.asarray(receipt.get("t"), dtype=float).reshape(-1)
    if rr.size != 1 or zz.size != 1 or tt.size != 1:
        raise ValueError("A4 expects exactly one ring per public production receipt")
    if abs(float(rr[0]) - spec.radius) > 2e-13:
        raise ValueError(f"production radius drifted for {spec.label}")
    if abs(float(zz[0]) - spec.z) > 2e-13 or abs(float(tt[0]) - spec.t) > 2e-13:
        raise ValueError(f"production z/t drifted for {spec.label}")

    backend = receipt.get("agent2_backend")
    if not isinstance(backend, Mapping):
        raise ValueError("production receipt lost Agent-2 backend identity")
    if backend.get("agent2_composite_pr") != AGENT2_PR or backend.get("agent2_composite_head") != AGENT2_HEAD:
        raise ValueError("production receipt A2 identity drifted")
    if backend.get("agent1_leading_pr") != AGENT1_PR or backend.get("agent1_leading_head") != AGENT1_HEAD:
        raise ValueError("production receipt A1 identity drifted")

    if float(receipt.get("pointwise_three_piece_closure_absolute_max")) > 5.0e-13:
        raise ValueError("production pointwise A+B+Q closure failed")
    if float(receipt.get("projected_three_piece_closure_absolute_max")) > 1.0e-12:
        raise ValueError("production projected A+B+Q closure failed")
    return {piece: _extract_receipt_vector(receipt, key) for piece, key in RECEIPT_KEYS.items()}


def _validate_a3_truth_boundary(boundary: Mapping[str, Any]) -> None:
    required_true = (
        "current_leading_plus_oscillation_through_RF40_power_law_consumed",
        "RF40_power_law_velocity_materialized",
        "current_RF40_power_law_mixed_nonlinear_mean_materialized",
        "current_RF40_power_law_quadratic_nonlinear_mean_materialized",
        "current_RF40_power_law_aggregate_nonlinear_mean_materialized",
    )
    required_false = (
        "velocity_after_RF40_power_law_materialized",
        "full_post_XR_RF40_current_lineage_materialized",
        "outer_global_leading_velocity_materialized",
        "radial_inverse_performed_in_this_increment",
        "pressure_gradient_included",
        "restricted_forcing_included",
        "complete_ns_defect",
        "scoped_nonlinear_mean_authorized_as_correction_target",
        "mean_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    )
    for key in required_true:
        if boundary.get(key) is not True:
            raise ValueError(f"A3 #1028 truth boundary drifted at {key}")
    for key in required_false:
        if boundary.get(key) is not False:
            raise ValueError(f"A3 #1028 truth boundary drifted at {key}")
    if float(boundary.get("final_normalized_momentum_gate")) != FINAL_NORMALIZED_MOMENTUM_GATE:
        raise ValueError("final momentum gate drifted")
    if float(boundary.get("final_normalized_divergence_gate")) != FINAL_NORMALIZED_DIVERGENCE_GATE:
        raise ValueError("final divergence gate drifted")


def validate_candidate_identity(field: object) -> object:
    field_type = type(field)
    if field_type.__module__ != AGENT2_MODULE or field_type.__name__ != AGENT2_CLASS:
        raise ValueError("unexpected A2 #1010 candidate runtime identity")
    source = inspect.getsourcefile(field_type)
    if source is None or _git_blob_sha1(source) != AGENT2_SOURCE_BLOB:
        raise ValueError("A2 #1010 source blob drifted")
    semantic = getattr(field, "semantic_sha256", None)
    if not _valid_sha256(semantic):
        raise ValueError("A2 #1010 semantic identity is malformed")
    leading = getattr(field, "leading_backend", None)
    if leading is None:
        raise ValueError("A2 #1010 candidate lost its leading backend")
    ltype = type(leading)
    if ltype.__module__ != AGENT1_MODULE or ltype.__name__ != AGENT1_CLASS:
        raise ValueError("unexpected A1 #1005 leading runtime identity")
    lsource = inspect.getsourcefile(ltype)
    if lsource is None or _git_blob_sha1(lsource) != AGENT1_SOURCE_BLOB:
        raise ValueError("A1 #1005 source blob drifted")
    return leading


def _stack_piece(level: Sequence[Mapping[str, np.ndarray]], piece: str) -> np.ndarray:
    return np.stack([np.asarray(item[piece], float) for item in level], axis=0)


def _rms(a: np.ndarray) -> float:
    x = np.asarray(a, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def _relative_rms(a: np.ndarray, b: np.ndarray) -> float:
    x = np.asarray(a, float)
    y = np.asarray(b, float)
    if x.shape != y.shape:
        raise ValueError("relative RMS requires matching shapes")
    return _rms(x - y) / max(_rms(x), _rms(y), np.finfo(float).tiny)


def _comparison_metrics(
    production: np.ndarray,
    independent: np.ndarray,
    radii: np.ndarray,
) -> dict[str, Any]:
    p = np.asarray(production, float)
    q = np.asarray(independent, float)
    r = np.asarray(radii, float)
    if p.shape != q.shape or p.ndim != 2 or p.shape[1] != 3 or r.shape != (p.shape[0],):
        raise ValueError("comparison arrays/radii have incompatible shapes")
    d = p - q
    weights = r / np.sum(r)
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
    delta = agg - mixed - q
    scale = max(_rms(agg), _rms(mixed), _rms(q), np.finfo(float).tiny)
    return _rms(delta) / scale


def _ab_closure(level: Sequence[Mapping[str, np.ndarray]]) -> float:
    A = _stack_piece(level, "A")
    B = _stack_piece(level, "B")
    mixed = _stack_piece(level, "mixed")
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
        }
    return out


def _axis_metrics(
    production_level: Sequence[Mapping[str, np.ndarray]],
    independent_level: Sequence[Mapping[str, np.ndarray]],
) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for piece in SCIENTIFIC_PIECES:
        p = _stack_piece(production_level, piece)
        q = _stack_piece(independent_level, piece)
        dmax = float(np.max(np.abs(p - q)))
        scale = max(float(np.max(np.abs(p))), float(np.max(np.abs(q))), POWER_LAW_ZERO_SUPPORT_ABS_GATE)
        out[piece] = {
            "absolute_max_error": dmax,
            "normalized_error": dmax / scale,
            "production_abs_max": float(np.max(np.abs(p))),
            "independent_abs_max": float(np.max(np.abs(q))),
        }
    return out


def audit_current_rf40_power_law_nonlinear_mean(
    field: object,
    production_receipts: Mapping[str, Mapping[str, Any]],
    a3_truth_boundary: Mapping[str, Any],
) -> dict[str, Any]:
    """Run the frozen implementation-distinct audit against public #1028 receipts."""
    leading = validate_candidate_identity(field)
    _validate_a3_truth_boundary(a3_truth_boundary)

    inner_specs = frozen_inner_ring_specs()
    axis_specs = frozen_axis_near_specs()
    power_specs = frozen_power_law_specs(field)
    all_specs = inner_specs + axis_specs + power_specs
    if set(production_receipts) != {s.label for s in all_specs}:
        raise ValueError("production receipt labels do not match frozen A4 probe set")

    production: dict[str, dict[str, np.ndarray]] = {
        spec.label: _validate_receipt(production_receipts[spec.label], spec)
        for spec in all_specs
    }

    # Spatial convergence: fixed finest angular quadrature.
    spatial_levels: dict[str, list[dict[str, np.ndarray]]] = {}
    for h in FROZEN_SPATIAL_STEPS:
        spatial_levels[f"h={h:.6g}"] = [
            independent_ring_means(
                field,
                leading,
                spec,
                spatial_step=h,
                angular_order=FROZEN_ANGULAR_ORDERS[-1],
            )
            for spec in inner_specs
        ]

    # Angular convergence: fixed finest spatial FD4 step.
    angular_levels: dict[str, list[dict[str, np.ndarray]]] = {}
    for ntheta in FROZEN_ANGULAR_ORDERS:
        angular_levels[f"n={ntheta}"] = [
            independent_ring_means(
                field,
                leading,
                spec,
                spatial_step=FROZEN_SPATIAL_STEPS[-1],
                angular_order=ntheta,
            )
            for spec in inner_specs
        ]

    fine_inner = angular_levels[f"n={FROZEN_ANGULAR_ORDERS[-1]}"]
    prod_inner = [production[s.label] for s in inner_specs]
    radii = np.asarray([s.radius for s in inner_specs], dtype=float)
    comparisons = {
        piece: _comparison_metrics(
            _stack_piece(prod_inner, piece), _stack_piece(fine_inner, piece), radii
        )
        for piece in PIECES
    }
    for piece in PIECES:
        comparisons[piece]["worst_ring_label"] = inner_specs[
            comparisons[piece]["worst_ring_index"]
        ].label

    axis_independent = [
        independent_ring_means(
            field,
            leading,
            spec,
            spatial_step=FROZEN_SPATIAL_STEPS[-1],
            angular_order=FROZEN_ANGULAR_ORDERS[-1],
        )
        for spec in axis_specs
    ]
    axis_production = [production[s.label] for s in axis_specs]

    power_independent = [
        independent_ring_means(
            field,
            leading,
            spec,
            spatial_step=FROZEN_SPATIAL_STEPS[-1],
            angular_order=FROZEN_ANGULAR_ORDERS[-1],
        )
        for spec in power_specs
    ]
    power_production = [production[s.label] for s in power_specs]
    power_abs: dict[str, Any] = {}
    for piece in SCIENTIFIC_PIECES:
        p = _stack_piece(power_production, piece)
        q = _stack_piece(power_independent, piece)
        power_abs[piece] = {
            "production_abs_max": float(np.max(np.abs(p))),
            "independent_abs_max": float(np.max(np.abs(q))),
            "difference_abs_max": float(np.max(np.abs(p - q))),
        }

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
            "power_stage_fractions": list(FROZEN_POWER_STAGE_FRACTIONS),
        },
        "inner_specs": [spec.__dict__ for spec in inner_specs],
        "axis_near_specs": [spec.__dict__ for spec in axis_specs],
        "power_law_specs": [spec.__dict__ for spec in power_specs],
        "fine_production_vs_independent": comparisons,
        "spatial_resolution": _resolution_report(spatial_levels),
        "angular_resolution": _resolution_report(angular_levels),
        "independent_piece_closure_relative_rms": _piece_closure(fine_inner),
        "independent_A_plus_B_closure_relative_rms": _ab_closure(fine_inner),
        "axis_near": _axis_metrics(axis_production, axis_independent),
        "power_law_zero_support": power_abs,
        "inner_production_quadratic_rms": _rms(_stack_piece(prod_inner, "quadratic")),
        "inner_independent_quadratic_rms": _rms(_stack_piece(fine_inner, "quadratic")),
        "inner_production_aggregate_rms": _rms(_stack_piece(prod_inner, "aggregate")),
        "inner_independent_aggregate_rms": _rms(_stack_piece(fine_inner, "aggregate")),
        "scientific_boundary": {
            "complete_ns_defect": False,
            "matched_pressure_materialized": False,
            "restricted_forcing_materialized": False,
            "correction_velocity_materialized": False,
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
    cm = float(metrics["coarse_to_medium_relative_rms"])
    mf = float(metrics["medium_to_fine_relative_rms"])
    if mf > FINE_PAIR_REL_RMS_GATE:
        return False
    if cm > RESOLUTION_NUMERICAL_FLOOR and mf > RESOLUTION_NUMERICAL_FLOOR:
        if mf > RESOLUTION_DEGRADATION_FACTOR * cm:
            return False
    return True


def scoped_gates_pass(report: Mapping[str, Any]) -> bool:
    try:
        comparisons = report["fine_production_vs_independent"]
        for piece in SCIENTIFIC_PIECES:
            m = comparisons[piece]
            if float(m["relative_rms"]) > FINE_REL_RMS_GATE:
                return False
            if float(m["relative_weighted_l2"]) > FINE_REL_WEIGHTED_L2_GATE:
                return False
            if float(m["relative_max"]) > FINE_REL_MAX_GATE:
                return False
        for family in ("spatial_resolution", "angular_resolution"):
            for piece in SCIENTIFIC_PIECES:
                if not _resolution_piece_pass(report[family]["pieces"][piece]):
                    return False
        if float(report["independent_piece_closure_relative_rms"]) > PIECE_CLOSURE_RELATIVE_GATE:
            return False
        if float(report["independent_A_plus_B_closure_relative_rms"]) > PIECE_CLOSURE_RELATIVE_GATE:
            return False
        for piece in SCIENTIFIC_PIECES:
            if float(report["axis_near"][piece]["normalized_error"]) > AXIS_NEAR_NORMALIZED_ERROR_GATE:
                return False
            power = report["power_law_zero_support"][piece]
            if max(
                float(power["production_abs_max"]),
                float(power["independent_abs_max"]),
                float(power["difference_abs_max"]),
            ) > POWER_LAW_ZERO_SUPPORT_ABS_GATE:
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
            "after_correction_residual_assessed",
            "heldout_normalized_ns_residual_assessed",
            "same_protocol_comparable_to_st006",
            "pde_validated",
        ):
            if boundary[key] is not False:
                return False
        if float(boundary["final_normalized_momentum_gate"]) != FINAL_NORMALIZED_MOMENTUM_GATE:
            return False
        if float(boundary["final_normalized_divergence_gate"]) != FINAL_NORMALIZED_DIVERGENCE_GATE:
            return False
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    return True


def enforce_scoped_gates(report: Mapping[str, Any]) -> None:
    if not scoped_gates_pass(report):
        raise AssertionError("preregistered K4-VAL-105 scoped nonlinear-mean gates failed")


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(audit_current_rf40_power_law_nonlinear_mean)
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "implementation_distinct_from_agent3_derivative_path": True,
        "implementation_distinct_from_agent3_angular_projector": True,
        "candidate_loaded_by_workflow_before_scientific_audit": True,
        "caller_tunable_threshold": False,
        "public_inputs": list(signature.parameters),
        "complete_ns_defect": False,
        "matched_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "correction_velocity_materialized": False,
        "heldout_normalized_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
