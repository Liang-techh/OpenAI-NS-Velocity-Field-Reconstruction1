"""Independent Agent-4 audit of the candidate-side Kokuno final bridge through Xi=110.

Agent 1 PR #940 materializes the current actual source-coordinate final bridge on
``100 <= X <= 110`` and exposes save/reload public ``values(X, eta)`` together
with production radial derivatives.  This module treats those production
derivatives strictly as quantities under audit.

The numerical reference is reconstructed only from the save/reloaded public
scalar value surface.  Interior derivatives use a seven-node asymmetric
degree-six local-polynomial stencil whose weights are solved from polynomial
moment equations.  Endpoint derivatives at X=100 and Xi=110 use separate
seven-node one-sided degree-six stencils, so the exact endpoint derivative
claims are not production-to-production replays.

This remains a source-coordinate prerequisite audit.  It does not construct a
global Cartesian velocity, matched pressure, restricted forcing, correction
velocity, or complete Navier--Stokes residual.  Repository-wide momentum
``1e-3`` and divergence ``1e-5`` gates remain frozen and unevaluated here.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import tempfile
from typing import Any, Callable, Mapping

import numpy as np


TASK = "K4-VAL-088"
SCHEMA = "kokuno-a4-actual-final-bridge-xi110-independent-audit-v1"
A1_FINAL_BRIDGE_PR = 940
A1_FINAL_BRIDGE_HEAD = "f430f8f97df3bfafce56bf094047c769ecec922d"
A1_FINAL_BRIDGE_SOURCE_BLOB = "0cbd49ad0726c1573f2e649c96a3b516d8df053e"

FROZEN_FINAL_GATES = {
    "normalized_momentum_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
}

INTERIOR_NODES = np.asarray(
    [-1.0, -0.73, -0.31, 0.0, 0.28, 0.67, 1.0], dtype=float
)
FORWARD_ENDPOINT_NODES = np.arange(7, dtype=float)
BACKWARD_ENDPOINT_NODES = -np.arange(7, dtype=float)


@dataclass(frozen=True)
class FrozenActualBridgeAuditProtocol:
    seed: int = 9173621
    random_per_region: int = 40
    interior_halfwidths: tuple[float, float, float] = (0.12, 0.06, 0.03)
    endpoint_steps: tuple[float, float, float] = (0.08, 0.04, 0.02)
    derivative_relative_rms_gate: float = 5.0e-3
    derivative_relative_max_gate: float = 2.0e-2
    endpoint_scale_normalized_rms_gate: float = 5.0e-3
    endpoint_scale_normalized_max_gate: float = 2.0e-2
    refinement_ratio_gate: float = 8.0
    refinement_floor: float = 2.0e-9
    value_identity_relative_gate: float = 5.0e-12
    derivative_chain_relative_gate: float = 5.0e-5
    endpoint_value_handoff_abs_gate: float = 2.0e-12
    xi_final_slope_abs_gate: float = 2.0e-4
    xi_zero_Ux_scale_normalized_gate: float = 5.0e-3
    value_rms_floor: float = 1.0e-10
    derivative_rms_floor: float = 1.0e-12


PROTOCOL = FrozenActualBridgeAuditProtocol()

TRUTH_BOUNDARY = {
    "actual_final_interpolation_100_to_Xi_materialized": True,
    "actual_G_i_at_Xi_materialized": True,
    "actual_ell_i_at_Xi_materialized": True,
    "actual_final_bridge_independently_audited": False,
    "x100_endpoint_radial_derivatives_independently_audited": False,
    "Xi_endpoint_radial_derivatives_independently_audited": False,
    "outer_global_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "correction_velocity_materialized": False,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "heldout_complete_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def public_contract() -> dict[str, Any]:
    return {
        "task": TASK,
        "schema": SCHEMA,
        "agent1_pr": A1_FINAL_BRIDGE_PR,
        "agent1_exact_head": A1_FINAL_BRIDGE_HEAD,
        "agent1_source_blob": A1_FINAL_BRIDGE_SOURCE_BLOB,
        "scientific_reference_surface": (
            "save/reloaded public values(X,eta) -> F_final_bridge/U_final_bridge/E_final_bridge"
        ),
        "production_surface_under_audit": (
            "radial_derivatives(X,eta) -> F_final_bridge_X/U_final_bridge_X/E_final_bridge_X"
        ),
        "independent_operator": (
            "degree-six local polynomial derivatives reconstructed from public values only; "
            "asymmetric seven-node interior stencil plus distinct one-sided endpoint stencils"
        ),
        "interior_nodes": INTERIOR_NODES.tolist(),
        "forward_endpoint_nodes": FORWARD_ENDPOINT_NODES.tolist(),
        "backward_endpoint_nodes": BACKWARD_ENDPOINT_NODES.tolist(),
        "protocol": asdict(PROTOCOL),
        "frozen_final_gates": dict(FROZEN_FINAL_GATES),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} returned non-finite values")
    return out


def _rms(value: Any) -> float:
    arr = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(arr * arr)))


def _relative_closure(lhs: Any, rhs: Any) -> float:
    left = np.asarray(lhs, dtype=float)
    right = np.asarray(rhs, dtype=float)
    return _rms(left - right) / max(_rms(left), _rms(right), 1.0e-30)


def first_derivative_weights(nodes: np.ndarray) -> np.ndarray:
    """Return d/ds weights at s=0 for the local interpolating polynomial."""
    s = np.asarray(nodes, dtype=float)
    if s.ndim != 1 or s.size < 3 or np.unique(s).size != s.size:
        raise ValueError("derivative nodes must be distinct one-dimensional values")
    V = np.vander(s, N=s.size, increasing=True)
    rhs = np.zeros(s.size, dtype=float)
    rhs[1] = 1.0
    weights = np.linalg.solve(V.T, rhs)
    if not np.all(np.isfinite(weights)):
        raise RuntimeError("failed to construct finite derivative weights")
    return weights


INTERIOR_WEIGHTS = first_derivative_weights(INTERIOR_NODES)
FORWARD_ENDPOINT_WEIGHTS = first_derivative_weights(FORWARD_ENDPOINT_NODES)
BACKWARD_ENDPOINT_WEIGHTS = first_derivative_weights(BACKWARD_ENDPOINT_NODES)


def _require_contract(field: Any) -> None:
    required = (
        "values",
        "radial_derivatives",
        "x100_handoff",
        "handoff_at_Xi",
        "save_configuration",
        "semantic_sha256",
        "X_b_ref",
        "X_axial_end",
        "X_angular_end",
        "X_i",
        "eta_interval",
        "C",
    )
    for name in required:
        if not hasattr(field, name):
            raise TypeError(f"actual final bridge missing required public surface {name!r}")
    for name in (
        "values",
        "radial_derivatives",
        "x100_handoff",
        "handoff_at_Xi",
        "save_configuration",
    ):
        if not callable(getattr(field, name)):
            raise TypeError(f"{name} must be callable")
    if abs(float(field.X_b_ref) - 100.0) > 0.0 or abs(float(field.X_i) - 110.0) > 0.0:
        raise ValueError("A4 actual-bridge audit is frozen to X=100 -> Xi=110")


def _public_scalar(field: Any, key: str) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    def evaluate(x: np.ndarray, eta: np.ndarray) -> np.ndarray:
        payload = field.values(x, eta)
        if key not in payload:
            raise KeyError(f"public values payload missing {key!r}")
        return _finite(payload[key], key)

    return evaluate


def polynomial_first(
    fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    x: np.ndarray,
    eta: np.ndarray,
    step: float,
    nodes: np.ndarray,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Differentiate one public scalar on one frozen polynomial micro-stencil."""
    if not (math.isfinite(step) and step > 0.0):
        raise ValueError("step must be finite and positive")
    xx = np.asarray(x, dtype=float)
    ee = np.asarray(eta, dtype=float)
    ww = first_derivative_weights(nodes) if weights is None else np.asarray(weights, dtype=float)
    out = np.zeros_like(xx, dtype=float)
    for node, weight in zip(np.asarray(nodes), ww, strict=True):
        out += weight * _finite(fn(xx + node * step, ee), "polynomial sample")
    return out / step


def _all_independent_derivatives(
    field: Any,
    x: np.ndarray,
    eta: np.ndarray,
    step: float,
    *,
    nodes: np.ndarray,
    weights: np.ndarray,
) -> dict[str, np.ndarray]:
    keys = ("F_final_bridge", "U_final_bridge", "E_final_bridge")
    accum = {key: np.zeros_like(np.asarray(x, dtype=float)) for key in keys}
    xx = np.asarray(x, dtype=float)
    ee = np.asarray(eta, dtype=float)
    for node, weight in zip(nodes, weights, strict=True):
        payload = field.values(xx + node * step, ee)
        for key in keys:
            accum[key] += weight * _finite(payload[key], f"polynomial {key}")
    return {
        "F_final_bridge_X": accum["F_final_bridge"] / step,
        "U_final_bridge_X": accum["U_final_bridge"] / step,
        "E_final_bridge_X": accum["E_final_bridge"] / step,
    }


def _match_metrics(production: Any, independent: Any) -> dict[str, float]:
    prod = np.asarray(production, dtype=float)
    indep = np.asarray(independent, dtype=float)
    err = indep - prod
    return {
        "absolute_rms": _rms(err),
        "absolute_max": float(np.max(np.abs(err))),
        "relative_rms": _rms(err) / max(_rms(prod), 1.0e-30),
        "relative_sampled_max": float(np.max(np.abs(err)))
        / max(float(np.max(np.abs(prod))), 1.0e-30),
    }


def _endpoint_metrics(
    production: Any, independent: Any, values: Any, Xscale: float
) -> dict[str, float]:
    prod = np.asarray(production, dtype=float)
    indep = np.asarray(independent, dtype=float)
    val = np.asarray(values, dtype=float)
    err = indep - prod
    rms_scale = max(_rms(prod), _rms(val) / Xscale, 1.0e-30)
    max_scale = max(
        float(np.max(np.abs(prod))),
        float(np.max(np.abs(val))) / Xscale,
        1.0e-30,
    )
    return {
        "absolute_rms": _rms(err),
        "absolute_max": float(np.max(np.abs(err))),
        "scale_normalized_rms": _rms(err) / rms_scale,
        "scale_normalized_max": float(np.max(np.abs(err))) / max_scale,
        "normalization_rms_scale": rms_scale,
        "normalization_max_scale": max_scale,
    }


def _stability(levels: list[np.ndarray]) -> dict[str, float | bool]:
    if len(levels) != 3:
        raise ValueError("exactly three independent resolution levels are required")
    coarse_medium = _rms(np.asarray(levels[0]) - np.asarray(levels[1]))
    medium_fine = _rms(np.asarray(levels[1]) - np.asarray(levels[2]))
    ratio = 1.0e300 if medium_fine == 0.0 else coarse_medium / medium_fine
    floor_hit = medium_fine <= PROTOCOL.refinement_floor
    return {
        "coarse_to_medium_rms": coarse_medium,
        "medium_to_fine_rms": medium_fine,
        "ratio": ratio,
        "numerical_floor_hit": floor_hit,
        "passed": bool(floor_hit or ratio >= PROTOCOL.refinement_ratio_gate),
    }


def _worst(
    x: np.ndarray,
    eta: np.ndarray,
    production: np.ndarray,
    independent: np.ndarray,
) -> dict[str, float]:
    err = np.abs(np.asarray(independent) - np.asarray(production))
    index = int(np.argmax(err))
    return {
        "X": float(np.asarray(x)[index]),
        "eta": float(np.asarray(eta)[index]),
        "production": float(np.asarray(production)[index]),
        "independent": float(np.asarray(independent)[index]),
        "absolute_error": float(err[index]),
    }


def _sample_points(field: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    eta_lo, eta_hi = (float(v) for v in field.eta_interval)
    margin_eta = max(0.08 * (eta_hi - eta_lo), 2.0e-3)
    if eta_lo + margin_eta >= eta_hi - margin_eta:
        raise ValueError("public eta interval is too narrow for frozen A4 probes")
    rng = np.random.default_rng(PROTOCOL.seed)
    hmax = max(PROTOCOL.interior_halfwidths)
    margin_x = hmax + 0.08

    intervals = (
        (float(field.X_b_ref) + margin_x, float(field.X_axial_end) - margin_x),
        (float(field.X_axial_end) + margin_x, float(field.X_angular_end) - margin_x),
        (float(field.X_angular_end) + margin_x, float(field.X_i) - margin_x),
    )
    xs: list[np.ndarray] = []
    es: list[np.ndarray] = []
    for lo, hi in intervals:
        if not lo < hi:
            raise ValueError("bridge region too narrow for frozen interior audit")
        xs.append(rng.uniform(lo, hi, PROTOCOL.random_per_region))
        es.append(
            rng.uniform(
                eta_lo + margin_eta,
                eta_hi - margin_eta,
                PROTOCOL.random_per_region,
            )
        )

    center_eta = np.asarray([0.0, -1.0e-8, 1.0e-8, 0.0, -1.0e-6, 1.0e-6])
    center_x = np.asarray(
        [
            100.7,
            101.5,
            0.5 * (field.X_axial_end + field.X_angular_end),
            field.X_angular_end + 0.6,
            106.5,
            109.0,
        ],
        dtype=float,
    )
    if np.any((center_eta < eta_lo) | (center_eta > eta_hi)):
        raise ValueError("public eta interval does not include frozen center probes")

    transition_x = np.asarray(
        [
            field.X_axial_end - 0.15,
            field.X_axial_end + 0.15,
            field.X_angular_end - 0.15,
            field.X_angular_end + 0.15,
        ],
        dtype=float,
    )
    transition_eta = np.asarray([-0.45, 0.37, -0.31, 0.49], dtype=float)
    transition_eta = np.clip(
        transition_eta, eta_lo + margin_eta, eta_hi - margin_eta
    )

    x = np.concatenate([*xs, center_x, transition_x])
    eta = np.concatenate([*es, center_eta, transition_eta])

    if np.min(x) - hmax < float(field.X_b_ref) or np.max(x) + hmax > float(field.X_i):
        raise ValueError("frozen interior probes do not leave room for coarse stencil")
    endpoint_eta = np.linspace(
        max(eta_lo + margin_eta, -0.6),
        min(eta_hi - margin_eta, 0.6),
        7,
    )
    return x, eta, endpoint_eta


def _mutation_checks(
    field: Any,
    x: np.ndarray,
    eta: np.ndarray,
    production: Mapping[str, np.ndarray],
    interior_fine: Mapping[str, np.ndarray],
    xi_eta: np.ndarray,
    xi_production: Mapping[str, np.ndarray],
    xi_values: Mapping[str, np.ndarray],
) -> dict[str, bool]:
    F_scaled = _match_metrics(
        0.99 * np.asarray(production["F_final_bridge_X"]),
        interior_fine["F_final_bridge_X"],
    )
    E_scaled = _match_metrics(
        0.90 * np.asarray(production["E_final_bridge_X"]),
        interior_fine["E_final_bridge_X"],
    )
    production_F_detected = not (
        F_scaled["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
        and F_scaled["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
    )
    production_E_detected = not (
        E_scaled["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
        and E_scaled["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
    )

    F = _public_scalar(field, "F_final_bridge")
    U = _public_scalar(field, "U_final_bridge")

    def F_plus_linear(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return F(xx, ee) + 1.0e-3 * np.asarray(xx, dtype=float)

    def U_plus_linear(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return U(xx, ee) + 1.0e-3 * np.asarray(xx, dtype=float)

    fine_h = PROTOCOL.interior_halfwidths[-1]
    F_mut = polynomial_first(
        F_plus_linear, x, eta, fine_h, INTERIOR_NODES, INTERIOR_WEIGHTS
    )
    U_mut = polynomial_first(
        U_plus_linear, x, eta, fine_h, INTERIOR_NODES, INTERIOR_WEIGHTS
    )
    public_F_detected = (
        _rms(F_mut - interior_fine["F_final_bridge_X"]) >= 5.0e-4
    )
    public_U_detected = (
        _rms(U_mut - interior_fine["U_final_bridge_X"]) >= 5.0e-4
    )

    Xi = float(field.X_i)

    def F_xi_tangent_mutation(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return F(xx, ee) + 1.0e-3 * (np.asarray(xx, dtype=float) - Xi)

    xi_x = np.full_like(xi_eta, Xi, dtype=float)
    xi_mut = polynomial_first(
        F_xi_tangent_mutation,
        xi_x,
        xi_eta,
        PROTOCOL.endpoint_steps[-1],
        BACKWARD_ENDPOINT_NODES,
        BACKWARD_ENDPOINT_WEIGHTS,
    )
    xi_mut_metrics = _endpoint_metrics(
        xi_production["F_final_bridge_X"],
        xi_mut,
        xi_values["F_final_bridge"],
        Xi,
    )
    xi_tangent_detected = not (
        xi_mut_metrics["scale_normalized_rms"]
        <= PROTOCOL.endpoint_scale_normalized_rms_gate
        and xi_mut_metrics["scale_normalized_max"]
        <= PROTOCOL.endpoint_scale_normalized_max_gate
    )

    return {
        "production_F_X_times_0p99_detected": bool(production_F_detected),
        "production_E_X_times_0p90_detected": bool(production_E_detected),
        "public_F_plus_1e_minus3_X_detected": bool(public_F_detected),
        "public_U_plus_1e_minus3_X_detected": bool(public_U_detected),
        "Xi_value_preserving_F_tangent_plus_1e_minus3_detected": bool(
            xi_tangent_detected
        ),
    }


def materialize_actual_final_bridge_audit(field: Any) -> dict[str, Any]:
    """Run the frozen implementation-distinct audit on one rebound bridge artifact."""
    _require_contract(field)
    x, eta, endpoint_eta = _sample_points(field)
    values = field.values(x, eta)
    production = field.radial_derivatives(x, eta)

    audited_keys = (
        "F_final_bridge_X",
        "U_final_bridge_X",
        "E_final_bridge_X",
    )
    for key in audited_keys:
        if key not in production:
            raise KeyError(f"production radial derivative payload missing {key!r}")

    interior_levels = [
        _all_independent_derivatives(
            field,
            x,
            eta,
            halfwidth,
            nodes=INTERIOR_NODES,
            weights=INTERIOR_WEIGHTS,
        )
        for halfwidth in PROTOCOL.interior_halfwidths
    ]
    interior_fine = interior_levels[-1]

    interior_metrics: dict[str, list[dict[str, float]]] = {}
    interior_stability: dict[str, dict[str, float | bool]] = {}
    interior_worst: dict[str, dict[str, float]] = {}
    interior_gates: dict[str, bool] = {}
    for key in audited_keys:
        prod = _finite(production[key], key)
        levels = [_finite(level[key], f"independent {key}") for level in interior_levels]
        interior_metrics[key] = [_match_metrics(prod, level) for level in levels]
        interior_stability[key] = _stability(levels)
        interior_worst[key] = _worst(x, eta, prod, levels[-1])
        fine = interior_metrics[key][-1]
        interior_gates[key] = bool(
            fine["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
            and fine["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
            and interior_stability[key]["passed"]
        )

    endpoint_results: dict[str, Any] = {}
    endpoint_payloads: dict[str, tuple[dict[str, np.ndarray], dict[str, np.ndarray]]] = {}
    for label, X0, nodes, weights in (
        (
            "X100",
            float(field.X_b_ref),
            FORWARD_ENDPOINT_NODES,
            FORWARD_ENDPOINT_WEIGHTS,
        ),
        (
            "Xi",
            float(field.X_i),
            BACKWARD_ENDPOINT_NODES,
            BACKWARD_ENDPOINT_WEIGHTS,
        ),
    ):
        xx = np.full_like(endpoint_eta, X0, dtype=float)
        vals = field.values(xx, endpoint_eta)
        prod = field.radial_derivatives(xx, endpoint_eta)
        levels = [
            _all_independent_derivatives(
                field,
                xx,
                endpoint_eta,
                step,
                nodes=nodes,
                weights=weights,
            )
            for step in PROTOCOL.endpoint_steps
        ]
        endpoint_payloads[label] = (vals, prod)
        metrics: dict[str, list[dict[str, float]]] = {}
        stability: dict[str, dict[str, float | bool]] = {}
        worst: dict[str, dict[str, float]] = {}
        gates: dict[str, bool] = {}
        value_key_for_derivative = {
            "F_final_bridge_X": "F_final_bridge",
            "U_final_bridge_X": "U_final_bridge",
            "E_final_bridge_X": "E_final_bridge",
        }
        for key in audited_keys:
            value_key = value_key_for_derivative[key]
            indep_levels = [level[key] for level in levels]
            metrics[key] = [
                _endpoint_metrics(prod[key], level, vals[value_key], X0)
                for level in indep_levels
            ]
            stability[key] = _stability(indep_levels)
            worst[key] = _worst(xx, endpoint_eta, prod[key], indep_levels[-1])
            fine = metrics[key][-1]
            gates[key] = bool(
                fine["scale_normalized_rms"]
                <= PROTOCOL.endpoint_scale_normalized_rms_gate
                and fine["scale_normalized_max"]
                <= PROTOCOL.endpoint_scale_normalized_max_gate
                and stability[key]["passed"]
            )
        endpoint_results[label] = {
            "metrics": metrics,
            "stability": stability,
            "worst_witness": worst,
            "gates": gates,
        }

    # Public value identity and independent derivative-chain identity.
    F = _finite(values["F_final_bridge"], "F_final_bridge")
    U = _finite(values["U_final_bridge"], "U_final_bridge")
    E = _finite(values["E_final_bridge"], "E_final_bridge")
    value_identity = _relative_closure(E, np.sqrt(2.0 * x) * F)
    derivative_chain_rhs = (
        F / np.sqrt(2.0 * x)
        + np.sqrt(2.0 * x) * interior_fine["F_final_bridge_X"]
    )
    derivative_chain = _relative_closure(
        interior_fine["E_final_bridge_X"], derivative_chain_rhs
    )

    # Exact public value handoffs, separated from independent endpoint derivatives.
    X100 = np.full_like(endpoint_eta, float(field.X_b_ref), dtype=float)
    Xi = np.full_like(endpoint_eta, float(field.X_i), dtype=float)
    x100_values = field.values(X100, endpoint_eta)
    x100_handoff = field.x100_handoff(endpoint_eta)
    xi_values = field.values(Xi, endpoint_eta)
    xi_handoff = field.handoff_at_Xi(endpoint_eta)

    x100_value_map = {
        "F_final_bridge": "F_fixed_kappa",
        "U_final_bridge": "U_fixed_kappa",
        "E_final_bridge": "E_fixed_kappa",
    }
    x100_handoff_abs_max = max(
        float(
            np.max(
                np.abs(
                    _finite(x100_values[child], child)
                    - _finite(x100_handoff[parent], parent)
                )
            )
        )
        for child, parent in x100_value_map.items()
    )
    xi_handoff_abs_max = max(
        float(
            np.max(
                np.abs(
                    _finite(xi_handoff[key], f"handoff {key}")
                    - _finite(xi_values[key], f"direct {key}")
                )
            )
        )
        for key in ("F_final_bridge", "U_final_bridge", "E_final_bridge")
    )
    xi_handoff_abs_max = max(
        xi_handoff_abs_max,
        float(np.max(np.abs(_finite(xi_handoff["G_i"], "G_i") - xi_values["U_final_bridge"]))),
        float(
            np.max(
                np.abs(
                    _finite(xi_handoff["ell_i"], "ell_i")
                    - np.log(float(field.C) * xi_values["E_final_bridge"])
                )
            )
        ),
    )

    # Reconstruct finest endpoint derivatives again for structural endpoint checks.
    xi_independent_derivatives = _all_independent_derivatives(
        field,
        Xi,
        endpoint_eta,
        PROTOCOL.endpoint_steps[-1],
        nodes=BACKWARD_ENDPOINT_NODES,
        weights=BACKWARD_ENDPOINT_WEIGHTS,
    )
    xi_slope = (
        float(field.X_i)
        * xi_independent_derivatives["F_final_bridge_X"]
        / _finite(xi_values["F_final_bridge"], "Xi F")
    )
    xi_slope_abs_max = float(np.max(np.abs(xi_slope + 0.4)))
    xi_Ux_scale = max(
        _rms(xi_values["U_final_bridge"]) / float(field.X_i), 1.0e-30
    )
    xi_Ux_scale_normalized = (
        _rms(xi_independent_derivatives["U_final_bridge_X"]) / xi_Ux_scale
    )

    value_rms = {
        "F_final_bridge": _rms(F),
        "U_final_bridge": _rms(U),
        "E_final_bridge": _rms(E),
    }
    derivative_rms = {key: _rms(production[key]) for key in audited_keys}
    nontrivial = bool(
        all(value >= PROTOCOL.value_rms_floor for value in value_rms.values())
        and derivative_rms["F_final_bridge_X"] >= PROTOCOL.derivative_rms_floor
        and derivative_rms["E_final_bridge_X"] >= PROTOCOL.derivative_rms_floor
    )

    mutations = _mutation_checks(
        field,
        x,
        eta,
        production,
        interior_fine,
        endpoint_eta,
        endpoint_payloads["Xi"][1],
        endpoint_payloads["Xi"][0],
    )

    gates = {
        "interior_derivative_surfaces": bool(all(interior_gates.values())),
        "X100_endpoint_derivative_surfaces": bool(
            all(endpoint_results["X100"]["gates"].values())
        ),
        "Xi_endpoint_derivative_surfaces": bool(
            all(endpoint_results["Xi"]["gates"].values())
        ),
        "E_equals_sqrt2X_F": bool(
            value_identity <= PROTOCOL.value_identity_relative_gate
        ),
        "independent_E_derivative_chain": bool(
            derivative_chain <= PROTOCOL.derivative_chain_relative_gate
        ),
        "X100_value_handoff": bool(
            x100_handoff_abs_max <= PROTOCOL.endpoint_value_handoff_abs_gate
        ),
        "Xi_value_handoff": bool(
            xi_handoff_abs_max <= PROTOCOL.endpoint_value_handoff_abs_gate
        ),
        "Xi_independent_final_F_slope": bool(
            xi_slope_abs_max <= PROTOCOL.xi_final_slope_abs_gate
        ),
        "Xi_independent_Ux_zero": bool(
            xi_Ux_scale_normalized <= PROTOCOL.xi_zero_Ux_scale_normalized_gate
        ),
        "nontrivial": nontrivial,
        "all_mutations_detected": bool(all(mutations.values())),
    }
    passed = bool(all(gates.values()))

    truth = dict(TRUTH_BOUNDARY)
    truth["actual_final_bridge_independently_audited"] = passed
    truth["x100_endpoint_radial_derivatives_independently_audited"] = bool(
        gates["X100_endpoint_derivative_surfaces"]
    )
    truth["Xi_endpoint_radial_derivatives_independently_audited"] = bool(
        gates["Xi_endpoint_derivative_surfaces"]
        and gates["Xi_independent_final_F_slope"]
        and gates["Xi_independent_Ux_zero"]
    )

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "artifact_semantic_sha256": str(field.semantic_sha256),
        "protocol": asdict(PROTOCOL),
        "operators": {
            "interior": {
                "name": "asymmetric_degree6_local_polynomial_first_derivative",
                "nodes": INTERIOR_NODES.tolist(),
                "weights": INTERIOR_WEIGHTS.tolist(),
            },
            "X100": {
                "name": "forward_degree6_endpoint_polynomial_first_derivative",
                "nodes": FORWARD_ENDPOINT_NODES.tolist(),
                "weights": FORWARD_ENDPOINT_WEIGHTS.tolist(),
            },
            "Xi": {
                "name": "backward_degree6_endpoint_polynomial_first_derivative",
                "nodes": BACKWARD_ENDPOINT_NODES.tolist(),
                "weights": BACKWARD_ENDPOINT_WEIGHTS.tolist(),
            },
        },
        "sample_count": int(np.asarray(x).size),
        "endpoint_sample_count": int(endpoint_eta.size),
        "interior_resolution_metrics": interior_metrics,
        "interior_resolution_stability": interior_stability,
        "interior_worst_witness": interior_worst,
        "endpoint_results": endpoint_results,
        "value_identity_relative": value_identity,
        "independent_derivative_chain_relative": derivative_chain,
        "X100_value_handoff_abs_max": x100_handoff_abs_max,
        "Xi_value_handoff_abs_max": xi_handoff_abs_max,
        "Xi_independent_final_slope_abs_max": xi_slope_abs_max,
        "Xi_independent_Ux_scale_normalized": xi_Ux_scale_normalized,
        "value_rms": value_rms,
        "derivative_rms": derivative_rms,
        "mutation_checks": mutations,
        "gates": gates,
        "passed": passed,
        "frozen_final_gates": dict(FROZEN_FINAL_GATES),
        "truth_boundary": truth,
    }
    receipt["receipt_sha256"] = _sha256(receipt)
    return receipt


def run_real_artifact_audit() -> dict[str, Any]:
    """Save/reload Agent-1 #940 and audit only rebound public surfaces."""
    from .kokuno_pa10_actual_final_bridge_xi110 import (
        KokunoPA10ActualFinalBridgeToXi110,
    )

    original = KokunoPA10ActualFinalBridgeToXi110()
    with tempfile.TemporaryDirectory(prefix="kokuno-a4-actual-xi110-") as tmp:
        path = Path(tmp) / "actual-final-bridge-xi110.json"
        original.save_configuration(path)
        rebound = KokunoPA10ActualFinalBridgeToXi110.load_configuration(path)
        if rebound.semantic_sha256 != original.semantic_sha256:
            raise AssertionError("save/reload changed actual final-bridge semantic identity")
        receipt = materialize_actual_final_bridge_audit(rebound)
    receipt["save_reload_exact_replay"] = True
    receipt.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = _sha256(receipt)
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    receipt = run_real_artifact_audit()
    rendered = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)
    return 0 if receipt["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
