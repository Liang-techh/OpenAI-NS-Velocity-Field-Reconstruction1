"""Independent Agent-4 audit of the autonomous reference-stress Xi extension.

Agent 1 PR #933 exposes a serializable public reference-stress artifact on
``100 <= X <= 110``.  The production artifact integrates ``p1`` and ``N_s``
with fixed Gauss--Legendre quadrature and exposes radial derivatives from the
public source ODEs.  This module treats those derivatives as quantities under
audit.

The numerical reference is reconstructed *only* from save/reloaded public
``values(X, eta)`` evaluations.  It uses a nonuniform degree-six local
polynomial interpolation derivative on seven Chebyshev--Lobatto nodes at three
frozen physical half-widths.  That path is implementation-distinct from the
Agent-1 source-ODE derivative and radial quadrature, and from the uniform FD6/
FD8 operators used by earlier Agent-4 audits.

This remains a source-coordinate prerequisite audit.  It does not execute the
actual candidate bridge on ``100 < X < 110`` and does not create a global
Cartesian velocity, matched pressure, restricted forcing, correction velocity,
or complete Navier--Stokes residual.  The repository-wide momentum ``1e-3``
and divergence ``1e-5`` gates are unchanged and unevaluated here.
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


TASK = "K4-VAL-087"
SCHEMA = "kokuno-a4-reference-stress-xi110-independent-audit-v1"
A1_REFERENCE_STRESS_PR = 933
A1_REFERENCE_STRESS_HEAD = "15791891fbbcfdaa614930791a843983afde2b4d"
A1_REFERENCE_STRESS_SOURCE_BLOB = "d534da9efed25df9d281bd0c78fc18ec9b706c09"

FROZEN_FINAL_GATES = {
    "normalized_momentum_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
}


@dataclass(frozen=True)
class FrozenReferenceStressXiAuditProtocol:
    seed: int = 9173611
    random_offgrid_count: int = 128
    x_min: float = 100.35
    x_max: float = 109.65
    x_halfwidths: tuple[float, float, float] = (0.24, 0.12, 0.06)
    derivative_relative_rms_gate: float = 5.0e-3
    derivative_relative_max_gate: float = 2.0e-2
    refinement_ratio_gate: float = 20.0
    refinement_floor: float = 2.0e-9
    algebraic_closure_relative_gate: float = 5.0e-12
    derivative_closure_relative_gate: float = 5.0e-8
    plateau_abs_gate: float = 2.0e-12
    xi_handoff_abs_gate: float = 2.0e-12
    value_rms_floor: float = 1.0e-10
    derivative_rms_floor: float = 1.0e-12


PROTOCOL = FrozenReferenceStressXiAuditProtocol()

# Seven nonuniform Chebyshev-Lobatto nodes.  We intentionally compute the
# interpolation derivative weights from the polynomial moment equations rather
# than hard-coding a standard uniform FD stencil.
LOCAL_POLY_NODES = np.asarray(
    [-1.0, -math.sqrt(3.0) / 2.0, -0.5, 0.0, 0.5, math.sqrt(3.0) / 2.0, 1.0],
    dtype=float,
)

TRUTH_BOUNDARY = {
    "autonomous_reference_stress_to_Xi_materialized": True,
    "autonomous_reference_stress_to_Xi_independently_audited": False,
    "actual_final_interpolation_100_to_Xi_materialized": False,
    "actual_G_i_at_Xi_materialized": False,
    "actual_ell_i_at_Xi_materialized": False,
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
        "agent1_pr": A1_REFERENCE_STRESS_PR,
        "agent1_exact_head": A1_REFERENCE_STRESS_HEAD,
        "agent1_source_blob": A1_REFERENCE_STRESS_SOURCE_BLOB,
        "scientific_reference_surface": (
            "save/reloaded public values(X,eta) -> "
            "F/U/E/Pi/p1/N_s/n_s/p2"
        ),
        "production_surface_under_audit": (
            "radial_derivatives(X,eta) -> "
            "E_X/Pi_X/p1_X/N_s_X/n_s_X"
        ),
        "independent_operator": (
            "seven-node nonuniform degree-six local-polynomial derivative "
            "reconstructed only from public save/reloaded scalar values"
        ),
        "local_polynomial_nodes": LOCAL_POLY_NODES.tolist(),
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


def _require_contract(field: Any) -> None:
    required = (
        "values",
        "radial_derivatives",
        "handoff_at_Xi",
        "save_configuration",
        "semantic_sha256",
        "X_b_ref",
        "X_i",
        "eta_interval",
        "h",
    )
    for name in required:
        if not hasattr(field, name):
            raise TypeError(f"Xi stress artifact missing required public surface {name!r}")
    for name in ("values", "radial_derivatives", "handoff_at_Xi", "save_configuration"):
        if not callable(getattr(field, name)):
            raise TypeError(f"{name} must be callable")
    if abs(float(field.X_b_ref) - 100.0) > 0.0 or abs(float(field.X_i) - 110.0) > 0.0:
        raise ValueError("A4 Xi audit is frozen to the public X=100 -> Xi=110 reference domain")


def local_polynomial_first_weights(nodes: np.ndarray = LOCAL_POLY_NODES) -> np.ndarray:
    """Weights for d/ds at s=0 of the degree-(n-1) interpolating polynomial."""
    s = np.asarray(nodes, dtype=float)
    if s.ndim != 1 or s.size < 3 or np.unique(s).size != s.size:
        raise ValueError("local-polynomial nodes must be distinct one-dimensional values")
    # V[j,k] = s_j^k.  Need sum_j w_j s_j^k = d/ds s^k |_{0},
    # which is 1 for k=1 and 0 for all other monomials through degree n-1.
    V = np.vander(s, N=s.size, increasing=True)
    rhs = np.zeros(s.size, dtype=float)
    rhs[1] = 1.0
    weights = np.linalg.solve(V.T, rhs)
    if not np.all(np.isfinite(weights)):
        raise RuntimeError("failed to construct finite local-polynomial derivative weights")
    return weights


LOCAL_POLY_WEIGHTS = local_polynomial_first_weights()


def local_polynomial_first(
    fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    x: np.ndarray,
    eta: np.ndarray,
    halfwidth: float,
) -> np.ndarray:
    """Differentiate one public scalar using the frozen nonuniform micro-stencil."""
    if not (math.isfinite(halfwidth) and halfwidth > 0.0):
        raise ValueError("halfwidth must be finite and positive")
    xx = np.asarray(x, dtype=float)
    ee = np.asarray(eta, dtype=float)
    out = np.zeros_like(xx, dtype=float)
    for node, weight in zip(LOCAL_POLY_NODES, LOCAL_POLY_WEIGHTS, strict=True):
        sampled = _finite(fn(xx + node * halfwidth, ee), "local-polynomial sample")
        out += weight * sampled
    return out / halfwidth


def _values_evaluator(field: Any, key: str) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    def evaluate(x: np.ndarray, eta: np.ndarray) -> np.ndarray:
        payload = field.values(x, eta)
        if key not in payload:
            raise KeyError(f"public values payload is missing {key!r}")
        return _finite(payload[key], key)

    return evaluate


def _all_independent_derivatives(
    field: Any, x: np.ndarray, eta: np.ndarray, halfwidth: float
) -> dict[str, np.ndarray]:
    """Share public value evaluations across all audited scalar derivatives."""
    xx = np.asarray(x, dtype=float)
    ee = np.asarray(eta, dtype=float)
    keys = (
        "E_reference",
        "Pi_reference_autonomous",
        "p1_reference",
        "N_s_reference_autonomous",
        "n_s_reference_autonomous",
    )
    accum = {key: np.zeros_like(xx, dtype=float) for key in keys}
    for node, weight in zip(LOCAL_POLY_NODES, LOCAL_POLY_WEIGHTS, strict=True):
        payload = field.values(xx + node * halfwidth, ee)
        for key in keys:
            accum[key] += weight * _finite(payload[key], f"local-polynomial {key}")
    return {
        "E_reference_X": accum["E_reference"] / halfwidth,
        "Pi_reference_autonomous_X": accum["Pi_reference_autonomous"] / halfwidth,
        "p1_reference_X": accum["p1_reference"] / halfwidth,
        "N_s_reference_autonomous_X": accum["N_s_reference_autonomous"] / halfwidth,
        "n_s_reference_autonomous_X": accum["n_s_reference_autonomous"] / halfwidth,
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


def _stability(levels: list[np.ndarray]) -> dict[str, float | bool]:
    if len(levels) != 3:
        raise ValueError("exactly three independent resolution levels are required")
    coarse_medium = _rms(np.asarray(levels[0]) - np.asarray(levels[1]))
    medium_fine = _rms(np.asarray(levels[1]) - np.asarray(levels[2]))
    ratio = math.inf if medium_fine == 0.0 else coarse_medium / medium_fine
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
    if not (eta_lo < eta_hi):
        raise ValueError("invalid public eta interval")
    margin = max(0.08 * (eta_hi - eta_lo), 2.0e-3)
    if eta_lo + margin >= eta_hi - margin:
        raise ValueError("public eta interval is too small for frozen A4 probes")

    rng = np.random.default_rng(PROTOCOL.seed)
    x_random = rng.uniform(PROTOCOL.x_min, PROTOCOL.x_max, PROTOCOL.random_offgrid_count)
    eta_random = rng.uniform(
        eta_lo + margin,
        eta_hi - margin,
        PROTOCOL.random_offgrid_count,
    )

    # Explicit center/axis-near probes at distinct X values.
    center_eta = np.asarray([0.0, -1.0e-8, 1.0e-8, 0.0, -1.0e-6, 1.0e-6])
    center_x = np.asarray([101.1, 102.7, 104.2, 105.8, 107.3, 108.9])
    if np.any((center_eta < eta_lo) | (center_eta > eta_hi)):
        raise ValueError("public eta interval does not include required center probes")

    # Four explicit edge-near points still leave room for the coarse 0.24 stencil.
    edge_x = np.asarray([100.30, 100.48, 109.52, 109.70])
    edge_eta = np.asarray(
        [
            max(eta_lo + margin, -0.55),
            min(eta_hi - margin, 0.35),
            max(eta_lo + margin, -0.25),
            min(eta_hi - margin, 0.55),
        ]
    )

    x = np.concatenate([x_random, center_x, edge_x])
    eta = np.concatenate([eta_random, center_eta, edge_eta])
    endpoint_eta = np.linspace(max(eta_lo + margin, -0.6), min(eta_hi - margin, 0.6), 5)

    coarse = max(PROTOCOL.x_halfwidths)
    if np.min(x) - coarse < float(field.X_b_ref) or np.max(x) + coarse > float(field.X_i):
        raise ValueError("frozen A4 points do not leave room for the coarse local stencil")
    return x, eta, endpoint_eta


def _mutation_checks(
    field: Any,
    x: np.ndarray,
    eta: np.ndarray,
    production: Mapping[str, np.ndarray],
    independent_fine: Mapping[str, np.ndarray],
) -> dict[str, bool]:
    p1_scaled = _match_metrics(
        0.99 * np.asarray(production["p1_reference_X"]),
        independent_fine["p1_reference_X"],
    )
    Ns_scaled = _match_metrics(
        0.90 * np.asarray(production["N_s_reference_autonomous_X"]),
        independent_fine["N_s_reference_autonomous_X"],
    )
    prod_p1_detected = not (
        p1_scaled["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
        and p1_scaled["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
    )
    prod_Ns_detected = not (
        Ns_scaled["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
        and Ns_scaled["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
    )

    p1 = _values_evaluator(field, "p1_reference")
    Ns = _values_evaluator(field, "N_s_reference_autonomous")
    Pi = _values_evaluator(field, "Pi_reference_autonomous")

    def p1_plus_linear(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return p1(xx, ee) + 1.0e-3 * np.asarray(xx, dtype=float)

    def Ns_minus_linear(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return Ns(xx, ee) - 1.0e-3 * np.asarray(xx, dtype=float)

    def Pi_plus_linear(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return Pi(xx, ee) + 1.0e-3 * np.asarray(xx, dtype=float)

    fine = PROTOCOL.x_halfwidths[-1]
    p1_mut = local_polynomial_first(p1_plus_linear, x, eta, fine)
    Ns_mut = local_polynomial_first(Ns_minus_linear, x, eta, fine)
    Pi_mut = local_polynomial_first(Pi_plus_linear, x, eta, fine)
    public_p1_detected = _rms(p1_mut - independent_fine["p1_reference_X"]) >= 5.0e-4
    public_Ns_detected = (
        _rms(Ns_mut - independent_fine["N_s_reference_autonomous_X"]) >= 5.0e-4
    )
    public_Pi_detected = (
        _rms(Pi_mut - independent_fine["Pi_reference_autonomous_X"]) >= 5.0e-4
    )
    return {
        "production_p1_X_times_0p99_detected": bool(prod_p1_detected),
        "production_Ns_X_times_0p90_detected": bool(prod_Ns_detected),
        "public_p1_plus_1e_minus3_X_detected": bool(public_p1_detected),
        "public_Ns_minus_1e_minus3_X_detected": bool(public_Ns_detected),
        "public_Pi_plus_1e_minus3_X_detected": bool(public_Pi_detected),
    }


def materialize_reference_stress_xi_audit(field: Any) -> dict[str, Any]:
    """Run the frozen implementation-distinct audit on one public Xi artifact."""
    _require_contract(field)
    x, eta, endpoint_eta = _sample_points(field)
    values = field.values(x, eta)
    production = field.radial_derivatives(x, eta)

    audited_keys = (
        "E_reference_X",
        "Pi_reference_autonomous_X",
        "p1_reference_X",
        "N_s_reference_autonomous_X",
        "n_s_reference_autonomous_X",
    )
    for key in audited_keys:
        if key not in production:
            raise KeyError(f"production radial derivative payload missing {key!r}")

    levels = [
        _all_independent_derivatives(field, x, eta, halfwidth)
        for halfwidth in PROTOCOL.x_halfwidths
    ]
    finest = levels[-1]

    level_metrics: dict[str, list[dict[str, float]]] = {}
    stability: dict[str, dict[str, float | bool]] = {}
    worst: dict[str, dict[str, float]] = {}
    derivative_gates: dict[str, bool] = {}
    for key in audited_keys:
        prod = _finite(production[key], key)
        independent_levels = [_finite(level[key], f"independent {key}") for level in levels]
        level_metrics[key] = [
            _match_metrics(prod, independent) for independent in independent_levels
        ]
        stability[key] = _stability(independent_levels)
        worst[key] = _worst(x, eta, prod, independent_levels[-1])
        fine_metrics = level_metrics[key][-1]
        derivative_gates[key] = bool(
            fine_metrics["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
            and fine_metrics["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
            and stability[key]["passed"]
        )

    # Algebraic closure in public values and, independently, in A4 derivatives.
    h = float(field.h)
    L = 1.0 - 2.0 * h * np.asarray(eta, dtype=float) ** 2
    Ns = _finite(values["N_s_reference_autonomous"], "N_s_reference_autonomous")
    ns = _finite(values["n_s_reference_autonomous"], "n_s_reference_autonomous")
    algebraic_closure = _relative_closure(Ns, L * ns)
    derivative_closure = _relative_closure(
        finest["N_s_reference_autonomous_X"],
        L * finest["n_s_reference_autonomous_X"],
    )

    # Plateau F/U should remain exactly X-independent on this reference interval.
    x100 = np.full_like(x, float(field.X_b_ref), dtype=float)
    start = field.values(x100, eta)
    F = _finite(values["F_reference"], "F_reference")
    U = _finite(values["U_reference"], "U_reference")
    F_start = _finite(start["F_reference"], "F_reference@X100")
    U_start = _finite(start["U_reference"], "U_reference@X100")
    plateau_F_abs_max = float(np.max(np.abs(F - F_start)))
    plateau_U_abs_max = float(np.max(np.abs(U - U_start)))

    # Public Xi handoff must replay direct public values/derivatives at X_i.
    Xi = np.full_like(endpoint_eta, float(field.X_i), dtype=float)
    direct_values = field.values(Xi, endpoint_eta)
    direct_derivatives = field.radial_derivatives(Xi, endpoint_eta)
    handoff = field.handoff_at_Xi(endpoint_eta)
    handoff_keys = (
        "F_reference",
        "U_reference",
        "E_reference",
        "Pi_reference_autonomous",
        "p1_reference",
        "N_s_reference_autonomous",
        "n_s_reference_autonomous",
        "p2_reference_autonomous",
        "E_reference_X",
        "Pi_reference_autonomous_X",
        "p1_reference_X",
        "N_s_reference_autonomous_X",
        "n_s_reference_autonomous_X",
    )
    xi_handoff_abs_max = 0.0
    for key in handoff_keys:
        reference_payload = direct_derivatives if key.endswith("_X") else direct_values
        if key not in handoff or key not in reference_payload:
            raise KeyError(f"Xi public handoff missing required replay key {key!r}")
        xi_handoff_abs_max = max(
            xi_handoff_abs_max,
            float(
                np.max(
                    np.abs(
                        _finite(handoff[key], f"handoff {key}")
                        - _finite(reference_payload[key], f"direct {key}")
                    )
                )
            ),
        )

    # Nontriviality: scalar state and nonzero derivative content must remain visible.
    value_rms = {
        key: _rms(values[key])
        for key in (
            "E_reference",
            "Pi_reference_autonomous",
            "p1_reference",
            "N_s_reference_autonomous",
            "n_s_reference_autonomous",
            "p2_reference_autonomous",
        )
    }
    derivative_rms = {key: _rms(production[key]) for key in audited_keys}
    nontrivial = bool(
        all(v >= PROTOCOL.value_rms_floor for v in value_rms.values())
        and all(v >= PROTOCOL.derivative_rms_floor for v in derivative_rms.values())
    )

    mutations = _mutation_checks(field, x, eta, production, finest)

    gates = {
        "all_derivative_surfaces": bool(all(derivative_gates.values())),
        "algebraic_Ns_equals_L_ns": bool(
            algebraic_closure <= PROTOCOL.algebraic_closure_relative_gate
        ),
        "independent_derivative_Ns_equals_L_ns": bool(
            derivative_closure <= PROTOCOL.derivative_closure_relative_gate
        ),
        "reference_F_plateau": bool(plateau_F_abs_max <= PROTOCOL.plateau_abs_gate),
        "reference_U_plateau": bool(plateau_U_abs_max <= PROTOCOL.plateau_abs_gate),
        "Xi_public_handoff_replay": bool(
            xi_handoff_abs_max <= PROTOCOL.xi_handoff_abs_gate
        ),
        "nontrivial": nontrivial,
        "all_mutations_detected": bool(all(mutations.values())),
    }
    passed = bool(all(gates.values()))

    truth = dict(TRUTH_BOUNDARY)
    truth["autonomous_reference_stress_to_Xi_independently_audited"] = passed

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "artifact_semantic_sha256": str(field.semantic_sha256),
        "protocol": asdict(PROTOCOL),
        "operator": {
            "name": "nonuniform_degree6_local_polynomial_first_derivative",
            "nodes": LOCAL_POLY_NODES.tolist(),
            "weights": LOCAL_POLY_WEIGHTS.tolist(),
        },
        "sample_count": int(np.asarray(x).size),
        "resolution_metrics": level_metrics,
        "resolution_stability": stability,
        "worst_witness": worst,
        "algebraic_closure_relative": algebraic_closure,
        "independent_derivative_closure_relative": derivative_closure,
        "plateau_F_abs_max": plateau_F_abs_max,
        "plateau_U_abs_max": plateau_U_abs_max,
        "Xi_public_handoff_abs_max": xi_handoff_abs_max,
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
    """Save/reload Agent-1 #933 and audit only its rebound public surfaces."""
    from .kokuno_pa10_autonomous_reference_stress_xi110 import (
        KokunoPA10AutonomousReferenceStressToXi110,
    )

    original = KokunoPA10AutonomousReferenceStressToXi110()
    with tempfile.TemporaryDirectory(prefix="kokuno-a4-xi110-") as tmp:
        path = Path(tmp) / "reference-stress-xi110.json"
        original.save_configuration(path)
        rebound = KokunoPA10AutonomousReferenceStressToXi110.load_configuration(path)
        if rebound.semantic_sha256 != original.semantic_sha256:
            raise AssertionError("save/reload changed the Xi reference-stress semantic identity")
        receipt = materialize_reference_stress_xi_audit(rebound)
    receipt["save_reload_exact_replay"] = True
    # Rebind checksum after adding the replay fact.
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
