"""Independent Agent-4 audit of the candidate-side fixed-kappa PA.10 continuation.

Agent 1 PR #925 exposes a serializable candidate-side continuation on
``X_1 <= X <= 100`` with public profile values ``F,U,E`` and separate analytic
radial derivatives.  This module treats those derivative surfaces as quantities
under audit.  The numerical reference is reconstructed only from save/reloaded
public ``values(X, eta)`` evaluations with a centered eighth-order finite-
difference operator at three frozen resolutions.

This is a source-coordinate profile audit, not complete Navier--Stokes
validation.  The current route still lacks the final ``100<X<110`` bridge,
global Cartesian leading velocity, matched Cartesian pressure, preregistered
restricted forcing and a materialized correction velocity.  The repository-wide
``1e-3`` momentum and ``1e-5`` divergence gates remain unchanged and unevaluated.
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

TASK = "K4-VAL-086"
SCHEMA = "kokuno-a4-fixed-kappa-x100-independent-audit-v1"
A1_FIXED_KAPPA_PR = 925
A1_FIXED_KAPPA_HEAD = "1f3dec711f8ce2052276cd951b2f44d2849579bf"
A1_FIXED_KAPPA_SOURCE_BLOB = "818d965d82bd92505bd6d3e805e698e1e0da64d8"

FROZEN_FINAL_GATES = {
    "normalized_momentum_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
}


@dataclass(frozen=True)
class FrozenFixedKappaX100AuditProtocol:
    seed: int = 9173601
    random_offgrid_count: int = 128
    x_fraction_min: float = 0.05
    x_fraction_max: float = 0.95
    relative_x_steps: tuple[float, float, float] = (1.6e-3, 8.0e-4, 4.0e-4)
    derivative_relative_rms_gate: float = 5.0e-3
    derivative_relative_max_gate: float = 2.0e-2
    refinement_ratio_gate: float = 20.0
    refinement_floor: float = 2.0e-9
    value_identity_relative_gate: float = 5.0e-12
    derivative_identity_relative_gate: float = 5.0e-5
    endpoint_composition_abs_gate: float = 2.0e-12
    value_rms_floor: float = 1.0e-10
    derivative_rms_floor: float = 1.0e-12


PROTOCOL = FrozenFixedKappaX100AuditProtocol()

TRUTH_BOUNDARY = {
    "candidate_side_fixed_kappa_continuation_materialized": True,
    "candidate_side_fixed_kappa_continuation_independently_audited": False,
    "source_prepared_appendixA_Pi0_materialized": False,
    "source_prepared_reference_nsr_materialized": False,
    "source_exact_fixed_kappa_continuation_materialized": False,
    "final_interpolation_to_Xi_materialized": False,
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
        "agent1_pr": A1_FIXED_KAPPA_PR,
        "agent1_exact_head": A1_FIXED_KAPPA_HEAD,
        "agent1_source_blob": A1_FIXED_KAPPA_SOURCE_BLOB,
        "scientific_reference_surface": "save/reloaded values(X,eta) -> F,U,E",
        "production_surfaces_under_audit": [
            "radial_derivatives(X,eta) -> F_X,U_X,E_X"
        ],
        "independent_operator": (
            "centered eighth-order X derivative reconstructed only from public "
            "save/reloaded values(X,eta) evaluations"
        ),
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
        "endpoint_values",
        "save_configuration",
        "semantic_sha256",
        "X_1",
        "X_b_ref",
        "eta_interval",
        "kappa_0",
    )
    for name in required:
        if not hasattr(field, name):
            raise TypeError(f"fixed-kappa artifact missing required public surface {name!r}")
    for name in ("values", "radial_derivatives", "endpoint_values", "save_configuration"):
        if not callable(getattr(field, name)):
            raise TypeError(f"{name} must be callable")


def _field_evaluator(field: Any, key: str) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    def evaluate(x: np.ndarray, eta: np.ndarray) -> np.ndarray:
        payload = field.values(x, eta)
        if key not in payload:
            raise KeyError(f"public values payload is missing {key!r}")
        return _finite(payload[key], key)

    return evaluate


def centered_fd8_first(
    fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    x: np.ndarray,
    eta: np.ndarray,
    relative_step: float,
) -> np.ndarray:
    """Centered eighth-order first derivative with pointwise scale-aware step."""
    x = np.asarray(x, dtype=float)
    eta = np.asarray(eta, dtype=float)
    h = relative_step * np.maximum(1.0, np.abs(x))
    fm4 = _finite(fn(x - 4.0 * h, eta), "fd8 fm4")
    fm3 = _finite(fn(x - 3.0 * h, eta), "fd8 fm3")
    fm2 = _finite(fn(x - 2.0 * h, eta), "fd8 fm2")
    fm1 = _finite(fn(x - h, eta), "fd8 fm1")
    fp1 = _finite(fn(x + h, eta), "fd8 fp1")
    fp2 = _finite(fn(x + 2.0 * h, eta), "fd8 fp2")
    fp3 = _finite(fn(x + 3.0 * h, eta), "fd8 fp3")
    fp4 = _finite(fn(x + 4.0 * h, eta), "fd8 fp4")
    return (
        3.0 * fm4
        - 32.0 * fm3
        + 168.0 * fm2
        - 672.0 * fm1
        + 672.0 * fp1
        - 168.0 * fp2
        + 32.0 * fp3
        - 3.0 * fp4
    ) / (840.0 * h)


def _profile_fd8_derivatives(
    field: Any, x: np.ndarray, eta: np.ndarray, relative_step: float
) -> dict[str, np.ndarray]:
    """Differentiate F/U/E together so each public stencil evaluation is shared."""
    x = np.asarray(x, dtype=float)
    eta = np.asarray(eta, dtype=float)
    h = relative_step * np.maximum(1.0, np.abs(x))
    offsets = (-4.0, -3.0, -2.0, -1.0, 1.0, 2.0, 3.0, 4.0)
    coefficients = (3.0, -32.0, 168.0, -672.0, 672.0, -168.0, 32.0, -3.0)
    sampled = [field.values(x + offset * h, eta) for offset in offsets]
    out: dict[str, np.ndarray] = {}
    for key in ("F_fixed_kappa", "U_fixed_kappa", "E_fixed_kappa"):
        numerator = np.zeros_like(x, dtype=float)
        for coefficient, payload in zip(coefficients, sampled, strict=True):
            numerator += coefficient * _finite(payload[key], f"fd8 {key}")
        out[key + "_X"] = numerator / (840.0 * h)
    return out


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
    x1 = float(field.X_1)
    xb = float(field.X_b_ref)
    eta_lo, eta_hi = (float(v) for v in field.eta_interval)
    if not (math.isfinite(x1) and math.isfinite(xb) and 0.0 < x1 < xb):
        raise ValueError("invalid public fixed-kappa X interval")
    if not (math.isfinite(eta_lo) and math.isfinite(eta_hi) and eta_lo < eta_hi):
        raise ValueError("invalid public eta interval")

    span = xb - x1
    rng = np.random.default_rng(PROTOCOL.seed)
    frac = rng.uniform(
        PROTOCOL.x_fraction_min,
        PROTOCOL.x_fraction_max,
        PROTOCOL.random_offgrid_count,
    )
    x_random = x1 + frac * span
    eta_margin = 0.08 * (eta_hi - eta_lo)
    eta_random = rng.uniform(
        eta_lo + eta_margin,
        eta_hi - eta_margin,
        PROTOCOL.random_offgrid_count,
    )

    center_eta = np.asarray([0.0, -1.0e-8, 1.0e-8, 0.0, -1.0e-6, 1.0e-6])
    center_frac = np.asarray([0.12, 0.28, 0.44, 0.61, 0.77, 0.90])
    x_center = x1 + center_frac * span
    if np.any((center_eta < eta_lo) | (center_eta > eta_hi)):
        raise ValueError("public eta interval does not include required center probes")

    x = np.concatenate([x_random, x_center])
    eta = np.concatenate([eta_random, center_eta])
    endpoint_eta = np.linspace(max(eta_lo, -0.6), min(eta_hi, 0.6), 5)
    return x, eta, endpoint_eta


def _mutation_checks(
    field: Any,
    x: np.ndarray,
    eta: np.ndarray,
    production_F_X: np.ndarray,
    production_U_X: np.ndarray,
    independent_F_X: np.ndarray,
    independent_U_X: np.ndarray,
) -> dict[str, bool]:
    scaled_F = _match_metrics(0.99 * production_F_X, independent_F_X)
    scaled_U = _match_metrics(0.90 * production_U_X, independent_U_X)
    prod_F_detected = not (
        scaled_F["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
        and scaled_F["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
    )
    prod_U_detected = not (
        scaled_U["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
        and scaled_U["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
    )

    base_F = _field_evaluator(field, "F_fixed_kappa")
    base_U = _field_evaluator(field, "U_fixed_kappa")

    def F_plus_linear(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return base_F(xx, ee) + 1.0e-3 * np.asarray(xx, dtype=float)

    def U_minus_linear(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return base_U(xx, ee) - 1.0e-3 * np.asarray(xx, dtype=float)

    mutated_F_X = centered_fd8_first(F_plus_linear, x, eta, PROTOCOL.relative_x_steps[-1])
    mutated_U_X = centered_fd8_first(U_minus_linear, x, eta, PROTOCOL.relative_x_steps[-1])
    public_F_detected = _rms(mutated_F_X - independent_F_X) >= 5.0e-4
    public_U_detected = _rms(mutated_U_X - independent_U_X) >= 5.0e-4
    return {
        "production_F_X_times_0p99_detected": bool(prod_F_detected),
        "production_U_X_times_0p90_detected": bool(prod_U_detected),
        "public_F_plus_1e_minus3_X_detected": bool(public_F_detected),
        "public_U_minus_1e_minus3_X_detected": bool(public_U_detected),
    }


def materialize_fixed_kappa_x100_audit(field: Any) -> dict[str, Any]:
    """Run the frozen independent audit on one public fixed-kappa artifact."""
    _require_contract(field)
    x, eta, endpoint_eta = _sample_points(field)
    values = field.values(x, eta)
    F = _finite(values["F_fixed_kappa"], "F")
    U = _finite(values["U_fixed_kappa"], "U")
    E = _finite(values["E_fixed_kappa"], "E")

    production = field.radial_derivatives(x, eta)
    production_F_X = _finite(production["F_fixed_kappa_X"], "F_X")
    production_U_X = _finite(production["U_fixed_kappa_X"], "U_X")
    production_E_X = _finite(production["E_fixed_kappa_X"], "E_X")

    independent_levels = [
        _profile_fd8_derivatives(field, x, eta, h) for h in PROTOCOL.relative_x_steps
    ]
    independent_F_X = [level["F_fixed_kappa_X"] for level in independent_levels]
    independent_U_X = [level["U_fixed_kappa_X"] for level in independent_levels]
    independent_E_X = [level["E_fixed_kappa_X"] for level in independent_levels]

    F_metrics = [_match_metrics(production_F_X, level) for level in independent_F_X]
    U_metrics = [_match_metrics(production_U_X, level) for level in independent_U_X]
    E_metrics = [_match_metrics(production_E_X, level) for level in independent_E_X]
    F_stability = _stability(independent_F_X)
    U_stability = _stability(independent_U_X)
    E_stability = _stability(independent_E_X)

    value_identity = _relative_closure(E, np.sqrt(2.0 * x) * F)
    independent_E_from_F = F / np.sqrt(2.0 * x) + np.sqrt(2.0 * x) * independent_F_X[-1]
    derivative_identity = _relative_closure(independent_E_X[-1], independent_E_from_F)

    x1 = float(field.X_1)
    endpoint_x = np.full_like(endpoint_eta, x1)
    continued = field.values(endpoint_x, endpoint_eta)
    endpoint = field.endpoint_values(endpoint_eta)
    endpoint_F_abs = float(np.max(np.abs(
        _finite(continued["F_fixed_kappa"], "endpoint F")
        - _finite(endpoint["F_activation"], "activation F")
    )))
    endpoint_U_abs = float(np.max(np.abs(
        _finite(continued["U_fixed_kappa"], "endpoint U")
        - _finite(endpoint["U_activation"], "activation U")
    )))

    fine_metrics = (F_metrics[-1], U_metrics[-1], E_metrics[-1])
    derivative_match_passed = all(
        m["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
        and m["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
        for m in fine_metrics
    )
    stability_passed = bool(F_stability["passed"] and U_stability["passed"] and E_stability["passed"])
    identity_passed = bool(
        value_identity <= PROTOCOL.value_identity_relative_gate
        and derivative_identity <= PROTOCOL.derivative_identity_relative_gate
        and endpoint_F_abs <= PROTOCOL.endpoint_composition_abs_gate
        and endpoint_U_abs <= PROTOCOL.endpoint_composition_abs_gate
    )
    nontriviality_passed = bool(
        min(_rms(F), _rms(E)) >= PROTOCOL.value_rms_floor
        and _rms(U) >= PROTOCOL.value_rms_floor
        and min(_rms(production_F_X), _rms(production_U_X), _rms(production_E_X))
        >= PROTOCOL.derivative_rms_floor
    )
    mutations = _mutation_checks(
        field,
        x,
        eta,
        production_F_X,
        production_U_X,
        independent_F_X[-1],
        independent_U_X[-1],
    )
    mutation_passed = all(mutations.values())
    passed = bool(
        derivative_match_passed
        and stability_passed
        and identity_passed
        and nontriviality_passed
        and mutation_passed
    )

    truth = dict(TRUTH_BOUNDARY)
    truth["candidate_side_fixed_kappa_continuation_independently_audited"] = passed
    result = {
        "task": TASK,
        "schema": SCHEMA,
        "candidate_semantic_sha256": str(field.semantic_sha256),
        "protocol": asdict(PROTOCOL),
        "sample_count": int(x.size),
        "resolution_metrics": {
            "F_X": F_metrics,
            "U_X": U_metrics,
            "E_X": E_metrics,
        },
        "stability": {
            "F_X": F_stability,
            "U_X": U_stability,
            "E_X": E_stability,
        },
        "fine_worst_witness": {
            "F_X": _worst(x, eta, production_F_X, independent_F_X[-1]),
            "U_X": _worst(x, eta, production_U_X, independent_U_X[-1]),
            "E_X": _worst(x, eta, production_E_X, independent_E_X[-1]),
        },
        "profile_identities": {
            "E_equals_sqrt_2X_F_relative": value_identity,
            "independent_E_X_chain_relative": derivative_identity,
            "activation_endpoint_F_abs_max": endpoint_F_abs,
            "activation_endpoint_U_abs_max": endpoint_U_abs,
        },
        "nontriviality": {
            "F_rms": _rms(F),
            "U_rms": _rms(U),
            "E_rms": _rms(E),
            "F_X_rms": _rms(production_F_X),
            "U_X_rms": _rms(production_U_X),
            "E_X_rms": _rms(production_E_X),
        },
        "mutation_checks": mutations,
        "gate_map": {
            "fine_derivative_match": derivative_match_passed,
            "three_resolution_stability": stability_passed,
            "public_profile_identities": identity_passed,
            "nontriviality": nontriviality_passed,
            "mutation_firewall": mutation_passed,
        },
        "passed": passed,
        "frozen_final_gates": dict(FROZEN_FINAL_GATES),
        "truth_boundary": truth,
    }
    result["receipt_sha256"] = _sha256(result)
    return result


def materialize_real_fixed_kappa_x100_audit() -> dict[str, Any]:
    """Save/reload the real #925 artifact and run the frozen audit on the rebound object."""
    from .kokuno_pa10_fixed_kappa_x100 import KokunoPA10FixedKappaContinuationToX100

    original = KokunoPA10FixedKappaContinuationToX100()
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "fixed_kappa_x100.json"
        original.save_configuration(path)
        rebound = KokunoPA10FixedKappaContinuationToX100.load_configuration(path)
        probe_x = np.asarray([rebound.X_1, 0.5 * (rebound.X_1 + rebound.X_b_ref), rebound.X_b_ref])
        probe_eta = np.asarray([-0.2, 0.0, 0.2])
        original_values = original.values(probe_x, probe_eta)
        rebound_values = rebound.values(probe_x, probe_eta)
        replay_exact = all(
            np.array_equal(
                np.asarray(original_values[key]), np.asarray(rebound_values[key])
            )
            for key in ("F_fixed_kappa", "U_fixed_kappa", "E_fixed_kappa")
        )
        if not replay_exact or rebound.semantic_sha256 != original.semantic_sha256:
            raise RuntimeError("save/reload changed fixed-kappa public artifact identity")
        report = materialize_fixed_kappa_x100_audit(rebound)
        report["save_reload_exact_replay"] = True
        report["receipt_sha256"] = _sha256({k: v for k, v in report.items() if k != "receipt_sha256"})
        return report


def save_real_receipt(path: str | Path) -> dict[str, Any]:
    report = materialize_real_fixed_kappa_x100_audit()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    report = save_real_receipt(args.out)
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    _main()
