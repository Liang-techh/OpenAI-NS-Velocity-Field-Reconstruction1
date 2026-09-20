"""Independent Agent-4 audit of the autonomous PA.10 reference-stress artifact.

The upstream Agent-1 artifact #918 exposes a serializable source-coordinate
continuation of the axial stress primitive through ``X=100`` using the
repository-autonomous pressure seed.  Its public value surface returns

    N_s(X, eta), n_s(X, eta), p_2(X, eta)

and it separately exposes source-ODE radial derivatives.  This module treats
those derivatives as quantities under audit.  The numerical reference is
reconstructed only from save/reloaded public ``values(X, eta)`` evaluations by
an implementation-distinct centered sixth-order finite-difference operator at
three frozen resolutions.

This is deliberately *not* a full Navier--Stokes validation.  The upstream
stress still depends on a repository-autonomous axis-pressure seed and there is
no corrected/global Cartesian leading field, matched Cartesian pressure,
preregistered restricted forcing, or Agent-3 correction velocity on this seam.
The repository-wide momentum/divergence gates remain unchanged and unevaluated.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import tempfile
from typing import Any, Callable, Mapping

import numpy as np


SCHEMA = "kokuno-a4-autonomous-reference-ns-independent-audit-v1"
A1_AUTONOMOUS_REFERENCE_NS_PR = 918
A1_AUTONOMOUS_REFERENCE_NS_HEAD = "c5ccdb8fd790c62cfb760139e5ec108c26d88a81"
A1_AUTONOMOUS_REFERENCE_NS_SOURCE_BLOB = "44159b63c440ee9764261dad41fb3f9d3fec7b0d"

FROZEN_FINAL_GATES = {
    "normalized_momentum_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
}


@dataclass(frozen=True)
class FrozenAutonomousReferenceStressAuditProtocol:
    seed: int = 9173591
    random_offgrid_count: int = 128
    x_fraction_min: float = 0.03
    x_fraction_max: float = 0.97
    eta_min: float = -0.72
    eta_max: float = 0.72
    relative_x_steps: tuple[float, float, float] = (2.0e-3, 1.0e-3, 5.0e-4)
    derivative_relative_rms_gate: float = 5.0e-3
    derivative_relative_max_gate: float = 2.0e-2
    refinement_ratio_gate: float = 12.0
    refinement_floor: float = 2.0e-9
    algebraic_relative_gate: float = 5.0e-10
    x0_composition_abs_gate: float = 2.0e-10
    value_rms_floor: float = 1.0e-10
    derivative_rms_floor: float = 1.0e-12


PROTOCOL = FrozenAutonomousReferenceStressAuditProtocol()

TRUTH_BOUNDARY = {
    "autonomous_reference_stress_independently_audited": False,
    "repository_autonomous_pressure_reference_stress_materialized": True,
    "source_prepared_appendixA_Pi0_materialized": False,
    "source_prepared_reference_nsr_materialized": False,
    "matched_global_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "correction_velocity_materialized": False,
    "restricted_forcing_materialized": False,
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
        "schema": SCHEMA,
        "agent1_pr": A1_AUTONOMOUS_REFERENCE_NS_PR,
        "agent1_exact_head": A1_AUTONOMOUS_REFERENCE_NS_HEAD,
        "agent1_source_blob": A1_AUTONOMOUS_REFERENCE_NS_SOURCE_BLOB,
        "scientific_reference_surface": "values(X,eta) -> N_s,n_s,p_2",
        "production_surfaces_under_audit": [
            "radial_derivatives(X,eta) -> N_s_X,n_s_X",
        ],
        "independent_operator": (
            "centered sixth-order X derivative reconstructed only from public "
            "save/reloaded values(X,eta) evaluations"
        ),
        "protocol": asdict(PROTOCOL),
        "frozen_final_gates": dict(FROZEN_FINAL_GATES),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def _finite_array(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} returned non-finite values")
    return out


def _require_contract(stress: Any) -> None:
    for name in (
        "values",
        "radial_derivatives",
        "initial_values",
        "semantic_sha256",
        "X_0",
        "X_b_ref",
        "h",
    ):
        if not hasattr(stress, name):
            raise TypeError(f"stress artifact is missing required public surface {name!r}")
    for name in ("values", "radial_derivatives", "initial_values"):
        if not callable(getattr(stress, name)):
            raise TypeError(f"{name} must be callable")


def _rms(value: np.ndarray) -> float:
    arr = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(arr * arr)))


def _field_evaluator(stress: Any, field: str) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    def evaluate(x: np.ndarray, eta: np.ndarray) -> np.ndarray:
        payload = stress.values(x, eta)
        if field not in payload:
            raise KeyError(f"public values payload is missing {field!r}")
        return _finite_array(payload[field], field)

    return evaluate


def _fd6_centered_first(
    fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    x: np.ndarray,
    eta: np.ndarray,
    relative_step: float,
) -> np.ndarray:
    """Centered sixth-order radial derivative with pointwise scale-aware step."""
    x = np.asarray(x, dtype=float)
    eta = np.asarray(eta, dtype=float)
    h = relative_step * np.maximum(1.0, np.abs(x))
    fm3 = _finite_array(fn(x - 3.0 * h, eta), "fd6 fm3")
    fm2 = _finite_array(fn(x - 2.0 * h, eta), "fd6 fm2")
    fm1 = _finite_array(fn(x - h, eta), "fd6 fm1")
    fp1 = _finite_array(fn(x + h, eta), "fd6 fp1")
    fp2 = _finite_array(fn(x + 2.0 * h, eta), "fd6 fp2")
    fp3 = _finite_array(fn(x + 3.0 * h, eta), "fd6 fp3")
    return (-fm3 + 9.0 * fm2 - 45.0 * fm1 + 45.0 * fp1 - 9.0 * fp2 + fp3) / (
        60.0 * h
    )


def _match_metrics(production: np.ndarray, independent: np.ndarray) -> dict[str, float]:
    production = np.asarray(production, dtype=float)
    independent = np.asarray(independent, dtype=float)
    error = independent - production
    return {
        "absolute_rms": _rms(error),
        "absolute_max": float(np.max(np.abs(error))),
        "relative_rms": _rms(error) / max(_rms(production), 1.0e-30),
        "relative_sampled_max": float(np.max(np.abs(error)))
        / max(float(np.max(np.abs(production))), 1.0e-30),
    }


def _relative_closure(lhs: np.ndarray, rhs: np.ndarray) -> float:
    lhs = np.asarray(lhs, dtype=float)
    rhs = np.asarray(rhs, dtype=float)
    return _rms(lhs - rhs) / max(_rms(lhs), _rms(rhs), 1.0e-30)


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
    error = np.abs(np.asarray(independent) - np.asarray(production))
    idx = int(np.argmax(error))
    return {
        "X": float(np.asarray(x)[idx]),
        "eta": float(np.asarray(eta)[idx]),
        "production": float(np.asarray(production)[idx]),
        "independent": float(np.asarray(independent)[idx]),
        "absolute_error": float(error[idx]),
    }


def _sample_points(stress: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x0 = float(stress.X_0)
    xb = float(stress.X_b_ref)
    if not (math.isfinite(x0) and math.isfinite(xb) and 0.0 < x0 < xb):
        raise ValueError("invalid public X propagation interval")
    span = xb - x0
    rng = np.random.default_rng(PROTOCOL.seed)
    fractions = rng.uniform(
        PROTOCOL.x_fraction_min,
        PROTOCOL.x_fraction_max,
        PROTOCOL.random_offgrid_count,
    )
    x_random = x0 + fractions * span
    eta_random = rng.uniform(PROTOCOL.eta_min, PROTOCOL.eta_max, PROTOCOL.random_offgrid_count)

    edge_fractions = np.asarray([0.01, 0.02, 0.98, 0.99], dtype=float)
    x_edge = x0 + edge_fractions * span
    eta_edge = np.asarray([-0.65, 0.0, 0.35, 0.65], dtype=float)

    x = np.concatenate([x_random, x_edge])
    eta = np.concatenate([eta_random, eta_edge])
    exact_x0_eta = np.asarray([-0.70, -0.25, 0.0, 0.25, 0.70], dtype=float)
    return x, eta, exact_x0_eta


def _mutation_checks(
    stress: Any,
    x: np.ndarray,
    eta: np.ndarray,
    production_N_X: np.ndarray,
    production_n_X: np.ndarray,
    independent_N_X: np.ndarray,
    independent_n_X: np.ndarray,
) -> dict[str, bool]:
    scaled_N = _match_metrics(0.99 * production_N_X, independent_N_X)
    scaled_n = _match_metrics(0.90 * production_n_X, independent_n_X)
    prod_scale_N_detected = not (
        scaled_N["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
        and scaled_N["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
    )
    prod_scale_n_detected = not (
        scaled_n["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
        and scaled_n["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
    )

    base_N = _field_evaluator(stress, "N_s_reference_autonomous")
    base_n = _field_evaluator(stress, "n_s_reference_autonomous")

    def N_plus_linear(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return base_N(xx, ee) + 1.0e-3 * np.asarray(xx, dtype=float)

    def n_plus_linear(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return base_n(xx, ee) - 1.0e-3 * np.asarray(xx, dtype=float)

    mutated_N_X = _fd6_centered_first(
        N_plus_linear, x, eta, PROTOCOL.relative_x_steps[-1]
    )
    mutated_n_X = _fd6_centered_first(
        n_plus_linear, x, eta, PROTOCOL.relative_x_steps[-1]
    )
    public_N_mutation_detected = _rms(mutated_N_X - production_N_X) >= 5.0e-4
    public_n_mutation_detected = _rms(mutated_n_X - production_n_X) >= 5.0e-4

    return {
        "production_Ns_X_times_0p99_detected": bool(prod_scale_N_detected),
        "production_ns_X_times_0p90_detected": bool(prod_scale_n_detected),
        "public_Ns_plus_1e_minus3_X_detected": bool(public_N_mutation_detected),
        "public_ns_minus_1e_minus3_X_detected": bool(public_n_mutation_detected),
    }


def materialize_autonomous_reference_stress_audit(stress: Any) -> dict[str, Any]:
    """Run the frozen implementation-distinct audit on one public stress artifact."""
    _require_contract(stress)
    x, eta, exact_x0_eta = _sample_points(stress)

    values = stress.values(x, eta)
    N_value = _finite_array(values["N_s_reference_autonomous"], "N_s")
    n_value = _finite_array(values["n_s_reference_autonomous"], "n_s")
    p2_value = _finite_array(values["p2_reference_autonomous"], "p_2")

    production = stress.radial_derivatives(x, eta)
    production_N_X = _finite_array(
        production["N_s_reference_autonomous_X"], "N_s_X"
    )
    production_n_X = _finite_array(
        production["n_s_reference_autonomous_X"], "n_s_X"
    )

    N_eval = _field_evaluator(stress, "N_s_reference_autonomous")
    n_eval = _field_evaluator(stress, "n_s_reference_autonomous")
    independent_N_X = [
        _fd6_centered_first(N_eval, x, eta, step) for step in PROTOCOL.relative_x_steps
    ]
    independent_n_X = [
        _fd6_centered_first(n_eval, x, eta, step) for step in PROTOCOL.relative_x_steps
    ]

    N_metrics = [_match_metrics(production_N_X, level) for level in independent_N_X]
    n_metrics = [_match_metrics(production_n_X, level) for level in independent_n_X]
    N_stability = _stability(independent_N_X)
    n_stability = _stability(independent_n_X)

    L = 1.0 - 2.0 * float(stress.h) * eta * eta
    value_relation_error = _relative_closure(N_value, L * n_value)
    derivative_relation_error = _relative_closure(independent_N_X[-1], L * independent_n_X[-1])

    x0 = float(stress.X_0)
    x0_array = np.full_like(exact_x0_eta, x0)
    x0_values = stress.values(x0_array, exact_x0_eta)
    x0_initial = stress.initial_values(exact_x0_eta)
    x0_N_abs = float(
        np.max(
            np.abs(
                _finite_array(x0_values["N_s_reference_autonomous"], "X0 N_s")
                - _finite_array(x0_initial["N_s_initial"], "X0 N_s initial")
            )
        )
    )
    x0_n_abs = float(
        np.max(
            np.abs(
                _finite_array(x0_values["n_s_reference_autonomous"], "X0 n_s")
                - _finite_array(x0_initial["n_s_initial"], "X0 n_s initial")
            )
        )
    )

    fine_N = N_metrics[-1]
    fine_n = n_metrics[-1]
    mutations = _mutation_checks(
        stress,
        x,
        eta,
        production_N_X,
        production_n_X,
        independent_N_X[-1],
        independent_n_X[-1],
    )

    checks = {
        "N_derivative_relative_rms": fine_N["relative_rms"]
        <= PROTOCOL.derivative_relative_rms_gate,
        "N_derivative_relative_sampled_max": fine_N["relative_sampled_max"]
        <= PROTOCOL.derivative_relative_max_gate,
        "n_derivative_relative_rms": fine_n["relative_rms"]
        <= PROTOCOL.derivative_relative_rms_gate,
        "n_derivative_relative_sampled_max": fine_n["relative_sampled_max"]
        <= PROTOCOL.derivative_relative_max_gate,
        "N_resolution_stability": bool(N_stability["passed"]),
        "n_resolution_stability": bool(n_stability["passed"]),
        "public_N_equals_L_n": value_relation_error <= PROTOCOL.algebraic_relative_gate,
        "independent_N_X_equals_L_n_X": derivative_relation_error
        <= PROTOCOL.algebraic_relative_gate,
        "X0_N_composition": x0_N_abs <= PROTOCOL.x0_composition_abs_gate,
        "X0_n_composition": x0_n_abs <= PROTOCOL.x0_composition_abs_gate,
        "N_nontrivial": _rms(N_value) >= PROTOCOL.value_rms_floor,
        "n_nontrivial": _rms(n_value) >= PROTOCOL.value_rms_floor,
        "p2_nontrivial": _rms(p2_value) >= PROTOCOL.value_rms_floor,
        "N_derivative_nontrivial": _rms(production_N_X) >= PROTOCOL.derivative_rms_floor,
        "n_derivative_nontrivial": _rms(production_n_X) >= PROTOCOL.derivative_rms_floor,
        **mutations,
    }
    passed = bool(all(checks.values()))

    truth = dict(TRUTH_BOUNDARY)
    truth["autonomous_reference_stress_independently_audited"] = passed

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "agent1_pr": A1_AUTONOMOUS_REFERENCE_NS_PR,
        "agent1_exact_head": A1_AUTONOMOUS_REFERENCE_NS_HEAD,
        "agent1_source_blob": A1_AUTONOMOUS_REFERENCE_NS_SOURCE_BLOB,
        "stress_semantic_sha256": str(stress.semantic_sha256),
        "protocol": asdict(PROTOCOL),
        "sample_counts": {
            "fresh_random_offgrid": PROTOCOL.random_offgrid_count,
            "edge_near": 4,
            "exact_X0": len(exact_x0_eta),
        },
        "production_rms": {
            "N_s": _rms(N_value),
            "n_s": _rms(n_value),
            "p_2": _rms(p2_value),
            "N_s_X": _rms(production_N_X),
            "n_s_X": _rms(production_n_X),
        },
        "N_three_resolution_match": N_metrics,
        "n_three_resolution_match": n_metrics,
        "N_resolution_stability": N_stability,
        "n_resolution_stability": n_stability,
        "public_value_relation_relative_error": value_relation_error,
        "independent_derivative_relation_relative_error": derivative_relation_error,
        "X0_composition_max_abs": {"N_s": x0_N_abs, "n_s": x0_n_abs},
        "worst_N_derivative_match": _worst(
            x, eta, production_N_X, independent_N_X[-1]
        ),
        "worst_n_derivative_match": _worst(
            x, eta, production_n_X, independent_n_X[-1]
        ),
        "mutation_checks": mutations,
        "checks": checks,
        "passed": passed,
        "frozen_final_gates": dict(FROZEN_FINAL_GATES),
        "truth_boundary": truth,
        "limitations": [
            "The upstream stress uses a repository-autonomous admissible axis-pressure seed, not the source-prepared Appendix-A Pi_0.",
            "This audit is source-coordinate stress consistency only; it does not form Cartesian velocity or grad(p).",
            "Global/outer leading velocity, matched pressure, preregistered restricted forcing, and correction velocity remain absent.",
            "No complete NS residual or whole-domain volume-L2 is evaluated; no result is same-protocol comparable to ST006.",
            "The final normalized momentum 1e-3 and divergence 1e-5 gates remain unchanged and unevaluated here.",
        ],
    }
    payload["receipt_sha256"] = _sha256(payload)
    return payload


def save_real_artifact_audit(path: str | Path) -> dict[str, Any]:
    """Save/reload the exact upstream artifact, audit it, and persist the receipt."""
    from .kokuno_pa10_autonomous_reference_ns_x100 import (
        KokunoPA10AutonomousReferenceNsToX100,
    )

    original = KokunoPA10AutonomousReferenceNsToX100()
    with tempfile.TemporaryDirectory() as tmp:
        config_path = Path(tmp) / "autonomous_reference_ns.json"
        original.save_configuration(config_path)
        rebound = KokunoPA10AutonomousReferenceNsToX100.load_configuration(config_path)
    if rebound.semantic_sha256 != original.semantic_sha256:
        raise AssertionError("save/reload changed autonomous reference-stress semantic identity")

    payload = materialize_autonomous_reference_stress_audit(rebound)
    payload["save_reload_semantic_identity_exact"] = True
    payload["receipt_sha256"] = _sha256(
        {key: value for key, value in payload.items() if key != "receipt_sha256"}
    )
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = save_real_artifact_audit(args.output)
    print("passed=", payload["passed"])
    print("receipt_sha256=", payload["receipt_sha256"])
    print("worst_N_derivative_match=", payload["worst_N_derivative_match"])
    print("worst_n_derivative_match=", payload["worst_n_derivative_match"])
    if not payload["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    _main()
