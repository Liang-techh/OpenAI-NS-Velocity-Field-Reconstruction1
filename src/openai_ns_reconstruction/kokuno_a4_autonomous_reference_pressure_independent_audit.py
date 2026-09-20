"""Independent Agent-4 audit of the autonomous PA.10 reference-pressure artifact.

The upstream Agent-1 artifact exposes an autonomous repository seed

    Pi0_seed(eta) = -(5/2) P_*^2 (1 + eta^2)^-2

and composes it with the source-determined reference pressure increment

    C_p(X, eta) = Pi_r(X, eta) - Pi_0(eta)

to obtain the executable scalar

    Pi_ref_seed(X, eta) = Pi0_seed(eta) + C_p(X, eta).

This module audits only that public scalar surface.  Independent derivatives are
reconstructed from ``pressure(X, eta)`` by centered sixth-order finite
differences.  The production ``radial_derivative``, ``eta_derivative`` and
``axis_pressure_eta`` implementations are treated as quantities under audit and
are never used to construct the numerical reference.

This is intentionally not a matched Cartesian pressure-gradient validation.
The Agent-1 seed is repository-autonomous rather than the source-prepared
Appendix-A axis datum, and no global leading join / restricted forcing /
complete NS residual exists on this seam.
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


SCHEMA = "kokuno-a4-autonomous-reference-pressure-independent-audit-v1"
A1_AUTONOMOUS_PRESSURE_PR = 913
A1_AUTONOMOUS_PRESSURE_HEAD = "a5903ce594bf1c1f775d6fbba8a5eeca226622fa"
A1_AUTONOMOUS_PRESSURE_SOURCE_BLOB = "68dddc4f6ad1fb5bb6955498099b380240be6750"

FROZEN_FINAL_GATES = {
    "normalized_momentum_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
}


@dataclass(frozen=True)
class FrozenAutonomousPressureAuditProtocol:
    seed: int = 9173581
    random_offgrid_count: int = 128
    x_min: float = 0.08
    x_max: float = 99.0
    eta_min: float = -0.70
    eta_max: float = 0.70
    x_steps: tuple[float, float, float] = (4.0e-3, 2.0e-3, 1.0e-3)
    eta_steps: tuple[float, float, float] = (4.0e-4, 2.0e-4, 1.0e-4)
    radial_relative_rms_gate: float = 2.0e-4
    radial_relative_max_gate: float = 1.0e-3
    eta_relative_rms_gate: float = 2.0e-3
    eta_relative_max_gate: float = 1.0e-2
    axis_eta_relative_rms_gate: float = 2.0e-5
    axis_eta_relative_max_gate: float = 1.0e-4
    refinement_ratio_gate: float = 20.0
    refinement_floor: float = 1.0e-10
    axis_composition_abs_gate: float = 1.0e-12
    pressure_rms_floor: float = 1.0e-10
    radial_gradient_rms_floor: float = 1.0e-10


PROTOCOL = FrozenAutonomousPressureAuditProtocol()

TRUTH_BOUNDARY = {
    "autonomous_reference_pressure_independently_audited": False,
    "repository_autonomous_axis_pressure_seed_materialized": True,
    "source_prepared_appendixA_Pi0_materialized": False,
    "source_exact_absolute_reference_pressure_materialized": False,
    "matched_global_pressure_materialized": False,
    "cartesian_pressure_gradient_assessed": False,
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
        "agent1_pr": A1_AUTONOMOUS_PRESSURE_PR,
        "agent1_exact_head": A1_AUTONOMOUS_PRESSURE_HEAD,
        "agent1_source_blob": A1_AUTONOMOUS_PRESSURE_SOURCE_BLOB,
        "scientific_input_surface": "pressure(X,eta)",
        "production_surfaces_under_audit": [
            "radial_derivative(X,eta)",
            "eta_derivative(X,eta)",
            "axis_pressure_eta(eta)",
        ],
        "independent_operator": (
            "centered Cartesian/source-coordinate sixth-order finite difference "
            "from only public scalar pressure/axis_pressure evaluations"
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


def _require_contract(pressure: Any) -> None:
    required = (
        "pressure",
        "axis_pressure",
        "radial_derivative",
        "eta_derivative",
        "axis_pressure_eta",
        "semantic_sha256",
    )
    for name in required:
        if not hasattr(pressure, name):
            raise TypeError(f"pressure artifact is missing required public surface {name!r}")
    for name in required[:-1]:
        if not callable(getattr(pressure, name)):
            raise TypeError(f"{name} must be callable")


def _fd6_centered_first(
    fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    x: np.ndarray,
    eta: np.ndarray,
    h: float,
    *,
    axis: str,
) -> np.ndarray:
    """Centered sixth-order first derivative using only scalar evaluations."""
    if axis == "x":
        fm3 = _finite_array(fn(x - 3.0 * h, eta), "pressure")
        fm2 = _finite_array(fn(x - 2.0 * h, eta), "pressure")
        fm1 = _finite_array(fn(x - h, eta), "pressure")
        fp1 = _finite_array(fn(x + h, eta), "pressure")
        fp2 = _finite_array(fn(x + 2.0 * h, eta), "pressure")
        fp3 = _finite_array(fn(x + 3.0 * h, eta), "pressure")
    elif axis == "eta":
        fm3 = _finite_array(fn(x, eta - 3.0 * h), "pressure")
        fm2 = _finite_array(fn(x, eta - 2.0 * h), "pressure")
        fm1 = _finite_array(fn(x, eta - h), "pressure")
        fp1 = _finite_array(fn(x, eta + h), "pressure")
        fp2 = _finite_array(fn(x, eta + 2.0 * h), "pressure")
        fp3 = _finite_array(fn(x, eta + 3.0 * h), "pressure")
    else:
        raise ValueError(f"unsupported axis {axis!r}")
    return (-fm3 + 9.0 * fm2 - 45.0 * fm1 + 45.0 * fp1 - 9.0 * fp2 + fp3) / (
        60.0 * h
    )


def _fd6_axis_eta(
    fn: Callable[[np.ndarray], np.ndarray],
    eta: np.ndarray,
    h: float,
) -> np.ndarray:
    fm3 = _finite_array(fn(eta - 3.0 * h), "axis_pressure")
    fm2 = _finite_array(fn(eta - 2.0 * h), "axis_pressure")
    fm1 = _finite_array(fn(eta - h), "axis_pressure")
    fp1 = _finite_array(fn(eta + h), "axis_pressure")
    fp2 = _finite_array(fn(eta + 2.0 * h), "axis_pressure")
    fp3 = _finite_array(fn(eta + 3.0 * h), "axis_pressure")
    return (-fm3 + 9.0 * fm2 - 45.0 * fm1 + 45.0 * fp1 - 9.0 * fp2 + fp3) / (
        60.0 * h
    )


def _rms(value: np.ndarray) -> float:
    arr = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(arr * arr)))


def _match_metrics(production: np.ndarray, independent: np.ndarray) -> dict[str, float]:
    production = np.asarray(production, dtype=float)
    independent = np.asarray(independent, dtype=float)
    error = independent - production
    scale_rms = max(_rms(production), 1.0e-30)
    scale_max = max(float(np.max(np.abs(production))), 1.0e-30)
    return {
        "absolute_rms": _rms(error),
        "absolute_max": float(np.max(np.abs(error))),
        "relative_rms": _rms(error) / scale_rms,
        "relative_sampled_max": float(np.max(np.abs(error))) / scale_max,
    }


def _stability(levels: list[np.ndarray]) -> dict[str, float | bool]:
    if len(levels) != 3:
        raise ValueError("exactly three resolution levels are required")
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
    idx = int(np.argmax(err))
    return {
        "X": float(np.asarray(x)[idx]),
        "eta": float(np.asarray(eta)[idx]),
        "production": float(np.asarray(production)[idx]),
        "independent": float(np.asarray(independent)[idx]),
        "absolute_error": float(err[idx]),
    }


def _sample_points() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(PROTOCOL.seed)
    x_random = rng.uniform(PROTOCOL.x_min, PROTOCOL.x_max, PROTOCOL.random_offgrid_count)
    eta_random = rng.uniform(
        PROTOCOL.eta_min, PROTOCOL.eta_max, PROTOCOL.random_offgrid_count
    )
    x_near = np.asarray([0.02, 0.04, 0.06], dtype=float)
    eta_near = np.asarray([-0.60, 0.0, 0.60], dtype=float)
    axis_eta = np.asarray([-0.70, -0.35, 0.0, 0.35, 0.70], dtype=float)
    return (
        np.concatenate([x_random, x_near]),
        np.concatenate([eta_random, eta_near]),
        x_near,
        eta_near,
        axis_eta,
    )


def _mutation_checks(
    pressure: Any,
    x: np.ndarray,
    eta: np.ndarray,
    axis_eta: np.ndarray,
    production_x: np.ndarray,
    production_eta: np.ndarray,
    independent_x_fine: np.ndarray,
    independent_eta_fine: np.ndarray,
) -> dict[str, bool]:
    radial_scaled = _match_metrics(0.99 * production_x, independent_x_fine)
    radial_scale_detected = not (
        radial_scaled["relative_rms"] <= PROTOCOL.radial_relative_rms_gate
        and radial_scaled["relative_sampled_max"] <= PROTOCOL.radial_relative_max_gate
    )

    eta_scaled = _match_metrics(0.90 * production_eta, independent_eta_fine)
    eta_scale_detected = not (
        eta_scaled["relative_rms"] <= PROTOCOL.eta_relative_rms_gate
        and eta_scaled["relative_sampled_max"] <= PROTOCOL.eta_relative_max_gate
    )

    def plus_x(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return _finite_array(pressure.pressure(xx, ee), "pressure") + 1.0e-3 * np.asarray(
            xx, dtype=float
        )

    mutated_x = _fd6_centered_first(
        plus_x, x, eta, PROTOCOL.x_steps[-1], axis="x"
    )
    pressure_plus_x_detected = _rms(mutated_x - production_x) >= 5.0e-4

    def plus_eta(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return _finite_array(pressure.pressure(xx, ee), "pressure") + 1.0e-3 * np.asarray(
            ee, dtype=float
        )

    mutated_eta = _fd6_centered_first(
        plus_eta, x, eta, PROTOCOL.eta_steps[-1], axis="eta"
    )
    pressure_plus_eta_detected = _rms(mutated_eta - production_eta) >= 5.0e-4

    shifted_axis_value = _finite_array(
        pressure.pressure(np.zeros_like(axis_eta), axis_eta), "axis pressure"
    ) + 1.0e-3
    axis_seed_value = _finite_array(pressure.axis_pressure(axis_eta), "axis_pressure")
    constant_axis_shift_detected = (
        float(np.max(np.abs(shifted_axis_value - axis_seed_value)))
        > PROTOCOL.axis_composition_abs_gate
    )

    return {
        "production_radial_derivative_times_0p99_detected": bool(radial_scale_detected),
        "production_eta_derivative_times_0p90_detected": bool(eta_scale_detected),
        "pressure_plus_1e_minus3_X_detected": bool(pressure_plus_x_detected),
        "pressure_plus_1e_minus3_eta_detected": bool(pressure_plus_eta_detected),
        "axis_scalar_plus_1e_minus3_detected": bool(constant_axis_shift_detected),
    }


def materialize_autonomous_reference_pressure_audit(pressure: Any) -> dict[str, Any]:
    """Run the frozen implementation-distinct audit on one public pressure artifact."""
    _require_contract(pressure)
    x, eta, x_near, eta_near, axis_eta = _sample_points()

    pressure_value = _finite_array(pressure.pressure(x, eta), "pressure")
    production_x = _finite_array(pressure.radial_derivative(x, eta), "radial_derivative")
    production_eta = _finite_array(pressure.eta_derivative(x, eta), "eta_derivative")
    production_axis_eta = _finite_array(
        pressure.axis_pressure_eta(axis_eta), "axis_pressure_eta"
    )

    independent_x = [
        _fd6_centered_first(pressure.pressure, x, eta, h, axis="x")
        for h in PROTOCOL.x_steps
    ]
    independent_eta = [
        _fd6_centered_first(pressure.pressure, x, eta, h, axis="eta")
        for h in PROTOCOL.eta_steps
    ]
    independent_axis_eta = [
        _fd6_axis_eta(pressure.axis_pressure, axis_eta, h) for h in PROTOCOL.eta_steps
    ]

    radial_metrics = [_match_metrics(production_x, level) for level in independent_x]
    eta_metrics = [_match_metrics(production_eta, level) for level in independent_eta]
    axis_eta_metrics = [
        _match_metrics(production_axis_eta, level) for level in independent_axis_eta
    ]
    radial_stability = _stability(independent_x)
    eta_stability = _stability(independent_eta)
    axis_eta_stability = _stability(independent_axis_eta)

    axis_pressure_from_total = _finite_array(
        pressure.pressure(np.zeros_like(axis_eta), axis_eta), "pressure at X=0"
    )
    axis_pressure_seed = _finite_array(pressure.axis_pressure(axis_eta), "axis_pressure")
    axis_composition_max_abs = float(
        np.max(np.abs(axis_pressure_from_total - axis_pressure_seed))
    )

    near_values = _finite_array(
        pressure.pressure(x_near, eta_near), "near-axis pressure"
    )
    x100_values = _finite_array(
        pressure.pressure(np.full_like(axis_eta, 100.0), axis_eta), "X100 pressure"
    )

    fine_radial = radial_metrics[-1]
    fine_eta = eta_metrics[-1]
    fine_axis_eta = axis_eta_metrics[-1]

    mutations = _mutation_checks(
        pressure,
        x,
        eta,
        axis_eta,
        production_x,
        production_eta,
        independent_x[-1],
        independent_eta[-1],
    )

    checks = {
        "radial_relative_rms": fine_radial["relative_rms"]
        <= PROTOCOL.radial_relative_rms_gate,
        "radial_relative_sampled_max": fine_radial["relative_sampled_max"]
        <= PROTOCOL.radial_relative_max_gate,
        "eta_relative_rms": fine_eta["relative_rms"] <= PROTOCOL.eta_relative_rms_gate,
        "eta_relative_sampled_max": fine_eta["relative_sampled_max"]
        <= PROTOCOL.eta_relative_max_gate,
        "axis_eta_relative_rms": fine_axis_eta["relative_rms"]
        <= PROTOCOL.axis_eta_relative_rms_gate,
        "axis_eta_relative_sampled_max": fine_axis_eta["relative_sampled_max"]
        <= PROTOCOL.axis_eta_relative_max_gate,
        "radial_resolution_stability": bool(radial_stability["passed"]),
        "eta_resolution_stability": bool(eta_stability["passed"]),
        "axis_eta_resolution_stability": bool(axis_eta_stability["passed"]),
        "exact_axis_composition": axis_composition_max_abs
        <= PROTOCOL.axis_composition_abs_gate,
        "pressure_nontrivial": _rms(pressure_value) >= PROTOCOL.pressure_rms_floor,
        "radial_gradient_nontrivial": _rms(production_x)
        >= PROTOCOL.radial_gradient_rms_floor,
        "near_axis_finite": bool(np.all(np.isfinite(near_values))),
        "X100_finite": bool(np.all(np.isfinite(x100_values))),
        **mutations,
    }
    passed = bool(all(checks.values()))

    truth = dict(TRUTH_BOUNDARY)
    truth["autonomous_reference_pressure_independently_audited"] = passed

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "agent1_pr": A1_AUTONOMOUS_PRESSURE_PR,
        "agent1_exact_head": A1_AUTONOMOUS_PRESSURE_HEAD,
        "agent1_source_blob": A1_AUTONOMOUS_PRESSURE_SOURCE_BLOB,
        "pressure_semantic_sha256": str(pressure.semantic_sha256),
        "protocol": asdict(PROTOCOL),
        "sample_counts": {
            "fresh_random_offgrid": PROTOCOL.random_offgrid_count,
            "axis_near": len(x_near),
            "exact_axis": len(axis_eta),
        },
        "production_rms": {
            "pressure": _rms(pressure_value),
            "radial_derivative": _rms(production_x),
            "eta_derivative": _rms(production_eta),
            "axis_pressure_eta": _rms(production_axis_eta),
        },
        "radial_three_resolution_match": radial_metrics,
        "eta_three_resolution_match": eta_metrics,
        "axis_eta_three_resolution_match": axis_eta_metrics,
        "radial_resolution_stability": radial_stability,
        "eta_resolution_stability": eta_stability,
        "axis_eta_resolution_stability": axis_eta_stability,
        "exact_axis_composition_max_abs": axis_composition_max_abs,
        "X100_pressure_range": [
            float(np.min(x100_values)),
            float(np.max(x100_values)),
        ],
        "worst_radial_match": _worst(
            x, eta, production_x, independent_x[-1]
        ),
        "worst_eta_match": _worst(
            x, eta, production_eta, independent_eta[-1]
        ),
        "mutation_checks": mutations,
        "checks": checks,
        "passed": passed,
        "frozen_final_gates": dict(FROZEN_FINAL_GATES),
        "truth_boundary": truth,
        "limitations": [
            "Pi0_seed is repository-autonomous; it is not the source-prepared Appendix-A Pi_0.",
            "The audit is in source (X,eta) coordinates and does not form Cartesian grad(p).",
            "Matched/global pressure, restricted forcing, global leading velocity and correction velocity are absent.",
            "No complete NS residual is evaluated and no result is same-protocol comparable to ST006.",
            "The final momentum 1e-3 and divergence 1e-5 gates remain unchanged and unevaluated here.",
        ],
    }
    payload["receipt_sha256"] = _sha256(payload)
    return payload


def save_real_artifact_audit(path: str | Path) -> dict[str, Any]:
    """Save/reload the exact upstream artifact, run the frozen audit, and persist it."""
    from .kokuno_pa10_autonomous_axis_pressure_seed_x100 import (
        KokunoPA10AutonomousAxisPressureSeedToX100,
    )

    original = KokunoPA10AutonomousAxisPressureSeedToX100()
    with tempfile.TemporaryDirectory() as tmp:
        config_path = Path(tmp) / "autonomous_reference_pressure.json"
        original.save_configuration(config_path)
        rebound = KokunoPA10AutonomousAxisPressureSeedToX100.load_configuration(config_path)
    if rebound.semantic_sha256 != original.semantic_sha256:
        raise AssertionError("save/reload changed autonomous reference-pressure semantic identity")

    payload = materialize_autonomous_reference_pressure_audit(rebound)
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
    print("worst_radial_match=", payload["worst_radial_match"])
    print("worst_eta_match=", payload["worst_eta_match"])
    if not payload["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    _main()
