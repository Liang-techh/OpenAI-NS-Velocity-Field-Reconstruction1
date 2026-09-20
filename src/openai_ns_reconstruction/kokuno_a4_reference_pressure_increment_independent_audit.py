"""Independent Agent-4 audit of the PA.10 reference pressure increment.

This module audits the new Agent-1 public surface

    C_p(X, eta) = Pi_r(X, eta) - Pi_0(eta)

without using Agent-1's production derivative implementations as the numerical
reference.  The reference derivatives are reconstructed from only the public
``pressure_increment(X, eta)`` callable using centered two-point differences
with Richardson extrapolation.

This is deliberately a scoped precursor to a future pressure-gradient audit.
``Pi_0(eta)`` is still missing, so ``C_p`` is not absolute/matched pressure and
this module cannot form the complete Cartesian NS residual.
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


SCHEMA = "kokuno-a4-reference-pressure-increment-independent-audit-v1"
A1_REFERENCE_PRESSURE_PR = 907
A1_REFERENCE_PRESSURE_HEAD = "3f43931cec1fc762678a5a8134adfceaf7edfa81"
A1_REFERENCE_PRESSURE_SOURCE_BLOB = "4704fed4271da2c6070110a6935d2f33f7932f56"

FROZEN_FINAL_GATES = {
    "normalized_momentum_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
}


@dataclass(frozen=True)
class FrozenPressureAuditProtocol:
    seed: int = 9173571
    random_offgrid_count: int = 128
    x_min: float = 0.08
    x_max: float = 99.0
    eta_min: float = -0.70
    eta_max: float = 0.70
    x_steps: tuple[float, float, float] = (8.0e-3, 4.0e-3, 2.0e-3)
    eta_steps: tuple[float, float, float] = (8.0e-4, 4.0e-4, 2.0e-4)
    axis_increment_abs_gate: float = 1.0e-12
    radial_relative_rms_gate: float = 2.0e-4
    radial_relative_max_gate: float = 1.0e-3
    eta_relative_rms_gate: float = 2.0e-3
    eta_relative_max_gate: float = 1.0e-2
    refinement_ratio_gate: float = 6.0
    refinement_floor: float = 1.0e-10
    pressure_rms_floor: float = 1.0e-10
    radial_gradient_rms_floor: float = 1.0e-10


PROTOCOL = FrozenPressureAuditProtocol()

TRUTH_BOUNDARY = {
    "reference_pressure_increment_independently_audited": False,
    "absolute_axis_pressure_Pi0_materialized": False,
    "absolute_reference_pressure_materialized": False,
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
        "agent1_pr": A1_REFERENCE_PRESSURE_PR,
        "agent1_exact_head": A1_REFERENCE_PRESSURE_HEAD,
        "agent1_source_blob": A1_REFERENCE_PRESSURE_SOURCE_BLOB,
        "scientific_input_surface": "pressure_increment(X,eta)",
        "production_surfaces_under_audit": ["radial_derivative(X,eta)", "eta_derivative(X,eta)"],
        "independent_operator": "centered-two-point + Richardson extrapolation D_R=(4D(h/2)-D(h))/3",
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
    for name in ("pressure_increment", "radial_derivative", "eta_derivative", "semantic_sha256"):
        if not hasattr(pressure, name):
            raise TypeError(f"pressure artifact is missing required public surface {name!r}")
    if not callable(pressure.pressure_increment):
        raise TypeError("pressure_increment must be callable")
    if not callable(pressure.radial_derivative) or not callable(pressure.eta_derivative):
        raise TypeError("production derivative surfaces must be callable")


def _centered_first(
    fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    x: np.ndarray,
    eta: np.ndarray,
    h: float,
    *,
    axis: str,
) -> np.ndarray:
    if axis == "x":
        return (_finite_array(fn(x + h, eta), "pressure_increment") - _finite_array(fn(x - h, eta), "pressure_increment")) / (2.0 * h)
    if axis == "eta":
        return (_finite_array(fn(x, eta + h), "pressure_increment") - _finite_array(fn(x, eta - h), "pressure_increment")) / (2.0 * h)
    raise ValueError(f"unsupported axis {axis!r}")


def _richardson_first(
    fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    x: np.ndarray,
    eta: np.ndarray,
    h: float,
    *,
    axis: str,
) -> np.ndarray:
    coarse = _centered_first(fn, x, eta, h, axis=axis)
    half = _centered_first(fn, x, eta, h / 2.0, axis=axis)
    return (4.0 * half - coarse) / 3.0


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
    coarse_medium = _rms(levels[0] - levels[1])
    medium_fine = _rms(levels[1] - levels[2])
    ratio = math.inf if medium_fine == 0.0 else coarse_medium / medium_fine
    floor_hit = medium_fine <= PROTOCOL.refinement_floor
    passed = bool(floor_hit or ratio >= PROTOCOL.refinement_ratio_gate)
    return {
        "coarse_to_medium_rms": coarse_medium,
        "medium_to_fine_rms": medium_fine,
        "ratio": ratio,
        "numerical_floor_hit": floor_hit,
        "passed": passed,
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
        "X": float(x[idx]),
        "eta": float(eta[idx]),
        "production": float(np.asarray(production)[idx]),
        "independent": float(np.asarray(independent)[idx]),
        "absolute_error": float(err[idx]),
    }


def _sample_points() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(PROTOCOL.seed)
    x = rng.uniform(PROTOCOL.x_min, PROTOCOL.x_max, PROTOCOL.random_offgrid_count)
    eta = rng.uniform(PROTOCOL.eta_min, PROTOCOL.eta_max, PROTOCOL.random_offgrid_count)
    # Explicit near-axis probes are appended but remain safely inside centered X stencils.
    x_near = np.asarray([0.02, 0.04, 0.06], dtype=float)
    eta_near = np.asarray([-0.60, 0.0, 0.60], dtype=float)
    return np.concatenate([x, x_near]), np.concatenate([eta, eta_near]), x_near, eta_near


def _mutation_checks(
    pressure: Any,
    x: np.ndarray,
    eta: np.ndarray,
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

    def mutated_increment(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return _finite_array(pressure.pressure_increment(xx, ee), "pressure_increment") + 1.0e-3 * np.asarray(xx, dtype=float)

    mutated_x = _richardson_first(mutated_increment, x, eta, PROTOCOL.x_steps[-1], axis="x")
    linear_x_detected = _rms(mutated_x - production_x) >= 5.0e-4
    return {
        "production_radial_derivative_times_0p99_detected": bool(radial_scale_detected),
        "production_eta_derivative_times_0p90_detected": bool(eta_scale_detected),
        "pressure_increment_plus_1e_minus3_X_detected": bool(linear_x_detected),
    }


def materialize_reference_pressure_increment_audit(pressure: Any) -> dict[str, Any]:
    """Run the frozen independent audit on a pressure-increment artifact.

    There is intentionally no caller-supplied protocol, threshold, seed, or
    resolution knob on this scientific entry point.
    """
    _require_contract(pressure)
    x, eta, x_near, eta_near = _sample_points()
    pressure_value = _finite_array(pressure.pressure_increment(x, eta), "pressure_increment")
    production_x = _finite_array(pressure.radial_derivative(x, eta), "radial_derivative")
    production_eta = _finite_array(pressure.eta_derivative(x, eta), "eta_derivative")

    independent_x = [
        _richardson_first(pressure.pressure_increment, x, eta, h, axis="x")
        for h in PROTOCOL.x_steps
    ]
    independent_eta = [
        _richardson_first(pressure.pressure_increment, x, eta, h, axis="eta")
        for h in PROTOCOL.eta_steps
    ]
    radial_metrics = [_match_metrics(production_x, value) for value in independent_x]
    eta_metrics = [_match_metrics(production_eta, value) for value in independent_eta]
    radial_stability = _stability(independent_x)
    eta_stability = _stability(independent_eta)

    axis_eta = np.asarray([-0.70, -0.35, 0.0, 0.35, 0.70], dtype=float)
    axis_values = _finite_array(
        pressure.pressure_increment(np.zeros_like(axis_eta), axis_eta), "axis pressure_increment"
    )
    boundary_values = _finite_array(
        pressure.pressure_increment(np.full_like(axis_eta, 100.0), axis_eta), "X100 pressure_increment"
    )
    near_values = _finite_array(pressure.pressure_increment(x_near, eta_near), "near-axis pressure_increment")

    fine_radial = radial_metrics[-1]
    fine_eta = eta_metrics[-1]
    mutations = _mutation_checks(
        pressure,
        x,
        eta,
        production_x,
        production_eta,
        independent_x[-1],
        independent_eta[-1],
    )

    checks = {
        "axis_increment": float(np.max(np.abs(axis_values))) <= PROTOCOL.axis_increment_abs_gate,
        "radial_relative_rms": fine_radial["relative_rms"] <= PROTOCOL.radial_relative_rms_gate,
        "radial_relative_sampled_max": fine_radial["relative_sampled_max"] <= PROTOCOL.radial_relative_max_gate,
        "eta_relative_rms": fine_eta["relative_rms"] <= PROTOCOL.eta_relative_rms_gate,
        "eta_relative_sampled_max": fine_eta["relative_sampled_max"] <= PROTOCOL.eta_relative_max_gate,
        "radial_resolution_stability": bool(radial_stability["passed"]),
        "eta_resolution_stability": bool(eta_stability["passed"]),
        "pressure_nontrivial": _rms(pressure_value) >= PROTOCOL.pressure_rms_floor,
        "radial_gradient_nontrivial": _rms(production_x) >= PROTOCOL.radial_gradient_rms_floor,
        "near_axis_finite": bool(np.all(np.isfinite(near_values))),
        "X100_finite": bool(np.all(np.isfinite(boundary_values))),
        **mutations,
    }
    passed = bool(all(checks.values()))

    truth = dict(TRUTH_BOUNDARY)
    truth["reference_pressure_increment_independently_audited"] = passed
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "agent1_exact_head": A1_REFERENCE_PRESSURE_HEAD,
        "agent1_source_blob": A1_REFERENCE_PRESSURE_SOURCE_BLOB,
        "pressure_semantic_sha256": str(pressure.semantic_sha256),
        "protocol": asdict(PROTOCOL),
        "sample_counts": {
            "fresh_random_offgrid": PROTOCOL.random_offgrid_count,
            "axis_near": len(x_near),
            "exact_axis": len(axis_eta),
        },
        "production_rms": {
            "pressure_increment": _rms(pressure_value),
            "radial_derivative": _rms(production_x),
            "eta_derivative": _rms(production_eta),
        },
        "radial_three_resolution_match": radial_metrics,
        "eta_three_resolution_match": eta_metrics,
        "radial_resolution_stability": radial_stability,
        "eta_resolution_stability": eta_stability,
        "axis_increment_max_abs": float(np.max(np.abs(axis_values))),
        "X100_pressure_increment_range": [float(np.min(boundary_values)), float(np.max(boundary_values))],
        "worst_radial_match": _worst(x, eta, production_x, independent_x[-1]),
        "worst_eta_match": _worst(x, eta, production_eta, independent_eta[-1]),
        "mutation_checks": mutations,
        "checks": checks,
        "passed": passed,
        "frozen_final_gates": dict(FROZEN_FINAL_GATES),
        "truth_boundary": truth,
        "limitations": [
            "C_p is Pi_r-Pi_0, not absolute pressure; Pi_0(eta) is not materialized.",
            "The audit is in source (X,eta) coordinates and does not form Cartesian grad(p).",
            "No velocity, forcing, correction, complete NS residual, or canonical 24/48/96 volume-L2 is assessed here.",
            "Passing this scoped audit cannot set pde_validated=true.",
        ],
    }
    payload["receipt_sha256"] = _sha256(payload)
    return payload


def save_audit_receipt(path: str | Path, payload: Mapping[str, Any]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(dict(payload), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def _main() -> None:
    import argparse

    from .kokuno_pa10_reference_pressure_increment_x100 import (
        KokunoPA10ReferencePressureIncrementToX100,
    )

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    # Scientific path audits a serialized/rebound artifact, not an in-memory
    # construction object with hidden state.
    original = KokunoPA10ReferencePressureIncrementToX100()
    with tempfile.TemporaryDirectory() as tmp:
        config_path = Path(tmp) / "reference-pressure-increment.json"
        original.save_configuration(config_path)
        rebound = KokunoPA10ReferencePressureIncrementToX100.load_configuration(config_path)
        if rebound.semantic_sha256 != original.semantic_sha256:
            raise RuntimeError("save/reload changed pressure semantic identity")
        payload = materialize_reference_pressure_increment_audit(rebound)

    save_audit_receipt(args.output, payload)
    print(json.dumps({"passed": payload["passed"], "receipt_sha256": payload["receipt_sha256"]}, sort_keys=True))
    if not payload["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    _main()
