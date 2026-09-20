"""Independent Agent-4 audit of the public PA.10 physical-center mapping.

This validator is downstream of the public
``KokunoPA10PhysicalCenterProfileContract`` but deliberately does not use its
partial-fraction ``zeta_primitive``, ``phi_star`` helper, analytic
``incompressibility_defect`` helper, or any Agent-1 receipt as a numerical
oracle.  The real-axis phi_* primitive is recomputed by composite Simpson
quadrature of the public zeta_* value path, and public first derivatives are
cross-checked with an independent FD4 implementation.

The audit is conditional on the upstream source-native axis/Phi0 public inputs.
It validates the new physical X/F/U/M/v0 assembly only.  It also records the
current configured-C normalization obstruction separately from formula
correctness: a sampled real-axis phi_* value above C already blocks the source
condition C >= sup_Omega |phi_*|, although this is still numerical evidence and
not a continuum complex-domain proof.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_physical_center_profile_contract import (
    KokunoPA10PhysicalCenterProfileContract,
)


SCHEMA = "kokuno-agent4-pa10-physical-center-independent-audit-v1"
FRESH_SEED = 9173451
FRESH_SAMPLE_COUNT = 2048
FD_SAMPLE_COUNT = 96
SIMPSON_LEVELS = (64, 128, 256)
FD_STEPS = (0.004, 0.002, 0.001)
MAPPING_RELATIVE_MAX_GATE = 2.0e-9
PHI_FINE_RELATIVE_MAX_GATE = 2.0e-9
DERIVATIVE_FINE_RELATIVE_RMS_GATE = 3.0e-6
DERIVATIVE_FINE_RELATIVE_MAX_GATE = 3.0e-5
INCOMPRESSIBILITY_FINE_RELATIVE_RMS_GATE = 3.0e-6
NORMALIZATION_NUMERICAL_MARGIN_MULTIPLIER = 50.0

_TRUTH_BOUNDARY = {
    "independent_agent4_physical_center_mapping_audit_executable": True,
    "independent_agent4_real_axis_phi_quadrature_executable": True,
    "independent_agent4_public_value_fd4_derivative_audit_executable": True,
    "source_native_physical_center_only": True,
    "source_complex_C_normalization_certified": False,
    "configured_C_real_axis_source_condition_admitted": False,
    "source_center_is_final_corrected_fixed_point": False,
    "fixed_point_correction_materialized": False,
    "global_leading_profile_reconstructed": False,
    "cartesian_spacetime_velocity_materialized": False,
    "matched_global_pressure_available": False,
    "restricted_forcing_composite_available": False,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "heldout_ns_momentum_residual_assessed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _relative_max(observed: Any, expected: Any) -> float:
    a = np.asarray(observed, dtype=float)
    b = np.asarray(expected, dtype=float)
    scale = np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
    return float(np.max(np.abs(a - b) / scale))


def _relative_rms(observed: Any, expected: Any) -> float:
    a = np.asarray(observed, dtype=float)
    b = np.asarray(expected, dtype=float)
    numerator = math.sqrt(float(np.mean(np.square(a - b))))
    denominator = max(math.sqrt(float(np.mean(np.square(b)))), 1.0)
    return numerator / denominator


def _simpson_zeta_primitive(
    profile: Any, eta: Any, intervals: int
) -> np.ndarray:
    """Independent real-axis integral int_0^eta zeta_*(w) dw."""
    n = int(intervals)
    if n <= 0 or n % 2:
        raise ValueError("Simpson intervals must be a positive even integer")
    eta_arr = np.asarray(eta, dtype=float)
    flat = eta_arr.reshape(-1)
    s = np.linspace(0.0, 1.0, n + 1, dtype=float)
    nodes = s[:, None] * flat[None, :]
    zeta = np.asarray(
        profile.axis_profiles.values(np.zeros_like(nodes), nodes)["zeta_star"],
        dtype=float,
    )
    weights = np.ones(n + 1, dtype=float)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    integral = flat * np.sum(weights[:, None] * zeta, axis=0) / (3.0 * n)
    return integral.reshape(eta_arr.shape)


def _independent_phi_star(profile: Any, eta: Any, intervals: int) -> np.ndarray:
    primitive = _simpson_zeta_primitive(profile, eta, intervals)
    return np.exp(float(profile.Lambda) * primitive)


def _independent_mapping(
    profile: Any, X: Any, eta: Any, *, intervals: int = SIMPSON_LEVELS[-1]
) -> dict[str, np.ndarray]:
    X_arr, eta_arr = np.broadcast_arrays(
        np.asarray(X, dtype=float), np.asarray(eta, dtype=float)
    )
    Y = float(profile.Lambda) * X_arr
    center = profile.center_profiles.values(Y, eta_arr)
    axis = profile.axis_profiles.values(Y, eta_arr)
    axis_d = profile.axis_profiles.derivatives(Y, eta_arr)
    phi_star = _independent_phi_star(profile, eta_arr, intervals)
    g = phi_star / float(profile.C)
    F = g * np.asarray(center["Phi_0"], dtype=float)
    U = np.asarray(axis["U_star"], dtype=float) + np.asarray(
        center["u_0"], dtype=float
    ) / float(profile.Lambda)
    slope = np.asarray(axis["u_0_Y"], dtype=float)
    slope_eta = np.asarray(axis_d["u_0_Y_eta"], dtype=float)
    M_over_X = np.asarray(axis["U_star"], dtype=float) + 0.5 * X_arr * slope
    M_eta_over_X = np.asarray(axis["U_star_eta"], dtype=float) + 0.5 * X_arr * slope_eta
    D = float(profile.axis_profiles.domain.D)
    numerator_over_X = (
        2.0 * eta_arr * U
        - 2.0 * D * eta_arr * M_over_X
        - np.asarray(axis["d"], dtype=float) * M_eta_over_X
    )
    v0 = numerator_over_X / np.asarray(axis["L"], dtype=float)
    return {
        "X": X_arr,
        "Y": Y,
        "phi_star": phi_star,
        "g": g,
        "F_0": F,
        "E_0": np.sqrt(2.0 * X_arr) * F,
        "U_0": U,
        "M_0": X_arr * M_over_X,
        "M_0_over_X": M_over_X,
        "v_0": v0,
        "V_0": X_arr * v0,
    }


def _fd4(profile: Any, X: np.ndarray, eta: np.ndarray, step: float, field: str, axis: str) -> np.ndarray:
    h = float(step)
    if axis == "X":
        def value(offset: float) -> np.ndarray:
            return np.asarray(profile.values(X + offset, eta)[field], dtype=float)
    elif axis == "eta":
        def value(offset: float) -> np.ndarray:
            return np.asarray(profile.values(X, eta + offset)[field], dtype=float)
    else:
        raise ValueError("axis must be X or eta")
    return (-value(2*h) + 8*value(h) - 8*value(-h) + value(-2*h)) / (12*h)


def audit_physical_center_profile(
    profile: Any | None = None,
    *,
    seed: int = FRESH_SEED,
    sample_count: int = FRESH_SAMPLE_COUNT,
) -> dict[str, Any]:
    if profile is None:
        profile = KokunoPA10PhysicalCenterProfileContract()
    rng = np.random.default_rng(int(seed))
    xmax = float(profile.source_X_interval[1])

    random_X = rng.uniform(0.0, xmax, int(sample_count))
    random_eta = rng.uniform(-0.95, 0.95, int(sample_count))
    probe_X = np.asarray([0.0, 1e-12, 1e-8, xmax, 0.0, xmax], dtype=float)
    probe_eta = np.asarray([-0.95, -1e-12, 0.0, 1e-12, 0.95, 0.5], dtype=float)
    X = np.concatenate([random_X, probe_X])
    eta = np.concatenate([random_eta, probe_eta])

    public = profile.values(X, eta)
    independent = _independent_mapping(profile, X, eta)
    mapping_fields = ("Y", "phi_star", "g", "F_0", "E_0", "U_0", "M_0", "M_0_over_X", "v_0", "V_0")
    mapping_errors = {
        name: _relative_max(public[name], independent[name]) for name in mapping_fields
    }
    max_mapping_error = max(mapping_errors.values())

    quad_eta = rng.uniform(-0.95, 0.95, 128)
    phi_levels = [
        _independent_phi_star(profile, quad_eta, n) for n in SIMPSON_LEVELS
    ]
    public_phi = np.asarray(profile.values(np.zeros_like(quad_eta), quad_eta)["phi_star"], dtype=float)
    phi_relative_max = [_relative_max(phi, public_phi) for phi in phi_levels]
    phi_successive_differences = [
        _relative_max(phi_levels[1], phi_levels[0]),
        _relative_max(phi_levels[2], phi_levels[1]),
    ]

    fd_X = rng.uniform(0.02, xmax - 0.02, FD_SAMPLE_COUNT)
    fd_eta = rng.uniform(-0.85, 0.85, FD_SAMPLE_COUNT)
    analytic = profile.derivatives(fd_X, fd_eta)
    derivative_pairs = (
        ("F_0", "F_0_X", "X"),
        ("U_0", "U_0_X", "X"),
        ("v_0", "v_0_X", "X"),
        ("F_0", "F_0_eta", "eta"),
        ("U_0", "U_0_eta", "eta"),
    )
    derivative_rms: dict[str, list[float]] = {}
    derivative_max: dict[str, list[float]] = {}
    for value_name, derivative_name, derivative_axis in derivative_pairs:
        key = f"{derivative_name}_from_{value_name}"
        derivative_rms[key] = []
        derivative_max[key] = []
        target = np.asarray(analytic[derivative_name], dtype=float)
        for step in FD_STEPS:
            numeric = _fd4(profile, fd_X, fd_eta, step, value_name, derivative_axis)
            derivative_rms[key].append(_relative_rms(numeric, target))
            derivative_max[key].append(_relative_max(numeric, target))
    max_fine_derivative_rms = max(values[-1] for values in derivative_rms.values())
    max_fine_derivative_max = max(values[-1] for values in derivative_max.values())

    incompressibility_rms: list[float] = []
    for step in FD_STEPS:
        V0_X = _fd4(profile, fd_X, fd_eta, step, "V_0", "X")
        U0_X = _fd4(profile, fd_X, fd_eta, step, "U_0", "X")
        U0_eta = _fd4(profile, fd_X, fd_eta, step, "U_0", "eta")
        values = profile.values(fd_X, fd_eta)
        axis = profile.axis_profiles.values(values["Y"], fd_eta)
        rhs = (
            2.0 * float(profile.axis_profiles.domain.A) * fd_eta * np.asarray(values["U_0"])
            - np.asarray(axis["d"]) * U0_eta
            + 2.0 * fd_eta * fd_X * U0_X
        ) / np.asarray(axis["L"])
        incompressibility_rms.append(_relative_rms(V0_X, rhs))

    axis_eta = np.asarray([-0.95, -1e-12, 0.0, 1e-12, 0.95])
    axis_values = profile.values(np.zeros_like(axis_eta), axis_eta)
    near_X = np.full_like(axis_eta, 1e-12)
    near_values = profile.values(near_X, axis_eta)
    axis_V0_zero = bool(np.all(np.asarray(axis_values["V_0"]) == 0.0))
    axis_E0_zero = bool(np.all(np.asarray(axis_values["E_0"]) == 0.0))
    near_E_scaling_error = _relative_max(
        np.asarray(near_values["E_0"]) / np.sqrt(2.0 * near_X),
        np.asarray(near_values["F_0"]),
    )

    norm_eta = np.linspace(-1.0, 1.0, 1025)
    phi_128 = _independent_phi_star(profile, norm_eta, 128)
    phi_256 = _independent_phi_star(profile, norm_eta, 256)
    normalization_refinement_delta = float(np.max(np.abs(phi_256 - phi_128)))
    real_phi_max = float(np.max(phi_256))
    normalization_margin = real_phi_max - float(profile.C)
    normalization_obstruction = bool(
        normalization_margin
        > NORMALIZATION_NUMERICAL_MARGIN_MULTIPLIER * normalization_refinement_delta
        + 1.0e-10
    )

    guards = {
        "fresh_physical_mapping_replay": max_mapping_error <= MAPPING_RELATIVE_MAX_GATE,
        "simpson_phi_finest_matches_public": phi_relative_max[-1] <= PHI_FINE_RELATIVE_MAX_GATE,
        "simpson_phi_refines": phi_successive_differences[1] <= phi_successive_differences[0] + 1e-15,
        "fd4_finest_relative_rms": max_fine_derivative_rms <= DERIVATIVE_FINE_RELATIVE_RMS_GATE,
        "fd4_finest_relative_max": max_fine_derivative_max <= DERIVATIVE_FINE_RELATIVE_MAX_GATE,
        "fd4_incompressibility_finest": incompressibility_rms[-1] <= INCOMPRESSIBILITY_FINE_RELATIVE_RMS_GATE,
        "fd4_incompressibility_refines": incompressibility_rms[-1] <= incompressibility_rms[0] + 1e-15,
        "axis_V0_exact_zero": axis_V0_zero,
        "axis_E0_exact_zero": axis_E0_zero,
        "near_axis_E_scaling": near_E_scaling_error <= 5e-12,
        "configured_C_real_axis_normalization_obstruction_detected": normalization_obstruction,
    }
    failed_guards = [name for name, passed in guards.items() if not passed]

    try:
        exact_head = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        exact_head = None

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "exact_head": exact_head,
        "seed": int(seed),
        "fresh_sample_count": int(sample_count),
        "fd_sample_count": FD_SAMPLE_COUNT,
        "simpson_levels": list(SIMPSON_LEVELS),
        "fd_steps": list(FD_STEPS),
        "gates": {
            "mapping_relative_max": MAPPING_RELATIVE_MAX_GATE,
            "phi_fine_relative_max": PHI_FINE_RELATIVE_MAX_GATE,
            "derivative_fine_relative_rms": DERIVATIVE_FINE_RELATIVE_RMS_GATE,
            "derivative_fine_relative_max": DERIVATIVE_FINE_RELATIVE_MAX_GATE,
            "incompressibility_fine_relative_rms": INCOMPRESSIBILITY_FINE_RELATIVE_RMS_GATE,
            "final_momentum_max_L2": 1.0e-3,
            "final_divergence_max_L2": 1.0e-5,
        },
        "mapping_relative_errors": mapping_errors,
        "maximum_mapping_relative_error": max_mapping_error,
        "phi_relative_max_by_simpson_level": dict(zip(map(str, SIMPSON_LEVELS), phi_relative_max)),
        "phi_successive_relative_differences": phi_successive_differences,
        "derivative_relative_rms_by_field": derivative_rms,
        "derivative_relative_max_by_field": derivative_max,
        "maximum_finest_derivative_relative_rms": max_fine_derivative_rms,
        "maximum_finest_derivative_relative_max": max_fine_derivative_max,
        "incompressibility_relative_rms_by_fd_level": dict(zip(map(str, FD_STEPS), incompressibility_rms)),
        "axis_V0_exact_zero": axis_V0_zero,
        "axis_E0_exact_zero": axis_E0_zero,
        "near_axis_E_scaling_relative_error": near_E_scaling_error,
        "configured_C": float(profile.C),
        "independent_real_axis_phi_star_sampled_max": real_phi_max,
        "normalization_refinement_delta": normalization_refinement_delta,
        "normalization_margin_phi_minus_C": normalization_margin,
        "configured_C_real_axis_normalization_obstruction_detected": normalization_obstruction,
        "mapping_formula_independently_consistent": all(
            guards[name]
            for name in guards
            if name != "configured_C_real_axis_normalization_obstruction_detected"
        ),
        "physical_center_profile_independently_admitted": False,
        "guards": guards,
        "failed_guards": failed_guards,
        "audit_passed": not failed_guards,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
    digest_payload = dict(report)
    report["receipt_sha256"] = hashlib.sha256(
        _canonical_json(digest_payload).encode("utf-8")
    ).hexdigest()
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=FRESH_SEED)
    parser.add_argument("--samples", type=int, default=FRESH_SAMPLE_COUNT)
    args = parser.parse_args()
    report = audit_physical_center_profile(seed=args.seed, sample_count=args.samples)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not report["audit_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
