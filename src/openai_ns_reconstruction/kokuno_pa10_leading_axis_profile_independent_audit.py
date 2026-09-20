"""Independent Agent-4 audit of the public PA.10 source-native leading-axis profile.

This validator is intentionally downstream of the candidate-facing public API in
``kokuno_pa10_leading_axis_profile_contract``. It does not call the Agent-1
``axis_state`` implementation or any private Fraction/product/rounding helper.
Instead it reconstructs the exported formulas from the serialized scalar
configuration and compares them with fresh public API evaluations.

The audit also differentiates the public value path with an independent centered
FD4 operator on three fixed step sizes. This is a profile-contract audit only:
it is not a Cartesian velocity, pressure/forcing composite, or Navier-Stokes
residual validation.
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

from .kokuno_pa10_leading_axis_profile_contract import (
    KokunoPA10LeadingAxisProfileContract,
)


SCHEMA = "kokuno-agent4-pa10-leading-axis-profile-independent-audit-v1"
FRESH_SEED = 9173441
FRESH_SAMPLE_COUNT = 8192
FD_SAMPLE_COUNT = 96
FD_STEPS = (0.004, 0.002, 0.001)
VALUE_RELATIVE_MAX_GATE = 5.0e-12
DERIVATIVE_RELATIVE_MAX_GATE = 5.0e-11
FD_FINE_RELATIVE_RMS_GATE = 2.0e-8
FD_REFINEMENT_RATIO_GATE = 12.0
AXIS_SCALING_RELATIVE_GATE = 5.0e-12

_TRUTH_BOUNDARY = {
    "independent_agent4_profile_formula_audit_executable": True,
    "independent_agent4_public_value_path_fd4_derivative_audit_executable": True,
    "source_native_profile_only": True,
    "repository_X_identified_with_source_Y": False,
    "global_leading_profile_reconstructed": False,
    "complete_kokuno_composite_velocity": False,
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
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1.0)))


def _relative_rms(observed: Any, expected: Any) -> float:
    a = np.asarray(observed, dtype=float)
    b = np.asarray(expected, dtype=float)
    numerator = math.sqrt(float(np.mean(np.square(a - b))))
    denominator = max(math.sqrt(float(np.mean(np.square(b)))), 1.0)
    return numerator / denominator


def _independent_parameters(configuration: Mapping[str, Any]) -> dict[str, float]:
    outer = dict(configuration["outer_schedule"])
    pressure = dict(configuration["pressure_binding"])
    axis = dict(configuration["source_axis_domain"])

    h = float(outer["h"])
    M_d = float(outer["M_d"])
    log_p_star_margin = float(outer["log_p_star_margin"])
    pressure_margin = float(pressure["pressure_margin"])
    j0 = float(pressure["j0"])
    sigma_star = float(axis["sigma_star"])
    enlarged_real_margin = float(axis["enlarged_real_margin"])

    T_d = math.exp(M_d) + 10.0
    log_P_star = T_d + log_p_star_margin
    P_star = math.exp(log_P_star)
    required_pressure_scale = math.sqrt(2.5) * P_star
    selected_pressure_scale = math.nextafter(
        pressure_margin * required_pressure_scale, math.inf
    )

    return {
        "h": h,
        "A": 0.5 + h,
        "D": 0.5 - h,
        "j0": j0,
        "sigma_star": sigma_star,
        "pressure_square": selected_pressure_scale * selected_pressure_scale,
        "enlarged_real_margin": enlarged_real_margin,
    }


def _independent_values(
    Y: Any, eta: Any, parameters: Mapping[str, float]
) -> dict[str, np.ndarray]:
    Y_arr, eta_arr = np.broadcast_arrays(
        np.asarray(Y, dtype=float), np.asarray(eta, dtype=float)
    )
    h = float(parameters["h"])
    A = float(parameters["A"])
    D = float(parameters["D"])
    j0 = float(parameters["j0"])
    sigma2 = float(parameters["sigma_star"]) ** 2
    p2 = float(parameters["pressure_square"])

    d = 1.0 - eta_arr * eta_arr
    L = 1.0 - 2.0 * h * eta_arr * eta_arr
    U = 4.0 * eta_arr + j0
    H = D * eta_arr + d * U
    W = 1.0 - 4.0 * d - 2.0 * D * eta_arr * U
    den = 1.0 + eta_arr * eta_arr
    Pi = -p2 / den**2
    Pi_eta = 4.0 * p2 * eta_arr / den**3
    B = 1.0 - 2.0 * eta_arr * U
    Z = -A * B * U - 4.0 * H - d * Pi_eta + 4.0 * A * eta_arr * Pi
    Q = H * H + sigma2
    chi = H * H / Q
    zeta = -L * H / Q
    slope = -Z / (2.0 * L)

    return {
        "d": d,
        "L": L,
        "U_star": U,
        "H_star": H,
        "W_star": W,
        "Z_star": Z,
        "Pi_0": Pi,
        "Pi_0_eta": Pi_eta,
        "chi": chi,
        "zeta_star": zeta,
        "u_0": Y_arr * slope,
        "u_0_Y": slope,
    }


def _independent_derivatives(
    Y: Any, eta: Any, parameters: Mapping[str, float]
) -> dict[str, np.ndarray]:
    Y_arr, eta_arr = np.broadcast_arrays(
        np.asarray(Y, dtype=float), np.asarray(eta, dtype=float)
    )
    h = float(parameters["h"])
    A = float(parameters["A"])
    D = float(parameters["D"])
    j0 = float(parameters["j0"])
    sigma2 = float(parameters["sigma_star"]) ** 2
    p2 = float(parameters["pressure_square"])

    d = 1.0 - eta_arr * eta_arr
    L = 1.0 - 2.0 * h * eta_arr * eta_arr
    U = 4.0 * eta_arr + j0
    H = D * eta_arr + d * U
    den = 1.0 + eta_arr * eta_arr
    Pi = -p2 / den**2
    Pi_eta = 4.0 * p2 * eta_arr / den**3
    B = 1.0 - 2.0 * eta_arr * U
    Z = -A * B * U - 4.0 * H - d * Pi_eta + 4.0 * A * eta_arr * Pi

    d_eta = -2.0 * eta_arr
    L_eta = -4.0 * h * eta_arr
    U_eta = np.full_like(eta_arr, 4.0)
    H_eta = D - 2.0 * eta_arr * U + 4.0 * d
    W_eta = -2.0 * D * j0 + 16.0 * h * eta_arr
    Pi_eta_eta = 4.0 * p2 * (1.0 - 5.0 * eta_arr * eta_arr) / den**4
    B_eta = -2.0 * U - 8.0 * eta_arr
    Z_eta = (
        -A * (B_eta * U + 4.0 * B)
        - 4.0 * H_eta
        - d_eta * Pi_eta
        - d * Pi_eta_eta
        + 4.0 * A * (Pi + eta_arr * Pi_eta)
    )

    Q = H * H + sigma2
    chi_eta = 2.0 * H * H_eta * sigma2 / Q**2
    zeta_eta = (
        -(L_eta * H + L * H_eta) * Q + 2.0 * L * H * H * H_eta
    ) / Q**2
    slope_eta = -(Z_eta * L - Z * L_eta) / (2.0 * L * L)

    return {
        "d_eta": d_eta,
        "L_eta": L_eta,
        "U_star_eta": U_eta,
        "H_star_eta": H_eta,
        "W_star_eta": W_eta,
        "Z_star_eta": Z_eta,
        "Pi_0_eta": Pi_eta,
        "Pi_0_eta_eta": Pi_eta_eta,
        "chi_eta": chi_eta,
        "zeta_star_eta": zeta_eta,
        "u_0_Y": -Z / (2.0 * L),
        "u_0_eta": Y_arr * slope_eta,
        "u_0_Y_eta": slope_eta,
    }


def _fd4_value_derivative(
    contract: Any, Y: np.ndarray, eta: np.ndarray, step: float, field: str
) -> np.ndarray:
    h = float(step)
    return (
        -np.asarray(contract.values(Y, eta + 2.0 * h)[field], dtype=float)
        + 8.0 * np.asarray(contract.values(Y, eta + h)[field], dtype=float)
        - 8.0 * np.asarray(contract.values(Y, eta - h)[field], dtype=float)
        + np.asarray(contract.values(Y, eta - 2.0 * h)[field], dtype=float)
    ) / (12.0 * h)


def _fresh_fd_cloud(
    rng: np.random.Generator,
    parameters: Mapping[str, float],
    count: int,
) -> tuple[np.ndarray, np.ndarray]:
    ys: list[float] = []
    etas: list[float] = []
    while len(ys) < count:
        eta = float(rng.uniform(-0.9, 0.9))
        Y = float(rng.uniform(0.25, 3.85))
        H = float(_independent_values(Y, eta, parameters)["H_star"])
        if abs(H) >= 0.5:
            ys.append(Y)
            etas.append(eta)
    return np.asarray(ys), np.asarray(etas)


def audit_leading_axis_profile(
    contract: Any | None = None,
    *,
    seed: int = FRESH_SEED,
    sample_count: int = FRESH_SAMPLE_COUNT,
) -> dict[str, Any]:
    """Run the independent value/derivative/off-grid audit and return a receipt."""
    if contract is None:
        contract = KokunoPA10LeadingAxisProfileContract()
    configuration = contract.configuration()
    parameters = _independent_parameters(configuration)
    rng = np.random.default_rng(int(seed))

    E = 1.0 + float(parameters["enlarged_real_margin"])
    random_Y = rng.uniform(0.0, 4.1, int(sample_count))
    random_eta = rng.uniform(-E, E, int(sample_count))
    probe_Y = np.asarray([0.0, 1.0e-12, 1.0e-8, 4.1, 0.0, 4.1])
    probe_eta = np.asarray([-E, -1.0e-12, 0.0, 1.0e-12, E, 1.0])
    Y = np.concatenate([random_Y, probe_Y])
    eta = np.concatenate([random_eta, probe_eta])

    public_values = contract.values(Y, eta)
    public_derivatives = contract.derivatives(Y, eta)
    expected_values = _independent_values(Y, eta, parameters)
    expected_derivatives = _independent_derivatives(Y, eta, parameters)

    value_fields = (
        "d", "L", "U_star", "H_star", "W_star", "Z_star", "Pi_0",
        "Pi_0_eta", "chi", "zeta_star", "u_0", "u_0_Y",
    )
    derivative_fields = (
        "d_eta", "L_eta", "U_star_eta", "H_star_eta", "W_star_eta",
        "Z_star_eta", "Pi_0_eta_eta", "chi_eta", "zeta_star_eta",
        "u_0_Y", "u_0_eta", "u_0_Y_eta",
    )
    value_errors = {
        name: _relative_max(public_values[name], expected_values[name])
        for name in value_fields
    }
    derivative_errors = {
        name: _relative_max(public_derivatives[name], expected_derivatives[name])
        for name in derivative_fields
    }

    axis_eta = np.asarray([-E, -1.0, -1.0e-12, 0.0, 1.0e-12, 1.0, E])
    axis_values = contract.values(np.zeros_like(axis_eta), axis_eta)
    near_values = contract.values(np.full_like(axis_eta, 1.0e-12), axis_eta)
    axis_exact_zero = bool(np.all(np.asarray(axis_values["u_0"]) == 0.0))
    axis_scaling_error = _relative_max(
        np.asarray(near_values["u_0"]) / 1.0e-12,
        np.asarray(near_values["u_0_Y"]),
    )
    nontrivial_u0Y_rms = math.sqrt(
        float(np.mean(np.square(np.asarray(public_values["u_0_Y"], dtype=float))))
    )

    identity = (
        np.asarray(public_values["H_star"])
        * np.asarray(public_values["zeta_star"])
        / np.asarray(public_values["L"])
        + np.asarray(public_values["chi"])
    )
    identity_max = float(np.max(np.abs(identity)))

    fd_Y, fd_eta = _fresh_fd_cloud(rng, parameters, FD_SAMPLE_COUNT)
    fd_fields = {
        "Z_star": "Z_star_eta",
        "Pi_0": "Pi_0_eta",
        "chi": "chi_eta",
        "zeta_star": "zeta_star_eta",
        "u_0": "u_0_eta",
    }
    fd_expected = _independent_derivatives(fd_Y, fd_eta, parameters)
    fd_errors: dict[str, list[float]] = {name: [] for name in fd_fields}
    for step in FD_STEPS:
        for value_name, derivative_name in fd_fields.items():
            numeric = _fd4_value_derivative(contract, fd_Y, fd_eta, step, value_name)
            fd_errors[value_name].append(
                _relative_rms(numeric, fd_expected[derivative_name])
            )

    fd_refinement_ratios = {
        name: [
            errors[0] / max(errors[1], np.finfo(float).tiny),
            errors[1] / max(errors[2], np.finfo(float).tiny),
        ]
        for name, errors in fd_errors.items()
    }
    min_fd_refinement_ratio = min(
        ratio
        for ratios in fd_refinement_ratios.values()
        for ratio in ratios
    )
    max_fd_fine_relative_rms = max(errors[-1] for errors in fd_errors.values())

    max_value_error = max(value_errors.values())
    max_derivative_error = max(derivative_errors.values())
    guards = {
        "fresh_value_formula_replay": max_value_error <= VALUE_RELATIVE_MAX_GATE,
        "fresh_analytic_derivative_replay": (
            max_derivative_error <= DERIVATIVE_RELATIVE_MAX_GATE
        ),
        "axis_u0_exact_zero": axis_exact_zero,
        "axis_near_linear_scaling": axis_scaling_error <= AXIS_SCALING_RELATIVE_GATE,
        "u0_nontrivial": nontrivial_u0Y_rms > 1.0e-8,
        "H_zeta_over_L_plus_chi_identity": identity_max <= 1.0e-12,
        "fd4_three_level_refinement": (
            min_fd_refinement_ratio >= FD_REFINEMENT_RATIO_GATE
        ),
        "fd4_finest_relative_rms": (
            max_fd_fine_relative_rms <= FD_FINE_RELATIVE_RMS_GATE
        ),
    }
    failed_guards = [name for name, passed in guards.items() if not passed]

    try:
        exact_head = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:
        exact_head = None

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "exact_head": exact_head,
        "seed": int(seed),
        "fresh_sample_count": int(sample_count),
        "axis_and_boundary_probe_count": int(probe_Y.size + axis_eta.size),
        "fd_sample_count": FD_SAMPLE_COUNT,
        "fd_steps": list(FD_STEPS),
        "gates": {
            "value_relative_max": VALUE_RELATIVE_MAX_GATE,
            "derivative_relative_max": DERIVATIVE_RELATIVE_MAX_GATE,
            "fd_fine_relative_rms": FD_FINE_RELATIVE_RMS_GATE,
            "fd_refinement_ratio": FD_REFINEMENT_RATIO_GATE,
            "axis_scaling_relative": AXIS_SCALING_RELATIVE_GATE,
            "final_momentum_max_L2": 1.0e-3,
            "final_divergence_max_L2": 1.0e-5,
        },
        "max_value_relative_error": max_value_error,
        "value_relative_errors": value_errors,
        "max_derivative_relative_error": max_derivative_error,
        "derivative_relative_errors": derivative_errors,
        "axis_u0_exact_zero": axis_exact_zero,
        "axis_near_scaling_relative_error": axis_scaling_error,
        "u0Y_rms": nontrivial_u0Y_rms,
        "H_zeta_over_L_plus_chi_abs_max": identity_max,
        "fd4_relative_rms_by_field": fd_errors,
        "fd4_refinement_ratios_by_field": fd_refinement_ratios,
        "minimum_fd4_refinement_ratio": min_fd_refinement_ratio,
        "maximum_finest_fd4_relative_rms": max_fd_fine_relative_rms,
        "guards": guards,
        "failed_guards": failed_guards,
        "independent_profile_contract_passed": not failed_guards,
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
    report = audit_leading_axis_profile(seed=args.seed, sample_count=args.samples)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not report["independent_profile_contract_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
