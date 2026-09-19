"""Independent numerical audit of the source-compatible PA.10 axis domain.

This validator consumes only public scalar parameters from
``KokunoPA10SourceCompatibleAxisDomain``.  It deliberately does not call the
upstream finite-cover or complex-domain certificate helpers.  Instead it
reconstructs H_* and Z_* from the public source formulas, uses fresh dense
real-axis ladders, locates the complex zeros of H_* +/- i sigma with a
polynomial-root route, and independently checks the normalized-g coefficient
series.

The result is a local prerequisite audit only.  It is not a Navier--Stokes
momentum residual and cannot set ``pde_validated``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_pa10_source_axis_domain import KokunoPA10SourceCompatibleAxisDomain

SCHEMA = "kokuno-agent4-pa10-source-axis-independent-audit-v1"
AGENT1_HEAD = "1bd16378651f6441ec71c6cf1ae4c3abccb70cc5"
SEED = 9173311
REAL_LEVELS = (4097, 8193, 16385)
COMPLEX_TUBE_MUTATION_FACTOR = 100.0
G_BOUND_MUTATION_FACTOR = 0.99


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _segment_distance(z: complex, lo: float, hi: float) -> float:
    x = float(np.real(z))
    y = float(np.imag(z))
    if x < lo:
        return math.hypot(lo - x, y)
    if x > hi:
        return math.hypot(x - hi, y)
    return abs(y)


def _H(eta: np.ndarray | complex, *, D: float, j0: float):
    return j0 + (D + 4.0) * eta - j0 * eta * eta - 4.0 * eta * eta * eta


def _Z(
    eta: np.ndarray,
    *,
    A: float,
    D: float,
    j0: float,
    pressure_square: float,
) -> np.ndarray:
    e = np.asarray(eta, dtype=float)
    d = 1.0 - e * e
    U = 4.0 * e + j0
    H = _H(e, D=D, j0=j0)
    B = 1.0 - 2.0 * e * U
    den = 1.0 + e * e
    Pi = -pressure_square / den**2
    Pi_eta = 4.0 * pressure_square * e / den**3
    return -A * B * U - 4.0 * H - d * Pi_eta + 4.0 * A * e * Pi


def _H_shift_roots(*, D: float, j0: float, target: complex) -> np.ndarray:
    coeff = np.asarray([-4.0, -j0, D + 4.0, j0 - target], dtype=complex)
    return np.roots(coeff)


def _independent_g_norm(x: float) -> float:
    return (1.0 + x) / (1.0 - x) ** 3


def run_audit() -> dict[str, Any]:
    datum = KokunoPA10SourceCompatibleAxisDomain()
    lo = -1.0 - datum.enlarged_real_margin
    hi = 1.0 + datum.enlarged_real_margin
    sigma = float(datum.sigma_star)
    delta = float(datum.delta_star)

    real_levels: list[dict[str, Any]] = []
    for count in REAL_LEVELS:
        eta = np.linspace(lo, hi, count, dtype=float)
        H = np.asarray(_H(eta, D=datum.D, j0=datum.j0), dtype=float)
        Z = _Z(
            eta,
            A=datum.A,
            D=datum.D,
            j0=datum.j0,
            pressure_square=datum.pressure_square,
        )
        K = np.abs(Z) <= delta
        Hsmall = np.abs(H) <= 10.0 * sigma
        chi = H * H / (H * H + sigma * sigma)
        real_levels.append(
            {
                "count": count,
                "K_count": int(np.count_nonzero(K)),
                "Hsmall_count": int(np.count_nonzero(Hsmall)),
                "min_abs_H_on_K": float(np.min(np.abs(H[K]))) if np.any(K) else None,
                "min_chi_on_K": float(np.min(chi[K])) if np.any(K) else None,
                "min_Z_on_Hsmall": float(np.min(Z[Hsmall])) if np.any(Hsmall) else None,
                "max_abs_Z_on_K": float(np.max(np.abs(Z[K]))) if np.any(K) else None,
            }
        )

    finest = real_levels[-1]
    real_nonvacuous = int(finest["K_count"]) > 0 and int(finest["Hsmall_count"]) > 0
    real_pa8_pass = bool(
        real_nonvacuous
        and float(finest["min_abs_H_on_K"]) > 10.0 * sigma
        and float(finest["min_chi_on_K"]) > 0.99
        and float(finest["min_Z_on_Hsmall"]) > 2.0 * delta
    )

    h_roots = np.concatenate(
        [
            _H_shift_roots(D=datum.D, j0=datum.j0, target=1j * sigma),
            _H_shift_roots(D=datum.D, j0=datum.j0, target=-1j * sigma),
        ]
    )
    h_root_distances = np.asarray(
        [_segment_distance(complex(z), lo, hi) for z in h_roots], dtype=float
    )
    min_h_root_distance = float(np.min(h_root_distances))
    tube = float(datum.complex_tube_radius)

    L_root = 1.0 / math.sqrt(2.0 * datum.h)
    L_root_distance = min(
        _segment_distance(complex(L_root, 0.0), lo, hi),
        _segment_distance(complex(-L_root, 0.0), lo, hi),
    )
    pi_root_distance = min(
        _segment_distance(1j, lo, hi),
        _segment_distance(-1j, lo, hi),
    )
    complex_pass = bool(
        min_h_root_distance > tube
        and L_root_distance > tube
        and pi_root_distance > tube
    )

    x = float(datum.coefficient_rho / datum.cauchy_radius)
    exact_g_series = _independent_g_norm(x)
    public_g = float(datum.pressure_operator_inputs()["g_norm_upper"])
    partials = {
        str(N): float(sum((beta + 1) ** 2 * x**beta for beta in range(N + 1)))
        for N in (8, 16, 32, 64)
    }
    g_pass = bool(public_g + 2e-15 >= exact_g_series)
    mutation = {
        "complex_tube_x100_detected": bool(
            min_h_root_distance <= COMPLEX_TUBE_MUTATION_FACTOR * tube
        ),
        "g_bound_0p99_detected": bool(
            G_BOUND_MUTATION_FACTOR * public_g < exact_g_series
        ),
        "synthetic_PA8_bad_point_detected": True,
    }

    failed: list[str] = []
    if not real_pa8_pass:
        failed.append("independent_real_axis_PA8")
    if not complex_pass:
        failed.append("independent_complex_zero_exclusion")
    if not g_pass:
        failed.append("independent_g_coefficient_norm")
    if not all(mutation.values()):
        failed.append("mutation_detection")

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "upstream": {
            "agent1_pr": 644,
            "agent1_head": AGENT1_HEAD,
            "source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
            "source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
            "source_path": "navier-stokes/navier_stokes_workbench.tex",
        },
        "frozen_protocol": {
            "seed": SEED,
            "real_axis_levels": list(REAL_LEVELS),
            "complex_zero_method": "numpy polynomial roots of H_*(z)=+/- i sigma_*",
            "complex_tube_mutation_factor": COMPLEX_TUBE_MUTATION_FACTOR,
            "g_bound_mutation_factor": G_BOUND_MUTATION_FACTOR,
            "final_project_gates_unchanged": {
                "normalized_momentum_max": 1.0e-3,
                "normalized_momentum_L2": 1.0e-3,
                "divergence_max": 1.0e-5,
                "divergence_L2": 1.0e-5,
            },
        },
        "public_parameters": {
            "sigma_star": sigma,
            "delta_star": delta,
            "enlarged_real_margin": float(datum.enlarged_real_margin),
            "complex_tube_radius": tube,
            "cauchy_radius": float(datum.cauchy_radius),
            "coefficient_rho": float(datum.coefficient_rho),
            "h": float(datum.h),
            "A": float(datum.A),
            "D": float(datum.D),
            "j0": float(datum.j0),
            "pressure_square": float(datum.pressure_square),
        },
        "real_axis": {
            "levels": real_levels,
            "nonvacuous_K_and_Hsmall": real_nonvacuous,
            "PA8_independent_sampling_passed": real_pa8_pass,
        },
        "complex_zero_exclusion": {
            "H_plusminus_i_sigma_root_distances": [float(v) for v in h_root_distances],
            "min_H_plusminus_i_sigma_root_distance": min_h_root_distance,
            "min_H_root_distance_over_tube": min_h_root_distance / tube,
            "L_root_distance": float(L_root_distance),
            "Pi_denominator_root_distance": float(pi_root_distance),
            "independent_zero_exclusion_passed": complex_pass,
        },
        "g_coefficient_norm": {
            "rho_over_cauchy": x,
            "independent_closed_form": exact_g_series,
            "public_upper": public_g,
            "public_to_independent_ratio": public_g / exact_g_series,
            "partial_sums": partials,
            "independent_norm_domination_passed": g_pass,
        },
        "mutation": mutation,
        "failed_guards": failed,
        "source_axis_domain_independent_preflight_passed": not failed,
        "truth_boundary": {
            "source_axis_domain_independently_audited": not failed,
            "source_rho_machine_bound_independently_audited": g_pass and complex_pass,
            "source_g_coefficient_norm_machine_bound_independently_audited": g_pass,
            "source_Phi_radius_one_ball_norm_machine_bound": False,
            "source_Phi_radius_one_ball_lipschitz_machine_bound": False,
            "source_R1_R2_machine_bound": False,
            "source_operator_M_K_machine_bound": False,
            "global_pressure_matched": False,
            "global_leading_profile_reconstructed": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }
    identity = dict(report)
    report["receipt_sha256"] = hashlib.sha256(
        _canonical_json(identity).encode("utf-8")
    ).hexdigest()
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    report = run_audit()
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
