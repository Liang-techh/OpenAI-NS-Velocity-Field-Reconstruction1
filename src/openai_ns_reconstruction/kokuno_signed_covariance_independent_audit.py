"""Independent Agent-4 audit of Agent-2 #472 signed covariance inverse.

The production object is treated as a public value/differential contract.  This
audit rebuilds its two-sign reference inverse with 80-digit Decimal arithmetic,
checks the public directional derivative against a separate FD5 derivative of
the high-precision reference, exercises the semantic sigma swap, and calibrates
mutations/negative controls.  It does not validate a complete-curl physical
velocity or the full Navier-Stokes equation.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_source_signed_covariance_pair import KokunoSourceSignedCovariancePair

SCHEMA = "kokuno-agent4-signed-covariance-independent-audit-v1"
TASK_ID = "KOKUNO-A4-SIGNED-COVARIANCE-INDEPENDENT-AUDIT-029"
BASE_PR = 472
BASE_HEAD = "2bfb4c6f3a9314eea19b7dc0aba4fc11f2877978"
SEED = 9173181
DECIMAL_PRECISION = 80
FD5_STEPS = (0.08, 0.04, 0.02)
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5
FROZEN_GUARDS = {
    "value_relative_error_max": 5.0e-13,
    "finest_directional_derivative_relative_rms_max": 1.0e-6,
    "directional_derivative_refinement_ratio_min": 8.0,
    "semantic_sigma_swap_relative_max": 5.0e-13,
    "wrong_q_sign_mutation_relative_rms_min": 0.10,
    "drop_pulse_derivative_mutation_relative_rms_min": 0.03,
    "cone_negative_rejection_fraction_min": 1.0,
    "direction_gap_negative_rejection_fraction_min": 1.0,
}


def _d(x: float) -> Decimal:
    return Decimal(repr(float(x)))


def _decimal_solution(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    A, u, hp, hm, TN, TK = (_d(x) for x in np.asarray(values, float))
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PRECISION
        p = -TN / A
        q = TK / u
        yp = (p - q) / (Decimal(2) * hp)
        ym = (p + q) / (Decimal(2) * hm)
        if yp <= 0 or ym <= 0:
            raise ValueError("reference path left the positive covariance cone")
        y = np.asarray((float(yp), float(ym)))
        a = np.asarray((float(ctx.sqrt(yp)), float(ctx.sqrt(ym))))
    return y, a


def _fd5_reference(values: np.ndarray, tangent: np.ndarray, h: float) -> np.ndarray:
    def f(s: float) -> np.ndarray:
        y, a = _decimal_solution(values + s * tangent)
        return np.concatenate((y, a))

    return (f(-2 * h) - 8 * f(-h) + 8 * f(h) - f(2 * h)) / (12 * h)


def _rms(x: Any) -> float:
    arr = np.asarray(x, float)
    return float(np.sqrt(np.mean(arr * arr)))


def _relative_rms(err: Any, ref: Any) -> float:
    return _rms(err) / max(_rms(ref), 1.0e-30)


def _relative_max(err: Any, ref: Any) -> float:
    e = np.asarray(err, float)
    r = np.asarray(ref, float)
    return float(np.max(np.abs(e) / np.maximum(np.abs(r), 1.0e-30)))


def _samples() -> list[tuple[np.ndarray, np.ndarray]]:
    rng = np.random.default_rng(SEED)
    out: list[tuple[np.ndarray, np.ndarray]] = []
    while len(out) < 48:
        A = 10.0 ** rng.uniform(-1.0, 1.0)
        u = 10.0 ** rng.uniform(-1.0, 1.0)
        hp = 10.0 ** rng.uniform(-1.0, 1.0)
        hm = 10.0 ** rng.uniform(-1.0, 1.0)
        p = 10.0 ** rng.uniform(-0.3, 0.3)
        q = rng.choice((-1.0, 1.0)) * rng.uniform(0.25, 0.65) * p
        TN = -A * p
        TK = u * q
        values = np.asarray((A, u, hp, hm, TN, TK), float)
        tangent = np.asarray(
            (0.07 * A, -0.05 * u, 0.09 * hp, -0.08 * hm, 0.04 * TN, -0.06 * TK),
            float,
        )
        try:
            _decimal_solution(values - 0.16 * tangent)
            _decimal_solution(values + 0.16 * tangent)
        except ValueError:
            continue
        out.append((values, tangent))
    return out


def run_audit() -> dict[str, Any]:
    pair = KokunoSourceSignedCovariancePair()
    samples = _samples()

    value_errors = []
    derivative_errors = {h: [] for h in FD5_STEPS}
    swap_errors = []
    wrong_q_mutation = []
    dropped_h_mutation = []
    cone_rejects = 0
    gap_rejects = 0

    for values, tangent in samples:
        A, u, hp, hm, TN, TK = values
        public = pair.solve_reference_amplitudes(A, u, hp, hm, TN, TK, direction_gap_eta=0.2)
        y_ref, a_ref = _decimal_solution(values)
        public_vec = np.concatenate((public["squared_amplitudes"], public["amplitudes"]))
        ref_vec = np.concatenate((y_ref, a_ref))
        value_errors.append(_relative_max(public_vec - ref_vec, ref_vec))

        dd = pair.directional_derivative(
            A,
            u,
            hp,
            hm,
            TN,
            TK,
            dA_c=tangent[0],
            du_star=tangent[1],
            dh_plus=tangent[2],
            dh_minus=tangent[3],
            dT_N=tangent[4],
            dT_K=tangent[5],
        )
        public_derivative = np.concatenate(
            (dd["d_squared_amplitudes"], dd["d_amplitudes"])
        )
        for h in FD5_STEPS:
            fd = _fd5_reference(values, tangent, h)
            derivative_errors[h].append(_relative_rms(public_derivative - fd, fd))

        swapped = pair.solve_reference_amplitudes(A, u, hm, hp, TN, -TK)
        swapped_vec = np.concatenate(
            (swapped["squared_amplitudes"][::-1], swapped["amplitudes"][::-1])
        )
        swap_errors.append(_relative_max(swapped_vec - public_vec, public_vec))

        wrong_values = values.copy()
        wrong_values[5] = -wrong_values[5]
        wrong_y, wrong_a = _decimal_solution(wrong_values)
        wrong_q_mutation.append(
            _relative_rms(public_vec - np.concatenate((wrong_y, wrong_a)), public_vec)
        )

        no_h_tangent = tangent.copy()
        no_h_tangent[2:4] = 0.0
        mutated_fd = _fd5_reference(values, no_h_tangent, FD5_STEPS[-1])
        dropped_h_mutation.append(_relative_rms(public_derivative - mutated_fd, public_derivative))

        p = -TN / A
        try:
            pair.solve_reference_amplitudes(A, u, hp, hm, TN, u * 1.05 * p)
        except ValueError:
            cone_rejects += 1
        try:
            pair.solve_reference_amplitudes(
                A, u, hp, hm, TN, u * 0.90 * p, direction_gap_eta=0.2
            )
        except ValueError:
            gap_rejects += 1

    derivative_rms = {
        str(h): _rms(np.asarray(derivative_errors[h])) for h in FD5_STEPS
    }
    ratios = (
        derivative_rms[str(FD5_STEPS[0])] / derivative_rms[str(FD5_STEPS[1])],
        derivative_rms[str(FD5_STEPS[1])] / derivative_rms[str(FD5_STEPS[2])],
    )
    metrics = {
        "sample_count": len(samples),
        "value_relative_error_max": float(max(value_errors)),
        "directional_derivative_relative_rms_by_step": derivative_rms,
        "directional_derivative_refinement_ratios": [float(x) for x in ratios],
        "semantic_sigma_swap_relative_max": float(max(swap_errors)),
        "wrong_q_sign_mutation_relative_rms_min": float(min(wrong_q_mutation)),
        "drop_pulse_derivative_mutation_relative_rms_min": float(min(dropped_h_mutation)),
        "cone_negative_rejection_fraction": cone_rejects / len(samples),
        "direction_gap_negative_rejection_fraction": gap_rejects / len(samples),
    }
    g = FROZEN_GUARDS
    passed = bool(
        metrics["value_relative_error_max"] <= g["value_relative_error_max"]
        and derivative_rms[str(FD5_STEPS[-1])]
        <= g["finest_directional_derivative_relative_rms_max"]
        and min(ratios) >= g["directional_derivative_refinement_ratio_min"]
        and metrics["semantic_sigma_swap_relative_max"]
        <= g["semantic_sigma_swap_relative_max"]
        and metrics["wrong_q_sign_mutation_relative_rms_min"]
        >= g["wrong_q_sign_mutation_relative_rms_min"]
        and metrics["drop_pulse_derivative_mutation_relative_rms_min"]
        >= g["drop_pulse_derivative_mutation_relative_rms_min"]
        and metrics["cone_negative_rejection_fraction"]
        >= g["cone_negative_rejection_fraction_min"]
        and metrics["direction_gap_negative_rejection_fraction"]
        >= g["direction_gap_negative_rejection_fraction_min"]
    )
    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_pr": BASE_PR,
        "base_head": BASE_HEAD,
        "seed": SEED,
        "decimal_precision": DECIMAL_PRECISION,
        "fd5_steps": list(FD5_STEPS),
        "frozen_guards": FROZEN_GUARDS,
        "metrics": metrics,
        "local_structural_preflight_passed": passed,
        "formal_project_gates": {
            "normalized_momentum_max_l2": FORMAL_MOMENTUM_GATE,
            "normalized_divergence_max_l2": FORMAL_DIVERGENCE_GATE,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
        "truth_boundary": {
            "reference_signed_covariance_inverse_independently_audited": passed,
            "actual_positive_order_background_bound": False,
            "actual_source_h_sigma_pulse_integrals_bound": False,
            "actual_auxiliary_torus_signed_mode_family_bound": False,
            "public_source_bound_xyz_t_oscillatory_velocity_ready": False,
            "genuinely_independent_second_covariance_column_ready": False,
            "correction_ready": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
        "limitations": [
            "This audits only the displayed two-sign reference covariance inverse and its slow-coordinate derivative contract.",
            "No complete-curl physical oscillatory velocity, pressure, forcing, or Navier-Stokes residual is evaluated here.",
            "The 1e-3 momentum and 1e-5 divergence gates remain unchanged and unassessed for the Kokuno composite.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = run_audit()
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if report["local_structural_preflight_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
