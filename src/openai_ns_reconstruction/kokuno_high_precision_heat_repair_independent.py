"""Independent high-precision audit of the Kokuno I2 heat repair.

Agent 1's construction-side implementation in ``kokuno_high_precision_heat_repair``
rebuilds the autonomous three-bump moment tensor with arbitrary-precision
*tanh-sinh* quadrature before solving the signed-log target.  This module does
not call that quadrature, its tensor builder, or its internal moment map.

Instead, it consumes only the public coefficient receipt returned by
``KokunoHighPrecisionHierarchicalHeatRepair.solve`` and reconstructs the same
continuous three-bump moments through the exact change of variables

    s = tanh(y),   y in R,

for which the compact bump

    exp(1 - 1/(1-s^2))

becomes ``exp(-sinh(y)^2)`` and ``ds = sech(y)^2 dy``.  The transformed
integrands are then integrated by a symmetric arbitrary-precision trapezoidal
rule on the real line.  This is a different discretization and a different
coordinate transform from the construction tanh-sinh path.

The audit is intentionally local.  Passing it can certify the continuous
three-row I2 moment seam for this autonomous bump realization; it cannot certify
heat/core matching, a global Kokuno velocity, or the full-domain Navier-Stokes
gate.
"""

from __future__ import annotations

from functools import lru_cache
import argparse
import json
import math
from pathlib import Path
from typing import Any

import mpmath as mp

from .kokuno_high_precision_heat_repair import (
    KokunoHighPrecisionHierarchicalHeatRepair,
)


SCHEMA = "kokuno-agent4-high-precision-heat-repair-independent-v1"
CHANNELS = ("C_p", "S", "I_sub")
ETA_PROBES = (-0.73, -0.31, 0.17, 0.59)
STEP_DENOMINATORS = (256, 512, 1024)
_BUMP_CENTERS = ("0.96", "1.24", "1.53")
_BUMP_HALF_WIDTH = "0.105"
_EXTRA_AUDIT_DIGITS = 96
_TAIL_EXTRA_DIGITS = 64
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5


@lru_cache(maxsize=16)
def _independent_base_moment_strings(
    lambda_text: str,
    precision_digits: int,
    step_denominator: int,
) -> tuple[tuple[str, str, str, str, str], ...]:
    """Integrate the five per-bump base moments with an independent y-grid.

    Returned entries are ``(L_Cp/f, L_S/f, L_I, Q_Cp, Q_S)``.  The first two
    linear entries are multiplied by ``f(eta)`` only after this eta-independent
    integration.  Cross-bump quadratic terms are exactly zero because the
    declared compact supports are disjoint.
    """

    digits = int(precision_digits)
    denominator = int(step_denominator)
    if digits < 64 or digits > 1100:
        raise ValueError("precision_digits must lie in [64,1100]")
    if denominator < 32 or denominator & (denominator - 1):
        raise ValueError("step_denominator must be a power of two >=32")

    lam_float = float(lambda_text)
    if not math.isfinite(lam_float) or not (0.0 < lam_float <= 0.5):
        raise ValueError("lambda_outer must lie in (0,0.5]")

    with mp.workdps(digits + 24):
        lam = mp.mpf(lambda_text)
        half_width = mp.mpf(_BUMP_HALF_WIDTH)
        centers = tuple(mp.mpf(value) for value in _BUMP_CENTERS)
        h = mp.mpf(1) / denominator

        # beta(tanh(y)) = exp(-sinh(y)^2).  Choose the finite symmetric window
        # so the omitted linear-bump tail is far below the requested precision.
        tail_digits = digits + _TAIL_EXTRA_DIGITS
        y_limit = mp.asinh(mp.sqrt(mp.mpf(tail_digits) * mp.log(10)))
        count = int(mp.ceil(y_limit / h))

        sums = [[mp.mpf(0) for _ in range(5)] for _ in centers]
        for k in range(count + 1):
            y_abs = k * h
            signed_y = (mp.mpf(0),) if k == 0 else (y_abs, -y_abs)
            for y in signed_y:
                sinh_y = mp.sinh(y)
                cosh_y = mp.cosh(y)
                s = mp.tanh(y)
                beta = mp.exp(-(sinh_y * sinh_y))
                jacobian_weight = half_width * h / (cosh_y * cosh_y)

                for index, center in enumerate(centers):
                    x = center + half_width * s
                    beta2 = beta * beta
                    values = (
                        mp.power(x, -mp.mpf("1.5") - lam) * beta,
                        -mp.power(x, -mp.mpf("0.5") - lam) * beta,
                        mp.sqrt(2 * x) * beta,
                        mp.mpf("0.5") * beta2 / x,
                        -mp.mpf("0.5") * beta2,
                    )
                    for row, value in enumerate(values):
                        sums[index][row] += value * jacobian_weight

        return tuple(
            tuple(
                mp.nstr(value, n=digits + 8, strip_zeros=False)
                for value in row
            )
            for row in sums
        )


def _max_relative_moment_delta(
    coarse: tuple[tuple[str, str, str, str, str], ...],
    fine: tuple[tuple[str, str, str, str, str], ...],
    *,
    precision_digits: int,
) -> mp.mpf:
    with mp.workdps(precision_digits):
        worst = mp.mpf(0)
        for coarse_row, fine_row in zip(coarse, fine):
            for coarse_value, fine_value in zip(coarse_row, fine_row):
                a = mp.mpf(coarse_value)
                b = mp.mpf(fine_value)
                scale = max(abs(b), mp.mpf(1))
                worst = max(worst, abs(a - b) / scale)
        return +worst


def _evaluate_receipt(
    repair: KokunoHighPrecisionHierarchicalHeatRepair,
    eta: float,
    base_moments: tuple[tuple[str, str, str, str, str], ...],
    *,
    precision_digits: int,
) -> dict[str, Any]:
    solution = repair.solve(float(eta))
    with mp.workdps(precision_digits):
        coefficients = [mp.mpf(value) for value in solution.coefficients]
        targets = [mp.mpf(value) for value in solution.target]
        f_eta = mp.mpf(repr(float(repair.outer_schedule.source_f(float(eta)))))
        base = [[mp.mpf(value) for value in row] for row in base_moments]

        achieved = [mp.mpf(0), mp.mpf(0), mp.mpf(0)]
        for column, coefficient in enumerate(coefficients):
            cp_linear, s_linear, i_linear, cp_quad, s_quad = base[column]
            achieved[0] += f_eta * cp_linear * coefficient + cp_quad * coefficient**2
            achieved[1] += f_eta * s_linear * coefficient + s_quad * coefficient**2
            achieved[2] += i_linear * coefficient

        residual = [achieved[i] - targets[i] for i in range(3)]
        relative = [
            mp.mpf(0) if targets[i] == 0 and residual[i] == 0
            else mp.inf if targets[i] == 0
            else abs(residual[i]) / abs(targets[i])
            for i in range(3)
        ]

        mutated = list(coefficients)
        mutated[0] *= mp.mpf("1.001")
        mutated_i = mp.fsum(base[i][2] * mutated[i] for i in range(3))
        mutated_i_relative = (
            mp.inf if targets[2] == 0
            else abs(mutated_i - targets[2]) / abs(targets[2])
        )

        def log10_or_none(value: mp.mpf) -> float | None:
            if value == 0:
                return None
            if not mp.isfinite(value):
                return math.inf
            return float(mp.log10(abs(value)))

        return {
            "eta": float(eta),
            "solution_sha256": solution.sha256,
            "solution_precision_digits": int(solution.precision_digits),
            "construction_max_relative_residual": float(solution.max_relative_residual),
            "target_log10_abs": [log10_or_none(value) for value in targets],
            "coefficient_log10_abs": [log10_or_none(value) for value in coefficients],
            "independent_residual_log10_abs": [
                log10_or_none(value) for value in residual
            ],
            "independent_relative_error": [mp.nstr(value, 20) for value in relative],
            "independent_log10_relative_error": [
                log10_or_none(value) for value in relative
            ],
            "mutation": {
                "operation": "multiply coefficient[0] by 1.001",
                "I_sub_relative_error": mp.nstr(mutated_i_relative, 20),
            },
        }


def build_report() -> dict[str, Any]:
    repair = KokunoHighPrecisionHierarchicalHeatRepair()
    solutions = [repair.solve(float(eta)) for eta in ETA_PROBES]
    audit_digits = max(solution.precision_digits for solution in solutions) + _EXTRA_AUDIT_DIGITS
    lambda_text = repr(float(repair.outer_schedule.lambda_outer))

    tensors = [
        _independent_base_moment_strings(lambda_text, audit_digits, denominator)
        for denominator in STEP_DENOMINATORS
    ]
    refinement = [
        _max_relative_moment_delta(
            tensors[index], tensors[index + 1], precision_digits=audit_digits
        )
        for index in range(len(tensors) - 1)
    ]

    cases = [
        _evaluate_receipt(
            repair,
            eta,
            tensors[-1],
            precision_digits=audit_digits,
        )
        for eta in ETA_PROBES
    ]

    max_relative_by_channel = []
    for channel in range(3):
        max_relative_by_channel.append(
            max(mp.mpf(case["independent_relative_error"][channel]) for case in cases)
        )
    min_mutation_i = min(
        mp.mpf(case["mutation"]["I_sub_relative_error"]) for case in cases
    )
    finest_refinement = refinement[-1]

    local_guards = {
        "independent_moment_refinement_le_1e-200": finest_refinement <= mp.mpf("1e-200"),
        "Cp_relative_error_le_1e-20": max_relative_by_channel[0] <= mp.mpf("1e-20"),
        "S_relative_error_le_1e-20": max_relative_by_channel[1] <= mp.mpf("1e-20"),
        "I_sub_relative_error_le_1e-20": max_relative_by_channel[2] <= mp.mpf("1e-20"),
        "coefficient_mutation_detected_in_I_sub_ge_1e-6": min_mutation_i >= mp.mpf("1e-6"),
    }
    local_pass = all(local_guards.values())

    return {
        "schema": SCHEMA,
        "task": "KOKUNO-A4-HIGHPREC-ALTQUAD-AUDIT-015",
        "construction_contract": {
            "consumed": "public KokunoHighPrecisionHierarchicalHeatRepair.solve coefficient receipt",
            "not_used": [
                "construction _base_moment_strings",
                "construction tanh-sinh quadrature",
                "construction internal moment map/Jacobian",
                "training tensors/loss",
                "pressure fit",
                "forcing fit",
                "free f=R",
            ],
        },
        "independent_operator": {
            "transform": "s=tanh(y), beta(tanh(y))=exp(-sinh(y)^2)",
            "quadrature": "symmetric arbitrary-precision trapezoidal rule on y in R",
            "step_denominators": list(STEP_DENOMINATORS),
            "eta_probes": list(ETA_PROBES),
            "audit_precision_digits": int(audit_digits),
            "tail_extra_digits": _TAIL_EXTRA_DIGITS,
        },
        "cases": cases,
        "summary": {
            "moment_refinement_relative": [mp.nstr(value, 20) for value in refinement],
            "finest_moment_refinement_relative": mp.nstr(finest_refinement, 20),
            "max_relative_error_by_channel": [
                mp.nstr(value, 20) for value in max_relative_by_channel
            ],
            "minimum_mutated_I_sub_relative_error": mp.nstr(min_mutation_i, 20),
            "local_continuous_moment_closure_certified": bool(local_pass),
        },
        "local_guards": local_guards,
        "local_audit_completed": bool(local_pass),
        "routing": {
            "agent1": (
                "If this independent alternate-quadrature audit passes, the new high-precision "
                "three-bump moment tensor is no longer the active local precision blocker; proceed "
                "to applying the repair and global core-to-heat leading assembly. If it fails, keep "
                "the heat repair unpromoted and inspect the reported channel/refinement mismatch."
            ),
            "agent2": (
                "No PDE promotion: actual-source public oscillatory velocity remains required "
                "before a full composite gate."
            ),
            "agent3": (
                "No correction-cycle promotion: a genuinely independent second covariance "
                "direction remains required before rerunning the finite cycle."
            ),
        },
        "formal_project_gates": {
            "normalized_momentum_max": FORMAL_MOMENTUM_GATE,
            "divergence_max": FORMAL_DIVERGENCE_GATE,
            "assessed_here": False,
            "pde_validated": False,
        },
        "truth_boundary": {
            "local_continuous_three_bump_moment_seam_assessed": True,
            "heat_compensation_completed": False,
            "core_to_heat_matching_completed": False,
            "global_leading_profile_reconstructed": False,
            "complete_kokuno_composite_velocity": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Independently audit the high-precision Kokuno heat repair"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report["summary"], indent=2, sort_keys=True))
    return 0 if report["local_audit_completed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
