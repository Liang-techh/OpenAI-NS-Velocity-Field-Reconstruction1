"""Independent audit of the hierarchical Kokuno I2 heat repair.

This validator intentionally does *not* call the construction-side Gauss--Legendre
moment map, its Jacobian, or its nonlinear map.  It consumes only the public
hierarchical coefficient solution and reconstructs the three moment rows with a
separate composite-Simpson implementation of the autonomous C-infinity bumps.

The point of the audit is precision, not another algebraic replay.  The target
channels can differ by more than 400 decades.  Therefore a Decimal Newton solve
can close a moment map to many digits only if the moment coefficients themselves
are known accurately enough.  This module measures what happens when the same
public coefficient vector is evaluated through an independently discretized
continuous moment integral.

This is still a local I2 structural diagnostic.  It is not a full-domain
Navier--Stokes residual test and never promotes ``pde_validated``.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_hierarchical_heat_repair import KokunoHierarchicalHeatRepair


SCHEMA = "kokuno-agent4-hierarchical-heat-repair-independent-v1"
CHANNELS = ("C_p", "S", "I_sub")
RESOLUTIONS = (1024, 2048, 4096)
ETA_PROBES = (-0.6, -0.2, 0.2, 0.6)

# Deliberately duplicated from the documented autonomous three-bump realization
# rather than imported from the construction module.  The validator owns its
# quadrature, basis evaluation and moment assembly.
_BUMP_CENTERS = (0.96, 1.24, 1.53)
_BUMP_HALF_WIDTH = 0.105

FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5


@dataclass(frozen=True)
class IndependentMomentTensors:
    resolution: int
    linear: np.ndarray
    quadratic: np.ndarray


def _simpson_weights(intervals: int) -> np.ndarray:
    if intervals <= 0 or intervals % 2:
        raise ValueError("Simpson interval count must be positive and even")
    weights = np.ones(intervals + 1, dtype=float)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    return weights


def _independent_bump(x: np.ndarray, center: float) -> np.ndarray:
    """Independent evaluation of exp(1-1/(1-s^2)) on one compact support."""

    values = np.asarray(x, dtype=float)
    s = (values - float(center)) / _BUMP_HALF_WIDTH
    out = np.zeros_like(values)
    active = np.abs(s) < 1.0
    if np.any(active):
        sa = s[active]
        out[active] = np.exp(1.0 - 1.0 / (1.0 - sa * sa))
    return out


def independent_moment_tensors(
    *, lambda_outer: float, f_eta: float, resolution: int
) -> IndependentMomentTensors:
    """Reconstruct the three-row map with independent composite Simpson rules.

    Each disjoint bump support is integrated separately, so support endpoints are
    represented exactly by the independent grid and cross-bump quadratic terms
    are exactly zero from the declared disjoint supports.
    """

    lam = float(lambda_outer)
    f_value = float(f_eta)
    intervals = int(resolution)
    if not (0.0 < lam <= 0.5):
        raise ValueError("lambda_outer must lie in (0,0.5]")
    if not math.isfinite(f_value) or f_value <= 0.0:
        raise ValueError("f_eta must be positive and finite")
    if intervals % 2 or intervals < 128:
        raise ValueError("resolution must be an even integer >=128")

    linear = np.zeros((3, 3), dtype=float)
    quadratic = np.zeros((3, 3, 3), dtype=float)
    base_weights = _simpson_weights(intervals)

    for column, center in enumerate(_BUMP_CENTERS):
        left = center - _BUMP_HALF_WIDTH
        right = center + _BUMP_HALF_WIDTH
        x = np.linspace(left, right, intervals + 1, dtype=float)
        h = (right - left) / intervals
        weights = base_weights * (h / 3.0)
        beta = _independent_bump(x, center)

        linear[0, column] = np.sum(
            weights * f_value * x ** (-1.5 - lam) * beta
        )
        linear[1, column] = np.sum(
            weights * (-f_value) * x ** (-0.5 - lam) * beta
        )
        linear[2, column] = np.sum(weights * np.sqrt(2.0 * x) * beta)

        beta2 = beta * beta
        quadratic[0, column, column] = np.sum(weights * 0.5 * beta2 / x)
        quadratic[1, column, column] = np.sum(weights * -0.5 * beta2)

    return IndependentMomentTensors(
        resolution=intervals,
        linear=linear,
        quadratic=quadratic,
    )


def _decimal(value: float) -> Decimal:
    if not math.isfinite(float(value)):
        raise ValueError("non-finite tensor entry")
    return Decimal(repr(float(value)))


def _map_decimal(
    coefficients: tuple[Decimal, Decimal, Decimal],
    tensors: IndependentMomentTensors,
    precision_digits: int,
) -> tuple[Decimal, Decimal, Decimal]:
    with localcontext() as ctx:
        ctx.prec = int(precision_digits)
        out: list[Decimal] = []
        for row in range(3):
            value = sum(
                (
                    _decimal(tensors.linear[row, column]) * coefficients[column]
                    for column in range(3)
                ),
                Decimal(0),
            )
            for i in range(3):
                for j in range(3):
                    q = tensors.quadratic[row, i, j]
                    if q != 0.0:
                        value += _decimal(q) * coefficients[i] * coefficients[j]
            out.append(+value)
    return tuple(out)  # type: ignore[return-value]


def _log10_abs(value: Decimal) -> float | None:
    if value == 0:
        return None
    with localcontext() as ctx:
        ctx.prec = max(80, len(value.as_tuple().digits) + 20)
        return float(abs(value).log10())


def _relative_error(actual: Decimal, target: Decimal) -> Decimal:
    if target == 0:
        return Decimal(0) if actual == 0 else Decimal("Infinity")
    return abs(actual - target) / abs(target)


def _log10_relative(actual: Decimal, target: Decimal) -> float | None:
    relative = _relative_error(actual, target)
    if relative == 0:
        return None
    if not relative.is_finite():
        return math.inf
    return float(relative.log10())


def _tensor_refinement_delta(
    coarse: IndependentMomentTensors, fine: IndependentMomentTensors
) -> float:
    return float(
        max(
            np.max(np.abs(fine.linear - coarse.linear)),
            np.max(np.abs(fine.quadratic - coarse.quadratic)),
        )
    )


def _eta_case(repair: KokunoHierarchicalHeatRepair, eta: float) -> dict[str, Any]:
    solution = repair.solve(float(eta))
    coefficients = solution.coefficient_decimals()
    targets = tuple(Decimal(value) for value in solution.target)
    f_eta = 1.0 / (1.0 + float(eta) ** 2)

    tensors = [
        independent_moment_tensors(
            lambda_outer=repair.outer_schedule.lambda_outer,
            f_eta=f_eta,
            resolution=resolution,
        )
        for resolution in RESOLUTIONS
    ]
    evaluations: list[dict[str, Any]] = []
    for tensor in tensors:
        achieved = _map_decimal(coefficients, tensor, solution.precision_digits)
        residual = tuple(achieved[i] - targets[i] for i in range(3))
        evaluations.append(
            {
                "resolution": tensor.resolution,
                "achieved": [str(value) for value in achieved],
                "residual": [str(value) for value in residual],
                "log10_abs_residual": [_log10_abs(value) for value in residual],
                "log10_relative_error": [
                    _log10_relative(achieved[i], targets[i]) for i in range(3)
                ],
            }
        )

    finest = tensors[-1]
    finest_achieved = _map_decimal(coefficients, finest, solution.precision_digits)
    finest_relative = [
        _relative_error(finest_achieved[i], targets[i]) for i in range(3)
    ]

    mutated = list(coefficients)
    mutated[0] *= Decimal("1.001")
    mutated_achieved = _map_decimal(
        tuple(mutated), finest, solution.precision_digits  # type: ignore[arg-type]
    )
    mutated_i_relative = _relative_error(mutated_achieved[2], targets[2])

    return {
        "eta": float(eta),
        "solution_sha256": solution.sha256,
        "precision_digits": solution.precision_digits,
        "target_sign": list(solution.target_sign),
        "target_log10_abs": [
            None if sign == 0 else log_abs / math.log(10.0)
            for sign, log_abs in zip(solution.target_sign, solution.target_log_abs)
        ],
        "coefficient_log10_abs": [
            _log10_abs(value) for value in coefficients
        ],
        "resolution_evaluations": evaluations,
        "tensor_refinement_max_abs": [
            _tensor_refinement_delta(tensors[i], tensors[i + 1])
            for i in range(len(tensors) - 1)
        ],
        "finest_relative_error": [
            str(value) if value.is_finite() else "Infinity"
            for value in finest_relative
        ],
        "finest_log10_relative_error": [
            None if value == 0 else float(value.log10())
            for value in finest_relative
        ],
        "mutation": {
            "operation": "multiply coefficient[0] by 1.001",
            "I_sub_relative_error": str(mutated_i_relative),
        },
    }


def build_report() -> dict[str, Any]:
    repair = KokunoHierarchicalHeatRepair()
    cases = [_eta_case(repair, eta) for eta in ETA_PROBES]

    max_tensor_refinement = max(
        delta for case in cases for delta in case["tensor_refinement_max_abs"]
    )
    max_i_relative = max(
        float(case["finest_relative_error"][2]) for case in cases
    )
    min_cp_log_rel = min(
        case["finest_log10_relative_error"][0] for case in cases
    )
    min_s_log_rel = min(
        case["finest_log10_relative_error"][1] for case in cases
    )
    min_mutation_i_relative = min(
        float(case["mutation"]["I_sub_relative_error"]) for case in cases
    )

    local_guards = {
        # These are audit guards, not replacements for the project PDE gates.
        "independent_tensor_refinement_max_abs_le_1e-10": (
            max_tensor_refinement <= 1.0e-10
        ),
        "dominant_I_sub_channel_relative_error_le_1e-10": (
            max_i_relative <= 1.0e-10
        ),
        # Negative-control expectation: if the construction moment tensor was
        # only float64-accurate, the tiny two channels should *not* magically
        # remain certified under an independent moment integration.
        "Cp_independent_mismatch_exceeds_target_by_50_decades": (
            min_cp_log_rel >= 50.0
        ),
        "S_independent_mismatch_exceeds_target_by_50_decades": (
            min_s_log_rel >= 50.0
        ),
        "coefficient_mutation_detected_in_I_sub_ge_1e-6": (
            min_mutation_i_relative >= 1.0e-6
        ),
    }

    return {
        "schema": SCHEMA,
        "task": "KOKUNO-A4-HIERARCHICAL-HEAT-REPAIR-INDEPENDENT-014",
        "construction_contract": {
            "consumed": "public KokunoHierarchicalHeatRepair.solve coefficient receipt",
            "not_used": [
                "construction Gauss-Legendre quadrature",
                "construction coefficient_jacobian",
                "construction nonlinear moment map",
                "training tensors/loss",
                "pressure fit",
                "forcing fit",
                "free f=R",
            ],
        },
        "independent_operator": {
            "basis": "separately coded autonomous C-infinity three-bump basis",
            "quadrature": "per-support composite Simpson",
            "resolutions": list(RESOLUTIONS),
            "eta_probes": list(ETA_PROBES),
            "moment_evaluation": "Decimal contraction of independently integrated float64 moment tensors",
        },
        "cases": cases,
        "summary": {
            "max_tensor_refinement_abs": max_tensor_refinement,
            "max_I_sub_relative_error": max_i_relative,
            "minimum_Cp_log10_relative_error": min_cp_log_rel,
            "minimum_S_log10_relative_error": min_s_log_rel,
            "minimum_mutated_I_sub_relative_error": min_mutation_i_relative,
            "independent_continuous_moment_closure": False,
            "frozen_discrete_map_only": True,
        },
        "local_guards": local_guards,
        "local_audit_completed": all(local_guards.values()),
        "routing": {
            "agent1": (
                "Rebuild the three-bump moment tensor itself at precision commensurate "
                "with the signed-log target, or use a hierarchical formulation that "
                "does not require cancellation of float64 moment coefficients. Do not "
                "promote the current frozen-discrete closure to source heat compensation."
            ),
            "agent2": "No change: actual-source public oscillatory velocity is still pending.",
            "agent3": "No change: finite correction cycle remains blocked pending an independent second covariance direction.",
        },
        "formal_project_gates": {
            "normalized_momentum_max": FORMAL_MOMENTUM_GATE,
            "divergence_max": FORMAL_DIVERGENCE_GATE,
            "assessed_here": False,
            "pde_validated": False,
        },
        "truth_boundary": {
            "full_domain_candidate_available": False,
            "formal_full_domain_pde_gate_assessed": False,
            "continuous_source_moment_compensation_certified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Independently audit the hierarchical Kokuno heat repair"
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
    if not report["local_audit_completed"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
