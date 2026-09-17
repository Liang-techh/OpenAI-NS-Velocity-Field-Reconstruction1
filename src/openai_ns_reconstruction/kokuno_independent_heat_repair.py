"""Independent numerical preflight for Agent-1's Kokuno heat-repair primitive.

This audit is deliberately separate from the construction quadrature.  The
production primitive uses fixed Gauss--Legendre quadrature and an analytic
coefficient Jacobian; here the same public pointwise profile/correction API is
integrated on independently refined uniform grids with composite Simpson
quadrature, and the coefficient Jacobian/profile derivative are reconstructed
by centered finite differences.

This is a local source-map preflight only.  It is not a Navier--Stokes residual
validation because Agent 1 has not yet supplied the actual heat discrepancy or
a globally matched core-to-heat velocity/pressure candidate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import simpson

from .kokuno_heat_discrepancy_repair import KokunoHeatDiscrepancyRepair

TASK = "KOKUNO-A4-INDEPENDENT-HEAT-REPAIR-AUDIT-007"
UPSTREAM_HEAD = "d2b20eeb52e1914e2aaa99be9ef6437d819a3e53"
SEED = 9_173_051
LEVELS = (1025, 2049, 4097)
X_STAR_INTERVAL = (0.80, 1.70)
JACOBIAN_STEP = 1.0e-6
PARAMETER_MUTATION = 2.0e-3
RESIDUAL_REFERENCE = 1.0e-3
DIVERGENCE_REFERENCE = 1.0e-5


def _independent_increment(
    repair: KokunoHeatDiscrepancyRepair,
    coefficients: np.ndarray,
    f_eta: float,
    points: int,
) -> np.ndarray:
    if points < 65 or points % 2 == 0:
        raise ValueError("points must be an odd integer >= 65")
    x = np.linspace(X_STAR_INTERVAL[0], X_STAR_INTERVAL[1], points)
    base = np.asarray(repair.background_over_e_star(x, f_eta), dtype=float)
    correction = np.asarray(repair.correction_over_e_star(x, coefficients), dtype=float)
    delta_square = 2.0 * base * correction + correction * correction
    integrands = np.vstack(
        (
            0.5 * delta_square / x,
            -0.5 * delta_square,
            np.sqrt(2.0 * x) * correction,
        )
    )
    return np.asarray(simpson(integrands, x=x, axis=1), dtype=float)


def _independent_jacobian(
    repair: KokunoHeatDiscrepancyRepair,
    coefficients: np.ndarray,
    f_eta: float,
    points: int,
    step: float = JACOBIAN_STEP,
) -> np.ndarray:
    jac = np.empty((3, 3), dtype=float)
    for column in range(3):
        plus = coefficients.copy()
        minus = coefficients.copy()
        plus[column] += step
        minus[column] -= step
        jac[:, column] = (
            _independent_increment(repair, plus, f_eta, points)
            - _independent_increment(repair, minus, f_eta, points)
        ) / (2.0 * step)
    return jac


def _norm_relative(error: np.ndarray, reference: np.ndarray) -> float:
    denominator = max(float(np.linalg.norm(reference)), 1.0e-15)
    return float(np.linalg.norm(error) / denominator)


def _profile_derivative_audit(
    repair: KokunoHeatDiscrepancyRepair,
    coefficients: np.ndarray,
) -> dict[str, Any]:
    X_star = 2.3
    e_star = 0.8
    x_star = np.asarray([0.93, 0.96, 1.21, 1.24, 1.50, 1.53], dtype=float)
    X = X_star * x_star
    analytic = np.asarray(
        repair.physical_profile_correction(
            X, X_star=X_star, e_star=e_star, coefficients=coefficients
        )["delta_E_X"],
        dtype=float,
    )
    rows: list[dict[str, float]] = []
    for step in (4.0e-4, 2.0e-4, 1.0e-4):
        plus = np.asarray(
            repair.physical_profile_correction(
                X + step, X_star=X_star, e_star=e_star, coefficients=coefficients
            )["delta_E"],
            dtype=float,
        )
        minus = np.asarray(
            repair.physical_profile_correction(
                X - step, X_star=X_star, e_star=e_star, coefficients=coefficients
            )["delta_E"],
            dtype=float,
        )
        finite_difference = (plus - minus) / (2.0 * step)
        error = finite_difference - analytic
        rows.append(
            {
                "step": step,
                "rms_error": float(np.sqrt(np.mean(error * error))),
                "max_abs_error": float(np.max(np.abs(error))),
                "relative_error": _norm_relative(error, analytic),
            }
        )
    return {"rows": rows, "finest": rows[-1]}


def run_audit(
    *,
    seed: int = SEED,
    levels: tuple[int, int, int] = LEVELS,
    sample_count: int = 6,
) -> dict[str, Any]:
    if len(levels) != 3 or any(n < 65 or n % 2 == 0 for n in levels):
        raise ValueError("levels must contain three odd integers >= 65")
    if not (levels[0] < levels[1] < levels[2]):
        raise ValueError("levels must be strictly increasing")
    if sample_count < 3:
        raise ValueError("sample_count must be >= 3")

    repair = KokunoHeatDiscrepancyRepair()
    rng = np.random.default_rng(seed)
    samples: list[dict[str, Any]] = []
    finest_rel_errors: list[float] = []
    refinement_rel_errors: list[float] = []
    jacobian_rel_errors: list[float] = []
    mutation_abs_errors: list[float] = []
    increment_norms: list[float] = []

    for index in range(sample_count):
        eta = float(rng.uniform(-0.72, 0.72))
        f_eta = float(repair.source_f(eta))
        coefficients = rng.uniform(-0.006, 0.006, 3)
        # Avoid an accidentally tiny manufactured target while staying far from
        # the declared +/-0.05 production bound.
        if np.linalg.norm(coefficients) < 0.002:
            coefficients[0] += 0.003

        public = np.asarray(
            repair.normalized_increment(coefficients, f_eta=f_eta), dtype=float
        )
        independent = [
            _independent_increment(repair, coefficients, f_eta, n) for n in levels
        ]
        finest_error = public - independent[-1]
        finest_rel = _norm_relative(finest_error, independent[-1])
        refine_rel = _norm_relative(independent[-1] - independent[-2], independent[-1])

        public_jacobian = np.asarray(
            repair.coefficient_jacobian(coefficients, f_eta=f_eta), dtype=float
        )
        independent_jacobian = _independent_jacobian(
            repair, coefficients, f_eta, levels[-1]
        )
        jac_rel = _norm_relative(public_jacobian - independent_jacobian, independent_jacobian)

        mutated = coefficients.copy()
        mutated[0] += PARAMETER_MUTATION
        mutated_public = np.asarray(
            repair.normalized_increment(mutated, f_eta=f_eta), dtype=float
        )
        mutation_abs = float(np.linalg.norm(mutated_public - independent[-1]))

        finest_rel_errors.append(finest_rel)
        refinement_rel_errors.append(refine_rel)
        jacobian_rel_errors.append(jac_rel)
        mutation_abs_errors.append(mutation_abs)
        increment_norms.append(float(np.linalg.norm(independent[-1])))
        samples.append(
            {
                "index": index,
                "eta": eta,
                "f_eta": f_eta,
                "coefficients": coefficients.tolist(),
                "public_increment": public.tolist(),
                "independent_by_level": [row.tolist() for row in independent],
                "finest_public_vs_independent_relative_error": finest_rel,
                "medium_to_fine_relative_change": refine_rel,
                "analytic_vs_independent_jacobian_relative_error": jac_rel,
                "parameter_mutation_abs_detection": mutation_abs,
            }
        )

    derivative = _profile_derivative_audit(
        repair, np.asarray(samples[0]["coefficients"], dtype=float)
    )
    max_public_rel = float(max(finest_rel_errors))
    max_refinement_rel = float(max(refinement_rel_errors))
    max_jacobian_rel = float(max(jacobian_rel_errors))
    min_mutation_abs = float(min(mutation_abs_errors))
    min_increment_norm = float(min(increment_norms))
    derivative_rel = float(derivative["finest"]["relative_error"])

    structural_preflight_passed = bool(
        max_public_rel <= 2.0e-5
        and max_refinement_rel <= 2.0e-5
        and max_jacobian_rel <= 2.0e-4
        and derivative_rel <= 2.0e-4
        and min_mutation_abs >= 1.0e-6
        and min_increment_norm >= 1.0e-7
    )

    return {
        "task": TASK,
        "upstream_agent1_head": UPSTREAM_HEAD,
        "seed": seed,
        "levels": list(levels),
        "sample_count": sample_count,
        "independent_operator": "uniform-grid composite Simpson + centered parameter/profile finite differences",
        "construction_operator": "upstream fixed Gauss-Legendre quadrature + analytic Jacobian",
        "samples": samples,
        "profile_derivative_audit": derivative,
        "summary": {
            "max_finest_public_vs_independent_relative_error": max_public_rel,
            "max_medium_to_fine_relative_change": max_refinement_rel,
            "max_analytic_vs_independent_jacobian_relative_error": max_jacobian_rel,
            "finest_profile_derivative_relative_error": derivative_rel,
            "min_parameter_mutation_abs_detection": min_mutation_abs,
            "min_independent_increment_norm": min_increment_norm,
            "structural_preflight_passed": structural_preflight_passed,
        },
        "fixed_references": {
            "normalized_ns_residual": RESIDUAL_REFERENCE,
            "divergence": DIVERGENCE_REFERENCE,
            "changed": False,
        },
        "truth_boundary": {
            "actual_heat_discrepancy_supplied": False,
            "global_leading_profile_reconstructed": False,
            "complete_kokuno_composite_velocity": False,
            "formal_full_domain_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e-3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/kokuno_agent4/independent_heat_repair_audit.json"),
    )
    args = parser.parse_args()
    report = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["summary"], indent=2, sort_keys=True))
    if not report["summary"]["structural_preflight_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
