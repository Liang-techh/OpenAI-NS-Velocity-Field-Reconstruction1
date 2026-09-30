"""Independent checks for the Section 7 algebraic closure nuclei.

These checks replay the source equations from scratch, perturb the supplied
discrepancies/energy target for finite-difference derivative checks, and
exercise unreachable or degenerate inputs.  They intentionally use
manufactured source-shaped weights; no actual outer profile is bound here.
"""

from __future__ import annotations

import json
import math
from decimal import Decimal
from pathlib import Path
import sys
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
LOCAL = Path(__file__).resolve().parent
if str(LOCAL) not in sys.path:
    sys.path.insert(0, str(LOCAL))

from lei_ren_part1_paper_outer_closure import (  # noqa: E402
    ANGULAR_EQUATION,
    AXIAL_ENERGY_EQUATION,
    AXIAL_LINEAR_EQUATION,
    SOURCE,
    SOURCE_SECTIONS,
    SOURCE_VERSION,
    angular_equation_residuals,
    paper_axial_weights,
    paper_bump_weights,
    assemble_stable_paper_axial_affine,
    stable_paper_axial_weights,
    solve_angular_bumps_from_paper_bump,
    solve_angular_bumps,
    solve_axial_affine,
    solve_axial_energy_root,
)


REPORT_PATH = LOCAL / "lei_ren_part1_paper_outer_closure_checks.json"
FLOAT_RHS_CUTOFF = 1.0e-12


def _relative(observed: float, expected: float) -> float:
    return abs(float(observed) - float(expected)) / max(abs(float(expected)), 1.0e-30)


def _angular_data() -> dict[str, Any]:
    bump = paper_bump_weights(mu=0.01, ell=3.0 / 20.0, quadrature_order=128)
    mu = bump.mu
    A_mu, B_mu, D_mu = bump.A_mu, bump.B_mu, bump.D_mu
    d1_true, d2_true = 0.02, -0.01
    r = A_mu * (
        math.exp(-3.0 * (1.0 - mu)) * d1_true
        + math.exp(-(1.0 - mu)) * d2_true
    )
    s_H = B_mu * (
        math.exp(3.0 * (1.0 + 2.0 * mu)) * d1_true
        + math.exp(1.0 + 2.0 * mu) * d2_true
    ) + 0.5 * D_mu * (
        math.exp(3.0 * (1.0 + 2.0 * mu)) * d1_true * d1_true
        + math.exp(1.0 + 2.0 * mu) * d2_true * d2_true
    )
    return {
        "mu": mu,
        "A_mu": A_mu,
        "B_mu": B_mu,
        "D_mu": D_mu,
        "r": r,
        "s_H": s_H,
        "d1_true": d1_true,
        "d2_true": d2_true,
        "ell": bump.ell,
        "quadrature_order": bump.quadrature_order,
        "normalization_integral": bump.normalization_integral,
        "beta_inf": bump.beta_inf,
        "beta_l2_squared": bump.beta_l2_squared,
    }


def _run_angular_checks() -> dict[str, Any]:
    data = _angular_data()
    solution = solve_angular_bumps_from_paper_bump(
        mu=data["mu"],
        ell=data["ell"],
        quadrature_order=data["quadrature_order"],
        r=data["r"],
        s_H=data["s_H"],
        require_positive=True,
    )
    bump64 = paper_bump_weights(mu=data["mu"], ell=data["ell"], quadrature_order=64)
    bump256 = paper_bump_weights(mu=data["mu"], ell=data["ell"], quadrature_order=256)
    quadrature_refinement = {
        "A_mu_abs_difference_64_256": abs(bump64.A_mu - bump256.A_mu),
        "B_mu_abs_difference_64_256": abs(bump64.B_mu - bump256.B_mu),
        "D_mu_abs_difference_64_256": abs(bump64.D_mu - bump256.D_mu),
        "normalization_abs_difference_64_256": abs(
            bump64.normalization_integral - bump256.normalization_integral
        ),
    }
    residuals = angular_equation_residuals(solution)
    # Replay (7.21), rather than the normalized equations used by the solver.
    r_replayed = data["A_mu"] * (
        math.exp(-3.0 * (1.0 - data["mu"])) * solution.d1
        + math.exp(-(1.0 - data["mu"])) * solution.d2
    )
    s_replayed = data["B_mu"] * (
        math.exp(3.0 * (1.0 + 2.0 * data["mu"])) * solution.d1
        + math.exp(1.0 + 2.0 * data["mu"]) * solution.d2
    ) + 0.5 * data["D_mu"] * (
        math.exp(3.0 * (1.0 + 2.0 * data["mu"])) * solution.d1 * solution.d1
        + math.exp(1.0 + 2.0 * data["mu"]) * solution.d2 * solution.d2
    )
    h_r = 1.0e-8
    h_s = 1.0e-7
    plus_r = solve_angular_bumps(
        mu=data["mu"], A_mu=data["A_mu"], B_mu=data["B_mu"], D_mu=data["D_mu"],
        r=data["r"] + h_r, s_H=data["s_H"], beta_inf=data["beta_inf"], require_positive=True,
    )
    minus_r = solve_angular_bumps(
        mu=data["mu"], A_mu=data["A_mu"], B_mu=data["B_mu"], D_mu=data["D_mu"],
        r=data["r"] - h_r, s_H=data["s_H"], beta_inf=data["beta_inf"], require_positive=True,
    )
    plus_s = solve_angular_bumps(
        mu=data["mu"], A_mu=data["A_mu"], B_mu=data["B_mu"], D_mu=data["D_mu"],
        r=data["r"], s_H=data["s_H"] + h_s, beta_inf=data["beta_inf"], require_positive=True,
    )
    minus_s = solve_angular_bumps(
        mu=data["mu"], A_mu=data["A_mu"], B_mu=data["B_mu"], D_mu=data["D_mu"],
        r=data["r"], s_H=data["s_H"] - h_s, beta_inf=data["beta_inf"], require_positive=True,
    )
    derivative_r = solution.implicit_derivative(dr=1.0)
    derivative_s = solution.implicit_derivative(ds_H=1.0)
    fd_r = np.array([
        (plus_r.d1 - minus_r.d1) / (2.0 * h_r),
        (plus_r.d2 - minus_r.d2) / (2.0 * h_r),
    ])
    fd_s = np.array([
        (plus_s.d1 - minus_s.d1) / (2.0 * h_s),
        (plus_s.d2 - minus_s.d2) / (2.0 * h_s),
    ])
    analytic_r = np.array([derivative_r.dd1, derivative_r.dd2])
    analytic_s = np.array([derivative_s.dd1, derivative_s.dd2])

    unreachable_rejected = False
    unreachable_error = ""
    # sigma=-1/2 is below the finite minimum of the quadratic for this
    # manufactured rho and therefore has negative discriminant.
    bad_s = -0.5 * data["B_mu"] * math.exp(3.0 * (1.0 + 2.0 * data["mu"]))
    try:
        solve_angular_bumps(
            mu=data["mu"], A_mu=data["A_mu"], B_mu=data["B_mu"], D_mu=data["D_mu"],
            r=0.0, s_H=bad_s,
        )
    except ValueError as exc:
        unreachable_rejected = True
        unreachable_error = str(exc)

    positivity_rejected = False
    positivity_error = ""
    # This exact data set has a negative small-branch multiplier relative to
    # the independently measured bump norm.
    bad_d1, bad_d2 = -0.2, 0.0
    bad_r = data["A_mu"] * (
        math.exp(-3.0 * (1.0 - data["mu"])) * bad_d1
        + math.exp(-(1.0 - data["mu"])) * bad_d2
    )
    bad_s_H = data["B_mu"] * math.exp(3.0 * (1.0 + 2.0 * data["mu"])) * bad_d1 + 0.5 * data["D_mu"] * math.exp(3.0 * (1.0 + 2.0 * data["mu"])) * bad_d1 * bad_d1
    try:
        solve_angular_bumps(
            mu=data["mu"], A_mu=data["A_mu"], B_mu=data["B_mu"], D_mu=data["D_mu"],
            r=bad_r, s_H=bad_s_H, beta_inf=data["beta_inf"], require_positive=True,
        )
    except ValueError as exc:
        positivity_rejected = True
        positivity_error = str(exc)

    degenerate_rejected = False
    degenerate_error = ""
    try:
        solve_angular_bumps(
            mu=data["mu"], A_mu=0.0, B_mu=data["B_mu"], D_mu=data["D_mu"],
            r=data["r"], s_H=data["s_H"],
        )
    except ValueError as exc:
        degenerate_rejected = True
        degenerate_error = str(exc)

    derivative_error_r = float(np.max(np.abs(fd_r - analytic_r)))
    derivative_error_s = float(np.max(np.abs(fd_s - analytic_s)))
    replay_error = max(abs(r_replayed - data["r"]), abs(s_replayed - data["s_H"]))
    normalized_error = max(abs(value) for value in residuals)
    passed = (
        replay_error < 1.0e-12
        and normalized_error < 1.0e-12
        and derivative_error_r < 1.0e-8
        and derivative_error_s < 1.0e-8
        and solution.positive is True
        # The compact C-infinity endpoint produces a finite high-order
        # quadrature tail; the observed 64-to-256 differences are retained
        # rather than hidden behind an exact-integration claim.
        and max(quadrature_refinement.values()) < 1.0e-9
        and unreachable_rejected
        and positivity_rejected
        and degenerate_rejected
    )
    return {
        "passed": bool(passed),
        "inputs": data,
        "paper_bump": {
            "source_definition": "(4.1), (4.13)",
            "support": [-data["ell"], data["ell"]],
            "quadrature_refinement_64_to_256": quadrature_refinement,
        },
        "solution": {
            "d1": solution.d1,
            "d2": solution.d2,
            "branch": solution.branch,
            "discriminant": solution.discriminant,
            "linear_determinant": solution.linear_determinant,
            "scalar_derivative": solution.scalar_derivative,
            "jacobian_condition": solution.jacobian_condition,
            "positive": solution.positive,
            "residual_equations": list(solution.residual_equations),
        },
        "equation_replay": {
            "r_replayed": r_replayed,
            "s_H_replayed": s_replayed,
            "max_absolute_error": replay_error,
            "normalized_residual_max": normalized_error,
        },
        "implicit_derivative": {
            "fd_dr": fd_r.tolist(),
            "analytic_dr": analytic_r.tolist(),
            "max_absolute_error_dr": derivative_error_r,
            "fd_ds_H": fd_s.tolist(),
            "analytic_ds_H": analytic_s.tolist(),
            "max_absolute_error_ds_H": derivative_error_s,
        },
        "rejection_checks": {
            "unreachable_rejected": unreachable_rejected,
            "unreachable_error": unreachable_error,
            "positivity_rejected": positivity_rejected,
            "positivity_error": positivity_error,
            "degenerate_rejected": degenerate_rejected,
            "degenerate_error": degenerate_error,
        },
    }


def _stable_decimal_fixture(mu: float) -> tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]:
    """Fixed Decimal source-row data with differences proportional to mu."""

    mu_decimal = Decimal(str(mu))
    base0 = Decimal("0.0002000000000000000000000000000000000001")
    pulse0 = Decimal("-0.0000170000000000000000000000000000000001")
    base1 = base0 - Decimal("1.25") * mu_decimal
    pulse1 = pulse0 + Decimal("0.75") * mu_decimal
    return (base0, base1), (pulse0, pulse1)


def _run_stable_axial_checks() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    all_passed = True
    for mu in (1.0e-8, 1.0e-28):
        base, pulse = _stable_decimal_fixture(mu)
        assembly = assemble_stable_paper_axial_affine(
            mu=mu,
            base=base,
            pulse=pulse,
            ell=3.0 / 20.0,
            quadrature_order=128,
        )
        weights = assembly.weights
        c = np.asarray(assembly.affine.coefficients(1.07), dtype=float)
        scaled_residual = np.asarray(weights.matrix, dtype=float) @ c
        scaled_residual += np.asarray(assembly.scaled_base, dtype=float)
        scaled_residual += 1.07 * np.asarray(assembly.scaled_pulse, dtype=float)
        expected_scaled_base = np.array([0.0002, 1.25], dtype=float)
        expected_scaled_pulse = np.array([-0.000017, -0.75], dtype=float)
        rhs_base_error = float(np.max(np.abs(np.asarray(assembly.scaled_base) - expected_scaled_base)))
        rhs_pulse_error = float(np.max(np.abs(np.asarray(assembly.scaled_pulse) - expected_scaled_pulse)))
        row = {
            "mu": mu,
            "base_decimal": [str(value) for value in base],
            "pulse_decimal": [str(value) for value in pulse],
            "lambda_values": list(weights.lambda_values),
            "matrix": [list(values) for values in weights.matrix],
            "normalized_determinant": weights.normalized_determinant,
            "source_determinant_estimate": weights.source_determinant_estimate,
            "source_determinant_over_mu": weights.source_determinant_over_mu,
            "condition": weights.condition,
            "scaled_base": list(assembly.scaled_base),
            "scaled_pulse": list(assembly.scaled_pulse),
            "rhs_exact": assembly.rhs_exact,
            "scaled_residual_at_a_1_07": scaled_residual.tolist(),
            "scaled_residual_max": float(np.max(np.abs(scaled_residual))),
            "rhs_base_fixture_error": rhs_base_error,
            "rhs_pulse_fixture_error": rhs_pulse_error,
            "gamma_offset_centers": list(weights.gamma_offset_centers),
            "gamma_offset_supports": [list(value) for value in weights.gamma_offset_supports],
            "absolute_geometry_supported": weights.absolute_geometry_supported,
            "absolute_gamma_centers": None if weights.gamma_centers is None else list(weights.gamma_centers),
            "absolute_gamma_supports": None if weights.gamma_supports is None else [list(value) for value in weights.gamma_supports],
        }
        if mu >= FLOAT_RHS_CUTOFF:
            direct = paper_axial_weights(mu=mu, ell=weights.ell, quadrature_order=128)
            direct_matrix = np.asarray(direct.matrix, dtype=float)
            stable_matrix = np.asarray(weights.matrix, dtype=float)
            direct_source_residual = direct_matrix @ c
            direct_source_residual += np.asarray([float(value) for value in base])
            direct_source_residual += 1.07 * np.asarray([float(value) for value in pulse])
            row["direct_matrix_row0_error"] = float(np.max(np.abs(stable_matrix[0] - direct_matrix[0])))
            row["direct_divided_row_replay_error"] = float(
                np.max(np.abs(stable_matrix[1] - (direct_matrix[0] - direct_matrix[1]) / mu))
            )
            row["direct_source_residual_max"] = float(np.max(np.abs(direct_source_residual)))
            row["source_determinant_relation_error"] = abs(
                np.linalg.det(direct_matrix) + mu * weights.normalized_determinant
            )
            row["source_rows_replay"] = "ordinary float rows available"
            row_passed = (
                row["direct_source_residual_max"] < 1.0e-12
                and row["source_determinant_relation_error"] < 1.0e-14
            )
        else:
            row["direct_matrix_row0_error"] = None
            row["direct_divided_row_replay_error"] = None
            row["direct_source_residual_max"] = None
            row["source_determinant_relation_error"] = None
            row["source_rows_replay"] = "scaled residual only; absolute float supports intentionally unavailable"
            row_passed = not weights.absolute_geometry_supported
        row_passed = row_passed and row["scaled_residual_max"] < 1.0e-12
        row_passed = row_passed and weights.condition < 1.0e12
        row_passed = row_passed and rhs_base_error < 1.0e-10 and rhs_pulse_error < 1.0e-10
        row["passed"] = bool(row_passed)
        rows.append(row)
        all_passed = all_passed and row_passed

    rejected_float_rhs = False
    float_rhs_error = ""
    try:
        assemble_stable_paper_axial_affine(
            mu=1.0e-28,
            base=[0.0002, 0.0002],
            pulse=[-0.000017, -0.000017],
            quadrature_order=128,
        )
    except ValueError as exc:
        rejected_float_rhs = True
        float_rhs_error = str(exc)
    all_passed = all_passed and rejected_float_rhs
    return {
        "passed": bool(all_passed),
        "rows": rows,
        "float_rhs_rejected_below_cutoff": rejected_float_rhs,
        "float_rhs_error": float_rhs_error,
        "cutoff": FLOAT_RHS_CUTOFF,
        "fixture": "fixed Decimal source-row targets; no velocity or residual fitting",
    }


def _run_axial_checks() -> dict[str, Any]:
    # The matrix and Gram weights are actual translated beta_ell quadratures
    # from (7.31)--(7.34).  Base and pulse are manufactured target vectors
    # because the source outer profile is not bound here.
    source_weights = paper_axial_weights(mu=0.01, ell=3.0 / 20.0, quadrature_order=128)
    source_weights64 = paper_axial_weights(mu=0.01, ell=3.0 / 20.0, quadrature_order=64)
    source_weights256 = paper_axial_weights(mu=0.01, ell=3.0 / 20.0, quadrature_order=256)
    mu = source_weights.mu
    matrix = np.asarray(source_weights.matrix, dtype=float)
    matrix_refinement = float(
        np.max(
            np.abs(
                np.asarray(source_weights64.matrix)
                - np.asarray(source_weights256.matrix)
            )
        )
    )
    gram_refinement = float(
        np.max(
            np.abs(
                np.asarray(source_weights64.K_bump)
                - np.asarray(source_weights256.K_bump)
            )
        )
    )
    base = np.array([2.0e-4, -1.5e-4])
    pulse = np.array([1.7e-5, -2.4e-5])
    affine = solve_axial_affine(matrix=matrix, base=base, pulse=pulse)
    affine_rows = []
    for amplitude in (0.9, 1.05, 1.2):
        residual = affine.residual(amplitude)
        affine_rows.append({
            "a": amplitude,
            "c": list(affine.coefficients(amplitude)),
            "residual": list(residual),
            "max_absolute_residual": max(abs(value) for value in residual),
        })
    K_p = 0.2445
    K_bump = np.asarray(source_weights.K_bump, dtype=float)
    a_true = 1.07
    c_true = np.asarray(affine.coefficients(a_true))
    energy_target = K_p * a_true * a_true + mu * float(np.sum(K_bump * c_true * c_true))
    energy = solve_axial_energy_root(
        affine,
        K_p=K_p,
        K_bump=K_bump,
        mu=mu,
        energy_target=energy_target,
    )
    qa, qb, qc = energy.quadratic_coefficients
    independent_energy = K_p * energy.a_p * energy.a_p + mu * float(
        np.sum(K_bump * np.asarray(energy.c_at_root) ** 2)
    ) - energy_target
    h = 1.0e-7
    plus = solve_axial_energy_root(
        affine, K_p=K_p, K_bump=K_bump, mu=mu, energy_target=energy_target + h
    )
    minus = solve_axial_energy_root(
        affine, K_p=K_p, K_bump=K_bump, mu=mu, energy_target=energy_target - h
    )
    fd_target = (plus.a_p - minus.a_p) / (2.0 * h)
    implicit_target = energy.implicit_derivative(d_energy_target=1.0)
    derivative_error = abs(fd_target - implicit_target)

    unreachable_rejected = False
    unreachable_error = ""
    try:
        solve_axial_energy_root(
            affine, K_p=K_p, K_bump=K_bump, mu=mu, energy_target=0.0
        )
    except ValueError as exc:
        unreachable_rejected = True
        unreachable_error = str(exc)
    degenerate_rejected = False
    degenerate_error = ""
    try:
        solve_axial_affine(matrix=np.zeros((2, 2)), base=base, pulse=pulse)
    except ValueError as exc:
        degenerate_rejected = True
        degenerate_error = str(exc)
    energy_degenerate_rejected = False
    energy_degenerate_error = ""
    try:
        solve_axial_energy_root(
            affine, K_p=0.0, K_bump=np.zeros(2), mu=mu, energy_target=energy_target
        )
    except ValueError as exc:
        energy_degenerate_rejected = True
        energy_degenerate_error = str(exc)

    affine_error = max(row["max_absolute_residual"] for row in affine_rows)
    passed = (
        affine_error < 1.0e-15
        and abs(independent_energy) < 1.0e-15
        and abs(energy.residual()) < 1.0e-15
        and derivative_error < 1.0e-8
        and 0.9 <= energy.a_p <= 1.2
        and matrix_refinement < 1.0e-12
        and gram_refinement < 1.0e-20
        and unreachable_rejected
        and degenerate_rejected
        and energy_degenerate_rejected
    )
    return {
        "passed": bool(passed),
        "inputs": {
            "mu": mu,
            "matrix": matrix.tolist(),
            "base": base.tolist(),
            "pulse": pulse.tolist(),
            "K_p": K_p,
            "K_bump": K_bump.tolist(),
            "energy_target": energy_target,
            "ell": source_weights.ell,
            "quadrature_order": source_weights.quadrature_order,
            "lambda_values": list(source_weights.lambda_values),
            "gamma_centers": list(source_weights.gamma_centers),
            "gamma_supports": [list(value) for value in source_weights.gamma_supports],
            "beta_l2_squared": source_weights.beta_l2_squared,
            "matrix_refinement_64_to_256_max_abs": matrix_refinement,
            "K_bump_refinement_64_to_256_max_abs": gram_refinement,
        },
        "affine": {
            "u": list(affine.u),
            "v": list(affine.v),
            "determinant": affine.determinant,
            "condition": affine.condition,
            "source_sign": affine.source_sign,
            "rows": affine_rows,
            "max_absolute_residual": affine_error,
        },
        "energy": {
            "a_true": a_true,
            "a_p": energy.a_p,
            "c_at_root": list(energy.c_at_root),
            "quadratic_coefficients": [qa, qb, qc],
            "roots": list(energy.roots),
            "endpoint_values": list(energy.endpoint_values),
            "derivative_at_root": energy.derivative_at_root,
            "independent_equation_residual": independent_energy,
            "solver_residual": energy.residual(),
        },
        "implicit_derivative": {
            "finite_difference_d_a_d_energy_target": fd_target,
            "analytic_d_a_d_energy_target": implicit_target,
            "absolute_error": derivative_error,
        },
        "rejection_checks": {
            "unreachable_rejected": unreachable_rejected,
            "unreachable_error": unreachable_error,
            "affine_degenerate_rejected": degenerate_rejected,
            "affine_degenerate_error": degenerate_error,
            "energy_degenerate_rejected": energy_degenerate_rejected,
            "energy_degenerate_error": energy_degenerate_error,
        },
    }


def main() -> int:
    angular = _run_angular_checks()
    axial = _run_axial_checks()
    stable_axial = _run_stable_axial_checks()
    report = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_sections": SOURCE_SECTIONS,
        "source_definitions": {
            "normalized_bump": "(4.1), with outer choice ell=3/20 in (4.13)",
            "gamma_end_bumps": "(4.12), (7.28)--(7.29), centers t=13/mu-3 and 13/mu-1",
        },
        "equations": {
            "angular": ANGULAR_EQUATION,
            "axial_linear": AXIAL_LINEAR_EQUATION,
            "axial_energy": AXIAL_ENERGY_EQUATION,
        },
        "scope": {
            "actual_outer_closed": False,
            "actual_source_outer_moment_data_bound": False,
            "waiting_length_root_bound": False,
            "paper_outer_profile_claim": False,
            "checks_are_manufactured_math_nucleus_only": True,
            "open_dependencies": [
                "Bind A_mu, B_mu, D_mu, r(Z), and s_H(Z) from the actual source outer profile.",
                "Bind (7.31) matrix/base/pulse and (7.34) Gram/target data from the same profile.",
                "Resolve the source waiting-length root and outer interval joins before any global claim.",
            ],
        },
        "angular": angular,
        "axial": axial,
        "stable_axial": stable_axial,
        "checks_passed": bool(angular["passed"] and axial["passed"] and stable_axial["passed"]),
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "report": str(REPORT_PATH),
        "checks_passed": report["checks_passed"],
        "angular_passed": angular["passed"],
        "axial_passed": axial["passed"],
        "stable_axial_passed": stable_axial["passed"],
        "angular_residual": angular["equation_replay"]["max_absolute_error"],
        "angular_derivative_error": max(
            angular["implicit_derivative"]["max_absolute_error_dr"],
            angular["implicit_derivative"]["max_absolute_error_ds_H"],
        ),
        "axial_affine_residual": axial["affine"]["max_absolute_residual"],
        "axial_energy_residual": axial["energy"]["independent_equation_residual"],
        "axial_derivative_error": axial["implicit_derivative"]["absolute_error"],
    }, indent=2))
    return 0 if report["checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
