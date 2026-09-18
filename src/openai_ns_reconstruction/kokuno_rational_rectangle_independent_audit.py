"""Independent Agent-4 audit of Kokuno rational rectangle separation.

This module deliberately does not reuse Agent-2's Fraction-based repeated
``J_g`` torus transform when certifying a rectangle family.  It embeds a fresh
caller-supplied rational family on an integer torus lattice, applies ``J_g`` by
binary modular matrix exponentiation, and recomputes the full finite separation
certificate from that representation.

The audit is structural only.  It does not materialize the still-missing actual
positive-order/background multi-band oscillatory velocity and it does not run or
relax the formal full-domain Navier--Stokes acceptance gates.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import random
from typing import Any

from .kokuno_rational_rectangle_separation import KokunoRationalRectangleSeparation

SCHEMA = "kokuno-agent4-rational-rectangle-independent-audit-v1"
SEED = 9173161
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5

# Frozen before execution.  These are local structural-audit guards, not PDE
# acceptance thresholds.
FROZEN_GUARDS = {
    "source_delta_max_mismatches": 0,
    "exact_minimum_mismatches": 0,
    "C_J_relative_error_max": 5.0e-14,
    "C_v_relative_error_max": 5.0e-15,
    "collision_mutation_misses": 0,
    "oversized_radius_mutation_misses": 0,
}

J_G = ((3, 1), (1, 5))


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _matrix_multiply(
    left: tuple[tuple[int, int], tuple[int, int]],
    right: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[tuple[int, int], tuple[int, int]]:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def _matrix_multiply_mod(
    left: tuple[tuple[int, int], tuple[int, int]],
    right: tuple[tuple[int, int], tuple[int, int]],
    modulus: int,
) -> tuple[tuple[int, int], tuple[int, int]]:
    raw = _matrix_multiply(left, right)
    return (
        (raw[0][0] % modulus, raw[0][1] % modulus),
        (raw[1][0] % modulus, raw[1][1] % modulus),
    )


def _matrix_power(power: int) -> tuple[tuple[int, int], tuple[int, int]]:
    result = ((1, 0), (0, 1))
    base = J_G
    exponent = int(power)
    while exponent:
        if exponent & 1:
            result = _matrix_multiply(result, base)
        base = _matrix_multiply(base, base)
        exponent //= 2
    return result


def _matrix_power_mod(
    power: int, modulus: int
) -> tuple[tuple[int, int], tuple[int, int]]:
    result = ((1, 0), (0, 1))
    base = ((3 % modulus, 1 % modulus), (1 % modulus, 5 % modulus))
    exponent = int(power)
    while exponent:
        if exponent & 1:
            result = _matrix_multiply_mod(result, base, modulus)
        base = _matrix_multiply_mod(base, base, modulus)
        exponent //= 2
    return result


def _apply_power_lattice(
    center: tuple[int, int], power: int, modulus: int
) -> tuple[int, int]:
    matrix = _matrix_power_mod(power, modulus)
    return (
        (matrix[0][0] * center[0] + matrix[0][1] * center[1]) % modulus,
        (matrix[1][0] * center[0] + matrix[1][1] * center[1]) % modulus,
    )


def _wrapped_integer_distance(delta: int, modulus: int) -> int:
    value = delta % modulus
    return min(value, modulus - value)


def _distance_squared_numerator(
    left: tuple[int, int], right: tuple[int, int], modulus: int
) -> int:
    dx = _wrapped_integer_distance(left[0] - right[0], modulus)
    dy = _wrapped_integer_distance(left[1] - right[1], modulus)
    return dx * dx + dy * dy


def _independent_minimum(
    centers: tuple[tuple[int, int], ...], modulus: int, delta_max: int
) -> tuple[Fraction, tuple[int, int, int]]:
    minimum: int | None = None
    argmin: tuple[int, int, int] | None = None
    for mu, left in enumerate(centers):
        for nu, right in enumerate(centers):
            for delta in range(delta_max + 1):
                if delta == 0 and mu == nu:
                    continue
                transformed = _apply_power_lattice(right, delta, modulus)
                value = _distance_squared_numerator(left, transformed, modulus)
                if minimum is None or value < minimum:
                    minimum = value
                    argmin = (mu, nu, delta)
    if minimum is None or argmin is None:
        raise ValueError("independent audit needs at least one non-tautological constraint")
    return Fraction(minimum, modulus * modulus), argmin


def _independently_valid_after_append(
    selected: list[tuple[int, int]],
    candidate: tuple[int, int],
    modulus: int,
    delta_max: int,
) -> bool:
    trial = tuple(selected + [candidate])
    minimum, _ = _independent_minimum(trial, modulus, delta_max)
    return minimum > 0


def _fresh_lattice_centers(
    *, modulus: int, color_count: int, delta_max: int, seed: int
) -> tuple[tuple[int, int], ...]:
    points = [(a, b) for a in range(modulus) for b in range(modulus)]
    random.Random(seed).shuffle(points)
    selected: list[tuple[int, int]] = []
    for candidate in points:
        if _independently_valid_after_append(selected, candidate, modulus, delta_max):
            selected.append(candidate)
            if len(selected) == color_count:
                return tuple(selected)
    raise RuntimeError("fresh independent modular witness search exhausted the lattice")


def _to_fraction_centers(
    centers: tuple[tuple[int, int], ...], modulus: int
) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple(
        (Fraction(first, modulus), Fraction(second, modulus))
        for first, second in centers
    )


def _independent_source_delta_max(*, h: float, ell0: int) -> int:
    t_g = 4.0 + math.sqrt(2.0)
    bound = 1.0 + (
        4.0 * (1.0 + h) * math.log(2.0) + 8.0 / ell0
    ) / math.log(t_g)
    return int(math.ceil(bound))


def _independent_c_v() -> float:
    b_g = math.sqrt(2.0) - 1.0
    return 2.0 * math.sqrt(1.0 + b_g * b_g)


def _independent_c_j(delta_max: int) -> float:
    # Different from Agent 2's T_g**Delta path: construct integer J^Delta and
    # evaluate its spectral norm through the two-by-two A^T A eigenvalue formula.
    matrix = _matrix_power(delta_max)
    a00 = matrix[0][0] * matrix[0][0] + matrix[1][0] * matrix[1][0]
    a01 = matrix[0][0] * matrix[0][1] + matrix[1][0] * matrix[1][1]
    a11 = matrix[0][1] * matrix[0][1] + matrix[1][1] * matrix[1][1]
    trace = float(a00 + a11)
    determinant = float(a00 * a11 - a01 * a01)
    discriminant = max(0.0, trace * trace - 4.0 * determinant)
    largest = 0.5 * (trace + math.sqrt(discriminant))
    return math.sqrt(largest)


def _relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(abs(right), 1.0e-300)


def _public_rejects_collision(
    centers: tuple[tuple[int, int], ...], modulus: int, delta_max: int, r0: float
) -> bool:
    mutated = list(centers)
    mutated[1] = _apply_power_lattice(mutated[0], 1, modulus)
    exact_minimum, _ = _independent_minimum(tuple(mutated), modulus, delta_max)
    if exact_minimum != 0:
        return False
    try:
        KokunoRationalRectangleSeparation(
            centers=_to_fraction_centers(tuple(mutated), modulus),
            delta_max=delta_max,
            r0=r0,
        )
    except ValueError:
        return True
    return False


def _public_rejects_oversized_radius(
    centers: tuple[tuple[int, int], ...],
    modulus: int,
    delta_max: int,
    strict_upper: float,
) -> bool:
    try:
        KokunoRationalRectangleSeparation(
            centers=_to_fraction_centers(centers, modulus),
            delta_max=delta_max,
            r0=strict_upper * 1.000001,
        )
    except ValueError:
        return True
    return False


def run_audit(*, seed: int = SEED) -> dict[str, Any]:
    """Run the frozen fresh-sample structural audit against Agent-2 #439."""
    cases = (
        {"h": 0.0025, "ell0": 12, "color_count": 5, "denominator": 19, "safety": 0.37},
        {"h": 0.0065, "ell0": 19, "color_count": 6, "denominator": 23, "safety": 0.41},
        {"h": 0.0090, "ell0": 37, "color_count": 5, "denominator": 29, "safety": 0.47},
    )

    reports: list[dict[str, Any]] = []
    delta_mismatches = 0
    minimum_mismatches = 0
    collision_misses = 0
    radius_misses = 0
    worst_cj_error = 0.0
    worst_cv_error = 0.0

    for index, spec in enumerate(cases):
        h = float(spec["h"])
        ell0 = int(spec["ell0"])
        color_count = int(spec["color_count"])
        denominator = int(spec["denominator"])
        safety = float(spec["safety"])
        delta_max = _independent_source_delta_max(h=h, ell0=ell0)
        public_delta_max = KokunoRationalRectangleSeparation.source_delta_max(h=h, ell0=ell0)
        delta_mismatches += int(public_delta_max != delta_max)

        centers = _fresh_lattice_centers(
            modulus=denominator,
            color_count=color_count,
            delta_max=delta_max,
            seed=seed + 101 * index,
        )
        exact_minimum, argmin = _independent_minimum(centers, denominator, delta_max)
        if exact_minimum <= 0:
            raise RuntimeError("fresh independent center family unexpectedly collides")

        d_c = math.sqrt(float(exact_minimum))
        c_v = _independent_c_v()
        c_j = _independent_c_j(delta_max)
        injectivity_upper = 1.0 / (4.0 * c_v)
        separation_upper = d_c / (2.0 * c_v * (c_j + 1.0))
        strict_upper = min(injectivity_upper, separation_upper)
        r0 = safety * strict_upper

        certificate = KokunoRationalRectangleSeparation(
            centers=_to_fraction_centers(centers, denominator),
            delta_max=delta_max,
            r0=r0,
            center_binding="agent4_independent_fresh_modular_lattice",
            r0_binding="agent4_preregistered_fraction_of_independent_upper",
        )
        minimum_match = certificate.d_c_squared_exact == exact_minimum
        minimum_mismatches += int(not minimum_match)

        cj_error = _relative_error(certificate.C_J, c_j)
        cv_error = _relative_error(certificate.C_v, c_v)
        worst_cj_error = max(worst_cj_error, cj_error)
        worst_cv_error = max(worst_cv_error, cv_error)

        collision_rejected = _public_rejects_collision(
            centers, denominator, delta_max, max(r0 * 0.5, 1.0e-18)
        )
        oversized_radius_rejected = _public_rejects_oversized_radius(
            centers, denominator, delta_max, strict_upper
        )
        collision_misses += int(not collision_rejected)
        radius_misses += int(not oversized_radius_rejected)

        receipt = certificate.receipt()
        truth = receipt["truth_boundary"]
        reports.append(
            {
                "case": index,
                "h": h,
                "ell0": ell0,
                "color_count": color_count,
                "denominator": denominator,
                "safety": safety,
                "delta_max_independent": delta_max,
                "delta_max_public": public_delta_max,
                "centers_lattice": [list(pair) for pair in centers],
                "independent_d_c_squared": {
                    "numerator": exact_minimum.numerator,
                    "denominator": exact_minimum.denominator,
                },
                "independent_argmin": list(argmin),
                "public_exact_minimum_match": minimum_match,
                "C_J_independent": c_j,
                "C_J_public": certificate.C_J,
                "C_J_relative_error": cj_error,
                "C_v_independent": c_v,
                "C_v_public": certificate.C_v,
                "C_v_relative_error": cv_error,
                "strict_radius_upper_independent": strict_upper,
                "r0": r0,
                "injectivity_left_side": 4.0 * r0 * c_v,
                "separation_left_over_d_c": 2.0 * r0 * c_v * (c_j + 1.0) / d_c,
                "collision_mutation_detected": collision_rejected,
                "oversized_radius_mutation_detected": oversized_radius_rejected,
                "truth_boundary": {
                    "source_rectangle_centers_recovered": truth["source_rectangle_centers_recovered"],
                    "source_rectangle_radius_r0_recovered": truth["source_rectangle_radius_r0_recovered"],
                    "actual_positive_order_background_bound": truth["actual_positive_order_background_bound"],
                    "actual_auxiliary_torus_mode_family_bound": truth["actual_auxiliary_torus_mode_family_bound"],
                    "public_xyz_t_velocity_correction_materialized": truth["public_xyz_t_velocity_correction_materialized"],
                    "formal_full_domain_pde_gate_assessed": truth["formal_full_domain_pde_gate_assessed"],
                    "pde_validated": truth["pde_validated"],
                },
            }
        )

    truth_fail_closed = all(
        all(value is False for value in case["truth_boundary"].values())
        for case in reports
    )
    summary = {
        "source_delta_max_mismatches": delta_mismatches,
        "exact_minimum_mismatches": minimum_mismatches,
        "C_J_relative_error_max": worst_cj_error,
        "C_v_relative_error_max": worst_cv_error,
        "collision_mutation_misses": collision_misses,
        "oversized_radius_mutation_misses": radius_misses,
        "truth_boundary_fail_closed": truth_fail_closed,
    }
    passed = (
        delta_mismatches == FROZEN_GUARDS["source_delta_max_mismatches"]
        and minimum_mismatches == FROZEN_GUARDS["exact_minimum_mismatches"]
        and worst_cj_error <= FROZEN_GUARDS["C_J_relative_error_max"]
        and worst_cv_error <= FROZEN_GUARDS["C_v_relative_error_max"]
        and collision_misses == FROZEN_GUARDS["collision_mutation_misses"]
        and radius_misses == FROZEN_GUARDS["oversized_radius_mutation_misses"]
        and truth_fail_closed
    )
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "seed": seed,
        "independent_operator": "integer-torus lattice + binary modular J_g exponentiation + direct 2x2 spectral norm",
        "fresh_center_generation": "shuffled modular lattice; does not call Agent-2 autonomous_rational_witness",
        "frozen_local_guards": dict(FROZEN_GUARDS),
        "formal_project_gates_unchanged": {
            "normalized_momentum_max_L2": FORMAL_MOMENTUM_GATE,
            "divergence_max_L2": FORMAL_DIVERGENCE_GATE,
        },
        "cases": reports,
        "summary": summary,
        "local_structural_preflight_passed": passed,
        "formal_full_domain_pde_gate_assessed": False,
        "pde_validated": False,
        "limitations": [
            "No actual positive-order/background multi-band oscillatory velocity exists in this ancestry.",
            "No public leading+oscillatory+correction composite is evaluated here.",
            "These separation-certificate diagnostics are not Navier-Stokes residuals and are not comparable to ST006 momentum metrics.",
        ],
    }
    payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="artifacts/constrained/kokuno_agent4_rational_rectangle_independent_audit.json",
    )
    args = parser.parse_args()
    report = run_audit()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))
    print(f"local_structural_preflight_passed={report['local_structural_preflight_passed']}")
    if not report["local_structural_preflight_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
