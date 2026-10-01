"""Independent resolved checks for the uniform flatten error enclosure.

The fixture uses a small resolved schedule written here rather than building
the expensive source profile.  Its weight and beta are analytic, positive,
and monotone on ``[0,1]``.  Independent scalar quadrature checks the true
coefficient integrals and the value/first/second derivative differences at a
few points.  The proof of the reported uniform bounds comes from the
coefficient intervals and geometric tails, not from those point samples.
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_uniform_flatten_error import (
    monotone_coefficient_intervals_directed,
    stored_flatten_coefficients_directed,
    uniform_error_from_coefficients_directed,
    weighted_geometric_tail,
)


def sigma_reference(value: mp.mpf) -> mp.mpf:
    value = mp.mpf(value)
    if value <= 0:
        return mp.mpf(0)
    if value >= 1:
        return mp.mpf(1)
    phase = 1 / (value * value) - 1 / ((1 - value) * (1 - value))
    if phase >= 0:
        tiny = mp.exp(-phase)
        return tiny / (1 + tiny)
    tiny = mp.exp(phase)
    return 1 / (1 + tiny)


def resolved_weight_beta(value: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    """A resolved stored-schedule analogue with beta in [0,2]."""

    value = mp.mpf(value)
    beta = 2 * (1 - sigma_reference(value))
    weight = mp.exp(-mp.mpf("1.4") * value) * mp.power(2, beta - 3)
    return weight, beta


def independent_coefficient(beta: mp.mpf, index: int) -> mp.mpf:
    if index == 0:
        return mp.mpf(1)
    return (-1 if index % 2 else 1) * mp.rf(beta, index) / mp.factorial(index)


def interval_bounds(value):
    if hasattr(value, "_mpi_"):
        return tuple(mp.make_mpf(item) for item in value._mpi_)
    value = mp.mpf(value)
    return value, value


def independent_jet(value: mp.mpf, *, stored: bool, atoms=None) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    value = mp.mpf(value)
    q = 1 + value * value
    if stored:
        rows = atoms
        mass = mp.fsum(
            row[0] * mp.power(q, -row[1])
            for row in rows
        )
        first = mp.fsum(
            row[0]
            * mp.power(q, -row[1])
            * (-2 * row[1] * value / q)
            for row in rows
        )
        second = mp.fsum(
            row[0]
            * mp.power(q, -row[1])
            * (
                -2 * row[1] / q
                + 4 * row[1] * (row[1] + 1) * value * value / (q * q)
            )
            for row in rows
        )
        return mass, first, second

    def integrand(u: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
        weight, beta = resolved_weight_beta(u)
        density = weight * mp.power(q, -beta)
        first = density * (-2 * beta * value / q)
        second = density * (
            -2 * beta / q + 4 * beta * (beta + 1) * value * value / (q * q)
        )
        return density, first, second

    result = []
    for component in range(3):
        result.append(
            mp.quad(
                lambda u: integrand(u)[component],
                [0, mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75"), 1],
            )
        )
    return tuple(result)


def encode(value):
    if hasattr(value, "_mpi_"):
        return [mp.nstr(mp.make_mpf(item), 45) for item in value._mpi_]
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 45)
    if isinstance(value, tuple):
        return [encode(item) for item in value]
    if isinstance(value, list):
        return [encode(item) for item in value]
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    return value


def main() -> dict:
    with mp.workdps(100):
        interval_context = MPIntervalContext()
        interval_context.dps = 100

        def interval_source(value):
            point = interval_context.mpf(str(value))
            sigma = interval_context.mpf(str(sigma_reference(value)))
            beta = 2 * (1 - sigma)
            weight = interval_context.exp(
                -interval_context.mpf("1.4") * point
                + (beta - 3) * interval_context.ln(2)
            )
            return weight, beta

        coefficient_order = 80
        coefficient_data = monotone_coefficient_intervals_directed(
            interval_source,
            0,
            1,
            coefficient_order=coefficient_order,
            panels=2048,
            monotone_decreasing=True,
            context=interval_context,
        )

        nodes, weights = mp.gauss_quadrature(24, "legendre")
        atoms = []
        for node, weight in zip(nodes, weights):
            u = (node + 1) / 2
            mass_weight, beta = resolved_weight_beta(u)
            atoms.append((weight / 2 * mass_weight, beta))
        stored_coefficients, stored_mass = stored_flatten_coefficients_directed(
            atoms, coefficient_order, interval_context
        )
        result = uniform_error_from_coefficients_directed(
            coefficient_data["coefficients"],
            stored_coefficients,
            stored_mass,
            z_radius=mp.mpf(".8"),
            context=interval_context,
        )

        # Independent scalar quadrature checks representative low, middle,
        # and high retained coefficients.  The endpoint enclosure itself is
        # constructed for every coefficient through ``coefficient_order``.
        coefficient_checks = []
        checked_indices = (0, 1, 2, 8, 32, coefficient_order)
        for index in checked_indices:
            interval = coefficient_data["coefficients"][index]
            target = mp.quad(
                lambda u, index=index: resolved_weight_beta(u)[0]
                * independent_coefficient(resolved_weight_beta(u)[1], index),
                [0, mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75"), 1],
            )
            lower, upper = interval_bounds(interval)
            if not (lower <= target <= upper):
                raise AssertionError(
                    f"coefficient {index} escaped interval: {lower} <= {target} <= {upper}"
                )
            coefficient_checks.append(
                {
                    "index": index,
                    "interval_width": upper - lower,
                    "reference": target,
                }
            )

        # Independent scalar jets are a regression check only; the uniform
        # result is already established by the series/tail construction.
        sample_checks = []
        for value in (mp.mpf("-.8"), mp.mpf("-.4"), mp.mpf(0), mp.mpf(".4"), mp.mpf(".8")):
            true_jet = independent_jet(value, stored=False)
            finite_jet = independent_jet(value, stored=True, atoms=atoms)
            differences = [abs(a - b) for a, b in zip(true_jet, finite_jet)]
            names = ("value", "first_Z", "second_Z")
            for name, difference in zip(names, differences):
                if difference > result["uniform_error_upper_bounds"][name]:
                    raise AssertionError(
                        f"sampled {name} difference escaped uniform bound at Z={value}"
                    )
            sample_checks.append(
                {
                    "Z": value,
                    "true_jet": true_jet,
                    "stored_jet": finite_jet,
                    "absolute_difference": differences,
                }
            )

        # Check the analytic tail factors against a long independent partial
        # series.  This explicitly exercises value, first, and second tails.
        radius = mp.mpf(".8")
        tail = weighted_geometric_tail(radius, coefficient_order)
        first = coefficient_order + 1
        tail_checks = {}
        for name, component in (
            ("value", 0),
            ("first_Z", 1),
            ("second_Z", 2),
        ):
            terms = []
            for index in range(first, first + 500):
                if component == 0:
                    term = (index + 1) * radius ** (2 * index)
                elif component == 1:
                    term = 2 * index * (index + 1) * radius ** (2 * index - 1)
                else:
                    term = (
                        2
                        * index
                        * (2 * index - 1)
                        * (index + 1)
                        * radius ** (2 * index - 2)
                    )
                terms.append(term)
            partial = mp.fsum(terms)
            if partial > tail[name]:
                raise AssertionError(f"analytic {name} tail undercovers partial sum")
            tail_checks[name] = {"partial_500": partial, "analytic": tail[name]}

        report = {
            "all_checks_passed": True,
            "uniform_interval_enclosure": True,
            "directed_interval_arithmetic": True,
            "coefficient_integrals_independently_quadrature_checked": True,
            "point_samples_are_regression_only": True,
            "analytic_tail_coverage_checked": True,
            "coefficient_order": coefficient_order,
            "stored_gauss_order": 24,
            "true_coefficient_panels": coefficient_data["panels"],
            "stored_atom_count": len(atoms),
            "uniform_error_upper_bounds": result["uniform_error_upper_bounds"],
            "coefficient_series_upper_bounds": result["coefficient_series_upper_bounds"],
            "positive_mass_fallback_upper_bounds": result[
                "positive_mass_fallback_upper_bounds"
            ],
            "weighted_geometric_tail_factors": result[
                "weighted_geometric_tail_factors"
            ],
            "true_mass_interval": coefficient_data["mass_interval"],
            "stored_mass": stored_mass,
            "coefficient_checks": coefficient_checks,
            "sample_checks": sample_checks,
            "tail_checks": tail_checks,
            "parameter_scope": "independent resolved fixture parameters",
            "original_parameter_errors_enclosed": False,
            "core_error_enclosed": False,
            "five_defect_interval_closure": False,
        }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(encode(report), indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": True, "output": str(output)}))
    return report


if __name__ == "__main__":
    main()
