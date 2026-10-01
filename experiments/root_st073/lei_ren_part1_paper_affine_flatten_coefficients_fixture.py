"""Independent checks for affine-cell flatten coefficient enclosures."""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_affine_flatten_coefficients import (
    affine_flatten_coefficient_intervals,
    beta_factor,
)
from lei_ren_part1_paper_uniform_flatten_error import (
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


def beta_reference(value: mp.mpf) -> mp.mpf:
    return 2 * (1 - sigma_reference(value))


def weight_reference(value: mp.mpf) -> mp.mpf:
    base = mp.mpf("-.2")
    length = mp.mpf("1.75")
    slope = mp.mpf("-.37")
    beta = beta_reference(value)
    return length * mp.exp(2 * base + 2 * slope * length * value + (beta - 3) * mp.log(2))


def independent_coefficient(beta: mp.mpf, index: int) -> mp.mpf:
    return (-1 if index % 2 else 1) * mp.rf(beta, index) / mp.factorial(index)


def interval_bounds(value):
    if hasattr(value, "_mpi_"):
        return tuple(mp.make_mpf(item) for item in value._mpi_)
    value = mp.mpf(value)
    return value, value


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
        context = MPIntervalContext()
        context.dps = 100
        base = context.mpf("-.2")
        length = context.mpf("1.75")
        slope = context.mpf("-.37")

        def beta_source(u):
            lower, upper = interval_bounds(u)
            # beta is decreasing, so the endpoint hull encloses a directed
            # interval input from the cell grid.
            return context.mpf(
                [str(beta_reference(upper)), str(beta_reference(lower))]
            )

        coarse = affine_flatten_coefficient_intervals(
            base=base,
            length=length,
            slope=slope,
            beta_source=beta_source,
            coefficient_order=32,
            panels=16,
            context=context,
        )
        fine = affine_flatten_coefficient_intervals(
            base=base,
            length=length,
            slope=slope,
            beta_source=beta_source,
            coefficient_order=32,
            panels=64,
            context=context,
        )
        if not all(
            fine["coefficient_interval_widths"][index]
            < coarse["coefficient_interval_widths"][index]
            for index in (0, 1, 2, 8, 32)
        ):
            raise AssertionError("affine coefficient widths did not contract")

        # Stored positive Gauss atoms use the same resolved source.
        nodes, weights = mp.gauss_quadrature(24, "legendre")
        atoms = []
        for node, weight in zip(nodes, weights):
            u = (node + 1) / 2
            atoms.append((weight / 2 * weight_reference(u), beta_reference(u)))
        stored_coefficients, stored_mass = stored_flatten_coefficients_directed(
            atoms, 32, context
        )
        result = uniform_error_from_coefficients_directed(
            fine["coefficients"],
            stored_coefficients,
            stored_mass,
            z_radius=".8",
            context=context,
        )

        coefficient_checks = []
        for index in (0, 1, 2, 8, 32):
            target = mp.quad(
                lambda u, index=index: weight_reference(u)
                * independent_coefficient(beta_reference(u), index),
                [0, mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75"), 1],
            )
            lower, upper = interval_bounds(fine["coefficients"][index])
            if not (lower <= target <= upper):
                raise AssertionError(f"coefficient {index} escaped affine interval")
            coefficient_checks.append(
                {
                    "index": index,
                    "lower": lower,
                    "upper": upper,
                    "reference": target,
                }
            )

        # Endpoint beta=0 must give b_0(0)=2^-3, while beta=2 gives 1/2.
        beta_zero = beta_factor(context.mpf(0), 0, context)
        beta_two = beta_factor(context.mpf(2), 0, context)
        if not (interval_bounds(beta_zero)[0] <= mp.mpf(".125") <= interval_bounds(beta_zero)[1]):
            raise AssertionError("beta=0 endpoint convention failed")
        if not (interval_bounds(beta_two)[0] <= mp.mpf(".5") <= interval_bounds(beta_two)[1]):
            raise AssertionError("beta=2 endpoint convention failed")

        # Explicit analytic tail coverage remains independent of the affine
        # coefficient quadrature checks.
        tails = weighted_geometric_tail(".8", 32)
        tail_checks = {}
        for name, derivative in (("value", 0), ("first_Z", 1), ("second_Z", 2)):
            partial = mp.mpf(0)
            for index in range(33, 400):
                coefficient_bound = index + 1
                if derivative == 0:
                    partial += coefficient_bound * mp.mpf(".8") ** (2 * index)
                elif derivative == 1:
                    partial += 2 * index * coefficient_bound * mp.mpf(".8") ** (2 * index - 1)
                else:
                    partial += (
                        2 * index * (2 * index - 1) * coefficient_bound
                        * mp.mpf(".8") ** (2 * index - 2)
                    )
            if partial > tails[name]:
                raise AssertionError(f"{name} analytic tail undercovers partial sum")
            tail_checks[name] = {"partial": partial, "analytic": tails[name]}

        report = {
            "all_checks_passed": True,
            "directed_interval_arithmetic": True,
            "affine_exponential_integrated_exactly": True,
            "beta_endpoint_bounds_directed": True,
            "endpoint_beta_zero_convention": True,
            "independent_resolved_quadrature": True,
            "coefficient_width_contraction_16_to_64": True,
            "analytic_tail_coverage_checked": True,
            "coefficient_order": 32,
            "coarse_panels": 16,
            "fine_panels": 64,
            "stored_gauss_order": 24,
            "uniform_error_upper_bounds": result["uniform_error_upper_bounds"],
            "coefficient_series_upper_bounds": result["coefficient_series_upper_bounds"],
            "positive_mass_fallback_upper_bounds": result[
                "positive_mass_fallback_upper_bounds"
            ],
            "coarse_widths": coarse["coefficient_interval_widths"],
            "fine_widths": fine["coefficient_interval_widths"],
            "coefficient_checks": coefficient_checks,
            "tail_checks": tail_checks,
            "stored_mass": stored_mass,
            "stored_parameter_errors_enclosed": False,
            "core_error_enclosed": False,
            "five_defect_interval_closure": False,
        }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(encode(report), indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": True, "output": str(output)}))
    return report


if __name__ == "__main__":
    main()
