"""Checks for directed partial Section 10.2 bump weights."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

import mpmath as mp

from lei_ren_part1_paper_bump_integral_enclosures import (
    CENTERS,
    RADIUS,
    _record_interval,
)
from lei_ren_part1_paper_bump_partial_enclosures import (
    FULL_MODULE,
    PaperBumpPartialEnclosures,
)


def _endpoints(value: Any) -> tuple[mp.mpf, mp.mpf]:
    raw = getattr(value, "_mpi_", None)
    if raw is None:
        point = mp.mpf(value)
        return point, point
    return mp.make_mpf(raw[0]), mp.make_mpf(raw[1])


def _contains(interval: Any, value: Any) -> bool:
    lower, upper = _endpoints(interval)
    point = mp.mpf(value)
    return lower <= point <= upper


def _raw_bump(t: mp.mpf, multiplicity: int) -> mp.mpf:
    if abs(t) >= 1:
        return mp.mpf(0)
    return mp.exp(-mp.mpf(multiplicity) / (1 - t * t))


def _gauss_integral(
    endpoint: mp.mpf,
    center: mp.mpf,
    power: Fraction,
    multiplicity: int,
    nodes_weights: tuple[Any, Any],
    radius: mp.mpf,
) -> mp.mpf:
    if endpoint <= -1:
        return mp.mpf(0)
    endpoint = min(endpoint, mp.mpf(1))
    nodes, weights = nodes_weights
    exponent = mp.mpf(power.numerator) / power.denominator
    total = mp.mpf(0)
    cuts = [mp.mpf(-1)]
    if -1 < 0 < endpoint:
        cuts.append(mp.mpf(0))
    cuts.append(endpoint)
    for left, right in zip(cuts, cuts[1:]):
        if right <= left:
            continue
        midpoint = (left + right) / 2
        half = (right - left) / 2
        for node, weight in zip(nodes, weights):
            t = midpoint + half * node
            x = center + radius * t
            total += half * weight * x**exponent * _raw_bump(t, multiplicity)
    return total


def _direct_partial(
    x: mp.mpf,
    center: mp.mpf,
    power: Fraction,
    multiplicity: int,
    nodes_weights: tuple[Any, Any],
    radius: mp.mpf,
    normalizer: mp.mpf,
) -> mp.mpf:
    endpoint = (x - center) / radius
    raw = _gauss_integral(endpoint, center, power, multiplicity, nodes_weights, radius)
    if multiplicity == 1:
        return raw / normalizer
    return raw / (radius * normalizer * normalizer)


def _full_saved_equal(receipt: dict[str, Any], key: str, value: Any) -> bool:
    item = receipt["weight_records"][key]["weight_interval"]
    expected = (
        tuple(item["lower_exact_mpf_tuple"]),
        tuple(item["upper_exact_mpf_tuple"]),
    )
    raw = getattr(value, "_mpi_", None)
    return raw is not None and (tuple(raw[0]), tuple(raw[1])) == expected


def run(
    *,
    precision: int = 80,
    order: int = 12,
    tolerance: str = "1e-20",
    output: str | None = None,
) -> dict[str, Any]:
    with mp.workdps(precision + 35):
        partial = PaperBumpPartialEnclosures(
            precision=precision,
            order=order,
            tolerance=tolerance,
        )
        receipt = partial.receipt
        radius = mp.mpf(1) / 40
        rules = mp.gauss_quadrature(256, "legendre")
        normalizer = _gauss_integral(
            mp.mpf(1),
            mp.mpf(5) / 4,
            Fraction(0),
            1,
            rules,
            radius,
        )
        support_checks: dict[str, Any] = {}
        all_contains = True
        all_nonnegative = True
        all_nominal_monotone = True
        all_full_exact = True
        cases = []
        for center_text in CENTERS:
            center = Fraction(center_text)
            left = center - RADIUS
            right = center + RADIUS
            for power, multiplicity in ((Fraction(0), 1), (Fraction(1, 2), 2)):
                key_base = f"{center}|{power}|k{multiplicity}"
                points = (left, center, right)
                values = []
                nominal_values = []
                per_point = {}
                for point in points:
                    value = partial.weight(center, power, point, multiplicity)
                    # Preserve exact support endpoint semantics in the
                    # independent diagnostic.  Converting the rational left
                    # endpoint to an ordinary MP scalar can otherwise move
                    # it inward by one ulp and manufacture an astronomically
                    # tiny fake tail mass.
                    if point <= left:
                        nominal = mp.mpf(0)
                    else:
                        nominal = _direct_partial(
                            mp.mpf(point.numerator) / point.denominator,
                            mp.mpf(center.numerator) / center.denominator,
                            power,
                            multiplicity,
                            rules,
                            radius,
                            normalizer,
                        )
                    status = partial.integral_record(
                        center, power, point, multiplicity
                    )["status"]
                    values.append(value)
                    nominal_values.append(nominal)
                    contains = _contains(value, nominal)
                    all_contains = all_contains and contains
                    lower, _ = _endpoints(value)
                    nonnegative = lower >= 0
                    all_nonnegative = all_nonnegative and nonnegative
                    per_point[str(point)] = {
                        "interval": _record_interval(value),
                        "nominal_order256": mp.nstr(nominal, 60),
                        "contains_nominal_order256": contains,
                        "lower_nonnegative": nonnegative,
                        "status": status,
                    }
                monotone = all(
                    nominal_values[index] <= nominal_values[index + 1]
                    for index in range(len(nominal_values) - 1)
                )
                all_nominal_monotone = all_nominal_monotone and monotone
                full_key = f"{center}|{power}|k{multiplicity}"
                exact_full = _full_saved_equal(
                    receipt, full_key, values[-1]
                )
                all_full_exact = all_full_exact and exact_full
                support_checks[key_base] = {
                    "points": per_point,
                    "nominal_order256_monotone": monotone,
                    "right_endpoint_exact_saved_interval": exact_full,
                }
                cases.append(key_base)

        zero = partial.weight("1.25", "0", "1", 1)
        lower_empty_exact = getattr(zero, "_mpi_", None) == (
            (0, 0, 0, 0),
            (0, 0, 0, 0),
        )
        report = {
            "module_sha256": hashlib.sha256(Path(__file__).with_name("lei_ren_part1_paper_bump_partial_enclosures.py").read_bytes()).hexdigest(),
            "frozen_full_module": FULL_MODULE,
            "frozen_full_module_sha256": receipt.get("module_sha256"),
            "settings": {
                "precision": precision,
                "Taylor_order": order,
                "tolerance": tolerance,
            },
            "normalization_loaded_from_full_receipt": True,
            "normalization_interval": _record_interval(partial.normalization),
            "support_checks": support_checks,
            "cases": cases,
            "all_order256_nominals_contained": all_contains,
            "all_partial_intervals_nonnegative": all_nonnegative,
            "all_nominal_partial_masses_monotone": all_nominal_monotone,
            "all_right_endpoints_exact_saved_intervals": all_full_exact,
            "lower_empty_exact_zero": lower_empty_exact,
            "independent_quadrature_order": 256,
            "nominal_comparisons_are_diagnostic_only": True,
            "source_error_enclosed": False,
            "parameter_error_enclosed": False,
            "cone_certified": False,
        }
        destination = Path(output) if output else Path(__file__).with_suffix(".json")
        destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({
            "all_order256_nominals_contained": all_contains,
            "all_partial_intervals_nonnegative": all_nonnegative,
            "all_nominal_partial_masses_monotone": all_nominal_monotone,
            "all_right_endpoints_exact_saved_intervals": all_full_exact,
            "lower_empty_exact_zero": lower_empty_exact,
            "output": str(destination),
        }, indent=2), flush=True)
        return report


if __name__ == "__main__":
    run()
