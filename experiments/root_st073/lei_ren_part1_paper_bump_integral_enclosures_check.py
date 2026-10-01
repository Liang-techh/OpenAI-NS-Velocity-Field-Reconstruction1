"""Bounded checks for the directed Section 10.2 bump integrals."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

import mpmath as mp

from lei_ren_part1_paper_bump_integral_enclosures import (
    CENTERS,
    RADIUS,
    PaperBumpIntegralEnclosures,
    _record_interval,
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


def _encode(value: Any) -> Any:
    if getattr(value, "_mpi_", None) is not None:
        return _record_interval(value)
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 60)
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_encode(item) for item in value]
    return value


def _raw_bump(t: mp.mpf, multiplicity: int) -> mp.mpf:
    if abs(t) >= 1:
        return mp.mpf(0)
    return mp.exp(-mp.mpf(multiplicity) / (1 - t * t))


def _gauss_integral(
    center: mp.mpf,
    power: Fraction,
    multiplicity: int,
    order: int,
    radius: mp.mpf,
    nodes_weights: tuple[Any, Any] | None = None,
) -> mp.mpf:
    nodes, weights = (
        nodes_weights
        if nodes_weights is not None
        else mp.gauss_quadrature(order, "legendre")
    )
    exponent = mp.mpf(power.numerator) / power.denominator
    total = mp.mpf(0)
    # Splitting at zero avoids relying on a single rule's symmetry when the
    # endpoint-flat bump is evaluated at finite precision.
    for left, right in ((mp.mpf(-1), mp.mpf(0)), (mp.mpf(0), mp.mpf(1))):
        midpoint = (left + right) / 2
        half = (right - left) / 2
        for node, weight in zip(nodes, weights):
            t = midpoint + half * node
            x = center + radius * t
            total += half * weight * x**exponent * _raw_bump(t, multiplicity)
    return total


def _gauss_weight(
    center: mp.mpf,
    power: Fraction,
    multiplicity: int,
    order: int,
    radius: mp.mpf,
    normalizer: mp.mpf,
    nodes_weights: tuple[Any, Any] | None = None,
) -> mp.mpf:
    raw = _gauss_integral(
        center, power, multiplicity, order, radius, nodes_weights
    )
    if multiplicity == 1:
        return raw / normalizer
    return raw / (radius * normalizer * normalizer)


def run(
    *,
    precision: int = 80,
    order: int = 12,
    tolerance: str = "1e-20",
    initial_cells: int = 16,
    max_depth: int = 24,
    output: str | None = None,
) -> dict[str, Any]:
    with mp.workdps(precision + 35):
        enclosure = PaperBumpIntegralEnclosures(
            precision=precision,
            order=order,
            tolerance=tolerance,
            initial_cells=initial_cells,
            max_depth=max_depth,
        )
        centers = tuple(Fraction(value) for value in CENTERS)
        # Every nonzero scalar in the Section 10.2 linear and quadratic map is
        # represented once here.  Cross-support products are exact zeros and
        # need no integration call.
        requested = []
        for center in centers:
            requested.extend(
                (center, power, 1)
                for power in (Fraction(0), Fraction(1, 2), Fraction(3, 5), Fraction(1, 10), Fraction(-9, 10))
            )
            requested.extend(
                (center, power, 2)
                for power in (Fraction(0), Fraction(1, 2), Fraction(-1))
            )
        unique = []
        seen = set()
        for item in requested:
            if item not in seen:
                seen.add(item)
                unique.append(item)
        weight_records = {}
        for center, power, multiplicity in unique:
            value = enclosure.weight(center, power, multiplicity)
            key = f"{center}|{power}|k{multiplicity}"
            weight_records[key] = {
                "center": str(center),
                "power": str(power),
                "multiplicity": multiplicity,
                "weight_interval": _record_interval(value),
                "raw_integral": enclosure.integral_record(center, power, multiplicity),
            }

        # Independent nominal comparisons.  These are diagnostics only: an
        # inaccurate finite Gauss order is recorded as an escape rather than
        # being promoted to an error in the directed enclosure.
        radius = mp.mpf(1) / 40
        gauss_rules = {
            order: mp.gauss_quadrature(order, "legendre")
            for order in (128, 256)
        }
        normalizer128 = _gauss_integral(
            mp.mpf(5) / 4,
            Fraction(0),
            1,
            128,
            radius,
            gauss_rules[128],
        )
        normalizer256 = _gauss_integral(
            mp.mpf(5) / 4,
            Fraction(0),
            1,
            256,
            radius,
            gauss_rules[256],
        )
        quadrature_checks = {
            "normalization_order128": mp.nstr(normalizer128, 60),
            "normalization_order256": mp.nstr(normalizer256, 60),
            "normalization_order128_within": _contains(enclosure.normalization, normalizer128),
            "normalization_order256_within": _contains(enclosure.normalization, normalizer256),
            "weights": {},
        }
        for center, power, multiplicity in unique:
            key = f"{center}|{power}|k{multiplicity}"
            nominal128 = _gauss_weight(
                mp.mpf(center.numerator) / center.denominator,
                power,
                multiplicity,
                128,
                radius,
                normalizer128,
                gauss_rules[128],
            )
            nominal256 = _gauss_weight(
                mp.mpf(center.numerator) / center.denominator,
                power,
                multiplicity,
                256,
                radius,
                normalizer256,
                gauss_rules[256],
            )
            interval = enclosure.weight(center, power, multiplicity)
            quadrature_checks["weights"][key] = {
                "order128": mp.nstr(nominal128, 60),
                "order256": mp.nstr(nominal256, 60),
                "order128_within": _contains(interval, nominal128),
                "order256_within": _contains(interval, nominal256),
            }

        normalization_checks = {
            str(center): _contains(enclosure.weight(center, Fraction(0), 1), 1)
            for center in centers
        }
        map_data = enclosure.map_enclosure()
        report = {
            "module_sha256": hashlib.sha256(Path(__file__).with_name("lei_ren_part1_paper_bump_integral_enclosures.py").read_bytes()).hexdigest(),
            "check_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "settings": {
                "precision": precision,
                "Taylor_order": order,
                "tolerance": tolerance,
                "initial_cells": initial_cells,
                "max_depth": max_depth,
            },
            "normalization": _record_interval(enclosure.normalization),
            "normalization_contains_one_for_all_centers": normalization_checks,
            "weight_records": weight_records,
            "map_enclosure": map_data,
            "quadrature_comparison": quadrature_checks,
            "metadata": enclosure.metadata,
            "directed_integrals_certified": True,
            "quadrature_comparison_is_diagnostic_only": True,
            "source_error_enclosed": False,
            "parameter_error_enclosed": False,
            "moment_defect_or_cone_certificate": False,
        }
        destination = Path(output) if output else Path(__file__).with_suffix(".json")
        destination.write_text(json.dumps(_encode(report), indent=2) + "\n", encoding="utf-8")
        print(json.dumps(_encode({
            "normalization": report["normalization"],
            "normalization_contains_one_for_all_centers": normalization_checks,
            "quadrature_comparison": quadrature_checks,
            "output": str(destination),
        }), indent=2), flush=True)
        return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--precision", type=int, default=80)
    parser.add_argument("--order", type=int, default=12)
    parser.add_argument("--tolerance", default="1e-20")
    parser.add_argument("--initial-cells", type=int, default=16)
    parser.add_argument("--max-depth", type=int, default=24)
    parser.add_argument("--output")
    args = parser.parse_args()
    run(
        precision=args.precision,
        order=args.order,
        tolerance=args.tolerance,
        initial_cells=args.initial_cells,
        max_depth=args.max_depth,
        output=args.output,
    )


if __name__ == "__main__":
    main()
