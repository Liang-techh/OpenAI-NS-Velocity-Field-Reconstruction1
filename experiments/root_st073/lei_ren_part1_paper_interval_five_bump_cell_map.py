"""Directed cumulative five-bump enclosures on a radial cell.

The frozen partial map evaluates cumulative integrals at one exact rational
endpoint.  This companion encloses the same cumulative map for every
``x`` in an exact rational cell ``[left, right]``.  A cumulative primitive is
monotone because every supported integrand is nonnegative.  Consequently a
cell that reaches a bump support is safely enclosed by ``[0, upper(full)]``;
after the support it is enclosed by the saved full-support interval.  The
map keeps all five coefficient rows and the quadratic terms from the frozen
partial map.  It does not solve controls or certify a field or cone.
"""

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
    _fraction,
)
from lei_ren_part1_paper_interval_partial_five_bump_map import (
    COEFFICIENT_ORDER,
    CENTER_FRACTIONS,
    RADIUS_FRACTION,
    FULL_MODULE,
    FULL_RECEIPT,
    IntervalPartialFiveBumpMap,
    _as_vector,
    _coerce_jet,
    _restore_in_context,
)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


HERE = Path(__file__).resolve().parent
PARTIAL_MODULE = "lei_ren_part1_paper_interval_partial_five_bump_map.py"
PARTIAL_RECEIPT = "lei_ren_part1_paper_interval_partial_five_bump_map.json"


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _pack(value: Any) -> Any:
    if isinstance(value, IntervalTaylor):
        return [_pack(item) for item in value.coefficients]
    if isinstance(value, dict):
        return {key: _pack(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_pack(item) for item in value]
    return value


def _as_fraction(value: Any, name: str) -> Fraction:
    return _fraction(value, name)


def _interval_upper(value: Any) -> mp.mpf:
    return mp.make_mpf(tuple(value._mpi_[1]))


def _interval_lower(value: Any) -> mp.mpf:
    return mp.make_mpf(tuple(value._mpi_[0]))


def _overlap(left: Any, right: Any) -> bool:
    left_lo, left_hi = endpoints(left)
    right_lo, right_hi = endpoints(right)
    return left_lo <= right_hi and right_lo <= left_hi


def _point_inside(value: Any, interval: Any) -> bool:
    lo, hi = endpoints(interval)
    return lo <= value <= hi


class IntervalFiveBumpCellMap:
    """Cumulative directed map for every point in one radial cell."""

    coefficient_order = COEFFICIENT_ORDER
    centers = CENTER_FRACTIONS
    radius = RADIUS_FRACTION

    def __init__(self, parent_map: IntervalPartialFiveBumpMap):
        if not isinstance(parent_map, IntervalPartialFiveBumpMap):
            raise TypeError("parent_map must be IntervalPartialFiveBumpMap")
        self.parent = parent_map
        self.ctx = parent_map.ctx
        self.precision = parent_map.precision
        self.order = parent_map.order
        self.tolerance = parent_map.tolerance
        self.normalization = _restore_in_context(self.ctx, parent_map.normalization)
        self.full_module_sha256 = parent_map.full_module_sha256
        self.full_receipt_sha256 = parent_map.full_receipt_sha256
        self.partial_module_sha256 = _hash(HERE / PARTIAL_MODULE)
        self.partial_receipt_sha256 = (
            _hash(HERE / PARTIAL_RECEIPT) if (HERE / PARTIAL_RECEIPT).exists() else None
        )
        self._full = self.parent.weights_at(Fraction(2))
        nlo = _interval_lower(self.normalization)
        if nlo <= 0:
            raise ValueError("normalization lower endpoint must be positive")
        radius = self.ctx.mpf(RADIUS_FRACTION.numerator) / self.ctx.mpf(
            RADIUS_FRACTION.denominator
        )
        # The raw bump satisfies exp(-1/q) <= exp(-1).  For its derivative,
        # 2|t|q^-2 exp(-1/q) <= 8 exp(-2), q=1-t^2.  We use the lower
        # endpoint of the positive normalizer to retain a directed upper.
        beta_upper = self.ctx.mpf(
            [
                mp.mpf("0"),
                _interval_upper(
                    self.ctx.exp(-1) / (radius * self.ctx.mpf(nlo))
                ),
            ]
        )
        derivative_upper = _interval_upper(
            self.ctx.mpf(8)
            * self.ctx.exp(-2)
            / (radius * radius * self.ctx.mpf(nlo))
        )
        self.beta_upper = beta_upper
        derivative_positive = self.ctx.mpf([mp.mpf(0), derivative_upper])
        derivative_negative = -derivative_positive
        self.derivative_upper = self.ctx.mpf(
            [endpoints(derivative_negative)[0], endpoints(derivative_positive)[1]]
        )
        self.metadata = {
            "cumulative_cell_map": True,
            "coefficient_order": list(COEFFICIENT_ORDER),
            "exact_fraction_cells": True,
            "positive_monotone_primitive_envelope": True,
            "full_support_saved_intervals": True,
            "bump_value_upper": beta_upper,
            "bump_derivative_symmetric_upper": derivative_upper,
            "support_endpoint_values_exact_zero": True,
            "source_error_enclosed": False,
            "parameter_error_enclosed": False,
            "field_installed": False,
            "controls_solved": False,
            "cone_certified": False,
        }

    def _cell(self, left: Any, right: Any) -> tuple[Fraction, Fraction]:
        left_fraction = _as_fraction(left, "left")
        right_fraction = _as_fraction(right, "right")
        if left_fraction < 1 or right_fraction > 2:
            raise ValueError("cell endpoints must lie in [1, 2]")
        if left_fraction > right_fraction:
            raise ValueError("cell endpoints must be ordered")
        return left_fraction, right_fraction

    def _primitive(
        self,
        center: Fraction,
        power: Any,
        multiplicity: int,
        left: Fraction,
        right: Fraction,
    ) -> Any:
        """Enclose cumulative ``int_1^x y**power gamma**k dy`` on a cell."""

        support_left = center - RADIUS_FRACTION
        support_right = center + RADIUS_FRACTION
        if right <= support_left:
            return self.ctx.mpf(0)
        if left >= support_right:
            # The frozen full interval is exactly the required cumulative
            # value after support completion.
            return self.parent._weight(center, power, Fraction(2), multiplicity)
        full = self.parent._weight(center, power, Fraction(2), multiplicity)
        upper = _interval_upper(full)
        return self.ctx.mpf([mp.mpf(0), upper])

    def _point_weights(self, x: Fraction) -> dict[str, Any]:
        return self.parent.weights_at(x)

    def weights_on(self, left: Any, right: Any) -> dict[str, Any]:
        """Enclose cumulative map weights at every ``x`` in ``[left,right]``."""

        left_fraction, right_fraction = self._cell(left, right)
        if left_fraction == right_fraction:
            return self._point_weights(left_fraction)
        c = self.ctx
        zero = c.mpf(0)

        def w(index: int, power: Any, multiplicity: int = 1) -> Any:
            return self._primitive(
                CENTER_FRACTIONS[index], power, multiplicity, left_fraction, right_fraction
            )

        L = [[zero for _ in range(5)] for _ in range(5)]
        L[0][0] = w(0, 0)
        L[0][1] = w(2, 0)
        L[1][0] = w(0, "3/5")
        L[1][1] = w(2, "3/5")
        for index in range(3):
            column = index + 2
            L[2][column] = w(index, "1/2")
            # Negating after the positive envelope preserves the correct
            # interval orientation for the row -int(x**.1 gamma).
            L[3][column] = -w(index, "1/10")
            L[4][column] = w(index, "-9/10")
        return {
            "L": L,
            "fg": [w(0, "1/2", 2), w(2, "1/2", 2)],
            "gg": [w(0, 0, 2), w(2, 0, 2)],
            "ff": [w(index, 0, 2) for index in range(3)],
            "ff_over_x": [w(index, -1, 2) for index in range(3)],
        }

    def bump_values_on(self, left: Any, right: Any) -> dict[str, Any]:
        """Enclose beta and d beta/dx over every point in a radial cell."""

        left_fraction, right_fraction = self._cell(left, right)
        if left_fraction == right_fraction:
            return self.parent.bump_values(left_fraction)
        values = []
        derivatives = []
        statuses = []
        for center in CENTER_FRACTIONS:
            support_left = center - RADIUS_FRACTION
            support_right = center + RADIUS_FRACTION
            if right_fraction <= support_left or left_fraction >= support_right:
                values.append(self.ctx.mpf(0))
                derivatives.append(self.ctx.mpf(0))
                statuses.append("outside_or_endpoint_exact_zero")
                continue
            values.append(self.beta_upper)
            if right_fraction <= center:
                derivatives.append(
                    self.ctx.mpf([mp.mpf(0), _interval_upper(self.derivative_upper)])
                )
                status = "support_left_half_nonnegative_derivative"
            elif left_fraction >= center:
                positive = self.ctx.mpf([mp.mpf(0), _interval_upper(self.derivative_upper)])
                negative = -positive
                derivatives.append(
                    self.ctx.mpf([endpoints(negative)[0], mp.mpf(0)])
                )
                status = "support_right_half_nonpositive_derivative"
            else:
                derivatives.append(self.derivative_upper)
                status = "support_crossing_symmetric_derivative"
            statuses.append(status)
        return {
            "left": str(left_fraction),
            "right": str(right_fraction),
            "values": values,
            "derivatives": derivatives,
            "centers": [str(center) for center in CENTER_FRACTIONS],
            "support_status": statuses,
        }

    def apply(self, h: Any, invAm2: Any, left: Any, right: Any) -> list[Any]:
        """Return the five C1 rows ``L(x)h+Q(x,h)`` over a radial cell."""

        left_fraction, right_fraction = self._cell(left, right)
        if left_fraction == right_fraction:
            return self.parent.apply(h, invAm2, left_fraction)
        weights = self.weights_on(left_fraction, right_fraction)
        coefficients = tuple(
            _coerce_jet(self.ctx, item, f"h[{index}]")
            for index, item in enumerate(_as_vector(h, "h"))
        )
        inv = _coerce_jet(self.ctx, invAm2, "invAm2")
        c1, c2, x1, x2, x3 = coefficients
        zero = coefficients[0] * 0
        linear = [
            sum(
                (coefficients[column] * weights["L"][row][column] for column in range(5)),
                zero,
            )
            for row in range(5)
        ]
        quadratic = [
            zero,
            (c1 * x1) * weights["fg"][0] + (c2 * x3) * weights["fg"][1],
            zero,
            inv
            * sum(
                ((value * value) * weights["gg"][index] for index, value in enumerate((c1, c2))),
                zero,
            )
            - sum(
                ((value * value) * weights["ff"][index] for index, value in enumerate((x1, x2, x3))),
                zero,
            )
            / 2,
            sum(
                ((value * value) * weights["ff_over_x"][index] for index, value in enumerate((x1, x2, x3))),
                zero,
            )
            / 2,
        ]
        return [linear[row] + quadratic[row] for row in range(5)]


def _scalar_bump(t: mp.mpf) -> mp.mpf:
    if abs(t) >= 1:
        return mp.mpf(0)
    return mp.exp(-1 / (1 - t * t))


def _independent_cell_integral(center: Fraction, power: Fraction, multiplicity: int, left: Fraction, right: Fraction) -> mp.mpf:
    with mp.workdps(100):
        radius = mp.mpf(1) / 40
        center_mp = mp.mpf(center.numerator) / center.denominator
        left_mp = mp.mpf(left.numerator) / left.denominator
        right_mp = mp.mpf(right.numerator) / right.denominator
        norm = mp.quad(_scalar_bump, [-1, 0, 1])
        support_left = max(left_mp, center_mp - radius)
        support_right = min(right_mp, center_mp + radius)
        if support_right <= support_left:
            return mp.mpf(0)
        return mp.quad(
            lambda y: y ** (mp.mpf(power.numerator) / power.denominator)
            * (_scalar_bump((y - center_mp) / radius) / (radius * norm)) ** multiplicity,
            [support_left, (support_left + support_right) / 2, support_right],
        )


def _run_checks() -> dict[str, Any]:
    parent = IntervalPartialFiveBumpMap()
    cell = IntervalFiveBumpCellMap(parent)
    c = cell.ctx
    point_x = Fraction(149, 100)
    point_weights = cell.weights_on(point_x, point_x)
    parent_weights = parent.weights_at(point_x)
    point_weight_equalities = 0
    for name in ("L", "fg", "gg", "ff", "ff_over_x"):
        left = point_weights[name]
        right = parent_weights[name]
        if isinstance(left, list):
            lv = [v for row in left for v in (row if isinstance(row, list) else [row])]
            rv = [v for row in right for v in (row if isinstance(row, list) else [row])]
        else:
            lv, rv = left, right
        for a, b in zip(lv, rv):
            if endpoints(a) != endpoints(b):
                raise AssertionError("degenerate cell did not delegate point weights")
            point_weight_equalities += 1
    point_bumps = cell.bump_values_on(point_x, point_x)
    parent_bumps = parent.bump_values(point_x)
    point_bump_equalities = 0
    for a, b in zip(point_bumps["values"] + point_bumps["derivatives"], parent_bumps["values"] + parent_bumps["derivatives"]):
        if endpoints(a) != endpoints(b):
            raise AssertionError("degenerate cell did not delegate point bumps")
        point_bump_equalities += 1
    if not (_interval_lower(point_bumps["values"][1]) > 0 or _interval_upper(point_bumps["values"][1]) > 0):
        raise AssertionError("interior point bump value is zero")
    if endpoints(point_bumps["derivatives"][1]) == (mp.mpf(0), mp.mpf(0)):
        raise AssertionError("interior point bump derivative unexpectedly zero")

    left, right = Fraction(6, 5), Fraction(13, 10)
    cell_weights = cell.weights_on(left, right)
    primitive_specs = [
        (CENTER_FRACTIONS[0], Fraction(0), 1),
        (CENTER_FRACTIONS[0], Fraction(3, 5), 1),
        (CENTER_FRACTIONS[0], Fraction(1, 2), 2),
        (CENTER_FRACTIONS[0], Fraction(0), 2),
        (CENTER_FRACTIONS[0], Fraction(-1), 2),
    ]
    primitive_keys = [
        ("L", 0, 0),
        ("L", 1, 0),
        ("fg", 0, 0),
        ("gg", 0, 0),
        ("ff_over_x", 0, 0),
    ]
    independent_weight_checks = 0
    for (center, power, multiplicity), (name, row, column) in zip(primitive_specs, primitive_keys):
        target = _independent_cell_integral(center, power, multiplicity, left, right)
        if name == "L":
            interval = cell_weights[name][row][column]
        else:
            interval = cell_weights[name][row]
        if not _point_inside(target, interval):
            raise AssertionError(("support-crossing primitive escaped", name, row, column))
        independent_weight_checks += 1

    bump_packet = cell.bump_values_on(left, right)
    dense_bump_checks = 0
    with mp.workdps(100):
        radius = mp.mpf(1) / 40
        norm = mp.quad(_scalar_bump, [-1, 0, 1])
        for index, center in enumerate(CENTER_FRACTIONS):
            center_mp = mp.mpf(center.numerator) / center.denominator
            for step in range(17):
                x = mp.mpf(left.numerator) / left.denominator + (
                    mp.mpf(step) / 16
                ) * (
                    mp.mpf(right.numerator) / right.denominator
                    - mp.mpf(left.numerator) / left.denominator
                )
                t = (x - center_mp) / radius
                beta = _scalar_bump(t) / (radius * norm)
                derivative = (
                    beta * (-2 * t / (radius * (1 - t * t) ** 2))
                    if abs(t) < 1
                    else mp.mpf(0)
                )
                if not _point_inside(beta, bump_packet["values"][index]):
                    raise AssertionError(("cell beta escaped", index, step))
                if not _point_inside(derivative, bump_packet["derivatives"][index]):
                    raise AssertionError(("cell beta derivative escaped", index, step))
                dense_bump_checks += 2

    report = {
        "precision": cell.precision,
        "coefficient_order": list(COEFFICIENT_ORDER),
        "cell": {"left": str(left), "right": str(right)},
        "point_interior": str(point_x),
        "point_weight_equalities": point_weight_equalities,
        "point_bump_equalities": point_bump_equalities,
        "independent_support_crossing_weight_checks": independent_weight_checks,
        "independent_dense_bump_value_derivative_checks": dense_bump_checks,
        "positive_primitive_envelope": "[0, upper(full)] when support is reached",
        "beta_upper": cell.beta_upper,
        "derivative_upper": cell.derivative_upper,
        "full_module_sha256": cell.full_module_sha256,
        "full_receipt_sha256": cell.full_receipt_sha256,
        "partial_module_sha256": cell.partial_module_sha256,
        "partial_receipt_sha256": cell.partial_receipt_sha256,
        "directed_interval_arithmetic": True,
        "source_error_enclosed": False,
        "parameter_error_enclosed": False,
        "field_installed": False,
        "controls_solved": False,
        "cone_certified": False,
    }
    report["input_hashes"] = {
        name: _hash(HERE / name)
        for name in (
            Path(__file__).name,
            PARTIAL_MODULE,
            PARTIAL_RECEIPT,
            FULL_MODULE,
            FULL_RECEIPT,
        )
    }
    return report


def run() -> dict[str, Any]:
    report = _run_checks()
    output = HERE / "lei_ren_part1_paper_interval_five_bump_cell_map.json"
    output.write_text(json.dumps(encode(_pack(report)), indent=2) + "\n", encoding="utf-8")
    print(
        "Cell map checks:",
        report["point_weight_equalities"],
        "point weights;",
        report["independent_support_crossing_weight_checks"],
        "support-crossing primitives;",
        report["independent_dense_bump_value_derivative_checks"],
        "dense bump checks",
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
