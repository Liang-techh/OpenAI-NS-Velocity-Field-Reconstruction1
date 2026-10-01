"""Directed enclosures for the Section 10.2 fixed five-bump integrals.

The production :mod:`lei_ren_part1_paper_five_bump_map` stores finite
Gauss--Legendre data.  This module evaluates the same fixed bump family with
outward MP interval arithmetic.  On each interior ``t`` cell it integrates a
finite Taylor polynomial and bounds the first omitted derivative on the whole
cell.  The two endpoint strips are handled by a positive analytic tail bound;
the singular-looking formula for the bump is therefore never evaluated at
``|t|=1``.

For ``r=1/40`` and ``gamma_s(x)=b((x-s)/r)/(r*N)``, where
``b(t)=exp(-1/(1-t**2))`` and ``N=int_{-1}^1 b``,

``weight(s,p,1) = N**-1 int (s+r*t)**p b(t) dt``

and

``weight(s,p,2) = (r*N**2)**-1 int (s+r*t)**p b(t)**2 dt``.

Only fixed support integrals are enclosed here.  No source, parameter, cone,
or moment-defect error is inferred.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
import operator
from pathlib import Path
from typing import Any, Iterable, Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_interval_taylor import IntervalTaylor, integrate_symmetric


SOURCE_SECTION = "Part I Section 10.2, equations (10.7)--(10.10)"
CENTERS = ("1.25", "1.5", "1.75")
RADIUS = Fraction(1, 40)
DEFAULT_TAIL_WIDTH = Fraction(1, 128)
DEFAULT_INITIAL_CELLS = 32
DEFAULT_MAX_DEPTH = 24


def _index(value: Any, name: str) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return int(result)


def _fraction(value: Any, name: str = "value") -> Fraction:
    """Parse the small rational inputs used by the bump map."""

    if isinstance(value, Fraction):
        return value
    if isinstance(value, bool):
        return Fraction(int(value), 1)
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return Fraction(int(value[0]), int(value[1]))
    if isinstance(value, mp.mpf):
        return Fraction(str(value))
    text = str(value).strip()
    if not text:
        raise ValueError(f"{name} must be a rational value")
    try:
        return Fraction(text)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a rational value: {value!r}") from exc


def _interval_endpoints(value: Any) -> tuple[mp.mpf, mp.mpf]:
    raw = getattr(value, "_mpi_", None)
    if raw is None:
        point = mp.mpf(value)
        return point, point
    return mp.make_mpf(raw[0]), mp.make_mpf(raw[1])


def _record_interval(value: Any) -> dict[str, Any]:
    lower, upper = _interval_endpoints(value)
    return {
        "lower": mp.nstr(lower, 60),
        "upper": mp.nstr(upper, 60),
        "width": mp.nstr(upper - lower, 30),
        "lower_exact_mpf_tuple": list(lower._mpf_),
        "upper_exact_mpf_tuple": list(upper._mpf_),
    }


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


class _RationalMonomial:
    """Sparse rational representation of ``sum c*t**i*(1-t**2)**(-j)``."""

    def __init__(self, terms: Mapping[tuple[int, int], Fraction] | None = None):
        self.terms = {key: value for key, value in (terms or {}).items() if value}

    def __add__(self, other: "_RationalMonomial") -> "_RationalMonomial":
        result = dict(self.terms)
        for key, value in other.terms.items():
            result[key] = result.get(key, Fraction(0)) + value
            if not result[key]:
                del result[key]
        return _RationalMonomial(result)

    def __mul__(self, other: "_RationalMonomial") -> "_RationalMonomial":
        result: dict[tuple[int, int], Fraction] = {}
        for (i, j), left in self.terms.items():
            for (k, ell), right in other.terms.items():
                key = (i + k, j + ell)
                result[key] = result.get(key, Fraction(0)) + left * right
        return _RationalMonomial(result)

    def derivative(self) -> "_RationalMonomial":
        result: dict[tuple[int, int], Fraction] = {}
        for (i, j), coefficient in self.terms.items():
            if i:
                key = (i - 1, j)
                result[key] = result.get(key, Fraction(0)) + coefficient * i
            if j:
                # d(1-t**2)^(-j) = 2*j*t*(1-t**2)^(-j-1).
                key = (i + 1, j + 1)
                result[key] = result.get(key, Fraction(0)) + coefficient * (2 * j)
        return _RationalMonomial(result)

    def evaluate(self, ctx: MPIntervalContext, t: Any) -> Any:
        q = 1 - t * t
        result = ctx.mpf(0)
        for (power_t, power_q), coefficient in self.terms.items():
            scalar = ctx.mpf(coefficient.numerator) / ctx.mpf(coefficient.denominator)
            term = scalar * t**power_t
            if power_q:
                term *= q**(-power_q)
            result += term
        return result


class PaperBumpIntegralEnclosures:
    """Directed integral backend for the fixed ``beta_{1/40}`` family.

    ``order`` is the Taylor remainder order on each adaptive interior cell.
    ``tolerance`` is a requested aggregate Taylor remainder budget.  The
    returned intervals also include directed arithmetic and endpoint-tail
    overestimation.  The defaults are deliberately modest so a check can run
    in a bounded time; callers may increase precision or tighten the target.
    """

    coefficient_order = ("c1", "c2", "xi1", "xi2", "xi3")
    centers_as_strings = CENTERS
    radius_as_fraction = RADIUS
    finite_quadrature = False
    quadrature_enclosed = True

    def __init__(
        self,
        *,
        precision: int = 80,
        order: int = 16,
        tolerance: Any = "1e-24",
        tail_width: Any = DEFAULT_TAIL_WIDTH,
        initial_cells: int = DEFAULT_INITIAL_CELLS,
        max_depth: int = DEFAULT_MAX_DEPTH,
    ) -> None:
        self.precision = _index(precision, "precision")
        if self.precision < 60:
            raise ValueError("precision must be at least 60")
        self.order = _index(order, "order")
        if self.order < 4:
            raise ValueError("order must be at least 4")
        self.initial_cells = _index(initial_cells, "initial_cells")
        if self.initial_cells < 1:
            raise ValueError("initial_cells must be positive")
        self.max_depth = _index(max_depth, "max_depth")
        if self.max_depth < 1:
            raise ValueError("max_depth must be positive")
        self.ctx = MPIntervalContext()
        self.ctx.dps = self.precision
        self.radius = self._point(RADIUS)
        self.tail_width_fraction = _fraction(tail_width, "tail_width")
        if not (Fraction(0) < self.tail_width_fraction < Fraction(1, 2)):
            raise ValueError("tail_width must lie in (0,1/2)")
        self.tail_width = self._point(self.tail_width_fraction)
        self.tail_cut_fraction = Fraction(1) - self.tail_width_fraction
        self.tail_cut = self._point(self.tail_cut_fraction)
        self.tolerance = self._point(_fraction(tolerance, "tolerance"))
        if _interval_endpoints(self.tolerance)[0] <= 0:
            raise ValueError("tolerance must be positive")
        self.centers = tuple(self._point(_fraction(value, "center")) for value in CENTERS)
        self._derivative_cache: dict[int, list[_RationalMonomial]] = {}
        self._bump_point_cache: dict[tuple[int, Fraction], tuple[Any, ...]] = {}
        self._bump_cell_cache: dict[tuple[int, Fraction, Fraction], tuple[Any, ...]] = {}
        self._integral_cache: dict[tuple[Fraction, Fraction, int], Any] = {}
        self._integral_records: dict[tuple[Fraction, Fraction, int], dict[str, Any]] = {}
        # The normalizer is independent of the physical center when p=0;
        # use a positive representative center so the common support checks
        # remain valid.
        self._normalization, normalization_record = self._raw_integral(
            Fraction(5, 4), Fraction(0), 1
        )
        self._normalization_record = normalization_record
        self.metadata = {
            "source_section": SOURCE_SECTION,
            "centers": list(CENTERS),
            "radius": _fraction_text(RADIUS),
            "precision": self.precision,
            "Taylor_order": self.order,
            "tolerance_requested": mp.nstr(_interval_endpoints(self.tolerance)[1], 30),
            "tail_width": _fraction_text(self.tail_width_fraction),
            "tail_cut": _fraction_text(self.tail_cut_fraction),
            "initial_cells": self.initial_cells,
            "max_depth": self.max_depth,
            "directed_interval_arithmetic": True,
            "Taylor_remainder_enclosed": True,
            "endpoint_tail_enclosed": True,
            "normalization_interval": _record_interval(self._normalization),
            "finite_quadrature": False,
            "quadrature_enclosed": True,
            "source_error_enclosed": False,
            "parameter_error_enclosed": False,
            "cone_certified": False,
        }

    # ------------------------------------------------------------------
    # Interval and rational helpers.
    # ------------------------------------------------------------------
    def _point(self, value: Any) -> Any:
        rational = _fraction(value) if not isinstance(value, Fraction) else value
        return self.ctx.mpf(rational.numerator) / self.ctx.mpf(rational.denominator)

    def _interval(self, left: Fraction, right: Fraction) -> Any:
        return self.ctx.mpf([self._point(left), self._point(right)])

    def _fraction_scalar(self, value: Fraction) -> Any:
        return self.ctx.mpf(value.numerator) / self.ctx.mpf(value.denominator)

    def _sup_abs(self, value: Any) -> Any:
        lower, upper = _interval_endpoints(value)
        return max(abs(lower), abs(upper))

    @staticmethod
    def _falling(value: Fraction, order: int) -> Fraction:
        result = Fraction(1)
        for index in range(order):
            result *= value - index
        return result

    def _power(self, base: Any, exponent: Fraction) -> Any:
        if _interval_endpoints(base)[0] <= 0:
            raise ArithmeticError("fractional bump weight encountered nonpositive base")
        return base ** self._fraction_scalar(exponent)

    # ------------------------------------------------------------------
    # Raw bump derivatives and integrand derivatives.
    # ------------------------------------------------------------------
    def _bump_derivative_polynomials(self, multiplicity: int) -> list[_RationalMonomial]:
        cached = self._derivative_cache.get(multiplicity)
        if cached is not None:
            return cached
        if multiplicity not in (1, 2):
            raise ValueError("multiplicity must be 1 or 2")
        # b_k' = a_k b_k, a_k=-2*k*t/(1-t^2)^2.
        multiplier = _RationalMonomial({(1, 2): Fraction(-2 * multiplicity)})
        values = [_RationalMonomial({(0, 0): Fraction(1)})]
        for _ in range(self.order + 1):
            previous = values[-1]
            values.append(previous.derivative() + previous * multiplier)
        self._derivative_cache[multiplicity] = values
        return values

    def _bump_derivative(self, multiplicity: int, order: int, t: Any) -> Any:
        q = 1 - t * t
        if _interval_endpoints(q)[0] <= 0:
            raise ArithmeticError("bump derivative cell reaches an endpoint")
        raw = self.ctx.exp(-self.ctx.mpf(multiplicity) / q)
        rational_part = self._bump_derivative_polynomials(multiplicity)[order].evaluate(
            self.ctx, t
        )
        return rational_part * raw

    def _bump_derivatives_at_point(
        self, multiplicity: int, point: Fraction
    ) -> tuple[Any, ...]:
        key = (multiplicity, point)
        cached = self._bump_point_cache.get(key)
        if cached is not None:
            return cached
        value = tuple(
            self._bump_derivative(multiplicity, derivative_order, self._point(point))
            for derivative_order in range(self.order + 1)
        )
        self._bump_point_cache[key] = value
        return value

    def _bump_derivatives_on_cell(
        self, multiplicity: int, left: Fraction, right: Fraction
    ) -> tuple[Any, ...]:
        key = (multiplicity, left, right)
        cached = self._bump_cell_cache.get(key)
        if cached is not None:
            return cached
        cell = self._interval(left, right)
        value = tuple(
            self._bump_derivative(multiplicity, derivative_order, cell)
            for derivative_order in range(self.order + 1)
        )
        self._bump_cell_cache[key] = value
        return value

    def _weight_derivative(
        self, center: Fraction, power: Fraction, order: int, t: Any
    ) -> Any:
        x = self._point(center) + self.radius * t
        if _interval_endpoints(x)[0] <= 0:
            raise ArithmeticError("bump support has nonpositive physical radius")
        coefficient = self._falling(power, order) * RADIUS**order
        return self._fraction_scalar(coefficient) * self._power(x, power - order)

    def _integrand_derivative(
        self,
        center: Fraction,
        power: Fraction,
        multiplicity: int,
        order: int,
        t: Any,
        bump_values: tuple[Any, ...] | None = None,
    ) -> Any:
        if power == 0:
            # Avoid needless fractional-power interval work for the
            # normalizer and the unweighted quadratic rows.
            if bump_values is not None:
                return bump_values[order]
            return self._bump_derivative(multiplicity, order, t)
        result = self.ctx.mpf(0)
        for weight_order in range(order + 1):
            bump = (
                bump_values[order - weight_order]
                if bump_values is not None
                else self._bump_derivative(
                    multiplicity, order - weight_order, t
                )
            )
            result += (
                self.ctx.mpf(math.comb(order, weight_order))
                * self._weight_derivative(center, power, weight_order, t)
                * bump
            )
        return result

    # ------------------------------------------------------------------
    # Cell and endpoint enclosure.
    # ------------------------------------------------------------------
    def _tail_bound(
        self, center: Fraction, power: Fraction, multiplicity: int
    ) -> Any:
        """Bound both endpoint strips by monotonicity and positivity."""

        if multiplicity not in (1, 2):
            raise ValueError("multiplicity must be 1 or 2")
        total = self.ctx.mpf(0)
        for left, right in (
            (Fraction(-1), -self.tail_cut_fraction),
            (self.tail_cut_fraction, Fraction(1)),
        ):
            x_left = self._point(center) + self.radius * self._point(left)
            x_right = self._point(center) + self.radius * self._point(right)
            x_min = min(_interval_endpoints(x_left)[0], _interval_endpoints(x_right)[0])
            x_max = max(_interval_endpoints(x_left)[1], _interval_endpoints(x_right)[1])
            x_bound = x_max if power >= 0 else x_min
            weight_bound = self._power(self.ctx.mpf([x_min, x_max]), power)
            # 1-t^2 <= 2*(1-|t|) <= 2*tail_width on either endpoint strip.
            bump_bound = self.ctx.exp(
                -self.ctx.mpf(multiplicity) / (2 * self.tail_width)
            )
            width = self._point(right - left)
            # x_bound is retained above for readability and endpoint auditing;
            # weight_bound is the directed interval used in the product.
            del x_bound
            total += width * weight_bound * bump_bound
        return total

    def _cell_integral(
        self,
        center: Fraction,
        power: Fraction,
        multiplicity: int,
        left: Fraction,
        right: Fraction,
    ) -> tuple[Any, Any]:
        midpoint = (left + right) / 2
        half_width = (right - left) / 2
        midpoint_value = self._point(midpoint)
        half_value = self._point(half_width)
        midpoint_bumps = self._bump_derivatives_at_point(multiplicity, midpoint)
        cell_bumps = self._bump_derivatives_on_cell(multiplicity, left, right)
        coefficients = [
            self._integrand_derivative(
                center,
                power,
                multiplicity,
                derivative_order,
                midpoint_value,
                midpoint_bumps,
            )
            / self.ctx.mpf(math.factorial(derivative_order))
            for derivative_order in range(self.order)
        ]
        # The center jet has one extra dummy coefficient so
        # integrate_symmetric uses the scalar order-th derivative bound as its
        # Taylor remainder coefficient.
        center_jet = IntervalTaylor(
            self.ctx, coefficients + [self.ctx.mpf(0)]
        )
        cell = self._interval(left, right)
        derivative_bound = (
            abs(
                self._integrand_derivative(
                    center,
                    power,
                    multiplicity,
                    self.order,
                    cell,
                    cell_bumps,
                )
            )
            / self.ctx.mpf(math.factorial(self.order))
        )
        enclosed = integrate_symmetric(center_jet, derivative_bound, half_value)
        # ``derivative_bound`` is already the coefficient bound
        # |f^(N)|/N! expected by integrate_symmetric.  Do not divide by N!
        # a second time when auditing the same Taylor remainder.
        remainder_enclosure = derivative_bound * (
            2 * half_value ** (self.order + 1) / (self.order + 1)
        )
        return enclosed, remainder_enclosure

    def _adaptive_cell(
        self,
        center: Fraction,
        power: Fraction,
        multiplicity: int,
        left: Fraction,
        right: Fraction,
        budget: Any,
        depth: int,
        stats: dict[str, Any],
    ) -> Any:
        enclosed, remainder = self._cell_integral(
            center, power, multiplicity, left, right
        )
        # Compare against the lower endpoint of the allocated budget.  This
        # keeps the requested directed error budget conservative even when
        # the budget interval itself has rounding width.
        budget_upper = _interval_endpoints(budget)[0]
        remainder_upper = _interval_endpoints(remainder)[1]
        if remainder_upper <= budget_upper:
            stats["accepted_cells"] += 1
            stats["max_depth"] = max(stats["max_depth"], depth)
            stats["remainder_interval"] += remainder
            return enclosed
        if depth >= self.max_depth:
            raise RuntimeError(
                "bump Taylor cell did not meet tolerance before max_depth; "
                f"remainder={remainder_upper} budget={budget_upper}"
            )
        midpoint = (left + right) / 2
        stats["splits"] += 1
        half_budget = budget / 2
        return self._adaptive_cell(
            center,
            power,
            multiplicity,
            left,
            midpoint,
            half_budget,
            depth + 1,
            stats,
        ) + self._adaptive_cell(
            center,
            power,
            multiplicity,
            midpoint,
            right,
            half_budget,
            depth + 1,
            stats,
        )

    def _raw_integral(
        self, center: Fraction, power: Fraction, multiplicity: int
    ) -> tuple[Any, dict[str, Any]]:
        key = (center, power, multiplicity)
        cached = self._integral_cache.get(key)
        if cached is not None:
            return cached, self._integral_records[key]
        if multiplicity not in (1, 2):
            raise ValueError("multiplicity must be 1 or 2")
        if self._point(center) - self.radius <= 0:
            raise ValueError("center-radius must be positive")
        tail = self._tail_bound(center, power, multiplicity)
        tail_upper = _interval_endpoints(tail)[1]
        target_upper = _interval_endpoints(self.tolerance)[1]
        # The endpoint contribution is part of the enclosure, rather than a
        # Taylor error.  Require it to fit a fixed half-budget so the receipt
        # reports a meaningful requested accuracy.
        if tail_upper > target_upper / 2:
            raise RuntimeError(
                "endpoint tail bound exceeds half the requested tolerance; "
                "reduce tail_width or loosen tolerance"
            )
        cut = self.tail_cut_fraction
        interior_length = 2 * cut
        target_interval = self.tolerance / 2
        total = self.ctx.mpf(0)
        stats = {
            "accepted_cells": 0,
            "splits": 0,
            "max_depth": 0,
            "remainder_interval": self.ctx.mpf(0),
        }
        for index in range(self.initial_cells):
            left = -cut + interior_length * index / self.initial_cells
            right = -cut + interior_length * (index + 1) / self.initial_cells
            width_ratio = self._point(right - left) / self._point(interior_length)
            budget = target_interval * width_ratio
            total += self._adaptive_cell(
                center,
                power,
                multiplicity,
                left,
                right,
                budget,
                0,
                stats,
            )
        # Add the interval itself rather than reconstructing it from a rounded
        # endpoint; this preserves the outward endpoint produced by MPIV.
        total += tail
        record = {
            "center": _fraction_text(center),
            "power": _fraction_text(power),
            "multiplicity": multiplicity,
            "interval": _record_interval(total),
            "endpoint_tail_upper": mp.nstr(tail_upper, 30),
            "accepted_cells": stats["accepted_cells"],
            "splits": stats["splits"],
            "max_depth": stats["max_depth"],
            "Taylor_remainder_budget": mp.nstr(target_upper / 2, 30),
            "Taylor_remainder_interval": _record_interval(stats["remainder_interval"]),
            "Taylor_remainder_upper": mp.nstr(
                _interval_endpoints(stats["remainder_interval"])[1], 30
            ),
            "endpoint_tail_bound_formula": "2*tail_width*sup_x^p*exp(-multiplicity/(2*tail_width))",
        }
        self._integral_cache[key] = total
        self._integral_records[key] = record
        return total, record

    # ------------------------------------------------------------------
    # Public normalized weights and the Section 10.2 coefficient rows.
    # ------------------------------------------------------------------
    @property
    def normalization(self) -> Any:
        return self._normalization

    @property
    def beta_normalization(self) -> Any:
        return self._normalization

    def weight(self, center: Any, power: Any, multiplicity: int = 1) -> Any:
        """Return ``int x**power * gamma_center(x)**multiplicity dx``."""

        center_fraction = _fraction(center, "center")
        power_fraction = _fraction(power, "power")
        multiplicity = _index(multiplicity, "multiplicity")
        if multiplicity not in (1, 2):
            raise ValueError("multiplicity must be 1 or 2")
        raw, _ = self._raw_integral(center_fraction, power_fraction, multiplicity)
        if multiplicity == 1:
            return raw / self._normalization
        return raw / (self.radius * self._normalization * self._normalization)

    def integral_record(self, center: Any, power: Any, multiplicity: int = 1) -> dict[str, Any]:
        key = (_fraction(center, "center"), _fraction(power, "power"), _index(multiplicity, "multiplicity"))
        self.weight(center, power, multiplicity)
        return dict(self._integral_records[key])

    @staticmethod
    def _zero(ctx: MPIntervalContext) -> Any:
        return ctx.mpf(0)

    def map_enclosure(self, amplitude: Any = 1) -> dict[str, Any]:
        """Assemble directed linear/quadratic rows of equation (10.8)."""

        c = tuple(_fraction(value, "center") for value in CENTERS)
        zero = self._zero(self.ctx)
        matrix = [[zero for _ in range(5)] for _ in range(5)]
        matrix[0][0] = self.weight(c[0], Fraction(0), 1)
        matrix[0][1] = self.weight(c[2], Fraction(0), 1)
        matrix[1][0] = self.weight(c[0], Fraction(3, 5), 1)
        matrix[1][1] = self.weight(c[2], Fraction(3, 5), 1)
        for index, center in enumerate(c):
            column = 2 + index
            matrix[2][column] = self.weight(center, Fraction(1, 2), 1)
            matrix[3][column] = -self.weight(center, Fraction(1, 10), 1)
            matrix[4][column] = self.weight(center, Fraction(-9, 10), 1)

        def same_or_zero(left: Fraction, right: Fraction, power: Fraction) -> Any:
            if left != right:
                return zero
            return self.weight(left, power, 2)

        fg = [
            [same_or_zero(c[axial], c[angular], Fraction(1, 2)) for angular in range(3)]
            for axial in (0, 2)
        ]
        gg = [
            [same_or_zero(c[left], c[right], Fraction(0)) for right in (0, 2)]
            for left in (0, 2)
        ]
        ff = [
            [same_or_zero(c[left], c[right], Fraction(0)) for right in range(3)]
            for left in range(3)
        ]
        ff_over_x = [
            [same_or_zero(c[left], c[right], Fraction(-1)) for right in range(3)]
            for left in range(3)
        ]
        return {
            "linear_matrix": matrix,
            "quadratic_weights": {
                "fg_sqrt": fg,
                "gg": gg,
                "ff": ff,
                "ff_over_x": ff_over_x,
            },
            "normalization": self._normalization,
            "amplitude": amplitude,
            "source_section": SOURCE_SECTION,
            "control_order": self.coefficient_order,
            "supports_disjoint_cross_products_zero": True,
            "directed_interval_arithmetic": True,
            "quadrature_enclosed": True,
            "source_error_enclosed": False,
            "cone_certified": False,
        }


def _module_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


__all__ = [
    "PaperBumpIntegralEnclosures",
    "CENTERS",
    "RADIUS",
    "SOURCE_SECTION",
]
