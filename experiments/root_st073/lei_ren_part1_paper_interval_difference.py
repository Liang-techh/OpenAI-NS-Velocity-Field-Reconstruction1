"""Paired directed intervals for a nominal value and its perturbation.

An :class:`IntervalDifference` stores ``nominal`` and ``difference`` as
separate intervals.  It represents the pair ``(x0, dx)`` and propagates the
two components without subtracting a nearly equal nominal value.  The class
is intentionally tied to a caller-supplied ``MPIntervalContext``; it never
changes mpmath's global precision context.
"""

from __future__ import annotations

import operator

import mpmath as mp


def _endpoints(value):
    return tuple(mp.make_mpf(item) for item in value._mpi_)


def _interval(ctx, value):
    return value if hasattr(value, "_mpi_") else ctx.mpf(value)


def _contains_zero(value):
    lower, upper = _endpoints(value)
    return lower <= 0 <= upper


class IntervalDifference:
    """Directed nominal/perturbation pair in one interval context.

    ``nominal`` and ``difference`` are both interval enclosures.  A scalar
    operand is promoted to a nominal interval with an exact zero difference.
    """

    def __init__(self, ctx, nominal, difference=0):
        self.ctx = ctx
        self.nominal = _interval(ctx, nominal)
        self.difference = _interval(ctx, difference)

    @property
    def x0(self):
        return self.nominal

    @property
    def dx(self):
        return self.difference

    @property
    def value(self):
        """Return the directed interval for nominal plus perturbation."""

        return self.nominal + self.difference

    @classmethod
    def constant(cls, ctx, value):
        return cls(ctx, value, ctx.mpf(0))

    @classmethod
    def converter(cls, ctx):
        """Return a one-argument scalar-converter closure for recurrences."""

        def convert(value):
            if isinstance(value, cls):
                if value.ctx is not ctx:
                    raise ValueError("IntervalDifference contexts must match")
                return value
            return cls(ctx, value)

        return convert

    factory = converter
    scalar_converter = converter

    def _coerce(self, other):
        if isinstance(other, IntervalDifference):
            if other.ctx is not self.ctx:
                raise ValueError("IntervalDifference contexts must match")
            return other
        return IntervalDifference.constant(self.ctx, other)

    def __add__(self, other):
        other = self._coerce(other)
        return IntervalDifference(
            self.ctx,
            self.nominal + other.nominal,
            self.difference + other.difference,
        )

    __radd__ = __add__
    add = __add__

    def __neg__(self):
        return IntervalDifference(self.ctx, -self.nominal, -self.difference)

    def __sub__(self, other):
        other = self._coerce(other)
        return IntervalDifference(
            self.ctx,
            self.nominal - other.nominal,
            self.difference - other.difference,
        )

    sub = __sub__

    def __rsub__(self, other):
        return self._coerce(other) - self

    def __mul__(self, other):
        other = self._coerce(other)
        ctx = self.ctx
        nominal = self.nominal * other.nominal
        difference = (
            self.nominal * other.difference
            + other.nominal * self.difference
            + self.difference * other.difference
        )
        return IntervalDifference(ctx, nominal, difference)

    __rmul__ = __mul__
    mul = __mul__

    def reciprocal(self):
        """Return ``(x0+dx)^-1`` with a nominal and difference component."""

        if _contains_zero(self.nominal) or _contains_zero(self.nominal + self.difference):
            raise ZeroDivisionError("IntervalDifference reciprocal denominator crosses zero")
        nominal = self.ctx.mpf(1) / self.nominal
        difference = (self.ctx.mpf(0) - nominal * self.difference) / (
            self.nominal + self.difference
        )
        return IntervalDifference(self.ctx, nominal, difference)

    inverse = reciprocal
    inv = reciprocal

    def __truediv__(self, other):
        other = self._coerce(other)
        if _contains_zero(other.nominal) or _contains_zero(other.nominal + other.difference):
            raise ZeroDivisionError("IntervalDifference division denominator crosses zero")
        nominal = self.nominal / other.nominal
        difference = (
            self.difference - nominal * other.difference
        ) / (other.nominal + other.difference)
        return IntervalDifference(self.ctx, nominal, difference)

    div = __truediv__

    def __rtruediv__(self, other):
        return self._coerce(other) / self

    def __pow__(self, exponent):
        try:
            exponent = operator.index(exponent)
        except TypeError as exc:
            raise TypeError("IntervalDifference powers must be integers") from exc
        if exponent < 0:
            return self.reciprocal() ** (-exponent)
        result = IntervalDifference.constant(self.ctx, self.ctx.mpf(1))
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            exponent //= 2
            if exponent:
                base = base * base
        return result

    def exp(self):
        """Return the paired enclosure of ``exp(nominal + difference)``.

        The perturbation is evaluated with ``expm1`` so a zero difference
        remains exactly zero even when the nominal interval is wide.
        """

        nominal = self.ctx.exp(self.nominal)
        difference = nominal * self.ctx.expm1(self.difference)
        return IntervalDifference(self.ctx, nominal, difference)

    def log(self):
        """Return the paired enclosure of ``log(nominal + difference)``."""

        if not _strictly_positive(self.nominal) or not _strictly_positive(self.value):
            raise ValueError(
                "IntervalDifference log requires positive nominal and total"
            )
        nominal = self.ctx.log(self.nominal)
        difference = self.ctx.log1p(self.difference / self.nominal)
        return IntervalDifference(self.ctx, nominal, difference)

    def sqrt(self):
        """Return the paired enclosure of ``sqrt(nominal + difference)``."""

        if not _strictly_positive(self.nominal) or not _strictly_positive(self.value):
            raise ValueError(
                "IntervalDifference sqrt requires positive nominal and total"
            )
        nominal = self.ctx.sqrt(self.nominal)
        actual_root = self.ctx.sqrt(self.value)
        difference = self.difference / (actual_root + nominal)
        return IntervalDifference(self.ctx, nominal, difference)

    def __repr__(self):
        return (
            "IntervalDifference(" f"nominal={self.nominal!r}, "
            f"difference={self.difference!r})"
        )


def _strictly_positive(value):
    """Whether an interval has a strictly positive lower endpoint."""

    lower, _ = _endpoints(value)
    return lower > 0


def make_converter(ctx):
    """Convenience alias for :meth:`IntervalDifference.converter`."""

    return IntervalDifference.converter(ctx)


scalar_converter = make_converter


__all__ = ["IntervalDifference", "make_converter", "scalar_converter"]
