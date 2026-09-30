"""Small directed interval Taylor algebra for bounded local calculations.

The coefficients are ordinary Taylor coefficients, so coefficient ``k`` is
the ``k``-th derivative divided by ``k!``.  This module only supplies finite
interval algebra and a symmetric-integral remainder wrapper; it makes no
claim about the source construction or a global source certificate.
"""

from __future__ import annotations

import operator


def _endpoints(value):
    import mpmath as mp

    return tuple(mp.make_mpf(item) for item in value._mpi_)


def _order(value):
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError("Taylor order must be an integer") from exc
    if result < 0:
        raise ValueError("Taylor order must be nonnegative")
    return int(result)


def _interval(ctx, value):
    return value if hasattr(value, "_mpi_") else ctx.mpf(value)


class IntervalTaylor:
    """Finite Taylor coefficients carried as intervals in one context."""

    def __init__(self, ctx, coeffs):
        if isinstance(coeffs, IntervalTaylor):
            if coeffs.ctx is not ctx:
                raise ValueError("Taylor contexts must match")
            coeffs = coeffs.coefficients
        values = tuple(_interval(ctx, value) for value in coeffs)
        if not values:
            raise ValueError("at least one Taylor coefficient is required")
        self.ctx = ctx
        self.coefficients = values
        self.order = len(values) - 1

    @classmethod
    def constant(cls, ctx, value, order=0):
        return constant(ctx, value, order)

    @classmethod
    def variable(cls, ctx, value, order):
        return variable(ctx, value, order)

    def truncate(self, order):
        order = min(self.order, _order(order))
        return IntervalTaylor(self.ctx, self.coefficients[: order + 1])

    def _pair(self, other):
        if isinstance(other, IntervalTaylor):
            if other.ctx is not self.ctx:
                raise ValueError("Taylor contexts must match")
            order = min(self.order, other.order)
            return self.truncate(order), other.truncate(order)
        return self, constant(self.ctx, other, self.order)

    def __getitem__(self, index):
        return self.coefficients[index]

    def __add__(self, other):
        left, right = self._pair(other)
        return IntervalTaylor(
            self.ctx,
            [left[k] + right[k] for k in range(left.order + 1)],
        )

    __radd__ = __add__

    add = __add__

    def __neg__(self):
        return IntervalTaylor(self.ctx, [-value for value in self.coefficients])

    def __sub__(self, other):
        left, right = self._pair(other)
        return IntervalTaylor(
            self.ctx,
            [left[k] - right[k] for k in range(left.order + 1)],
        )

    sub = __sub__

    def __rsub__(self, other):
        return constant(self.ctx, other, self.order) - self

    def __mul__(self, other):
        left, right = self._pair(other)
        coefficients = []
        for degree in range(left.order + 1):
            value = self.ctx.mpf(0)
            for index in range(degree + 1):
                value += left[index] * right[degree - index]
            coefficients.append(value)
        return IntervalTaylor(self.ctx, coefficients)

    __rmul__ = __mul__

    mul = __mul__

    def reciprocal(self):
        """Return the reciprocal series through the retained order."""

        ctx = self.ctx
        result = [ctx.mpf(1) / self[0]]
        for degree in range(1, self.order + 1):
            convolution = ctx.mpf(0)
            for index in range(1, degree + 1):
                convolution += self[index] * result[degree - index]
            result.append(-convolution / self[0])
        return IntervalTaylor(ctx, result)

    inverse = reciprocal
    inv = reciprocal

    def __truediv__(self, other):
        left, right = self._pair(other)
        return left * right.reciprocal()

    div = __truediv__

    def __rtruediv__(self, other):
        return constant(self.ctx, other, self.order) / self

    def __pow__(self, exponent):
        try:
            exponent = operator.index(exponent)
        except TypeError as exc:
            raise TypeError("Taylor powers must be integers") from exc
        if exponent < 0:
            return self.reciprocal() ** (-exponent)
        result = constant(self.ctx, 1, self.order)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            exponent //= 2
            if exponent:
                base = base * base
        return result

    def exp(self):
        """Return the exponential series through the retained order."""

        ctx = self.ctx
        result = [ctx.exp(self[0])]
        for degree in range(1, self.order + 1):
            numerator = ctx.mpf(0)
            for index in range(1, degree + 1):
                numerator += ctx.mpf(index) * self[index] * result[degree - index]
            result.append(numerator / ctx.mpf(degree))
        return IntervalTaylor(ctx, result)

    def __repr__(self):
        return f"IntervalTaylor(order={self.order}, coefficients={self.coefficients!r})"


def constant(ctx, value, order=0):
    """Return an order-``order`` constant Taylor series."""

    order = _order(order)
    zero = ctx.mpf(0)
    return IntervalTaylor(ctx, [_interval(ctx, value)] + [zero] * order)


def variable(ctx, value, order):
    """Return a Taylor variable with interval C0 value and unit slope."""

    order = _order(order)
    zero = ctx.mpf(0)
    coefficients = [_interval(ctx, value), ctx.mpf(1)]
    coefficients.extend(zero for _ in range(max(0, order - 1)))
    return IntervalTaylor(ctx, coefficients[: order + 1])


def genericintegrate_symmetric(fcenter_jet, derivative_cell_jet, radius):
    """Enclose the integral of a local Taylor function over ``[-radius,radius]``.

    The common convention is that both jets have order ``N``: center terms
    ``0`` through ``N - 1`` are integrated exactly, while coefficient ``N``
    from the cell jet supplies the absolute remainder bound.  A constant cell
    jet is also accepted as a direct bound for that coefficient.  If the
    center jet stores only terms through ``N - 1`` and the cell jet has order
    ``N``, that representation is accepted as well.
    """

    if not isinstance(fcenter_jet, IntervalTaylor):
        raise TypeError("fcenter_jet must be an IntervalTaylor")
    ctx = fcenter_jet.ctx
    center_order = fcenter_jet.order
    if isinstance(derivative_cell_jet, IntervalTaylor):
        if derivative_cell_jet.ctx is not ctx:
            raise ValueError("Taylor contexts must match")
        cell = derivative_cell_jet
        if cell.order == center_order + 1:
            remainder_order = center_order + 1
            cell_coefficient = cell[remainder_order]
        elif cell.order >= center_order:
            remainder_order = center_order
            cell_coefficient = cell[remainder_order]
        elif cell.order == 0:
            remainder_order = center_order
            cell_coefficient = cell[0]
        else:
            raise ValueError("derivative cell jet has no requested remainder coefficient")
    else:
        remainder_order = center_order
        cell_coefficient = _interval(ctx, derivative_cell_jet)

    r = _interval(ctx, radius)
    rlo, _ = _endpoints(r)
    if rlo < 0:
        raise ValueError("radius must be nonnegative")

    total = ctx.mpf(0)
    for degree in range(remainder_order):
        if degree % 2:
            continue
        factor = ctx.mpf(2) * (r ** (degree + 1)) / ctx.mpf(degree + 1)
        total += fcenter_jet[degree] * factor

    # Use a positive interval upper endpoint as a point interval, then form
    # the symmetric error entirely through the supplied interval context.
    absolute = abs(cell_coefficient)
    supremum = ctx.mpf(_endpoints(absolute)[1])
    remainder = supremum * ctx.mpf(2) * (r ** (remainder_order + 1)) / ctx.mpf(remainder_order + 1)
    neg = ctx.mpf(0) - remainder
    error = ctx.mpf([_endpoints(neg)[0], _endpoints(remainder)[1]])
    return total + error


# Short integration name used by the radial caller.
integrate_symmetric = genericintegrate_symmetric


__all__ = [
    "IntervalTaylor",
    "constant",
    "variable",
    "genericintegrate_symmetric",
    "integrate_symmetric",
]
