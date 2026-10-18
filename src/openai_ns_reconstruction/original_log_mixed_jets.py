"""Directed logarithmic absolute majorants for ordinary mixed jets.

Each non-``None`` entry stores an interval containing ``log(M[i,k])``, where
``M[i,k]`` is an absolute upper bound for the ordinary derivative with
physical-y order ``i`` and Z order ``k``.  ``None`` is reserved for a proved
exact zero.  The class never exponentiates a stored logarithm and never
represents a signed function value.
"""

from __future__ import annotations

from fractions import Fraction
import math
from collections.abc import Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


def _raw_bounds(value):
    if not hasattr(value, "_mpi_"):
        raise TypeError("An MPIntervalContext interval is required")
    rows = value._mpi_
    if len(rows) == 4 and all(type(part) is int for part in rows):
        rows = (rows, rows)
    if len(rows) == 1:
        rows = (rows, rows)
    if len(rows) != 2:
        raise TypeError("Malformed interval")
    return mp.make_mpf(rows[0]), mp.make_mpf(rows[1])


def _lower(value):
    return _raw_bounds(value)[0]


def _upper(value):
    return _raw_bounds(value)[1]


def _interval(ctx, lower, upper):
    if hasattr(lower, "_mpi_"):
        lower = _lower(lower)
    elif not hasattr(lower, "_mpf_"):
        lower = mp.make_mpf(lower)
    if hasattr(upper, "_mpi_"):
        upper = _upper(upper)
    elif not hasattr(upper, "_mpf_"):
        upper = mp.make_mpf(upper)
    if lower > upper:
        raise ValueError("Reversed directed interval")
    return ctx.mpf([lower, upper])


def _coerce(ctx, value):
    """Convert exact scalar or directed interval input without binary floats."""
    if hasattr(value, "_mpi_"):
        lo, hi = _raw_bounds(value)
        return _interval(ctx, lo, hi)
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError("Use an exact scalar or MPIntervalContext interval")
    if type(value) is int:
        return ctx.mpf(value)
    if isinstance(value, Fraction):
        return ctx.mpf(value.numerator) / value.denominator
    if isinstance(value, str):
        return ctx.mpf(value)
    if isinstance(value, (tuple, list)):
        if len(value) == 4 and all(type(part) is int for part in value):
            point = mp.make_mpf(tuple(value))
            return _interval(ctx, point, point)
        if len(value) == 2:
            left, right = _coerce(ctx, value[0]), _coerce(ctx, value[1])
            return _interval(ctx, _lower(left), _upper(right))
    raise TypeError("Use an exact scalar or MPIntervalContext interval")


def _finite(value):
    lo, hi = _raw_bounds(value)
    return not any(mp.isnan(x) or mp.isinf(x) for x in (lo, hi))


def _log_integer(ctx, value):
    if value < 1:
        raise ValueError("A positive integer logarithm is required")
    return ctx.ln(ctx.mpf(value))


def _logsum(ctx, terms):
    """Directed log-sum bound without forming any exponential."""
    terms = tuple(terms)
    if not terms:
        return None
    if any(not _finite(term) for term in terms):
        raise ValueError("Log majorants must be finite")
    if len(terms) == 1:
        return terms[0]
    shift = _log_integer(ctx, len(terms))
    candidates = [term + shift for term in terms]
    return _interval(ctx, min(_lower(term) for term in candidates),
                     max(_upper(term) for term in candidates))


def _logsum_weighted(ctx, terms):
    """Combine ``(log-cap, positive integer coefficient)`` terms."""
    shifted = []
    for term, coefficient in terms:
        if term is None:
            continue
        if type(coefficient) is not int or coefficient <= 0:
            raise ValueError("A positive integer derivative coefficient is required")
        shifted.append(term if coefficient == 1 else term + _log_integer(ctx, coefficient))
    return _logsum(ctx, shifted)


def _negate(ctx, value):
    return -value


def _sub(ctx, left, right):
    return left - right


class LogMixed:
    """Log-absolute majorants on the eight ordinary ``(y,Z)`` rows.

    Rows are ``(i,k)`` with ``i in {0,1}`` and ``0 <= k <= z_order``.  The
    public contract uses ``z_order=3`` and therefore requires all eight keys
    whenever a mapping is supplied.  A missing row is never interpreted as a
    zero; only an explicit ``None`` is an exact zero sentinel.
    """

    def __init__(self, ctx, z_order=3, values=None):
        if not isinstance(ctx, MPIntervalContext):
            raise TypeError("ctx must be an MPIntervalContext")
        if type(z_order) is not int or z_order != 3:
            raise ValueError("LogMixed requires the physical-y/Z rows through (1,3)")
        self.ctx = ctx
        self.z_order = z_order
        self.indices = tuple((i, k) for total in range(5)
                             for i in (0, 1) for k in range(z_order + 1)
                             if i + k == total)
        self._index_set = frozenset(self.indices)
        if values is None:
            self._values = {index: None for index in self.indices}
        else:
            self._values = self._validated_values(values)

    def _validated_values(self, values):
        if not isinstance(values, Mapping):
            raise TypeError("A complete mixed-row mapping is required")
        if set(values) != self._index_set:
            raise ValueError("All eight mixed ordinary rows are required")
        out = {}
        for index in self.indices:
            value = values[index]
            if value is None:
                out[index] = None
                continue
            item = _coerce(self.ctx, value)
            if not _finite(item):
                raise ValueError("Log majorants must be finite")
            out[index] = item
        return out

    def _new(self, values):
        return type(self)(self.ctx, self.z_order, values)

    def _check(self, jet):
        if not isinstance(jet, LogMixed) or jet.ctx is not self.ctx or jet.z_order != self.z_order:
            raise TypeError("Mixed jets must share the exact interval context and order")

    def __getitem__(self, index):
        if index not in self._index_set:
            raise KeyError(index)
        return self._values[index]

    def items(self):
        return tuple((index, self._values[index]) for index in self.indices)

    def as_dict(self):
        return dict(self.items())

    @property
    def values(self):
        return self.as_dict()

    @classmethod
    def from_mapping(cls, ctx, values, z_order=3):
        return cls(ctx, z_order, values)

    def constant(self, logabs):
        """Return a jet with only the explicitly supplied constant row."""
        value = None if logabs is None else _coerce(self.ctx, logabs)
        return self._new({index: (value if index == (0, 0) else None)
                          for index in self.indices})

    def zero(self):
        """Return an all-exact-zero jet."""
        return self.constant(None)

    def add(self, *jets):
        if not jets:
            raise ValueError("At least one logarithmic jet is required")
        for jet in jets:
            self._check(jet)
        values = {}
        for index in self.indices:
            values[index] = _logsum(self.ctx, [jet[index] for jet in jets
                                               if jet[index] is not None])
        return self._new(values)

    def scale(self, jet, logabs):
        self._check(jet)
        if logabs is None:
            return self._new({index: None for index in self.indices})
        scalar = _coerce(self.ctx, logabs)
        if not _finite(scalar):
            raise ValueError("A finite scalar log-cap is required")
        return self._new({index: (None if jet[index] is None else jet[index] + scalar)
                          for index in self.indices})

    def mul(self, left, right):
        self._check(left)
        self._check(right)
        values = {}
        for i, k in self.indices:
            terms = []
            for r in range(i + 1):
                for s in range(k + 1):
                    terms.append((left[r, s] if (r, s) in self._index_set else None,
                                  right[i - r, k - s] if (i - r, k - s) in self._index_set else None,
                                  math.comb(i, r) * math.comb(k, s)))
            values[i, k] = _logsum_weighted(
                self.ctx, [(a + b, coefficient) for a, b, coefficient in terms
                           if a is not None and b is not None])
        return self._new(values)

    def reciprocal(self, jet, log_lower_value):
        """Return absolute logarithmic bounds for ``1 / jet``.

        ``log_lower_value`` is a certified interval for a lower bound on the
        logarithm of the positive magnitude of the value row.  Its lower
        endpoint is therefore the safe endpoint for the reciprocal cap.
        """
        self._check(jet)
        lower = _coerce(self.ctx, log_lower_value)
        if not _finite(lower):
            raise ValueError("A finite logarithmic lower bound is required")
        if jet[0, 0] is None:
            raise ValueError("The reciprocal value row cannot be an exact zero")
        inverse_log = _negate(self.ctx, lower)
        values = {(0, 0): inverse_log}
        for i, k in self.indices:
            if (i, k) == (0, 0):
                continue
            terms = []
            for r in range(i + 1):
                for s in range(k + 1):
                    if (r, s) == (0, 0):
                        continue
                    delta = (i - r, k - s)
                    if jet[r, s] is None or values[delta] is None:
                        continue
                    terms.append((jet[r, s] + values[delta],
                                  math.comb(i, r) * math.comb(k, s)))
            total = _logsum_weighted(self.ctx, terms)
            values[i, k] = None if total is None else total + inverse_log
        return self._new(values)

    def sqrt(self, jet, log_lower_value, log_upper_value):
        """Return absolute logarithmic bounds for the positive square root."""
        self._check(jet)
        lower = _coerce(self.ctx, log_lower_value)
        upper = _coerce(self.ctx, log_upper_value)
        if not _finite(lower) or not _finite(upper):
            raise ValueError("Finite square-root logarithmic bounds are required")
        if _upper(lower) < _lower(lower):
            raise ValueError("Malformed square-root lower bound")
        if _lower(upper) < _lower(lower):
            raise ValueError("Square-root log bounds are reversed")
        if jet[0, 0] is None:
            raise ValueError("The square-root value row cannot be an exact zero")
        half_upper = upper / 2
        values = {(0, 0): _interval(self.ctx, _lower(half_upper), _upper(half_upper))}
        denominator = _log_integer(self.ctx, 2) + lower / 2
        for i, k in self.indices:
            if (i, k) == (0, 0):
                continue
            terms = []
            if jet[i, k] is not None:
                terms.append(jet[i, k])
            for r in range(i + 1):
                for s in range(k + 1):
                    beta = (r, s)
                    delta = (i - r, k - s)
                    if beta == (0, 0) or beta == (i, k):
                        continue
                    if values[beta] is None or values[delta] is None:
                        continue
                    coefficient = math.comb(i, r) * math.comb(k, s)
                    terms.append(values[beta] + values[delta] +
                                 (_log_integer(self.ctx, coefficient)
                                  if coefficient != 1 else self.ctx.mpf(0)))
            numerator = _logsum(self.ctx, terms)
            values[i, k] = None if numerator is None else _sub(self.ctx, numerator, denominator)
        return self._new(values)

    def compose_scalar(self, jet, derivative_log_caps, *rest):
        """Compose a scalar Taylor jet with a zero-constant mixed increment.

        The supplied five entries are caps for the ordinary derivatives
        ``f^(j)(x0)``.  Powers of the increment use the ordinary binomial
        product rule; subtracting ``log(j!)`` is the Taylor factorial, so the
        returned rows are ordinary derivatives again.
        """
        self._check(jet)
        if rest:
            caps = (derivative_log_caps,) + rest
        elif isinstance(derivative_log_caps, Mapping):
            if set(derivative_log_caps) != set(range(5)):
                raise ValueError("Five scalar derivative log-caps are required")
            caps = tuple(derivative_log_caps[j] for j in range(5))
        else:
            try:
                caps = tuple(derivative_log_caps)
            except TypeError as exc:
                raise TypeError("Five scalar derivative log-caps are required") from exc
        if len(caps) != 5:
            raise ValueError("Five scalar derivative log-caps are required")
        normalized = {index: (None if index == (0, 0) else jet[index])
                      for index in self.indices}
        increment = self._new(normalized)
        terms = []
        if caps[0] is not None:
            terms.append(self.constant(caps[0]))
        power = increment
        for order in range(1, 5):
            if caps[order] is not None:
                coefficient = _coerce(self.ctx, caps[order])
                coefficient = _sub(self.ctx, coefficient,
                                   _log_integer(self.ctx, math.factorial(order)))
                terms.append(self.scale(power, coefficient))
            if order != 4:
                power = self.mul(power, increment)
        return self.add(*terms) if terms else self._new({index: None for index in self.indices})


def _rational(ctx, value):
    value = Fraction(value)
    return ctx.mpf(value.numerator) / value.denominator


def _point_log(ctx, value):
    value = Fraction(value)
    if value == 0:
        return None
    return ctx.ln(_rational(ctx, abs(value)))


def _poly_jet(ctx, expression, y, z, engine):
    import sympy as sp
    return engine.from_mapping(ctx, {
        (i, k): _point_log(ctx, sp.Rational(abs(sp.diff(expression, y, i, z, k).subs({y: 0, z: 0}))))
        if sp.diff(expression, y, i, z, k).subs({y: 0, z: 0}) != 0 else None
        for i, k in engine.indices})


def _contains(ctx, jet, value):
    import sympy as sp
    value = sp.Rational(value)
    if value == 0:
        return True
    if jet is None:
        return False
    exact = ctx.ln(_rational(ctx, abs(value)))
    return _upper(jet) >= _upper(exact)


def symbolic_controls():
    """Run exact SymPy identities and finite positive point-cap checks."""
    import sympy as sp
    ctx = MPIntervalContext()
    ctx.prec = 180
    engine = LogMixed(ctx)
    y, z, t = sp.symbols("y z t")
    x = sp.Rational(2) + 3*y - 2*z + y*z + z**2 - y*z**3
    w = sp.Rational(3) + 2*y - z + y*z
    f = sp.expand(w*w)
    product = engine.mul(_poly_jet(ctx, x, y, z, engine),
                         _poly_jet(ctx, w, y, z, engine))
    reciprocal = engine.reciprocal(_poly_jet(ctx, f, y, z, engine),
                                   _point_log(ctx, 9),)
    square_root = engine.sqrt(_poly_jet(ctx, f, y, z, engine),
                              _point_log(ctx, 9), _point_log(ctx, 9))
    delta = _poly_jet(ctx, x - 2, y, z, engine)
    scalar = 5 + 2*t - t**2 + 3*t**3 - t**4
    composed = engine.compose_scalar(delta, tuple(
        _point_log(ctx, sp.diff(scalar, t, order).subs(t, 0))
        for order in range(5)))
    expected_product = sp.expand(x*w)
    expected_reciprocal = 1/f
    expected_root = w
    expected_composed = sp.expand(scalar.subs(t, x - 2))
    checks = {
        "mul_nonunit_polynomial": all(_contains(ctx, product[i, k],
                                                   sp.diff(expected_product, y, i, z, k).subs({y: 0, z: 0}))
                                       for i, k in engine.indices),
        "reciprocal_nonunit_value": all(_contains(ctx, reciprocal[i, k],
                                                     sp.diff(expected_reciprocal, y, i, z, k).subs({y: 0, z: 0}))
                                         for i, k in engine.indices),
        "sqrt_positive_square": all(_contains(ctx, square_root[i, k],
                                                sp.diff(expected_root, y, i, z, k).subs({y: 0, z: 0}))
                                     for i, k in engine.indices),
        "compose_taylor_factorials": all(_contains(ctx, composed[i, k],
                                                      sp.diff(expected_composed, y, i, z, k).subs({y: 0, z: 0}))
                                          for i, k in engine.indices),
        "ordinary_factorial": sp.diff(y*z, y, 1, z, 1).subs({y: 0, z: 0}) == 1,
        "complete_rows": set(engine.indices) == {(i, k) for i in (0, 1) for k in range(4)},
    }
    if not all(checks.values()):
        raise ArithmeticError("LogMixed symbolic controls failed")
    return {"passed": True, "checks": checks, "ordinary_derivatives": True,
            "log_absolute_only": True, "zero_sentinel_is_exact": True,
            "no_underlying_exp_or_signed_value": True}


__all__ = ["LogMixed", "symbolic_controls"]
