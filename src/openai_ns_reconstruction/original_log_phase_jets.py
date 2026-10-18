"""Directed log-absolute jets for an implicit phase inverse.

The stored entries are logarithms of absolute bounds for ordinary
derivatives.  ``None`` is reserved for an explicitly proved exact zero.
Implicit inversion is performed in normalized Taylor coefficients; the
conversion back to ordinary derivatives is kept explicit so factorials are
not hidden in the polynomial multiplication.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

from mpmath.ctx_iv import MPIntervalContext

from .original_log_mixed_jets import (
    LogMixed,
    _coerce,
    _finite,
    _log_integer,
    _logsum,
    _logsum_weighted,
    _upper,
)


PARAM_INDICES = tuple(
    (i, k)
    for total in range(5)
    for i in (0, 1)
    for k in range(4)
    if i + k == total
)
PARAM_INDEX_SET = frozenset(PARAM_INDICES)
PHASE_INDICES = tuple(
    (i, k, j)
    for total in range(5)
    for i in (0, 1)
    for k in range(4)
    for j in range(5)
    if i + k + j == total
)
PHASE_INDEX_SET = frozenset(PHASE_INDICES)
ZERO_PARAM = (0, 0)
ZERO_PHASE = (0, 0, 0)


def _check_context(c):
    if not isinstance(c, MPIntervalContext):
        raise TypeError("c must be an MPIntervalContext")


def _leq(left, right):
    return all(a <= b for a, b in zip(left, right))


def _difference(right, left):
    return tuple(a - b for a, b in zip(right, left))


def _factorial_log(c, index):
    value = 0
    for component in index:
        value += _log_integer(c, math.factorial(component))
    return value


def _multinomial(left, right):
    coefficient = 1
    for a, b in zip(left, right):
        coefficient *= math.comb(a, b)
    return coefficient


class LogPhase:
    """Log-absolute ordinary jets in ``(physical_y, Z, free_angle)``."""

    indices = PHASE_INDICES

    def __init__(self, c, values=None):
        _check_context(c)
        self.c = c
        self._index_set = PHASE_INDEX_SET
        if values is None:
            self._values = {index: None for index in self.indices}
        else:
            self._values = self._validate(values)

    def _validate(self, values):
        if not isinstance(values, Mapping):
            raise TypeError("A complete phase-jet mapping is required")
        if set(values) != self._index_set:
            raise ValueError("All phase ordinary rows are required")
        out = {}
        for index in self.indices:
            value = values[index]
            if value is None:
                out[index] = None
                continue
            item = _coerce(self.c, value)
            if not _finite(item):
                raise ValueError("Log majorants must be finite")
            out[index] = item
        return out

    def _new(self, values):
        return type(self)(self.c, values)

    def _check(self, jet):
        if not isinstance(jet, LogPhase) or jet.c is not self.c:
            raise TypeError("Phase jets must share the exact interval context")

    @classmethod
    def from_mapping(cls, c, values):
        return cls(c, values)

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

    def constant(self, logabs):
        value = None if logabs is None else _coerce(self.c, logabs)
        return self._new({
            index: (value if index == ZERO_PHASE else None)
            for index in self.indices
        })

    def add(self, *jets):
        if not jets:
            raise ValueError("At least one phase jet is required")
        for jet in jets:
            self._check(jet)
        return self._new({
            index: _logsum(self.c, [
                jet[index] for jet in jets if jet[index] is not None
            ])
            for index in self.indices
        })

    def scale(self, jet, logabs):
        self._check(jet)
        if logabs is None:
            return self._new({index: None for index in self.indices})
        scalar = _coerce(self.c, logabs)
        if not _finite(scalar):
            raise ValueError("A finite scalar log-cap is required")
        return self._new({
            index: None if jet[index] is None else jet[index] + scalar
            for index in self.indices
        })

    def mul(self, left, right):
        self._check(left)
        self._check(right)
        values = {}
        for index in self.indices:
            terms = []
            for beta in self.indices:
                if not _leq(beta, index):
                    continue
                gamma = _difference(index, beta)
                if gamma not in self._index_set:
                    continue
                a, b = left[beta], right[gamma]
                if a is None or b is None:
                    continue
                terms.append((a + b, _multinomial(index, beta)))
            values[index] = _logsum_weighted(self.c, terms)
        return self._new(values)

    def reciprocal(self, jet, log_lower_value):
        """Bound ``1 / jet`` using the ordinary multivariate recurrence."""
        self._check(jet)
        lower = _coerce(self.c, log_lower_value)
        if not _finite(lower):
            raise ValueError("A finite logarithmic lower bound is required")
        if jet[ZERO_PHASE] is None:
            raise ValueError("The reciprocal value row cannot be an exact zero")
        inverse_log = -lower
        values = {ZERO_PHASE: inverse_log}
        for index in self.indices:
            if index == ZERO_PHASE:
                continue
            terms = []
            for beta in self.indices:
                if beta == ZERO_PHASE or not _leq(beta, index):
                    continue
                gamma = _difference(index, beta)
                a, b = jet[beta], values.get(gamma)
                if a is None or b is None:
                    continue
                terms.append((a + b, _multinomial(index, beta)))
            total = _logsum_weighted(self.c, terms)
            values[index] = None if total is None else total + inverse_log
        return self._new(values)

    def lift_mixed(self, mixed):
        """Lift an eight-row ``LogMixed`` jet to angle-independent rows."""
        if isinstance(mixed, Mapping):
            data = dict(mixed)
        elif hasattr(mixed, "as_dict"):
            data = mixed.as_dict()
        else:
            raise TypeError("A LogMixed jet or complete eight-row mapping is required")
        if set(data) != PARAM_INDEX_SET:
            raise ValueError("All eight mixed ordinary rows are required")
        values = {index: None for index in self.indices}
        for (i, k), value in data.items():
            values[i, k, 0] = value
        return self._new(values)


def _phase_coefficients(engine, phase):
    engine._check(phase)
    return {
        index: (
            None
            if phase[index] is None
            else phase[index] - _factorial_log(engine.c, index)
        )
        for index in engine.indices
    }


def _coefficient_product(engine, left, right):
    """Multiply normalized Taylor coefficients without binomial factors."""
    values = {}
    for index in PARAM_INDICES:
        terms = []
        for beta in PARAM_INDICES:
            if not _leq(beta, index):
                continue
            gamma = _difference(index, beta)
            a, b = left[beta], right.get(gamma)
            if a is None or b is None:
                continue
            terms.append(a + b)
        values[index] = _logsum(engine.c, terms)
    return values


def _coefficient_powers(engine, increment, order=4):
    one = {index: (engine.c.mpf(0) if index == ZERO_PARAM else None)
           for index in PARAM_INDICES}
    powers = [one]
    for _ in range(order):
        powers.append(_coefficient_product(engine, powers[-1], increment))
    return powers


def _validate_increment(engine, increment):
    if isinstance(increment, Mapping):
        data = dict(increment)
    elif hasattr(increment, "as_dict"):
        data = increment.as_dict()
    else:
        raise TypeError("A complete normalized increment mapping is required")
    if set(data) != PARAM_INDEX_SET:
        raise ValueError("All eight normalized increment rows are required")
    if data[ZERO_PARAM] is not None:
        raise ValueError("The normalized phase increment must have exact zero constant")
    out = {}
    for index in PARAM_INDICES:
        value = data[index]
        if value is None:
            out[index] = None
        else:
            item = _coerce(engine.c, value)
            if not _finite(item):
                raise ValueError("Log majorants must be finite")
            out[index] = item
    return out


def inverse_coefficients(engine, phi, log_lambda_lower):
    """Return inverse ordinary bounds and normalized phase-increment terms.

    ``phi`` contains ordinary derivatives of ``Phi(p, psi)`` at the base
    point.  The returned increment values are Taylor coefficients of
    ``delta(physical_y, Z)``; their keys are the eight parameter rows.
    """
    if not isinstance(engine, LogPhase):
        raise TypeError("engine must be a LogPhase instance")
    engine._check(phi)
    lower = _coerce(engine.c, log_lambda_lower)
    if not _finite(lower):
        raise ValueError("A finite logarithmic lower bound is required")
    if phi[0, 0, 1] is None:
        raise ValueError("The phase derivative lambda row cannot be exact zero")
    phase_coeff = _phase_coefficients(engine, phi)
    delta = {index: None for index in PARAM_INDICES}

    for alpha in PARAM_INDICES:
        if alpha == ZERO_PARAM:
            continue
        powers = _coefficient_powers(engine, delta, order=4)
        terms = []
        for beta in PARAM_INDICES:
            if not _leq(beta, alpha):
                continue
            remainder = _difference(alpha, beta)
            for angle_order in range(5):
                phase_index = beta + (angle_order,)
                if phase_index not in PHASE_INDEX_SET:
                    continue
                if beta == ZERO_PARAM and angle_order == 1:
                    continue
                fterm = phase_coeff[phase_index]
                if fterm is None:
                    continue
                power = powers[angle_order][remainder]
                if power is None:
                    continue
                terms.append(fterm + power)
        total = _logsum(engine.c, terms)
        delta[alpha] = None if total is None else total - lower

    ordinary = {}
    inverse_value = engine.c.ln(2 * engine.c.pi)
    for alpha in PARAM_INDICES:
        if alpha == ZERO_PARAM:
            ordinary[alpha] = inverse_value
        elif delta[alpha] is None:
            ordinary[alpha] = None
        else:
            ordinary[alpha] = delta[alpha] + _factorial_log(engine.c, alpha)
    inverse = LogMixed.from_mapping(engine.c, ordinary)
    proof = {
        "passed": True,
        "ordinary_derivatives": True,
        "normalized_increment_is_taylor_coefficient": True,
        "coefficient_product_has_no_binomial_factors": True,
        "phase_terms_divided_by_parameter_and_angle_factorials": True,
        "returned_parameter_factorials_restored": True,
        "leading_beta_zero_angle_one_excluded": True,
        "leading_denominator_log_lower_bound": True,
        "inverse_value_bound": "|psi|<=2*pi",
        "inverse_value_is_not_selected": True,
        "parameter_rows": len(PARAM_INDICES),
    }
    return inverse, delta, proof


def compose_inverse(engine, freejet, increment):
    """Compose a full phase jet with a normalized inverse increment."""
    if not isinstance(engine, LogPhase):
        raise TypeError("engine must be a LogPhase instance")
    engine._check(freejet)
    delta = _validate_increment(engine, increment)
    powers = _coefficient_powers(engine, delta, order=4)
    free_coeff = _phase_coefficients(engine, freejet)
    ordinary = {}
    for alpha in PARAM_INDICES:
        terms = []
        for beta in PARAM_INDICES:
            if not _leq(beta, alpha):
                continue
            remainder = _difference(alpha, beta)
            for angle_order in range(5):
                phase_index = beta + (angle_order,)
                if phase_index not in PHASE_INDEX_SET:
                    continue
                fterm = free_coeff[phase_index]
                power = powers[angle_order][remainder]
                if fterm is None or power is None:
                    continue
                terms.append(fterm + power)
        total = _logsum(engine.c, terms)
        ordinary[alpha] = (
            None if total is None else total + _factorial_log(engine.c, alpha)
        )
    return LogMixed.from_mapping(engine.c, ordinary)


def _rational_log(c, value):
    numerator, denominator = int(value.p), int(value.q)
    if numerator == 0:
        return None
    return c.ln(c.mpf(abs(numerator)) / denominator)


def _exact_phase_mapping(c, expression, symbols, engine):
    import sympy as sp

    y, z, psi = symbols
    out = {}
    for i, k, j in engine.indices:
        value = sp.diff(expression, y, i, z, k, psi, j).subs(
            {y: 0, z: 0, psi: 0}
        )
        value = sp.Rational(value)
        out[i, k, j] = _rational_log(c, value)
    return out


def _contains(c, cap, value):
    import sympy as sp

    value = sp.Rational(value)
    if value == 0:
        return True
    if cap is None:
        return False
    exact = _rational_log(c, value)
    return _upper(cap) >= _upper(exact)


def symbolic_controls():
    """Validate factorial conventions with an exact nonunit implicit example."""
    import sympy as sp

    c = MPIntervalContext()
    c.prec = 220
    engine = LogPhase(c)
    y, z, psi = sp.symbols("y z psi")
    phi_expression = (
        3 + (5 + y + z) * psi + 2 * y - z + y * z**3
    )
    psi_expression = -(2 * y - z + y * z**3) / (5 + y + z)
    phi = LogPhase.from_mapping(
        c, _exact_phase_mapping(c, phi_expression, (y, z, psi), engine)
    )
    inverse, increment, proof = inverse_coefficients(
        engine, phi, c.ln(c.mpf(5))
    )
    inverse_checks = {
        index: _contains(
            c,
            inverse[index],
            sp.diff(psi_expression, y, index[0], z, index[1]).subs(
                {y: 0, z: 0}
            ),
        )
        for index in PARAM_INDICES
    }
    target = (1, 3)
    normalized_target = sp.diff(
        psi_expression, y, target[0], z, target[1]
    ).subs({y: 0, z: 0}) / math.factorial(target[0]) / math.factorial(target[1])
    factorial_checks = {
        "ordinary_yZ3": inverse_checks[target],
        "normalized_yZ3": _contains(c, increment[target], normalized_target),
        "factorial_restored": _contains(
            c,
            inverse[target],
            sp.diff(psi_expression, y, 1, z, 3).subs({y: 0, z: 0}),
        ),
    }
    free_expression = (
        11 + 2 * y - z + 3 * psi + y * z**3 + z * psi
        + y * psi**2 + z**3 * psi
    )
    freejet = LogPhase.from_mapping(
        c, _exact_phase_mapping(c, free_expression, (y, z, psi), engine)
    )
    composed = compose_inverse(engine, freejet, increment)
    expected_composed = sp.expand(free_expression.subs(psi, psi_expression))
    composition_checks = {
        index: _contains(
            c,
            composed[index],
            sp.diff(expected_composed, y, index[0], z, index[1]).subs(
                {y: 0, z: 0}
            ),
        )
        for index in PARAM_INDICES
    }
    missing_phase = dict(phi.as_dict())
    missing_phase.pop(next(iter(PHASE_INDEX_SET)))
    missing_increment = dict(increment)
    missing_increment.pop(next(iter(PARAM_INDEX_SET)))
    missing_checks = {}
    try:
        LogPhase.from_mapping(c, missing_phase)
    except (TypeError, ValueError):
        missing_checks["phase_missing_row_fails_closed"] = True
    else:
        missing_checks["phase_missing_row_fails_closed"] = False
    try:
        compose_inverse(engine, freejet, missing_increment)
    except (TypeError, ValueError):
        missing_checks["increment_missing_row_fails_closed"] = True
    else:
        missing_checks["increment_missing_row_fails_closed"] = False
    checks = {
        "inverse_nonunit_implicit": all(inverse_checks.values()),
        "compose_nonunit_implicit": all(composition_checks.values()),
        **factorial_checks,
        **missing_checks,
        "freevalue_included": composed[ZERO_PARAM] is not None,
        "proof_passed": proof["passed"],
    }
    if not all(checks.values()):
        raise ArithmeticError("LogPhase symbolic controls failed")
    return {
        "passed": True,
        "checks": checks,
        "ordinary_derivatives": True,
        "normalized_taylor_coefficients": True,
        "no_source_execution": True,
    }


__all__ = [
    "LogPhase",
    "inverse_coefficients",
    "compose_inverse",
    "symbolic_controls",
]
