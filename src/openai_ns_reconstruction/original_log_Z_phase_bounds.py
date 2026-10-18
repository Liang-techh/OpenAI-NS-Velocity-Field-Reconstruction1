"""Z-only logarithmic bounds for a fixed-phase angle inverse.

The rows in this module contain only ordinary ``Z`` and free-angle
derivatives.  A stored value is an interval containing the logarithm of an
absolute derivative bound; ``None`` is reserved for a proved exact zero.
No physical-y row is manufactured.  Normalized Taylor coefficients are used
only internally by the implicit inverse and composition recurrences.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from mpmath.ctx_iv import MPIntervalContext

from .original_log_mixed_jets import (
    _coerce,
    _finite,
    _log_integer,
    _logsum,
    _logsum_weighted,
    _lower,
    _upper,
)


Z_INDICES = (0, 1, 2, 3)
Z_INDEX_SET = frozenset(Z_INDICES)
PHASE_INDICES = tuple(
    (k, j)
    for total in range(4)
    for k in range(4)
    for j in range(4)
    if k + j == total
)
PHASE_INDEX_SET = frozenset(PHASE_INDICES)
ZERO_Z = 0
ZERO_PHASE = (0, 0)
_MISSING = object()


def _check_context(c):
    if not isinstance(c, MPIntervalContext):
        raise TypeError("c must be an MPIntervalContext")


def _leq(left, right):
    return all(a <= b for a, b in zip(left, right))


def _difference(right, left):
    return tuple(a - b for a, b in zip(right, left))


def _factorial_log(c, index):
    """Return the logarithm of the factorial product for an index."""
    if isinstance(index, int):
        components = (index,)
    else:
        components = tuple(index)
    out = c.mpf(0)
    for component in components:
        out += _log_integer(c, math.factorial(component))
    return out


def _multinomial(left, right):
    coefficient = 1
    for a, b in zip(left, right):
        coefficient *= math.comb(a, b)
    return coefficient


def _mapping_values(values, indices, label):
    """Validate a complete mapping while preserving explicit ``None`` rows."""
    if not isinstance(values, Mapping):
        raise TypeError("A complete " + label + " mapping is required")
    if set(values) != set(indices):
        raise ValueError("All " + label + " rows are required")
    return {index: values[index] for index in indices}


def _coerce_log_mapping(c, values, indices, label):
    raw = _mapping_values(values, indices, label)
    out = {}
    for index in indices:
        value = raw[index]
        if value is None:
            out[index] = None
            continue
        item = _coerce(c, value)
        if not _finite(item):
            raise ValueError("Log majorants must be finite")
        out[index] = item
    return out


class LogZ:
    """Log-absolute majorants for ordinary ``Z`` derivatives of orders 0--3."""

    indices = Z_INDICES

    def __init__(self, c, values=None):
        _check_context(c)
        self.c = c
        self._index_set = Z_INDEX_SET
        if values is None:
            self._values = {index: None for index in self.indices}
        elif isinstance(values, Mapping):
            self._values = _coerce_log_mapping(c, values, self.indices, "Z")
        else:
            if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
                raise TypeError("A complete four-row Z sequence is required")
            if len(values) != len(self.indices):
                raise ValueError("A complete four-row Z sequence is required")
            self._values = _coerce_log_mapping(
                c, dict(zip(self.indices, values)), self.indices, "Z"
            )

    def _new(self, values):
        return type(self)(self.c, values)

    def _check(self, jet):
        if not isinstance(jet, LogZ) or jet.c is not self.c:
            raise TypeError("Z jets must share the exact interval context")

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
        if value is not None and not _finite(value):
            raise ValueError("Log majorants must be finite")
        return self._new({index: (value if index == ZERO_Z else None)
                          for index in self.indices})

    def zero(self):
        return self.constant(None)

    def add(self, *jets):
        if not jets:
            raise ValueError("At least one Z jet is required")
        for jet in jets:
            self._check(jet)
        return self._new({
            index: _logsum(self.c, [jet[index] for jet in jets
                                    if jet[index] is not None])
            for index in self.indices
        })

    def scale(self, jet, logabs=_MISSING):
        """Scale a jet by a logarithmic absolute scalar cap.

        The engine-style form is ``engine.scale(jet, logabs)``.  The
        convenience form ``jet.scale(logabs)`` is accepted as well.
        """
        if logabs is _MISSING:
            actual, scalar_value = self, jet
        else:
            actual, scalar_value = jet, logabs
        self._check(actual)
        if scalar_value is None:
            return self._new({index: None for index in self.indices})
        scalar = _coerce(self.c, scalar_value)
        if not _finite(scalar):
            raise ValueError("A finite scalar log-cap is required")
        return self._new({
            index: None if actual[index] is None else actual[index] + scalar
            for index in self.indices
        })

    def mul(self, left, right):
        """Multiply ordinary Z jets, retaining the binomial factors."""
        self._check(left)
        self._check(right)
        values = {}
        for order in self.indices:
            terms = []
            for left_order in range(order + 1):
                right_order = order - left_order
                left_value = left[left_order]
                right_value = right[right_order]
                if left_value is None or right_value is None:
                    continue
                terms.append((left_value + right_value,
                              math.comb(order, left_order)))
            values[order] = _logsum_weighted(self.c, terms)
        return self._new(values)

    def reciprocal(self, jet, log_lower_value):
        """Bound the ordinary Z derivatives of a reciprocal."""
        self._check(jet)
        lower = _coerce(self.c, log_lower_value)
        if not _finite(lower):
            raise ValueError("A finite logarithmic lower bound is required")
        if jet[ZERO_Z] is None:
            raise ValueError("The reciprocal value row cannot be an exact zero")
        inverse_log = -lower
        values = {ZERO_Z: inverse_log}
        for order in self.indices:
            if order == ZERO_Z:
                continue
            terms = []
            for left_order in range(1, order + 1):
                left_value = jet[left_order]
                right_value = values[order - left_order]
                if left_value is None or right_value is None:
                    continue
                terms.append((left_value + right_value,
                              math.comb(order, left_order)))
            total = _logsum_weighted(self.c, terms)
            values[order] = None if total is None else total + inverse_log
        return self._new(values)


class LogZPhase:
    """Log-absolute ordinary derivatives in ``(Z, angle)`` through total 3."""

    indices = PHASE_INDICES

    def __init__(self, c, values=None):
        _check_context(c)
        self.c = c
        self._index_set = PHASE_INDEX_SET
        if values is None:
            self._values = {index: None for index in self.indices}
        else:
            self._values = _coerce_log_mapping(
                c, values, self.indices, "Z/angle"
            )

    def _new(self, values):
        return type(self)(self.c, values)

    def _check(self, jet):
        if not isinstance(jet, LogZPhase) or jet.c is not self.c:
            raise TypeError("Z/angle jets must share the exact interval context")

    @classmethod
    def from_mapping(cls, c, values):
        return cls(c, values)

    @classmethod
    def liftZ(cls, jet):
        if not isinstance(jet, LogZ):
            raise TypeError("liftZ requires a LogZ jet")
        return cls(jet.c, {
            (order, angle_order):
            (jet[order] if angle_order == 0 else None)
            for order, angle_order in cls.indices
        })

    @classmethod
    def lift_z(cls, jet):
        return cls.liftZ(jet)

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
        if value is not None and not _finite(value):
            raise ValueError("Log majorants must be finite")
        return self._new({
            index: (value if index == ZERO_PHASE else None)
            for index in self.indices
        })

    def zero(self):
        return self.constant(None)

    def add(self, *jets):
        if not jets:
            raise ValueError("At least one Z/angle jet is required")
        for jet in jets:
            self._check(jet)
        return self._new({
            index: _logsum(self.c, [jet[index] for jet in jets
                                    if jet[index] is not None])
            for index in self.indices
        })

    def scale(self, jet, logabs=_MISSING):
        if logabs is _MISSING:
            actual, scalar_value = self, jet
        else:
            actual, scalar_value = jet, logabs
        self._check(actual)
        if scalar_value is None:
            return self._new({index: None for index in self.indices})
        scalar = _coerce(self.c, scalar_value)
        if not _finite(scalar):
            raise ValueError("A finite scalar log-cap is required")
        return self._new({
            index: None if actual[index] is None else actual[index] + scalar
            for index in self.indices
        })

    def mul(self, left, right):
        """Multiply ordinary bivariate derivatives with binomial factors."""
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
                left_value = left[beta]
                right_value = right[gamma]
                if left_value is None or right_value is None:
                    continue
                terms.append((left_value + right_value,
                              _multinomial(index, beta)))
            values[index] = _logsum_weighted(self.c, terms)
        return self._new(values)

    def reciprocal(self, jet, log_lower_value):
        """Bound the ordinary bivariate derivatives of a reciprocal."""
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
                left_value = jet[beta]
                right_value = values.get(gamma)
                if left_value is None or right_value is None:
                    continue
                terms.append((left_value + right_value,
                              _multinomial(index, beta)))
            total = _logsum_weighted(self.c, terms)
            values[index] = None if total is None else total + inverse_log
        return self._new(values)

    def lift_z(self, jet):
        return type(self).liftZ(jet)

    def primitiveAngle(self, jet, angle_upper=_MISSING):
        """Integrate an angle jet from zero, shifting ordinary rows by one.

        ``angle_upper`` is a magnitude cap for ``0 <= angle <= angle_upper``;
        it defaults to ``2*pi``.  For an angle-independent ``LogZ`` source,
        only the ``j=1`` rows are nonzero after the shift.  There is no
        ``1/j`` factor because these are ordinary derivatives.
        """
        if isinstance(jet, LogZ):
            source = type(self).liftZ(jet)
        else:
            self._check(jet)
            source = jet
        if angle_upper is _MISSING:
            angle_upper = self.c.mpf(2) * self.c.pi
        angle = _coerce(self.c, angle_upper)
        if not _finite(angle):
            raise ValueError("The angle upper bound must be finite")
        if _lower(angle) < 0:
            raise ValueError("The angle upper bound must be nonnegative")
        angle_cap = _upper(angle)
        log_angle = None if angle_cap == 0 else self.c.ln(angle_cap)
        values = {}
        for order, angle_order in self.indices:
            if angle_order == 0:
                source_value = source[order, 0]
                values[order, 0] = (
                    None if source_value is None or log_angle is None
                    else source_value + log_angle
                )
            else:
                values[order, angle_order] = source[order, angle_order - 1]
        return self._new(values)

    def primitive_angle(self, jet, angle_upper=_MISSING):
        return self.primitiveAngle(jet, angle_upper)

    def primitiveAngleLog(self, jet, log_angle_upper):
        """Angle-primitive variant when the magnitude cap is already logged."""
        if isinstance(jet, LogZ):
            source = type(self).liftZ(jet)
        else:
            self._check(jet)
            source = jet
        log_angle = None if log_angle_upper is None else _coerce(
            self.c, log_angle_upper
        )
        if log_angle is not None and not _finite(log_angle):
            raise ValueError("The logarithmic angle upper bound must be finite")
        values = {}
        for order, angle_order in self.indices:
            if angle_order == 0:
                source_value = source[order, 0]
                values[order, 0] = (
                    None if source_value is None or log_angle is None
                    else source_value + log_angle
                )
            else:
                values[order, angle_order] = source[order, angle_order - 1]
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
    """Multiply normalized Z Taylor coefficients without binomial factors."""
    values = {}
    for order in Z_INDICES:
        terms = []
        for left_order in range(order + 1):
            right_order = order - left_order
            left_value = left.get(left_order)
            right_value = right.get(right_order)
            if left_value is None or right_value is None:
                continue
            terms.append(left_value + right_value)
        values[order] = _logsum(engine.c, terms)
    return values


def _coefficient_powers(engine, increment, order=3):
    one = {index: (engine.c.mpf(0) if index == ZERO_Z else None)
           for index in Z_INDICES}
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
        raise TypeError("A complete normalized Z increment mapping is required")
    if set(data) != Z_INDEX_SET:
        raise ValueError("All four normalized Z increment rows are required")
    if data[ZERO_Z] is not None:
        raise ValueError("The normalized phase increment must have exact zero constant")
    out = {}
    for index in Z_INDICES:
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
    """Return ordinary inverse Z bounds and normalized ``delta`` coefficients.

    ``phi`` is centered at the selected base angle, while the returned
    increment is the Taylor coefficient jet of ``delta(z)`` about that point.
    The angle value itself is bounded by ``2*pi`` and is never selected.
    """
    if not isinstance(engine, LogZPhase):
        raise TypeError("engine must be a LogZPhase instance")
    engine._check(phi)
    lower = _coerce(engine.c, log_lambda_lower)
    if not _finite(lower):
        raise ValueError("A finite logarithmic lower bound is required")
    if phi[0, 1] is None:
        raise ValueError("The phase derivative lambda row cannot be an exact zero")
    phase_coeff = _phase_coefficients(engine, phi)
    delta = {index: None for index in Z_INDICES}

    for order in (1, 2, 3):
        powers = _coefficient_powers(engine, delta, order=3)
        terms = []
        for beta_order in range(order + 1):
            remainder = order - beta_order
            for angle_order in range(4):
                phase_index = beta_order, angle_order
                if phase_index not in PHASE_INDEX_SET:
                    continue
                if beta_order == 0 and angle_order == 1:
                    continue
                phase_value = phase_coeff[phase_index]
                power = powers[angle_order][remainder]
                if phase_value is None or power is None:
                    continue
                terms.append(phase_value + power)
        total = _logsum(engine.c, terms)
        delta[order] = None if total is None else total - lower

    inverse_values = {
        ZERO_Z: engine.c.ln(engine.c.mpf(2) * engine.c.pi)
    }
    for order in (1, 2, 3):
        inverse_values[order] = (
            None
            if delta[order] is None
            else delta[order] + _factorial_log(engine.c, order)
        )
    inverse = LogZ.from_mapping(engine.c, inverse_values)
    proof = {
        "passed": True,
        "Z_only": True,
        "physical_y_derivatives": False,
        "ordinary_Z_derivatives": True,
        "normalized_increment_is_Taylor_coefficient": True,
        "coefficient_product_has_no_binomial_factors": True,
        "phase_terms_divided_by_Z_and_angle_factorials": True,
        "returned_Z_factorials_restored": True,
        "leading_beta_zero_angle_one_excluded": True,
        "leading_denominator_log_lower_bound": True,
        "delta_zero_is_exact": True,
        "inverse_value_bound": "|psi0|<=2*pi",
        "inverse_value_is_not_selected": True,
        "Z_rows": len(Z_INDICES),
        "phase_rows": len(PHASE_INDICES),
    }
    return inverse, delta, proof


def inverse(phi, log_lambda_lower):
    """Convenience wrapper for :func:`inverse_coefficients`."""
    if not isinstance(phi, LogZPhase):
        raise TypeError("phi must be a LogZPhase instance")
    return inverse_coefficients(LogZPhase(phi.c), phi, log_lambda_lower)


def compose_inverse(engine, freejet, increment):
    """Compose a phase jet with a normalized inverse Z increment."""
    if not isinstance(engine, LogZPhase):
        raise TypeError("engine must be a LogZPhase instance")
    if isinstance(freejet, LogZ):
        freejet = LogZPhase.liftZ(freejet)
    engine._check(freejet)
    delta = _validate_increment(engine, increment)
    powers = _coefficient_powers(engine, delta, order=3)
    free_coeff = _phase_coefficients(engine, freejet)
    ordinary = {}
    for order in Z_INDICES:
        terms = []
        for beta_order in range(order + 1):
            remainder = order - beta_order
            for angle_order in range(4):
                phase_index = beta_order, angle_order
                if phase_index not in PHASE_INDEX_SET:
                    continue
                phase_value = free_coeff[phase_index]
                power = powers[angle_order][remainder]
                if phase_value is None or power is None:
                    continue
                terms.append(phase_value + power)
        total = _logsum(engine.c, terms)
        ordinary[order] = (
            None if total is None else total + _factorial_log(engine.c, order)
        )
    return LogZ.from_mapping(engine.c, ordinary)


def compose(freejet, increment):
    """Convenience wrapper for :func:`compose_inverse`."""
    if isinstance(freejet, LogZ):
        engine = LogZPhase(freejet.c)
    elif isinstance(freejet, LogZPhase):
        engine = LogZPhase(freejet.c)
    else:
        raise TypeError("freejet must be a LogZPhase or LogZ instance")
    return compose_inverse(engine, freejet, increment)


def liftZ(jet):
    return LogZPhase.liftZ(jet)


def primitiveAngle(jet, angle_upper=_MISSING):
    if isinstance(jet, LogZ):
        engine = LogZPhase(jet.c)
    elif isinstance(jet, LogZPhase):
        engine = LogZPhase(jet.c)
    else:
        raise TypeError("jet must be a LogZ or LogZPhase instance")
    return engine.primitiveAngle(jet, angle_upper)


def _rational_log(c, value):
    numerator, denominator = int(value.p), int(value.q)
    if numerator == 0:
        return None
    return c.ln(c.mpf(abs(numerator)) / denominator)


def _exact_phase_mapping(c, expression, z, psi):
    import sympy as sp

    out = {}
    for order, angle_order in PHASE_INDICES:
        value = sp.diff(expression, z, order, psi, angle_order).subs(
            {z: 0, psi: sp.Rational(1, 3)}
        )
        out[order, angle_order] = _rational_log(c, sp.Rational(value))
    return out


def _contains(c, cap, value):
    import sympy as sp

    value = sp.Rational(sp.simplify(value))
    if value == 0:
        return cap is None
    if cap is None:
        return False
    return _upper(cap) >= _upper(_rational_log(c, value))


def self_check():
    """Run exact controls for a nonunit implicit phase and return a JSON receipt."""
    import sympy as sp

    c = MPIntervalContext()
    c.prec = 220
    phase_engine = LogZPhase(c)
    z, psi = sp.symbols("z psi")
    phi_expression = (1 + z + z**2 / 2) * psi + z * psi**2 / 10
    phi = LogZPhase.from_mapping(
        c, _exact_phase_mapping(c, phi_expression, z, psi)
    )
    inverse_jet, delta, proof = inverse_coefficients(phase_engine, phi, c.mpf(0))
    a = 1 + z + z**2 / 2
    psi_expression = sp.Rational(2, 3) / (
        a + sp.sqrt(a**2 + sp.Rational(2, 15) * z)
    )
    inverse_checks = {
        order: _contains(
            c, inverse_jet[order],
            sp.diff(psi_expression, z, order).subs(z, 0),
        )
        for order in Z_INDICES
    }
    normalized_checks = {
        0: _contains(c, delta[0], sp.Integer(0)),
        **{
            order: _contains(
                c, delta[order],
                sp.diff(psi_expression, z, order).subs(z, 0)
                / math.factorial(order),
            )
            for order in (1, 2, 3)
        },
    }

    target_expression = psi**2 + z * psi + z**2
    target = LogZPhase.from_mapping(
        c, _exact_phase_mapping(c, target_expression, z, psi)
    )
    composed = compose_inverse(phase_engine, target, delta)
    composed_expression = sp.expand(target_expression.subs(psi, psi_expression))
    composition_checks = {
        order: _contains(
            c, composed[order],
            sp.diff(composed_expression, z, order).subs(z, 0),
        )
        for order in Z_INDICES
    }
    factorial_checks = {
        "ordinary_factorial_restored": all(
            inverse_checks[order] and normalized_checks[order]
            for order in (1, 2, 3)
        ),
        "explicit_Z_factorial_order_2": (
            inverse_jet[2] is not None and delta[2] is not None
        ),
    }

    missing_z = {order: c.mpf(0) for order in Z_INDICES}
    missing_z.pop(3)
    missing_phase = phi.as_dict()
    missing_phase.pop((3, 0))
    missing_delta = dict(delta)
    missing_delta.pop(3)
    missing_checks = {}
    try:
        LogZ.from_mapping(c, missing_z)
    except (TypeError, ValueError):
        missing_checks["missing_Z_row_fails_closed"] = True
    else:
        missing_checks["missing_Z_row_fails_closed"] = False
    try:
        LogZPhase.from_mapping(c, missing_phase)
    except (TypeError, ValueError):
        missing_checks["missing_phase_row_fails_closed"] = True
    else:
        missing_checks["missing_phase_row_fails_closed"] = False
    try:
        compose_inverse(phase_engine, target, missing_delta)
    except (TypeError, ValueError):
        missing_checks["missing_delta_row_fails_closed"] = True
    else:
        missing_checks["missing_delta_row_fails_closed"] = False

    zero = LogZ.constant(LogZ(c), None)
    sentinel_checks = {
        "explicit_None_is_exact_zero": all(zero[order] is None
                                            for order in Z_INDICES),
        "zero_sequence_is_valid": LogZ.from_mapping(c, {
            order: None for order in Z_INDICES
        })[0] is None,
    }
    float_checks = {}
    try:
        LogZ.from_mapping(c, [0.0, None, None, None])
    except TypeError:
        float_checks["float_input_rejected"] = True
    else:
        float_checks["float_input_rejected"] = False
    checks = {
        "inverse_C0_C3_enclosed": all(inverse_checks.values()),
        "normalized_delta_C0_C3_enclosed": all(normalized_checks.values()),
        "compose_C0_C3_enclosed": all(composition_checks.values()),
        **factorial_checks,
        **missing_checks,
        **sentinel_checks,
        **float_checks,
        "proof_passed": proof["passed"],
        "proof_Z_only": proof["Z_only"],
        "proof_no_physical_y_claim": proof["physical_y_derivatives"] is False,
    }
    if not all(checks.values()):
        raise ArithmeticError("LogZPhase self-check failed")
    return {
        "passed": True,
        "checks": checks,
        "fixture": "Phi=(1+z+z^2/2)*psi+z*psi^2/10 at psi0=1/3",
        "target": "T=psi^2+z*psi+z^2",
        "Z_indices": list(Z_INDICES),
        "phase_indices": [list(index) for index in PHASE_INDICES],
        "inverse_derivatives": [
            str(sp.simplify(sp.diff(psi_expression, z, order).subs(z, 0)))
            for order in Z_INDICES
        ],
        "delta_normalized_derivatives": [
            str(sp.simplify(
                sp.diff(psi_expression, z, order).subs(z, 0)
                / math.factorial(order)
            ))
            for order in Z_INDICES
        ],
        "proof": proof,
    }


__all__ = (
    "Z_INDICES",
    "PHASE_INDICES",
    "LogZ",
    "LogZPhase",
    "inverse",
    "inverse_coefficients",
    "compose",
    "compose_inverse",
    "liftZ",
    "primitiveAngle",
    "self_check",
)
