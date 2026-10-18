"""Pure directed finite-``N`` pressure and error algebra.

The public functions in this module consume *raw* logarithmic absolute
majorants.  A row is an interval containing the logarithm of an absolute
bound; ``None`` is the only exact-zero sentinel.  No saved packet is loaded,
and no producer, builder, recurrence, owner, checker, norm, or quadrature is
called here.

Rows use ordinary ``Z`` derivatives.  In particular, the product rule keeps
the binomial factors and never treats a row vector as a Taylor coefficient
vector.  All arithmetic stays in the interval logarithm domain.  The only
signed arithmetic below is the small geometric coefficient calculation used
as a receipt by :func:`D_difference`; it is performed directly on the caller's
intervals and is never selected as a function value.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from mpmath.ctx_iv import MPIntervalContext

from .original_directed_abs_log import abs_log_upper
from .original_log_mixed_jets import (
    _coerce as _coerce_interval,
    _finite as _finite_interval,
    _logsum as _raw_logsum,
)


I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")
D_COEFFICIENT_NAMES = (
    "minus_inv_L",
    "hZ_over_2L",
    "minus_half",
    "d_over_L",
    "two_h_Z2_over_L2",
    "two_Zd_over_L2",
)
_D_TERM_COEFFICIENT = dict(zip(I1_NAMES, D_COEFFICIENT_NAMES))


def safeiv(c: MPIntervalContext, value):
    """Box an exact scalar or directed interval in ``c``.

    Binary floats are deliberately rejected by the shared coercion routine;
    callers should pass exact integers, fractions, strings, MPF values, or
    ``MPIntervalContext`` intervals.
    """

    if not isinstance(c, MPIntervalContext):
        raise TypeError("c must be an MPIntervalContext")
    result = _coerce_interval(c, value)
    if not _finite_interval(result):
        raise ValueError("finite directed intervals are required")
    return result


def _safeiv(c, value):
    return safeiv(c, value)


def _validate_rows(c, rows, label, minimum, *, value_required=False):
    if (not isinstance(rows, Sequence)
            or isinstance(rows, (str, bytes))):
        raise TypeError(label + " must be a row sequence")
    if len(rows) < minimum:
        raise ValueError(label + " must contain at least " + str(minimum)
                         + " rows")
    out = []
    for index, row in enumerate(rows):
        # None is never inferred from an omitted row: it is accepted only when
        # the caller explicitly supplies the exact-zero sentinel.
        out.append(None if row is None else safeiv(c, row))
        if value_required and index == 0 and out[-1] is None:
            raise ValueError(label + " value row cannot be exact zero")
    return out


def _logsum(c, values):
    """Directed logarithmic sum of absolute terms, skipping exact zeros."""

    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise TypeError("logsum requires a term sequence")
    terms = [value for value in values if value is not None]
    if not terms:
        return None
    terms = [safeiv(c, value) for value in terms]
    return _raw_logsum(c, terms)


def logsum(c, *values):
    """Public wrapper for :func:`_logsum` accepting either form.

    ``logsum(c, [a, b])`` and ``logsum(c, a, b)`` are equivalent.
    """

    if len(values) == 1 and isinstance(values[0], Sequence) \
            and not isinstance(values[0], (str, bytes)):
        values = tuple(values[0])
    return _logsum(c, values)


def _log_product(*args):
    """Multiply two logarithmic caps without exponentiating them.

    The two-argument form is native to this module.  A compatibility
    three-argument form, ``_log_product(c, left, right)``, is accepted because
    the historical arithmetic helper carried an unused context parameter.
    """

    if len(args) == 2:
        left, right = args
    elif len(args) == 3:
        _, left, right = args
    else:
        raise TypeError("_log_product expects two logarithmic terms")
    if left is None or right is None:
        return None
    return left + right


def log_product(left, right):
    return _log_product(left, right)


def _abs_log(c, value):
    """Tuple-safe logarithmic absolute cap for a signed interval."""

    return abs_log_upper(c, safeiv(c, value))


def _signed_logs(c, values):
    return [_abs_log(c, value) for value in values]


def _log_scale(c, row, scalar):
    if row is None:
        return None
    return _log_product(row, _abs_log(c, scalar))


def _log_jet_product(c, left, right, order):
    """Ordinary derivative product of two logarithmic jets."""

    if len(left) < order + 1 or len(right) < order + 1:
        raise ValueError("incomplete ordinary derivative product")
    out = []
    for k in range(order + 1):
        terms = []
        for j in range(k + 1):
            term = _log_product(left[j], right[k - j])
            if term is not None:
                coefficient = math.comb(k, j)
                if coefficient != 1:
                    term = _log_product(
                        term, _abs_log(c, c.mpf(coefficient)))
                terms.append(term)
        out.append(_logsum(c, terms))
    return out


def _log_add_jets(c, *jets):
    if not jets:
        raise ValueError("at least one jet is required")
    order = min(len(jet) for jet in jets)
    return [_logsum(c, [jet[k] for jet in jets]) for k in range(order)]


def _signed_product(left, right):
    """C1 ordinary product for signed interval coefficient receipts."""

    if len(left) < 2 or len(right) < 2:
        raise ValueError("signed C1 product requires two rows")
    return [left[0] * right[0], left[1] * right[0] + left[0] * right[1]]


def _signed_scale(row, scalar):
    return [scalar * value for value in row]


def _signed_add(*rows):
    return [sum((row[k] for row in rows), rows[0][k] * 0)
            for k in range(2)]


def _z_signed_jet(c, z, order):
    z = safeiv(c, z)
    values = [z, c.mpf(1), c.mpf(0)]
    return values[:order + 1]


def _geometric_log_jets(c, z, delta, order):
    """Return log caps for ``Z``, ``h=1-delta`` and ``d=1-Z^2``."""

    z = safeiv(c, z)
    delta = safeiv(c, delta)
    one = c.mpf(1)
    zero = c.mpf(0)
    z_signed = _z_signed_jet(c, z, order)
    h_signed = [one - delta] + [zero] * order
    d_signed = [one - z * z, -2 * z]
    if order >= 2:
        d_signed.append(-2 * one)
    return (_signed_logs(c, z_signed),
            _signed_logs(c, h_signed),
            _signed_logs(c, d_signed[:order + 1]))


def pressure_Q(c, U, M, z, delta, inverseL):
    """Bound the ordinary ``Z0..Z2`` rows of the pressure operator.

    ``U`` and ``M`` are raw logarithmic absolute caps; four rows are required
    because the order-two row consumes ``M_ZZZ``.  ``inverseL`` is an
    independently supplied logarithmic cap for ``1/L`` and requires rows
    through order two (a five-row packet is accepted and extra rows ignored).
    The result is the direct ordinary-product majorant of

    ``Q=(2*Z*U-(1-delta)*Z*M-(1-Z**2)*M_Z)/L``.
    """

    u = _validate_rows(c, U, "U", 4)
    m = _validate_rows(c, M, "M", 4)
    inv_l = _validate_rows(c, inverseL, "inverseL", 3,
                           value_required=True)
    z_logs, h_logs, d_logs = _geometric_log_jets(c, z, delta, 2)
    h_z_logs = _log_jet_product(c, h_logs, z_logs, 2)

    two_z_u = _log_jet_product(c, z_logs, u[:3], 2)
    two_z_u = [_log_scale(c, value, c.mpf(2)) for value in two_z_u]
    h_z_m = _log_jet_product(c, h_z_logs, m[:3], 2)
    d_m_z = _log_jet_product(c, d_logs, m[1:4], 2)
    numerator = _log_add_jets(c, two_z_u, h_z_m, d_m_z)
    return _log_jet_product(c, numerator, inv_l[:3], 2)


def _require_defect_rows(c, defects, keys):
    if not isinstance(defects, Mapping):
        raise TypeError("defects must be a mapping")
    out = {}
    for key in keys:
        if key not in defects:
            raise KeyError("missing defect row: " + key)
        out[key] = _validate_rows(c, defects[key], key, 2)
    return out


def six_differences(c, dU4, dQ3, U4, Q3):
    """Bound the six C1 integrand differences by exact polynomial expansion.

    The expansions are formed symbolically at the level of product terms
    before each term is absolute-majorized.  Thus the quadratic rows retain
    ``2*base*correction + correction^2`` and the mixed rows retain all three
    Leibniz terms.
    """

    d_u = _validate_rows(c, dU4, "dU", 4)
    d_q = _validate_rows(c, dQ3, "dQ", 3)
    u = _validate_rows(c, U4, "U", 4)
    q = _validate_rows(c, Q3, "Q", 3)

    d_u_c1, d_u_z = d_u[:2], d_u[1:3]
    d_q_c1, d_q_z = d_q[:2], d_q[1:3]
    u_c1, u_z = u[:2], u[1:3]
    q_c1, q_z = q[:2], q[1:3]

    # Exact polynomial differences, then logarithmic absolute majorization:
    # (Q+dQ)^2-Q^2 = 2 Q dQ+dQ^2.
    q_base_correction = _log_jet_product(c, q_c1, d_q_c1, 1)
    q_correction_squared = _log_jet_product(c, d_q_c1, d_q_c1, 1)
    q_squared = [
        _logsum(c, (_log_scale(c, q_base_correction[k], c.mpf(2)),
                    q_correction_squared[k]))
        for k in range(2)
    ]

    # (U+dU)(Q+dQ)_Z-U Q_Z.
    u_q_z = [
        _logsum(c, (
            _log_jet_product(c, u_c1, d_q_z, 1)[k],
            _log_jet_product(c, d_u_c1, q_z, 1)[k],
            _log_jet_product(c, d_u_c1, d_q_z, 1)[k],
        ))
        for k in range(2)
    ]

    # (U+dU)^2-U^2 = 2 U dU+dU^2.
    u_base_correction = _log_jet_product(c, u_c1, d_u_c1, 1)
    u_correction_squared = _log_jet_product(c, d_u_c1, d_u_c1, 1)
    u_squared = [
        _logsum(c, (_log_scale(c, u_base_correction[k], c.mpf(2)),
                    u_correction_squared[k]))
        for k in range(2)
    ]

    # (U+dU)(U+dU)_Z-U U_Z.
    u_u_z = [
        _logsum(c, (
            _log_jet_product(c, u_c1, d_u_z, 1)[k],
            _log_jet_product(c, d_u_c1, u_z, 1)[k],
            _log_jet_product(c, d_u_c1, d_u_z, 1)[k],
        ))
        for k in range(2)
    ]

    return {
        "q": list(d_q_c1),
        "q_Z": list(d_q_z),
        "q_squared": q_squared,
        "U_q_Z": u_q_z,
        "U_squared": u_squared,
        "U_U_Z": u_u_z,
    }


def boundary_difference(c, defects, z, delta, inverseL):
    """Return the two endpoint ``delta B/W`` logarithmic caps.

    The returned rows include the endpoint pair factor ``2``.  Before that
    factor, the exact difference is majorized as

    ``deltaG=deltaQ/L+delta(Q^2)/2-2*Z^2*delta(U^2)/L^2``.
    """

    d = _require_defect_rows(c, defects, ("q", "q_squared", "U_squared"))
    inv_l = _validate_rows(c, inverseL, "inverseL", 2,
                           value_required=True)
    z_logs, _, _ = _geometric_log_jets(c, z, delta, 1)
    z_squared = _log_jet_product(c, z_logs, z_logs, 1)
    inv_l_squared = _log_jet_product(c, inv_l[:2], inv_l[:2], 1)

    first = _log_jet_product(c, inv_l[:2], d["q"][:2], 1)
    second = [_log_scale(c, row, c.mpf("0.5"))
              for row in d["q_squared"][:2]]
    third = _log_jet_product(c, z_squared, inv_l_squared, 1)
    third = [_log_scale(c, row, c.mpf(2)) for row in third]
    third = _log_jet_product(c, third, d["U_squared"][:2], 1)
    delta_g = [
        _logsum(c, (first[k], second[k], third[k])) for k in range(2)
    ]
    pair_log = _abs_log(c, c.mpf(2))
    return [None if row is None else row + pair_log for row in delta_g]


def _signed_geometric_coefficients(c, z, delta):
    """Build the six signed C1 geometric coefficient receipts."""

    z = safeiv(c, z)
    delta = safeiv(c, delta)
    one = c.mpf(1)
    zero = c.mpf(0)
    z_jet = [z, one]
    z_squared = _signed_product(z_jet, z_jet)
    h = [one - delta, zero]
    d = [one - z * z, -2 * z]
    L = one - delta * z * z
    if L.a <= 0:
        raise ValueError("L=1-delta*z^2 must be positive")
    inv_l0 = one / L
    inv_l1 = (2 * delta * z) / (L * L)
    inv_l = [inv_l0, inv_l1]
    inv_l_squared = _signed_product(inv_l, inv_l)

    return {
        "minus_inv_L": _signed_scale(inv_l, -one),
        "hZ_over_2L": _signed_scale(
            _signed_product(_signed_product(h, z_jet), inv_l),
            one / 2),
        "minus_half": [-one / 2, zero],
        "d_over_L": _signed_product(d, inv_l),
        "two_h_Z2_over_L2": _signed_scale(
            _signed_product(_signed_product(h, z_squared), inv_l_squared),
            2 * one),
        "two_Zd_over_L2": _signed_scale(
            _signed_product(_signed_product(z_jet, d), inv_l_squared),
            2 * one),
    }


def _coefficient_log_rows(c, z, delta, inverseL):
    """Build absolute log caps for the six coefficients from inverseL caps."""

    inv_l = _validate_rows(c, inverseL, "inverseL", 2,
                           value_required=True)[:2]
    z_logs, h_logs, d_logs = _geometric_log_jets(c, z, delta, 1)
    z_squared = _log_jet_product(c, z_logs, z_logs, 1)
    inv_l_squared = _log_jet_product(c, inv_l, inv_l, 1)
    return {
        "minus_inv_L": list(inv_l),
        "hZ_over_2L": [
            _log_scale(c, row, c.mpf("0.5"))
            for row in _log_jet_product(
                c, _log_jet_product(c, h_logs, z_logs, 1), inv_l, 1)
        ],
        "minus_half": [_abs_log(c, c.mpf("0.5")), None],
        "d_over_L": _log_jet_product(c, d_logs, inv_l, 1),
        "two_h_Z2_over_L2": [
            _log_scale(c, row, c.mpf(2))
            for row in _log_jet_product(
                c, _log_jet_product(c, h_logs, z_squared, 1),
                inv_l_squared, 1)
        ],
        "two_Zd_over_L2": [
            _log_scale(c, row, c.mpf(2))
            for row in _log_jet_product(
                c, _log_jet_product(c, z_logs, d_logs, 1),
                inv_l_squared, 1)
        ],
    }


def D_difference(c, defects, z, delta, inverseL):
    """Compose the six finite-``N`` ``I1`` errors into the ``D`` error.

    The coefficient rows are returned as signed interval C1 rows for the
    caller's full signed ``D`` reconstruction.  ``term_log_bounds`` and
    ``rows`` remain logarithmic absolute caps, and use the independently
    supplied inverse-``L`` cap throughout.
    """

    d = _require_defect_rows(c, defects, I1_NAMES)
    coefficient_rows = _signed_geometric_coefficients(c, z, delta)
    coefficient_log_rows = _coefficient_log_rows(c, z, delta, inverseL)
    term_log_bounds = {
        name: _log_jet_product(c, coefficient_log_rows[_D_TERM_COEFFICIENT[name]],
                               d[name][:2], 1)
        for name in I1_NAMES
    }
    rows = [
        _logsum(c, [term_log_bounds[name][k] for name in I1_NAMES])
        for k in range(2)
    ]
    source_formula = (
        "D=-I[q]/L+h*Z*I[q_Z]/(2L)-I[q_squared]/2+"
        "d*I[U_q_Z]/L+2*h*Z^2*I[U_squared]/L^2+"
        "2*Z*d*I[U_U_Z]/L^2"
    )
    return {
        "rows": rows,
        "coefficient_rows": coefficient_rows,
        "coefficient_log_rows": coefficient_log_rows,
        "term_log_bounds": term_log_bounds,
        "source_formula": source_formula,
        "formula": "deltaD=sum(C_i*deltaI1_i), ordinary C1 Leibniz products",
        "composed_from_six_I1": True,
    }


__all__ = [
    "I1_NAMES",
    "D_COEFFICIENT_NAMES",
    "safeiv",
    "logsum",
    "log_product",
    "pressure_Q",
    "six_differences",
    "boundary_difference",
    "D_difference",
]
