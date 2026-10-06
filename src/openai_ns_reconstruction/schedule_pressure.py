"""Constructive scalar inputs for the pinned outgoing schedule pressure.

The official Lean development defines the outgoing smooth step by

    sigma(x) = edge(1,x) / (edge(1,x) + edge(1,1-x)),
    edge(1,x) = 0                    for x <= 0,
                exp(-1/x^2)         for x > 0,

and then defines ``SchedulePressure.shapeExponent`` by

    a(y) = 1 - sigma((y - endpoint) / flattenLength).

``OutgoingTail.flattenLength`` currently depends on ``stepBound := Classical.choose``
from the existence of a global derivative bound for sigma.  This module makes
one admissible choice executable: ``S = 32``.  The elementary bound is obtained
by symmetry and, on ``0 < x <= 1/2``,

    sigma'(x)
      = 2 a b/(a+b)^2 * (x^-3 + (1-x)^-3)
      <= 2 e^4 exp(-1/x^2) * (x^-3 + 8)
      <= 32.

Indeed ``b >= e^-4``, ``exp(-1/x^2) <= e^-4``, and
``x^-3 exp(-1/x^2) <= 8 e^-4`` on that half interval.  Symmetry gives the
other half.  Therefore ``S=32`` is a genuine witness for the theorem-side
existential derivative bound, yielding the constructive flattening length

    L = 10 * (S + 1) * log(2) + 1.

This is a constructive replacement for the noncomputable choice, not a claim
that Lean's opaque ``Classical.choose`` is definitionally equal to 32.  The
actual ``SchedulePressure.axisPressure`` still requires the complete
``finalAngular`` outgoing schedule and is intentionally not synthesized here.
"""

from __future__ import annotations

import math


EXPLICIT_STEP_BOUND = 32.0
"""A proved admissible global upper bound for the pinned outgoing ``sigma'``."""

EXPLICIT_FLATTEN_LENGTH = 10.0 * (EXPLICIT_STEP_BOUND + 1.0) * math.log(2.0) + 1.0
"""Constructive ``flattenLength`` obtained by instantiating the bound with 32."""


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def flat_edge_unit(x: float) -> float:
    """Pinned ``FlatCutoff.edge 1 x`` evaluated robustly in floating point."""

    x = _finite(x, "x")
    if x <= 0.0:
        return 0.0
    inv = 1.0 / x
    exponent = -(inv * inv)
    if not math.isfinite(exponent):
        return 0.0
    return math.exp(exponent)


def outgoing_sigma(x: float) -> float:
    """The exact scalar formula of ``OutgoingSchedule.sigma``.

    A log-ratio form avoids dividing two simultaneously tiny exponentials in
    the transition interval while preserving the same mathematical function.
    """

    x = _finite(x, "x")
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0

    log_a = -(1.0 / x) ** 2
    log_b = -(1.0 / (1.0 - x)) ** 2
    delta = log_b - log_a
    if delta >= 0.0:
        e = math.exp(-delta)
        return e / (1.0 + e)
    e = math.exp(delta)
    return 1.0 / (1.0 + e)


def outgoing_sigma_derivative(x: float) -> float:
    """Analytic first derivative of the pinned outgoing ``sigma``.

    It is identically zero off the transition interval.  Inside ``(0,1)`` we
    evaluate ``2*sigma*(1-sigma)*(x^-3+(1-x)^-3)`` in logarithmic form to
    avoid the endpoint ``0 * infinity`` cancellation.
    """

    x = _finite(x, "x")
    if x <= 0.0 or x >= 1.0:
        return 0.0
    s = outgoing_sigma(x)
    if s <= 0.0 or s >= 1.0:
        return 0.0

    log_prefactor = math.log(2.0) + math.log(s) + math.log1p(-s)
    left = math.exp(log_prefactor - 3.0 * math.log(x))
    right = math.exp(log_prefactor - 3.0 * math.log1p(-x))
    return left + right


def schedule_shape_exponent(
    y: float,
    endpoint: float,
    *,
    flatten_length: float = EXPLICIT_FLATTEN_LENGTH,
) -> float:
    """Executable ``SchedulePressure.shapeExponent`` for a fixed endpoint.

    ``endpoint`` remains an upstream outgoing-schedule parameter.  The default
    flattening length uses the explicit derivative-bound witness above rather
    than an arbitrary hand-tuned scale.
    """

    y = _finite(y, "y")
    endpoint = _finite(endpoint, "endpoint")
    flatten_length = _finite(flatten_length, "flatten_length")
    if flatten_length <= 0.0:
        raise ValueError("flatten_length must be positive")
    return 1.0 - outgoing_sigma((y - endpoint) / flatten_length)


def flatten_end(endpoint: float, *, flatten_length: float = EXPLICIT_FLATTEN_LENGTH) -> float:
    """Constructive counterpart of ``TailData.flattenEnd``."""

    endpoint = _finite(endpoint, "endpoint")
    flatten_length = _finite(flatten_length, "flatten_length")
    if flatten_length <= 0.0:
        raise ValueError("flatten_length must be positive")
    return endpoint + flatten_length
