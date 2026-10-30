r"""Independent ordinary-``Z`` transport evaluation for a ``Q^2`` integral.

The public callable evaluates the convention-dependent expression

.. math::

   -\tfrac14\int_l^u M_Z Q_0^2\,d\xi
   +\tfrac12\int_l^u M Q_0(Q_1-d_ZQ_0)\,d\xi
   -\tfrac14(MQ_0^2)(u)u_Z
   +\tfrac14(MQ_0^2)(l)l_Z.

``M``, ``M_Z``, ``Q0``, ``Q1``, and ``dZ_Q0`` are caller-owned scalar
functions of ``xi``.  The ``*_Z`` names denote fixed-``xi`` ordinary-
``Z`` derivatives; ``l_Z`` and ``u_Z`` are the corresponding endpoint
velocities.  ``quadrature`` is also caller-owned and is called as
``quadrature(integrand, a, b)`` to approximate the oriented ``d xi`` integral.
Here ``M`` is already the caller's physical-measure factor (if one is part of
the convention), so ``quadrature`` must not add another physical weight.  An
extra weight would be counted twice and would not match the endpoint terms.
Quadrature accuracy and numerical reliability remain the caller's
responsibility; this module does not estimate an error or combine this result
with any prior uncertainty bound.

This is a calculation under the stated convention and callback contract.  It
does not claim a real ``F`` identity, source authentication, certified error,
or closure.  No source graph, saved packet, or producer module is consulted.
The endpoint order is kept exactly as supplied: equal endpoints and reverse
oriented intervals are valid, and no absolute value or endpoint sorting is
performed.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable


ScalarCallback = Callable[[float], float]
Quadrature = Callable[[ScalarCallback, float, float], float]


__all__ = [
    "QSquareLeibnizInputError",
    "QSquareLeibnizTerms",
    "qsquare_leibniz_transport",
    "compute_qsquare_leibniz_transport",
    "compute_qsquare_leibniz",
]


class QSquareLeibnizInputError(ValueError):
    """Raised for a non-callable input or a non-finite scalar value."""


@dataclass(frozen=True)
class QSquareLeibnizTerms:
    """The four signed terms and their arithmetic sum.

    ``term_mz_q0_squared`` and ``term_m_q0_q1_minus_dz_q0`` already include
    their coefficients ``-1/4`` and ``+1/2``.  The endpoint fields already
    include ``-1/4`` at ``u`` and ``+1/4`` at ``l``.  Nothing in this result is
    an interval, uncertainty, or certification object.
    """

    term_mz_q0_squared: float
    term_m_q0_q1_minus_dz_q0: float
    term_upper_endpoint: float
    term_lower_endpoint: float
    total: float

    @property
    def integral_mz_q0_squared(self) -> float:
        """Compatibility name for the signed first integral term."""

        return self.term_mz_q0_squared

    @property
    def integral_m_q0_q1_minus_dz_q0(self) -> float:
        """Compatibility name for the signed second integral term."""

        return self.term_m_q0_q1_minus_dz_q0

    @property
    def upper_endpoint(self) -> float:
        """Compatibility name for the signed ``u`` endpoint term."""

        return self.term_upper_endpoint

    @property
    def lower_endpoint(self) -> float:
        """Compatibility name for the signed ``l`` endpoint term."""

        return self.term_lower_endpoint

    def as_dict(self) -> dict[str, float]:
        """Return only the five computed numbers, without claims of error."""

        return {
            "term_mz_q0_squared": self.term_mz_q0_squared,
            "term_m_q0_q1_minus_dz_q0": self.term_m_q0_q1_minus_dz_q0,
            "term_upper_endpoint": self.term_upper_endpoint,
            "term_lower_endpoint": self.term_lower_endpoint,
            "total": self.total,
        }


def _finite_scalar(value: Any, label: str) -> float:
    """Convert a scalar to ``float`` and reject booleans and non-finite data."""

    if isinstance(value, bool):
        raise QSquareLeibnizInputError(f"{label} must be a finite scalar")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise QSquareLeibnizInputError(f"{label} must be a finite scalar") from exc
    if not math.isfinite(result):
        raise QSquareLeibnizInputError(f"{label} must be finite")
    return result


def _require_callable(value: Any, label: str) -> None:
    if not callable(value):
        raise QSquareLeibnizInputError(f"{label} must be callable")


def _at(callback: ScalarCallback, xi: Any, callback_name: str) -> float:
    """Evaluate one caller callback at a finite ``xi`` and check its scalar."""

    x = _finite_scalar(xi, "xi")
    return _finite_scalar(callback(x), f"{callback_name}(xi)")


def qsquare_leibniz_transport(
    M: ScalarCallback,
    M_Z: ScalarCallback,
    Q0: ScalarCallback,
    Q1: ScalarCallback,
    dZ_Q0: ScalarCallback,
    l: float,
    u: float,
    l_Z: float,
    u_Z: float,
    quadrature: Quadrature,
) -> QSquareLeibnizTerms:
    """Evaluate the signed ordinary-``Z`` transport terms for ``M Q0^2``.

    Parameters
    ----------
    M, M_Z, Q0, Q1, dZ_Q0:
        Scalar callbacks with the single argument ``xi``.  ``M_Z`` and
        ``dZ_Q0`` are supplied ordinary-``Z`` derivatives at fixed ``xi``;
        this function does not differentiate callbacks itself.
    l, u:
        The two ``xi`` endpoints.  They are passed to ``quadrature`` in this
        order, without sorting, so ``u < l`` gives an oriented integral.
    l_Z, u_Z:
        Scalar ``Z`` velocities of the endpoints ``l`` and ``u``.
    quadrature:
        Caller-owned callable with signature ``quadrature(f, a, b)`` that
        approximates the oriented ``d xi`` integral.  The physical-measure
        factor is already in ``M``; do not apply another weight in this
        callable.  Its numerical reliability remains the caller's
        responsibility.

    Returns
    -------
    QSquareLeibnizTerms
        Four signed contributions and their finite arithmetic sum.  No error
        estimate or previously certified uncertainty is inferred or changed.
    """

    callbacks = {
        "M": M,
        "M_Z": M_Z,
        "Q0": Q0,
        "Q1": Q1,
        "dZ_Q0": dZ_Q0,
    }
    for name, callback in callbacks.items():
        _require_callable(callback, name)
    _require_callable(quadrature, "quadrature")

    # Do not impose l <= u: orientation is part of the caller's convention.
    l_value = _finite_scalar(l, "l")
    u_value = _finite_scalar(u, "u")
    l_z_value = _finite_scalar(l_Z, "l_Z")
    u_z_value = _finite_scalar(u_Z, "u_Z")

    def first_integrand(xi: float) -> float:
        mz = _at(M_Z, xi, "M_Z")
        q0 = _at(Q0, xi, "Q0")
        return _finite_scalar(mz * q0 * q0, "M_Z Q0^2 integrand")

    def second_integrand(xi: float) -> float:
        m = _at(M, xi, "M")
        q0 = _at(Q0, xi, "Q0")
        q1 = _at(Q1, xi, "Q1")
        dz_q0 = _at(dZ_Q0, xi, "dZ_Q0")
        return _finite_scalar(
            m * q0 * (q1 - dz_q0),
            "M Q0 (Q1 - dZ_Q0) integrand",
        )

    first_raw = _finite_scalar(
        quadrature(first_integrand, l_value, u_value),
        "quadrature(M_Z Q0^2)",
    )
    second_raw = _finite_scalar(
        quadrature(second_integrand, l_value, u_value),
        "quadrature(M Q0 (Q1 - dZ_Q0))",
    )
    first = _finite_scalar(-0.25 * first_raw, "-1/4 quadrature(M_Z Q0^2)")
    second = _finite_scalar(
        0.5 * second_raw,
        "+1/2 quadrature(M Q0 (Q1 - dZ_Q0))",
    )

    def endpoint_product(xi: float, label: str) -> float:
        m = _at(M, xi, "M")
        q0 = _at(Q0, xi, "Q0")
        return _finite_scalar(m * q0 * q0, f"{label} M Q0^2")

    upper = _finite_scalar(
        -0.25 * endpoint_product(u_value, "u") * u_z_value,
        "upper endpoint term",
    )
    lower = _finite_scalar(
        0.25 * endpoint_product(l_value, "l") * l_z_value,
        "lower endpoint term",
    )
    total = _finite_scalar(math.fsum((first, second, upper, lower)), "total")
    return QSquareLeibnizTerms(first, second, upper, lower, total)


# A descriptive alias keeps the operation discoverable without introducing a
# second implementation or a second source of arithmetic.
compute_qsquare_leibniz_transport = qsquare_leibniz_transport
compute_qsquare_leibniz = qsquare_leibniz_transport
