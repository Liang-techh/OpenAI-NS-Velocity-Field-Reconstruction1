"""Interval-only partial reserved-minus-pulse coefficients.

The interval context is supplied by the caller.  This module does not create
a global numerical context and does not coerce scalar inputs into intervals.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class PartialReservedMinusPulseResult:
    """C0, baseline Z1, optional full Z1, and signed interval terms."""

    C0: Any
    Z1_baseline: Any
    Z1: Any | None
    sensitivity_coefficients: tuple[tuple[str, Any], ...]
    C0_terms: tuple[tuple[str, Any], ...]
    Z1_baseline_terms: tuple[tuple[str, Any], ...]
    Z1_terms: tuple[tuple[str, Any], ...] | None
    A: Any
    A_mu_Z: Any
    A_geometry_Z: Any


def _require_interval_context(iv: Any) -> None:
    if iv is None:
        raise TypeError("iv must be the caller-provided mpmath interval context")
    if not callable(getattr(iv, "mpf", None)):
        raise TypeError("iv must provide an mpf constructor")
    if not callable(getattr(iv, "exp", None)):
        raise TypeError("iv must provide interval exp")


def _is_interval(value: Any) -> bool:
    """Recognize an interval value without converting or evaluating it."""

    if value is None or isinstance(value, (bool, int, float, complex, str, bytes)):
        return False
    # mpmath's ivmpf exposes lower and upper endpoints as ``a`` and ``b``.
    return hasattr(value, "a") and hasattr(value, "b")


def _require_interval(name: str, value: Any) -> Any:
    if not _is_interval(value):
        raise TypeError(f"{name} must be an mpmath interval value")
    return value


def _require_optional_interval(name: str, value: Any) -> Any | None:
    if value is None:
        return None
    return _require_interval(name, value)


def _require_vector(name: str, values: Any, length: int) -> tuple[Any, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise TypeError(f"{name} must be a sequence of interval values")
    if len(values) != length:
        raise ValueError(f"{name} must contain exactly {length} interval values")
    return tuple(_require_interval(f"{name}[{index}]", value) for index, value in enumerate(values))


def _require_optional_vector(
    name: str,
    values: Any,
    length: int,
) -> tuple[Any | None, ...] | None:
    if values is None:
        return None
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise TypeError(f"{name} must be a sequence of interval values or None")
    if len(values) != length:
        raise ValueError(f"{name} must contain exactly {length} interval values")
    return tuple(
        None if value is None else _require_interval(f"{name}[{index}]", value)
        for index, value in enumerate(values)
    )


def _constant(iv: Any, text: str) -> Any:
    value = iv.mpf(text)
    if not _is_interval(value):
        raise TypeError("iv.mpf must return an mpmath interval value")
    return value


def _add_terms(iv: Any, terms: tuple[tuple[str, Any], ...]) -> Any:
    total = _constant(iv, "0")
    for _, term in terms:
        total = total + term
    return total


def evaluate_partial_reserved_minus_pulse(
    iv: Any,
    r_i: Sequence[Any],
    r_i_Z: Sequence[Any] | None,
    mu: Any,
    mu_Z: Any | None,
    theta_i: Sequence[Any],
    theta_i_Z: Sequence[Any],
    J: Any,
    J_Z: Any,
    P: Any,
    P_Z: Any,
    Z: Any,
) -> PartialReservedMinusPulseResult:
    """Evaluate partial reserved-minus-pulse C0 and ordinary Z1.

    ``r_i`` and ``r_i_Z`` are ``(r0, r1, r2)`` and
    ``(r0_Z, r1_Z, r2_Z)``.  ``theta_i`` and ``theta_i_Z`` are
    ``(theta0, theta1)`` and ``(theta0_Z, theta1_Z)``.  The four optional
    sensitivity inputs may be omitted independently: an omitted entry in
    ``r_i_Z`` is represented by ``None``, while ``mu_Z`` is an optional scalar
    interval.  When any of the four coefficients' inputs are absent, the
    result contains the baseline terms and coefficient intervals but ``Z1``
    is ``None``.
    """

    _require_interval_context(iv)

    r0, r1, r2 = _require_vector("r_i", r_i, 3)
    r_i_Z = _require_optional_vector("r_i_Z", r_i_Z, 3)
    theta0, theta1 = _require_vector("theta_i", theta_i, 2)
    theta0_Z, theta1_Z = _require_vector("theta_i_Z", theta_i_Z, 2)

    mu = _require_interval("mu", mu)
    mu_Z = _require_optional_interval("mu_Z", mu_Z)
    J = _require_interval("J", J)
    J_Z = _require_interval("J_Z", J_Z)
    P = _require_interval("P", P)
    P_Z = _require_interval("P_Z", P_Z)
    Z = _require_interval("Z", Z)

    one = _constant(iv, "1")
    two = _constant(iv, "2")
    five = _constant(iv, "5")
    seventeen = _constant(iv, "17")
    ten = _constant(iv, "10")

    denominator = one + Z ** 2
    exponent = -(seventeen / ten) - (five * mu / two)
    A = iv.exp(exponent) / denominator

    # Keep the two signed derivatives of A separate.  In particular, the
    # caller-visible mu contribution remains exactly -(5/2) * mu_Z * A.
    A_mu_Z = None if mu_Z is None else -(five / two) * mu_Z * A
    A_geometry_Z = -(two * Z / denominator) * A

    C0_terms = (
        ("minus_r0_A_theta0", -(r0 * A * theta0)),
        ("minus_r1_A_theta1", -(r1 * A * theta1)),
        ("plus_r2_J", r2 * J),
        ("plus_one_minus_r2_P", (one - r2) * P),
    )

    Z1_baseline_terms = (
        ("minus_r0_A_geometry_Z_theta0", -(r0 * A_geometry_Z * theta0)),
        ("minus_r1_A_geometry_Z_theta1", -(r1 * A_geometry_Z * theta1)),
        ("minus_r0_A_theta0_Z", -(r0 * A * theta0_Z)),
        ("minus_r1_A_theta1_Z", -(r1 * A * theta1_Z)),
        ("plus_r2_J_Z", r2 * J_Z),
        ("plus_one_minus_r2_P_Z", (one - r2) * P_Z),
    )

    sensitivity_coefficients = (
        ("r0_Z", -(A * theta0)),
        ("r1_Z", -(A * theta1)),
        ("r2_Z", J - P),
        ("mu_Z", -(five / two) * (C0_terms[0][1] + C0_terms[1][1])),
    )

    full_Z1_terms: tuple[tuple[str, Any], ...] | None = None
    full_Z1: Any | None = None
    if r_i_Z is not None and all(value is not None for value in r_i_Z) and mu_Z is not None:
        r0_Z, r1_Z, r2_Z = r_i_Z
        full_Z1_terms = (
            ("minus_r0_Z_A_theta0", -(r0_Z * A * theta0)),
            ("minus_r1_Z_A_theta1", -(r1_Z * A * theta1)),
            ("plus_r2_Z_times_J_minus_P", r2_Z * (J - P)),
            ("plus_mu_Z_times_minus_five_over_two_own", mu_Z * sensitivity_coefficients[3][1]),
            *Z1_baseline_terms,
        )
        full_Z1 = _add_terms(iv, full_Z1_terms)

    return PartialReservedMinusPulseResult(
        C0=_add_terms(iv, C0_terms),
        Z1_baseline=_add_terms(iv, Z1_baseline_terms),
        Z1=full_Z1,
        sensitivity_coefficients=sensitivity_coefficients,
        C0_terms=C0_terms,
        Z1_baseline_terms=Z1_baseline_terms,
        Z1_terms=full_Z1_terms,
        A=A,
        A_mu_Z=A_mu_Z,
        A_geometry_Z=A_geometry_Z,
    )


def partial_reserved_minus_pulse_c0_z1(
    iv: Any,
    r_i: Sequence[Any],
    r_i_Z: Sequence[Any] | None,
    mu: Any,
    mu_Z: Any | None,
    theta_i: Sequence[Any],
    theta_i_Z: Sequence[Any],
    J: Any,
    J_Z: Any,
    P: Any,
    P_Z: Any,
    Z: Any,
) -> tuple[Any, Any | None, tuple[tuple[str, Any], ...]]:
    """Return ``(C0, Z1, sensitivity_coefficients)``.

    ``Z1`` is ``None`` when one or more of ``r_i_Z`` and ``mu_Z`` is absent;
    unknown derivatives are never replaced with zero.
    """

    result = evaluate_partial_reserved_minus_pulse(
        iv,
        r_i,
        r_i_Z,
        mu,
        mu_Z,
        theta_i,
        theta_i_Z,
        J,
        J_Z,
        P,
        P_Z,
        Z,
    )
    return result.C0, result.Z1, result.sensitivity_coefficients


__all__ = [
    "PartialReservedMinusPulseResult",
    "evaluate_partial_reserved_minus_pulse",
    "partial_reserved_minus_pulse_c0_z1",
]
