"""Finite-interval cumulative moments from Lei--Ren Part I, Section 2.5.

Source pinned to Z. Lei and X. Ren, arXiv:2609.35406v1,
Section 2.5, equations (2.21)--(2.22):

* ``M_theta = 2 * integral(rho * F, d rho)``;
* ``M_z = integral(Uz, d rho)``;
* ``M_theta_z = 2 * integral(rho * F * Uz, d rho)``;
* ``M_z_theta = integral(Uz**2 - rho * F**2, d rho)``; and
* ``M_p = integral(F**2, d rho)``.

The implementation below evaluates those integrals only on the finite
interval ``0 <= rho <= Rmax`` with an n-point Gauss--Legendre rule.  ``Rmax``
and ``Z`` are broadcast together, and the supplied profile callables are
evaluated at array-valued ``(rho, Z)`` nodes.  This is a numerical interface
for finite supplied profiles.  It does not evaluate an infinite tail, certify
the heat exterior, establish a PDE solution, or establish scale recursion.
"""

from __future__ import annotations

import operator
from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTION = "Section 2.5, equations (2.21)--(2.22)"
MOMENT_NAMES = ("angular", "axial", "mixed", "quadratic", "pressure")


def _finite_array(value: Any, name: str) -> np.ndarray:
    """Convert a numeric input to float and require finite entries."""

    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real numeric value or array") from exc
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _validate_order(n: Any) -> int:
    """Validate a Gauss--Legendre order without accepting booleans."""

    if isinstance(n, (bool, np.bool_)):
        raise TypeError("n must be a positive integer")
    try:
        order = operator.index(n)
    except TypeError as exc:
        raise TypeError("n must be a positive integer") from exc
    if order <= 0:
        raise ValueError("n must be a positive integer")
    return int(order)


def _call_profile(function: Callable[..., Any], rho: np.ndarray,
                  z_nodes: np.ndarray, name: str) -> np.ndarray:
    """Evaluate a profile on nodes, with a scalar-callable fallback.

    Array-aware callables are preferred.  A fallback through ``np.vectorize``
    keeps the interface usable for scalar-only profiles such as the finite
    ``PaperCoreReference`` approximation, while retaining the same output
    shape and finite-value checks.
    """

    if not callable(function):
        raise TypeError(f"{name} must be callable as {name}(R, Z)")
    try:
        values = function(rho, z_nodes)
        array = np.asarray(values, dtype=float)
        array = np.broadcast_to(array, rho.shape)
    except (TypeError, ValueError, IndexError):
        try:
            scalar_function = np.vectorize(function, otypes=[float])
            array = np.asarray(scalar_function(rho, z_nodes), dtype=float)
        except (TypeError, ValueError, IndexError) as exc:
            raise ValueError(
                f"{name}(R, Z) must return values broadcastable to {rho.shape}"
            ) from exc
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name}(R, Z) returned non-finite values")
    return array


def _call_pressure(value: Any, rmax: np.ndarray, z: np.ndarray,
                   name: str) -> np.ndarray:
    """Evaluate an optional scalar/array or pressure callable.

    ``P0`` is commonly a function of ``Z`` and ``P`` may be a function of
    ``(Rmax, Z)``.  Both forms are accepted here.  Numeric values must be
    broadcastable to the common ``(Rmax, Z)`` shape.
    """

    if callable(value):
        try:
            candidate = value(rmax, z)
        except TypeError as two_argument_error:
            try:
                candidate = value(z)
            except (TypeError, ValueError, IndexError) as one_argument_error:
                raise ValueError(
                    f"{name} callable must accept (Rmax, Z) or (Z)"
                ) from one_argument_error
            del two_argument_error
    else:
        candidate = value
    try:
        array = np.asarray(candidate, dtype=float)
        array = np.broadcast_to(array, rmax.shape)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be broadcastable to shape {rmax.shape}") from exc
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def pressure_consistency_defect(P: Any, P0: Any, M_p: Any) -> np.ndarray:
    """Return the finite-interval pressure defect ``P - P0 - M_p``.

    This is the local consistency identity associated with the Part I
    convention ``P(R,Z) = P0(Z) + M_p(R,Z)``.  It only compares supplied
    values and a finite cumulative pressure moment; it does not infer a
    pressure tail or certify radial centrifugal balance.

    Source: arXiv:2609.35406v1, Section 2.5, equation (2.22).
    """

    p, p0, mp = np.broadcast_arrays(
        _finite_array(P, "P"), _finite_array(P0, "P0"), _finite_array(M_p, "M_p")
    )
    return np.asarray(p - p0 - mp, dtype=float)


def cumulative_five_moments(
    F: Callable[..., Any],
    Uz: Callable[..., Any],
    Rmax: Any,
    Z: Any,
    *,
    n: Any = 64,
    P0: Any | None = None,
    P: Any | None = None,
) -> dict[str, np.ndarray]:
    """Compute the five finite cumulative moments for supplied profiles.

    Parameters
    ----------
    F, Uz:
        Callables of ``(R, Z)``.  They may accept NumPy arrays directly or
        scalar values through the vectorized fallback.
    Rmax, Z:
        Finite, nonnegative terminal radius and finite axial coordinate(s).
        They are broadcast together; every returned array has that shape.
    n:
        Positive integer Gauss--Legendre order on ``[0, Rmax]``.
    P0, P:
        Optional pressure datum and terminal pressure.  Supplying both adds
        a ``pressure_defect`` entry equal to ``P - P0 - M_p``.  Each may be a
        numeric scalar/array or a callable accepting ``Z`` (or ``Rmax, Z``).

    Returns
    -------
    dict
        Arrays under ``angular``, ``axial``, ``mixed``, ``quadratic`` and
        ``pressure``.  If both pressure values are supplied, the dictionary
        also contains ``pressure_defect``.

    Notes
    -----
    The source formulas are equations (2.21)--(2.22) of
    ``arXiv:2609.35406v1``.  This helper is intentionally finite interval
    only: no infinite-tail estimate, heat-exterior matching, PDE validation,
    or scale-recursion certificate is implied.
    """

    order = _validate_order(n)
    rmax, z = np.broadcast_arrays(
        _finite_array(Rmax, "Rmax"), _finite_array(Z, "Z")
    )
    if np.any(rmax < 0):
        raise ValueError("Rmax must be nonnegative")

    # A zero-length interval has exactly zero moments.  Avoid evaluating a
    # potentially singular profile at R=0 in this degenerate case.
    if rmax.size == 0 or np.all(rmax == 0):
        zeros = {name: np.zeros(rmax.shape, dtype=float) for name in MOMENT_NAMES}
        if P0 is not None or P is not None:
            if P0 is None or P is None:
                raise ValueError("P0 and P must be supplied together")
            p0 = _call_pressure(P0, rmax, z, "P0")
            p = _call_pressure(P, rmax, z, "P")
            zeros["pressure_defect"] = pressure_consistency_defect(
                p, p0, zeros["pressure"]
            )
        return zeros

    nodes, weights = leggauss(order)
    rho = rmax[..., None] * (nodes + 1.0) / 2.0
    quadrature_weights = rmax[..., None] * weights / 2.0
    z_nodes = np.broadcast_to(z[..., None], rho.shape)
    f = _call_profile(F, rho, z_nodes, "F")
    uz = _call_profile(Uz, rho, z_nodes, "Uz")

    # Equations (2.22), in the order named by MOMENT_NAMES.
    values = {
        "angular": np.sum(quadrature_weights * (2.0 * rho * f), axis=-1),
        "axial": np.sum(quadrature_weights * uz, axis=-1),
        "mixed": np.sum(quadrature_weights * (2.0 * rho * f * uz), axis=-1),
        "quadratic": np.sum(
            quadrature_weights * (uz * uz - rho * f * f), axis=-1
        ),
        "pressure": np.sum(quadrature_weights * (f * f), axis=-1),
    }
    values = {key: np.asarray(value, dtype=float) for key, value in values.items()}

    if P0 is not None or P is not None:
        if P0 is None or P is None:
            raise ValueError("P0 and P must be supplied together")
        p0 = _call_pressure(P0, rmax, z, "P0")
        p = _call_pressure(P, rmax, z, "P")
        values["pressure_defect"] = pressure_consistency_defect(
            p, p0, values["pressure"]
        )
    return values


def five_moments(
    F: Callable[..., Any], Uz: Callable[..., Any], Rmax: Any, Z: Any,
    *, n: Any = 64, P0: Any | None = None, P: Any | None = None,
) -> dict[str, np.ndarray]:
    """Alias for :func:`cumulative_five_moments`.

    The alias keeps the short name convenient for finite profile diagnostics;
    the source and finite-interval limitations are those of Section 2.5 in
    ``arXiv:2609.35406v1``.
    """

    return cumulative_five_moments(F, Uz, Rmax, Z, n=n, P0=P0, P=P)


__all__ = [
    "MOMENT_NAMES",
    "SOURCE",
    "SOURCE_SECTION",
    "SOURCE_VERSION",
    "cumulative_five_moments",
    "five_moments",
    "pressure_consistency_defect",
]
