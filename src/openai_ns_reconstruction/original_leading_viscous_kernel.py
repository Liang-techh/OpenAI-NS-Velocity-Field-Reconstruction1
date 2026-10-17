"""Ordinary axial operators for source-bound leading viscous coefficients.

Inputs are ordinary derivatives of signed source functions, never positive
error caps.  The operators do not materialize an implicit width or amplitude.
"""
from __future__ import annotations

from typing import Any, Sequence

from mpmath.libmp import mpf_abs, mpf_cmp


def endpoint_coefficient(z: Any, delta: Any, source: Sequence[Any], *, sign: int = 1) -> list[Any]:
    """Return ordinary Z0/Z1 of sign * 2 Z B / (1-delta Z^2)."""
    if len(source) < 2 or sign not in (-1, 1):
        raise ValueError("Signed B and ordinary B_Z are required")
    b, bz = source[:2]
    denominator = 1 - delta * z * z
    return [sign * 2 * z * b / denominator,
            sign * (2 * (b + z * bz) / denominator
                    + 4 * delta * z * z * b / (denominator * denominator))]


def mean_source_value(z: Any, delta: Any, source: Sequence[Any]) -> Any:
    """Return C(B) from signed B/B_Z; no second derivative is invented."""
    if len(source) < 2:
        raise ValueError("Signed B and ordinary B_Z are required")
    return ((1 - delta) * z * source[0] + (1 - z * z) * source[1]) / (1 - delta * z * z)


def mean_source_coefficient(z: Any, delta: Any, source: Sequence[Any]) -> list[Any]:
    """Return ordinary Z0/Z1 of C(B), requiring the actual B_ZZ input."""
    if len(source) < 3:
        raise ValueError("C(B) ordinary Z1 requires signed B_ZZ; missing rows cannot be zero-padded")
    b, bz, bzz = source[:3]
    denominator = 1 - delta * z * z
    numerator = (1 - delta) * z * b + (1 - z * z) * bz
    numerator_z = (1 - delta) * b - (1 + delta) * z * bz + (1 - z * z) * bzz
    return [numerator / denominator,
            numerator_z / denominator + 2 * delta * z * numerator / (denominator * denominator)]


def decode_pair(ctx: Any, iv: Any, row: dict[str, Any]) -> tuple[Any, Any]:
    """Hydrate the saved approximation and exact directed enclosure separately."""
    if row.get("positive_error_bound_atom_only") or row.get("source_log_atom_only"):
        raise ValueError("A cap or logarithmic scale atom is not a signed coefficient")
    if not isinstance(row.get("approximation"), str):
        raise ValueError("A saved numerical approximation is required")
    endpoints = row.get("source_enclosure_exact_mpf_tuples")
    if not isinstance(endpoints, (list, tuple)) or len(endpoints) != 2:
        raise ValueError("Exact directed source endpoints are required; rounded bounds are insufficient")
    if any(not isinstance(v, (list, tuple)) or len(v) != 4 for v in endpoints):
        raise ValueError("Malformed exact mpf endpoints")
    point = ctx.mpf(row["approximation"])
    truth = iv.mpf([ctx.make_mpf(tuple(v)) for v in endpoints])
    return point, truth


def encode_pair(ctx: Any, iv: Any, point: Any, truth: Any) -> dict[str, Any]:
    """Serialize a new arithmetic approximation without promoting it to truth."""
    displayed = ctx.nstr(point, ctx.dps)
    difference = truth - iv.mpf(displayed)
    endpoints = [mpf_abs(value) for value in difference._mpi_]
    error = endpoints[0] if mpf_cmp(endpoints[0], endpoints[1]) >= 0 else endpoints[1]
    return {
        "approximation": displayed,
        "source_enclosure_exact_mpf_tuples": [list(v) for v in truth._mpi_],
        "serialized_approximation_absolute_error_upper_exact_mpf_tuple": list(error),
        "approximation_is_exact_source_value": False,
        "source_enclosure_is_not_a_selected_value": True,
        "ordinary_derivatives_not_Taylor_coefficients": True,
    }
