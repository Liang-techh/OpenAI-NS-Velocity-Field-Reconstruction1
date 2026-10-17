"""Pure ordinary radial kernels for frozen macro source modes.

For a signed ordinary source row ``B_j100`` at the native-two to ``R100``
port, the radial mode is

``B_j(y) = B_j100 * exp((q-j)*(y-L)) * A(y,Z)``.

The caller supplies a uniform truth row when the relative angular factor
``A`` is not one.  After the positive transport factor ``exp(-(L-y))`` is
included, the scalar rate is ``a=q-j+1`` and

``K_a(L) = -expm1(-a*L)/a`` for nonzero ``a`` and ``K_0(L)=L``.

This module only evaluates those kernels and sums ordinary ``Z0/Z1/Z2``
rows.  It does not load source reports, call owners or producers, materialize
a width, select an enclosure midpoint, or attach physical ``R100``, epsilon,
or ``hbB`` factors.  Both source identity and Z-independent length proofs are
explicit inputs to the public evaluator.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import mpmath as _mp
from mpmath.libmp import fzero, mpf_cmp


def _sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _scalar(value: Any) -> bool:
    return not _sequence(value) and not isinstance(value, Mapping)


def _interval(value: Any) -> bool:
    return hasattr(value, "_mpi_")


def _row(value: Any, name: str = "ordinary row") -> list[Any]:
    if not _sequence(value) or len(value) != 3:
        raise ValueError(f"{name} must contain ordinary Z0/Z1/Z2 rows")
    if any(item is None or not _scalar(item) for item in value):
        raise ValueError(f"{name} contains a missing or nested derivative")
    return list(value)


def _integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not _scalar(value):
        raise ValueError(f"{name} must be an integer radial rate")
    try:
        result = int(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{name} must be an integer radial rate") from None
    if value != result:
        raise ValueError(f"{name} must be an integer radial rate")
    return result


def _mode_parts(mode: Any, index: int) -> tuple[list[Any], int]:
    if isinstance(mode, Mapping):
        rows = mode.get("ordinary_rows", mode.get("rows"))
        rate = mode.get("rate", mode.get("a"))
    elif _sequence(mode) and len(mode) == 2 and _sequence(mode[0]):
        rows, rate = mode
    else:
        raise ValueError(
            "modes must contain {ordinary_rows: [B, B_Z, B_ZZ], rate: a}"
        )
    if rows is None or rate is None:
        raise ValueError(f"mode {index} requires ordinary_rows and rate")
    return _row(rows, f"mode {index}"), _integer(rate, f"mode {index} rate")


def _modes(value: Any) -> tuple[list[list[Any]], list[int]]:
    if isinstance(value, Mapping):
        value = value.get("modes", value.get("source_modes"))
    if not _sequence(value) or not value:
        raise ValueError("modes must contain at least one signed source mode")
    parts = [_mode_parts(mode, index) for index, mode in enumerate(value)]
    return [part[0] for part in parts], [part[1] for part in parts]


def _bounds(context: Any, value: Any) -> tuple[Any, Any]:
    if _interval(value):
        lower, upper = value._mpi_
        return _mp.make_mpf(lower), _mp.make_mpf(upper)
    point = context.mpf(value)
    return point, point


def _length(context: Any, value: Any, name: str = "loglength") -> Any:
    if not _scalar(value):
        raise ValueError(f"{name} must be a scalar radial length")
    lower, upper = _bounds(context, value)
    if upper < lower or lower < 0:
        raise ValueError(f"{name} must be a nonnegative directed length")
    return value


def _expm1(context: Any, value: Any) -> Any:
    function = getattr(context, "expm1", None)
    if function is not None:
        try:
            return function(value)
        except (AttributeError, NotImplementedError, TypeError):
            pass
    return context.exp(value) - 1


def radial_kernel(context: Any, rate: Any, length: Any) -> Any:
    """Evaluate the positive ``K_rate(length)`` kernel in a point or iv ctx."""

    rate = _integer(rate, "rate")
    length = _length(context, length, "length")
    if _bounds(context, length) == (0, 0):
        return length
    if rate == 0:
        return length
    return -_expm1(context, -rate * length) / rate


def directed_radial_kernel(context: Any, rate: Any, directedL: Any) -> Any:
    """Evaluate ``radial_kernel`` with an explicitly directed length."""

    return radial_kernel(context, rate, directedL)


K_a = radial_kernel


def weighted_source_rows(modes: Any, kernels: Sequence[Any]) -> list[Any]:
    """Sum signed ordinary rows against an already evaluated kernel list."""

    rows, _ = _modes(modes)
    if not _sequence(kernels) or len(kernels) != len(rows):
        raise ValueError("kernels must have one entry per source mode")
    if any(not _scalar(kernel) for kernel in kernels):
        raise ValueError("kernels must be scalar values")
    return [
        sum((row[order] * kernel for row, kernel in zip(rows, kernels)), rows[0][order] * 0)
        for order in range(3)
    ]


def evaluate_macro_modes(
    context: Any,
    modes: Any,
    loglength: Any,
    *,
    source_proved: Any,
    L_Z_independent: Any = True,
    z: Any | None = None,
    delta: Any | None = None,
) -> dict[str, Any]:
    """Evaluate frozen macro modes and return response rows plus kernels.

    ``modes`` is a list of mappings with ``ordinary_rows`` and integer
    ``rate`` (``a=q-j+1``).  ``context`` may be ``mp`` for a point result or
    ``MPIntervalContext`` for a directed truth result.  A missing ordinary
    ``Z2`` row, source proof, or length independence proof fails closed.  If
    both ``z`` and ``delta`` are supplied, ``C_delta(response)`` through
    ordinary ``Z1`` is added to the packet.
    """

    if source_proved is not True:
        raise ValueError("an explicit same-source proof is required")
    if L_Z_independent is not True:
        raise ValueError("an explicit Z-independent length proof is required")
    rows, rates = _modes(modes)
    loglength = _length(context, loglength)
    kernels = [radial_kernel(context, rate, loglength) for rate in rates]
    response = weighted_source_rows(modes, kernels)
    result: dict[str, Any] = {
        "response": response,
        "ordinary_rows": response,
        "W": response,
        "kernels": kernels,
        "rates": rates,
        "loglength": loglength,
        "certificate": {
            "passed": True,
            "source_proved": True,
            "L_Z_independent": True,
            "ordinary_rows_through_Z2": True,
            "signed_modes_retained": True,
            "quadrature_used": False,
            "width_materialized": False,
            "physical_R100_epsilon_hbB_factors_attached": False,
        },
    }
    if z is not None or delta is not None:
        if z is None or delta is None:
            raise ValueError("z and delta must be supplied together for C_delta")
        result["C_delta"] = apply_c_delta(z, delta, response)
        result["c_delta"] = result["C_delta"]
        result["z"] = z
        result["delta"] = delta
    return result


def weighted_response(
    context: Any,
    modes: Any,
    loglength: Any,
    *,
    source_proved: Any,
    L_Z_independent: Any = True,
) -> list[Any]:
    """Short row-only wrapper around :func:`evaluate_macro_modes`."""

    return evaluate_macro_modes(
        context,
        modes,
        loglength,
        source_proved=source_proved,
        L_Z_independent=L_Z_independent,
    )["response"]


def _denominator_zero(value: Any) -> bool:
    if _interval(value):
        lower, upper = value._mpi_
        return mpf_cmp(lower, fzero) <= 0 and mpf_cmp(upper, fzero) >= 0
    return value == 0


def apply_c_delta(z: Any, delta: Any, source_rows: Sequence[Any]) -> list[Any]:
    """Apply ``C_delta`` through ordinary ``Z1`` to a ``Z0/Z1/Z2`` row."""

    b, bz, bzz = _row(source_rows, "source_rows")
    if not _scalar(z) or not _scalar(delta):
        raise ValueError("z and delta must be Z-independent scalars")
    denominator = 1 - delta * z * z
    if _denominator_zero(denominator):
        raise ValueError("C_delta denominator crosses zero")
    numerator = (1 - delta) * z * b + (1 - z * z) * bz
    numerator_z = (1 - delta) * b - (1 + delta) * z * bz + (1 - z * z) * bzz
    return [
        numerator / denominator,
        numerator_z / denominator
        + 2 * delta * z * numerator / (denominator * denominator),
    ]


c_delta = apply_c_delta
C_delta = apply_c_delta


def macro_log_length(context: Any, R0: Any, outer_radius: Any = 100) -> Any:
    """Compute the actual positive ``log(outer_radius/R0)`` length."""

    R0 = context.mpf(R0)
    outer_radius = context.mpf(outer_radius)
    if R0 <= 0 or outer_radius <= 0 or outer_radius < R0:
        raise ValueError("require 0 < R0 <= outer_radius")
    return context.log(outer_radius / R0)


native2_to_R100_length = macro_log_length


__all__ = [
    "C_delta",
    "K_a",
    "apply_c_delta",
    "c_delta",
    "directed_radial_kernel",
    "evaluate_macro_modes",
    "macro_log_length",
    "native2_to_R100_length",
    "radial_kernel",
    "weighted_response",
    "weighted_source_rows",
]
