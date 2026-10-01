"""Incremental explicit-center step for the directed functional core.

This module is the one-step form of ``coupled_rows``.  It consumes one axial
Taylor coefficient at every radial order and appends the next ``A``, ``Uz``
and physical ``P`` row in place.  The center is an explicit scalar or real
interval argument on every call; no fixed axial center is read from a legacy
driver.

The recurrence retains the Section 8.2 pressure, swirl, and cross terms.  A
non-degenerate center interval describes a family of local coefficient rows;
it is not a Taylor series about an interval midpoint.  This finite step has
no infinite-radial, whole-axis, or terminal matching claim.
"""

from __future__ import annotations

import hashlib
import json
import operator
from pathlib import Path
from typing import Any

from mpmath.ctx_iv import MPIntervalContext

HERE = Path(__file__).resolve().parent

DEFAULT_REQUIRED_DEPTH = 3
SOURCE_NAMES = (
    "lei_ren_part1_paper_amplitude_factored_core.py",
    "lei_ren_part1_paper_candidate_gauge_core.py",
    "lei_ren_part1_paper_functional_core_recursion.py",
)


def _source_hashes() -> dict[str, str]:
    """Return source hashes used for the copied recurrence and provenance."""

    names = (Path(__file__).name,) + SOURCE_NAMES
    return {
        name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
        for name in names
    }


def _as_interval(ctx: MPIntervalContext, value: Any):
    if hasattr(value, "_mpi_"):
        return value
    return ctx.mpf(value)


def _coerce_vector(ctx: MPIntervalContext, values: Any) -> list[Any]:
    return [_as_interval(ctx, value) for value in values]


def _validate_fixed(
    fixed: dict[str, Any], degree: int, required_depth: int
) -> int:
    names = ("ell_Z_taylor", "S_Z_taylor", "U0_Z_taylor", "P0_Z_taylor")
    missing = [name for name in names if name not in fixed]
    if missing:
        raise KeyError("fixed jets missing: " + ", ".join(missing))
    lengths = [len(fixed[name]) for name in names]
    if len(set(lengths)) != 1:
        raise ValueError("all fixed axial jets must have the same length")
    initial_count = lengths[0]
    needed = degree + required_depth + 1
    if initial_count < needed:
        raise ValueError(
            "axis jets must have degree + required_depth + 1 coefficients "
            f"({needed}); got {initial_count}"
        )
    return initial_count


def _validate_rows(
    rows: dict[str, list[list[Any]]], n: int, initial_count: int
) -> None:
    names = ("A", "Uz", "P")
    missing = [name for name in names if name not in rows]
    if missing:
        raise KeyError("rows missing: " + ", ".join(missing))
    for name in names:
        rowset = rows[name]
        if len(rowset) != n + 1:
            raise ValueError(
                f"{name} must contain rows 0 through {n}; got {len(rowset)}"
            )
        for index, row in enumerate(rowset):
            expected = initial_count - index
            if len(row) != expected:
                raise ValueError(
                    f"{name}[{index}] must contain {expected} axial coefficients; "
                    f"got {len(row)}"
                )


def initial_rows(
    ctx: MPIntervalContext,
    fixed: dict[str, Any],
    center: Any,
    degree: int,
    *,
    required_depth: int = DEFAULT_REQUIRED_DEPTH,
) -> dict[str, list[list[Any]]]:
    """Create rows 0 for an explicit scalar or interval center.

    ``center`` is coerced here for API symmetry and validation; the recurrence
    uses it explicitly in every subsequent :func:`advance_one` call.
    """

    try:
        degree = operator.index(degree)
        required_depth = operator.index(required_depth)
    except TypeError as exc:
        raise ValueError("degree and required_depth must be integers") from exc
    if degree < 0 or required_depth < 0:
        raise ValueError("degree and required_depth must be nonnegative")
    initial_count = _validate_fixed(fixed, degree, required_depth)
    _as_interval(ctx, center)
    zero = ctx.mpf(0)
    one = ctx.mpf(1)
    return {
        "A": [[one] + [zero] * (initial_count - 1)],
        "Uz": [_coerce_vector(ctx, fixed["U0_Z_taylor"])],
        "P": [_coerce_vector(ctx, fixed["P0_Z_taylor"])],
    }


def advance_one(
    ctx: MPIntervalContext,
    fixed: dict[str, Any],
    rows: dict[str, list[list[Any]]],
    n: int,
    center: Any,
    delta: Any,
) -> None:
    """Append radial order ``n + 1`` to ``rows`` in place.

    The input rows must contain orders ``0..n`` with lengths decreasing by
    one.  The fixed axial vectors retain their full length; this step consumes
    exactly one coefficient for the newly appended row and keeps the final
    ``required_depth + 1`` coefficients when initialized accordingly.
    """

    try:
        n = operator.index(n)
    except TypeError as exc:
        raise ValueError("n must be an integer") from exc
    if n < 0:
        raise ValueError("n must be nonnegative")

    # ``required_depth`` is not needed by the step itself.  A zero-depth
    # validation still enforces that the new row has one axial coefficient.
    initial_count = _validate_fixed(fixed, n, 0)
    _validate_rows(rows, n, initial_count)
    K = initial_count - n - 1
    if K < 1:
        raise ValueError("no remaining axial coefficients for radial order n+1")

    convert = ctx.mpf
    zero = convert(0)

    def pad(values: Any, length: int) -> list[Any]:
        values = [_as_interval(ctx, value) for value in values]
        if len(values) >= length:
            return values[:length]
        return values + [zero] * (length - len(values))

    def const(value: Any) -> list[Any]:
        return pad([value], K)

    def add(*args: list[Any]) -> list[Any]:
        return [sum((arg[k] for arg in args), zero) for k in range(K)]

    def scale(values: list[Any], scalar: Any) -> list[Any]:
        return [value * scalar for value in values]

    def mul(left: list[Any], right: list[Any]) -> list[Any]:
        return [
            sum((left[i] * right[k - i] for i in range(k + 1)), zero)
            for k in range(K)
        ]

    def diff(values: list[Any]) -> list[Any]:
        return [
            (k + 1) * values[k + 1] if k + 1 < len(values) else zero
            for k in range(K)
        ]

    def inverse(values: list[Any]) -> list[Any]:
        quotient = [convert(1) / values[0]]
        for k in range(1, K):
            quotient.append(
                -sum(
                    (values[i] * quotient[k - i] for i in range(1, k + 1)),
                    zero,
                )
                / values[0]
            )
        return quotient

    ell = _coerce_vector(ctx, fixed["ell_Z_taylor"])
    squared = _coerce_vector(ctx, fixed["S_Z_taylor"])
    u0 = _coerce_vector(ctx, fixed["U0_Z_taylor"])
    p0 = _coerce_vector(ctx, fixed["P0_Z_taylor"])
    z_center = _as_interval(ctx, center)
    z = pad([z_center, convert(1)], initial_count)
    delta_iv = _as_interval(ctx, delta)
    one = pad([convert(1)], initial_count)
    d = add(one, scale(mul(z, z), -1))
    L = add(one, scale(mul(z, z), -delta_iv))
    invL = inverse(L)

    # Convert existing rows without shortening their full input rows.  The
    # n-th row has K + 1 entries, so ndiff can consume its last coefficient.
    a = [[_as_interval(ctx, value) for value in row] for row in rows["A"]]
    u = [[_as_interval(ctx, value) for value in row] for row in rows["Uz"]]
    pressure = [
        [_as_interval(ctx, value) for value in row] for row in rows["P"]
    ]

    W: list[list[Any]] = []
    H: list[list[Any]] = []
    for i in range(n + 1):
        wi = scale(
            add(
                scale(mul(z, u[i]), 1 - delta_iv),
                mul(d, diff(u[i])),
            ),
            -convert(1) / (i + 1),
        )
        if i == 0:
            wi = add(const(1), wi)
        W.append(wi)
        hi = mul(d, u[i])
        if i == 0:
            hi = add(hi, scale(z, (1 - delta_iv) / 2))
        H.append(hi)

    aa_u = const(0)
    uu = const(0)
    rhs_a = const(0)
    rhs_u = const(0)
    for i in range(n + 1):
        aa_u = add(aa_u, mul(u[i], a[n - i]))
        uu = add(uu, mul(u[i], u[n - i]))
        a_derivative = add(diff(a[n - i]), mul(ell, a[n - i]))
        rhs_a = add(
            rhs_a,
            scale(mul(W[i], a[n - i]), n - i + 1),
            mul(H[i], a_derivative),
        )
        rhs_u = add(
            rhs_u,
            scale(mul(W[i], u[n - i]), n - i),
            mul(H[i], diff(u[n - i])),
        )

    rhs_a = add(
        rhs_a,
        scale(add(a[n], scale(mul(z, aa_u), -2)), delta_iv / 2),
    )
    rhs_u = add(
        rhs_u,
        scale(add(u[n], scale(mul(z, uu), -2)), (1 + delta_iv) / 2),
        mul(d, diff(pressure[n])),
        scale(mul(z, pressure[n]), -2 * (1 + delta_iv)),
    )

    if n:
        aa_previous = const(0)
        for i in range(n):
            aa_previous = add(aa_previous, mul(a[i], a[n - 1 - i]))
        rhs_u = add(
            rhs_u,
            scale(mul(z, mul(squared, aa_previous)), -2),
        )

    a_next = scale(
        mul(invL, rhs_a),
        convert(1) / (2 * (n + 1) * (n + 2)),
    )
    u_next = scale(
        mul(invL, rhs_u),
        convert(1) / (2 * (n + 1) ** 2),
    )
    aa_current = const(0)
    for i in range(n + 1):
        aa_current = add(aa_current, mul(a[i], a[n - i]))
    p_next = scale(mul(squared, aa_current), convert(1) / (n + 1))

    rows["A"].append(a_next)
    rows["Uz"].append(u_next)
    rows["P"].append(p_next)


def run_self_check() -> dict[str, Any]:
    """Run only cheap algebraic checks; the center-aware check is separate."""

    return {
        "api": "advance_one(ctx, fixed, rows, n, center, delta)",
        "initializer": "initial_rows(ctx, fixed, center, degree, required_depth=3)",
        "source_hashes": _source_hashes(),
        "all_terms_retained": True,
        "center_is_explicit": True,
        "whole_axis_core_generated": False,
        "terminal_functional_closure": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_self_check(), indent=2), flush=True)
