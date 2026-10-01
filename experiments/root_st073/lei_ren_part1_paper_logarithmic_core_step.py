"""Direct r=Lambda R re-expression of the full Section 8.2 radial step.

Pressure is additionally scaled by epsilon to avoid a huge P0 coefficient.
Every frozen nonlinear term is retained; swirl is represented through the
scaled F0-squared jets and is never discarded. Seed generation is separate.
"""
from typing import Any
from mpmath.ctx_iv import MPIntervalContext
import operator
from lei_ren_part1_paper_functional_core_step import (
    _as_interval, _coerce_vector, _validate_fixed, _validate_rows,
)

def advance_scaled_one(
    ctx: MPIntervalContext,
    fixed: dict[str, Any],
    rows: dict[str, list[list[Any]]],
    n: int,
    center: Any,
    delta: Any,
    epsilon: Any,
) -> None:
    """Append directly normalized radial rows.

    A[n]=physical_A[n]/Lambda**n, Uz[n]=physical_Uz[n]/Lambda**n,
    P[n]=epsilon*physical_P[n]/Lambda**n, epsilon=1/Lambda.
    Fixed jets ell=epsilon*physical_ell and S=epsilon**2*F0**2.
    U0 is unchanged and P0 is epsilon*physical_P0. The ordinary axial
    Taylor convention and decreasing row lengths match the frozen step.
    This is radial profile generation, NOT temporal coefficient recursion.
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
    eps = _as_interval(ctx, epsilon)
    if not (eps > 0 and eps <= 1):
        raise ValueError('epsilon must lie in (0,1]')
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
        rhs_a = add(
            rhs_a,
            scale(mul(W[i], a[n - i]), eps * (n - i + 1)),
            mul(H[i], add(scale(diff(a[n - i]), eps), mul(ell, a[n - i]))),
        )
        rhs_u = add(
            rhs_u,
            scale(mul(W[i], u[n - i]), eps * (n - i)),
            scale(mul(H[i], diff(u[n - i])), eps),
        )

    rhs_a = add(
        rhs_a,
        scale(add(a[n], scale(mul(z, aa_u), -2)), eps * delta_iv / 2),
    )
    rhs_u = add(
        rhs_u,
        scale(add(u[n], scale(mul(z, uu), -2)), eps * (1 + delta_iv) / 2),
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

