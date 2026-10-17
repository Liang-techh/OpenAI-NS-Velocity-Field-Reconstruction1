"""Directed ordinary-derivative jets and their absolute majorants."""

import math
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sp


def _check(ctx, jet):
    if not isinstance(ctx, MPIntervalContext) or not isinstance(jet, (list, tuple)) or not jet:
        raise TypeError("Use an MPIntervalContext and a nonempty interval jet")
    if any(not hasattr(x, "_mpi_") for x in jet):
        raise TypeError("Every jet entry must be a directed interval")
    return len(jet)


def _lower(x):
    return mp.make_mpf(x._mpi_[0])


def _upper(x):
    return mp.make_mpf(x._mpi_[1])


def _abs(ctx, x):
    return ctx.fabs(x)


def _majorant(ctx, x):
    return ctx.mpf([0, _upper(_abs(ctx, x))])


def jet_add(ctx, f, g):
    n, m = _check(ctx, f), _check(ctx, g)
    if n != m:
        raise ValueError("Jet orders must agree")
    return [f[k] + g[k] for k in range(n)]


def jet_neg(ctx, f):
    _check(ctx, f)
    return [-x for x in f]


def jet_mul(ctx, f, g):
    n, m = _check(ctx, f), _check(ctx, g)
    if n != m:
        raise ValueError("Jet orders must agree")
    return [sum((math.comb(k, j) * f[j] * g[k - j] for j in range(k + 1)), ctx.mpf(0))
            for k in range(n)]


def jet_derivative(ctx, f):
    """Shift an ordinary jet down one order without inventing a top row."""
    n = _check(ctx, f)
    if n < 2:
        raise ValueError("A derivative shift requires at least two jet rows")
    return list(f[1:])


def jet_reciprocal(ctx, f):
    n = _check(ctx, f)
    if _lower(f[0]) <= 0:
        raise ValueError("Reciprocal requires a positive lower bound for f0")
    q = [ctx.mpf(1) / f[0]]
    for order in range(1, n):
        rest = sum((math.comb(order, k) * f[k] * q[order - k]
                    for k in range(1, order + 1)), ctx.mpf(0))
        q.append(-q[0] * rest)
    return q


def jet_div(ctx, f, g):
    _check(ctx, f); _check(ctx, g)
    return jet_mul(ctx, f, jet_reciprocal(ctx, g))


def jet_exp(ctx, g, base=None):
    n = _check(ctx, g)
    h = [ctx.exp(g[0]) if base is None else base]
    for order in range(1, n):
        h.append(sum((math.comb(order - 1, k - 1) * g[k] * h[order - k]
                      for k in range(1, order + 1)), ctx.mpf(0)))
    return h


def jet_log(ctx, f):
    n = _check(ctx, f)
    if _lower(f[0]) <= 0:
        raise ValueError("Log requires a positive lower bound for f0")
    g = [ctx.log(f[0])]
    for order in range(1, n):
        rest = sum((math.comb(order - 1, k) * f[k] * g[order - k]
                    for k in range(1, order)), ctx.mpf(0))
        g.append((f[order] - rest) / f[0])
    return g


def absolute_majorant(ctx, f):
    _check(ctx, f)
    return [_majorant(ctx, x) for x in f]


def bound_product(ctx, f, g):
    """Multiply nonnegative absolute-majorant jets with binomial factors."""
    n, m = _check(ctx, f), _check(ctx, g)
    if n != m:
        raise ValueError("Jet orders must agree")
    f = absolute_majorant(ctx, f); g = absolute_majorant(ctx, g)
    return [sum((math.comb(k, j) * f[j] * g[k - j] for j in range(k + 1)), ctx.mpf(0))
            for k in range(n)]


def bound_reciprocal(ctx, f, lower_f0=None):
    """Absolute reciprocal majorant using a certified positive lower f0."""
    n = _check(ctx, f)
    floor = f[0] if lower_f0 is None else lower_f0
    if _lower(floor) <= 0:
        raise ValueError("A positive lower f0 is required")
    a = absolute_majorant(ctx, f)
    q = [_majorant(ctx, ctx.mpf(1) / floor)]
    for order in range(1, n):
        rest = sum((math.comb(order, k) * a[k] * q[order - k]
                    for k in range(1, order + 1)), ctx.mpf(0))
        q.append(q[0] * rest)
    return q


def exp_bound(ctx_or_logjets, logjets=None, base_upper=None):
    """Absolute majorant for exp(logjets), with supplied positive base upper."""
    if isinstance(ctx_or_logjets, MPIntervalContext):
        ctx, logs = ctx_or_logjets, logjets
    else:
        logs, base_upper = ctx_or_logjets, logjets
        if not logs or not hasattr(logs[0], "ctx"):
            raise TypeError("Pass ctx explicitly when interval context is unavailable")
        ctx = logs[0].ctx
    n = _check(ctx, logs)
    if base_upper is None or _lower(base_upper) < 0:
        raise ValueError("A nonnegative exp base upper bound is required")
    a = absolute_majorant(ctx, logs)
    h = [_majorant(ctx, base_upper)]
    for order in range(1, n):
        h.append(sum((math.comb(order - 1, k - 1) * a[k] * h[order - k]
                      for k in range(1, order + 1)), ctx.mpf(0)))
    return h


def symbolic_controls():
    """Exact SymPy controls for ordinary, unfactorialized jet arithmetic."""
    x = sp.symbols("x")
    f = sp.Rational(3, 2) + 2*x - sp.Rational(5, 3)*x**2 + sp.Rational(7, 4)*x**3 - sp.Rational(2, 5)*x**4
    g = sp.Rational(5, 3) - x + sp.Rational(4, 7)*x**2 + sp.Rational(3, 2)*x**3 - sp.Rational(1, 3)*x**4
    d = lambda h: [sp.diff(h, x, k).subs(x, 0) for k in range(5)]
    F, G = d(f), d(g)
    product = [sum(sp.binomial(k, j)*F[j]*G[k-j] for j in range(k+1)) for k in range(5)]
    reciprocal = [1 / F[0]]
    logjet = [sp.log(F[0])]
    exp_zero = [sp.Integer(1)]
    for order in range(1, 5):
        reciprocal.append(-reciprocal[0] * sum(sp.binomial(order, k) * F[k] * reciprocal[order-k]
                                               for k in range(1, order+1)))
        logjet.append((F[order] - sum(sp.binomial(order-1, k) * F[k] * logjet[order-k]
                                        for k in range(1, order))) / F[0])
        exp_zero.append(sum(sp.binomial(order-1, k-1) * F[k] * exp_zero[order-k]
                            for k in range(1, order+1)))
    reciprocal_expected, log_expected, exp_expected = d(1/f), d(sp.log(f)), d(sp.exp(f-f.subs(x, 0)))
    checks = {
        "add": all(F[k] + G[k] == d(f+g)[k] for k in range(5)),
        "neg": all(-F[k] == d(-f)[k] for k in range(5)),
        "mul_binomial": product == d(f*g),
        "reciprocal_nonunit_f0": all(sp.simplify(reciprocal[k] - reciprocal_expected[k]) == 0 for k in range(5)),
        "log_recurrence": all(sp.simplify(logjet[k] - log_expected[k]) == 0 for k in range(5)),
        "exp_recurrence": all(sp.simplify(exp_zero[k] - exp_expected[k]) == 0 for k in range(5)),
        "derivative_shift": all(F[k+1] == d(sp.diff(f, x))[k] for k in range(4)),
        "ordinary_factorial": d(x**2)[2] == 2 and d(x**2)[1] == 0,
    }
    if not all(checks.values()):
        raise ArithmeticError("Ordinary jet symbolic controls failed")
    return {"passed": True, "order": 4, "factorial_normalized": False, "checks": checks,
            "nonunit_f0": str(F[0]), "no_source_claim": True}


__all__ = ["jet_add", "jet_neg", "jet_mul", "jet_reciprocal", "jet_div", "jet_exp",
           "jet_log", "jet_derivative", "absolute_majorant", "bound_product",
           "bound_reciprocal", "exp_bound", "symbolic_controls"]
