"""Numerically stable solver for the implicit Eq. (4.5) scale q.

For h >= 0, tau = 1-t > 0 and finite z, solve
    q - z^2 q^(-2h) = tau
for the unique positive q >= tau.

This module only implements the coordinate backbone. It does not identify
the paper profiles v0, F or U and makes no paper-exact or PDE-valid claim.
"""
from __future__ import annotations

import numpy as np


def solve_eq45_q(z, time, h, *, rtol=1e-12, atol=1e-14, max_iter=96):
    """Return the unique positive Eq. (4.5) scale q, vectorized over z/time.

    ``time`` must satisfy ``time < 1`` and ``h`` must be a finite scalar >= 0.
    For z=0 the exact solution q=1-time is returned. For h=0 the exact
    solution q=(1-time)+z^2 is returned. Otherwise monotone bisection is used.
    """
    z = np.asarray(z, dtype=float)
    time = np.asarray(time, dtype=float)
    try:
        z, time = np.broadcast_arrays(z, time)
    except ValueError as exc:
        raise ValueError("z and time must be broadcastable") from exc

    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(time)):
        raise ValueError("z and time must be finite")
    if np.any(time >= 1.0):
        raise ValueError("Eq. (4.5) q solver requires time < 1")
    if not np.isscalar(h) or not np.isfinite(h) or h < 0:
        raise ValueError("h must be a finite scalar >= 0")
    if not np.isfinite(rtol) or not np.isfinite(atol) or rtol < 0 or atol < 0:
        raise ValueError("rtol and atol must be finite and nonnegative")
    if not isinstance(max_iter, (int, np.integer)) or max_iter <= 0:
        raise ValueError("max_iter must be a positive integer")

    tau = 1.0 - time
    z2 = z * z

    if h == 0:
        return tau + z2

    # f(q)=q-z^2*q^(-2h)-tau is strictly increasing for q>0 when h>=0.
    lo = tau.copy()
    hi = np.maximum(tau + z2 + 1.0, 2.0 * tau)

    def f(q):
        return q - z2 * np.power(q, -2.0 * h) - tau

    f_hi = f(hi)
    for _ in range(max_iter):
        mask = f_hi < 0.0
        if not np.any(mask):
            break
        hi = np.where(mask, 2.0 * hi, hi)
        f_hi = f(hi)
    else:
        raise RuntimeError("failed to bracket Eq. (4.5) q root")

    # Exact axis values do not need iteration and avoid gratuitous roundoff.
    axis = z2 == 0.0
    lo = np.where(axis, tau, lo)
    hi = np.where(axis, tau, hi)

    for _ in range(max_iter):
        width = hi - lo
        tol = atol + rtol * np.maximum(np.abs(lo), np.abs(hi))
        active = width > tol
        if not np.any(active):
            break
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        go_right = (f_mid < 0.0) & active
        lo = np.where(go_right, mid, lo)
        hi = np.where((~go_right) & active, mid, hi)
    else:
        width = hi - lo
        tol = atol + rtol * np.maximum(np.abs(lo), np.abs(hi))
        if np.any(width > tol):
            raise RuntimeError("Eq. (4.5) q solver did not converge")

    q = 0.5 * (lo + hi)
    q = np.where(axis, tau, q)
    if not np.all(np.isfinite(q)) or np.any(q <= 0.0):
        raise RuntimeError("Eq. (4.5) q solver produced a nonpositive/nonfinite value")
    return q


def eq45_q_residual(q, z, time, h):
    """Evaluate q - z^2 q^(-2h) - (1-t) for diagnostics."""
    q = np.asarray(q, dtype=float)
    z = np.asarray(z, dtype=float)
    time = np.asarray(time, dtype=float)
    q, z, time = np.broadcast_arrays(q, z, time)
    if not np.all(np.isfinite(q)) or np.any(q <= 0):
        raise ValueError("q must be finite and positive")
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(time)):
        raise ValueError("z and time must be finite")
    if not np.isscalar(h) or not np.isfinite(h) or h < 0:
        raise ValueError("h must be a finite scalar >= 0")
    return q - z * z * np.power(q, -2.0 * h) - (1.0 - time)
