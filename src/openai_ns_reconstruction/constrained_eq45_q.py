"""Stable positive-root solver for the sourced Eq. (4.1)/(4.5) scale q.

For 0 < h < 1/2 and tau = 1-t > 0, solve

    q - z^2 q^(2h) = tau

for the unique physical root q >= tau.  The sourced similarity coordinates then
give D = 1/2-h and eta = z / q^D, so tau > 0 implies |eta| < 1 away from the
axis limit.

This module implements only the similarity-coordinate primitive.  It does not
identify the profiles v0, F or U, does not validate Navier--Stokes, and does not
claim recovery of the complete paper field.
"""
from __future__ import annotations

import numpy as np


def _validate_h(h):
    if not np.isscalar(h) or not np.isfinite(h) or not (0.0 < float(h) < 0.5):
        raise ValueError("h must be a finite scalar with 0 < h < 1/2")
    return float(h)


def solve_eq45_q(z, time, h, *, rtol=1e-12, atol=1e-14, max_iter=128):
    """Return the unique positive physical root of q-z^2*q^(2h)=1-time.

    Parameters are NumPy-broadcast over ``z`` and ``time``.  The implementation
    uses a sign-safe bracket followed by bisection.  On the symmetry axis the
    exact value ``q=1-time`` is returned.
    """
    h = _validate_h(h)
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
    if not np.isfinite(rtol) or not np.isfinite(atol) or rtol < 0 or atol < 0:
        raise ValueError("rtol and atol must be finite and nonnegative")
    if not isinstance(max_iter, (int, np.integer)) or max_iter <= 0:
        raise ValueError("max_iter must be a positive integer")

    tau = 1.0 - time
    z2 = z * z
    axis = z2 == 0.0

    # f(tau) = -z^2*tau^(2h) <= 0.  Because 2h<1, the linear q term
    # dominates at large q, so doubling eventually yields f(hi)>=0.
    lo = tau.copy()
    hi = np.maximum(2.0 * tau, tau + z2 + 1.0)

    def f(q):
        return q - z2 * np.power(q, 2.0 * h) - tau

    hi = np.where(axis, tau, hi)
    f_hi = f(hi)
    for _ in range(max_iter):
        mask = (~axis) & (f_hi < 0.0)
        if not np.any(mask):
            break
        hi = np.where(mask, 2.0 * hi, hi)
        f_hi = f(hi)
        if not np.all(np.isfinite(hi)) or not np.all(np.isfinite(f_hi)):
            raise RuntimeError("failed to bracket Eq. (4.5) q root finitely")
    else:
        raise RuntimeError("failed to bracket Eq. (4.5) q root")

    for _ in range(max_iter):
        width = hi - lo
        tol = atol + rtol * np.maximum(np.abs(lo), np.abs(hi))
        active = (~axis) & (width > tol)
        if not np.any(active):
            break
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        go_right = active & (f_mid < 0.0)
        lo = np.where(go_right, mid, lo)
        hi = np.where(active & (~go_right), mid, hi)
    else:
        width = hi - lo
        tol = atol + rtol * np.maximum(np.abs(lo), np.abs(hi))
        if np.any((~axis) & (width > tol)):
            raise RuntimeError("Eq. (4.5) q solver did not converge")

    q = np.where(axis, tau, 0.5 * (lo + hi))
    if not np.all(np.isfinite(q)) or np.any(q <= 0.0):
        raise RuntimeError("Eq. (4.5) q solver produced a nonpositive/nonfinite value")

    # Fail closed on the sourced physical branch: tau=q(1-eta^2)>0.
    D = 0.5 - h
    eta = z / np.power(q, D)
    eps = 32.0 * np.finfo(float).eps
    if np.any(np.abs(eta) > 1.0 + eps):
        raise RuntimeError("Eq. (4.5) q solver left the physical |eta|<=1 branch")
    return q


def eq45_q_residual(q, z, time, h):
    """Evaluate q-z^2*q^(2h)-(1-time) for diagnostics."""
    h = _validate_h(h)
    q = np.asarray(q, dtype=float)
    z = np.asarray(z, dtype=float)
    time = np.asarray(time, dtype=float)
    q, z, time = np.broadcast_arrays(q, z, time)
    if not np.all(np.isfinite(q)) or np.any(q <= 0):
        raise ValueError("q must be finite and positive")
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(time)):
        raise ValueError("z and time must be finite")
    return q - z * z * np.power(q, 2.0 * h) - (1.0 - time)
