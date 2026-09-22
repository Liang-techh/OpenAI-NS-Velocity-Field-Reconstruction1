"""Source-native implicit similarity coordinates; not a velocity-field solver.

OpenAI PDF (3.2): tau=q(1-eta^2), z=q^(1/2-h)*eta,
X=r^2/(2q).  The q used here depends on z and tau. It is NOT the
ST064--66 dimensionless scale ratio 2**(-k). This module uses the source
unit-viscosity coordinates; conversion from nu=.01 needs the separate
x/sqrt(nu) transformation before reusing this map in a physical evaluator.
"""
from __future__ import annotations
import numpy as np


def coordinates(r, z, tau, h: float = 0.005):
    """Return q, eta, X and first physical derivatives using safeguarded inversion.

    Bisection is performed in dimensionless q/max(tau,abs(z)**(1/D))
    on [1,4]. Keeping tau/q separately avoids falsely resolving 1-eta^2
    by cancellation when eta rounds to an endpoint. Derivatives are first
    derivatives only; no new velocity, pressure or PDE result is implied.
    """
    if not np.isfinite(h) or not 0 < h < .01:
        raise ValueError('Source regime requires 0<h<0.01')
    r, z, tau = np.broadcast_arrays(np.asarray(r, float), np.asarray(z, float), np.asarray(tau, float))
    if not (np.isfinite(r).all() and np.isfinite(z).all() and np.isfinite(tau).all()):
        raise ValueError('Finite inputs required')
    if np.any(r < 0) or np.any(tau <= 0):
        raise ValueError('Require r>=0 and tau>0; singular endpoint excluded')
    D = .5-h
    with np.errstate(over='raise', under='ignore', invalid='raise', divide='raise'):
        qz = np.abs(z)**(1/D)
        scale = np.maximum(tau, qz)
        c = (qz/scale)**(1-2*h)
        b = tau/scale
        lo, hi = np.ones_like(scale), np.full_like(scale, 4.)
        for _ in range(60):
            mid = (lo+hi)/2
            f = mid-c*mid**(2*h)-b
            lo = np.where(f < 0, mid, lo)
            hi = np.where(f >= 0, mid, hi)
        Q = (lo+hi)/2
        q = scale*Q
        eta = z/q**D
        eta_gap = tau/q
        L = 1-2*h*(1-eta_gap)
        qt = 1/L  # q_tau, not q_physical_time
        qzder = 2*z*q**(2*h)/L
        X = r*r/(2*q)
        out = dict(q=q, eta=eta, X=X, one_minus_eta_squared=eta_gap,
                   q_tau=qt, q_z=qzder, eta_tau=-D*eta*qt/q,
                   eta_z=q**(-D)-D*eta*qzder/q,
                   X_r=r/q, X_z=-X*qzder/q, X_tau=-X*qt/q,
                   inversion_scaled_error=Q-c*Q**(2*h)-b)
    if not all(np.isfinite(v).all() for v in out.values()):
        raise ValueError('Inputs exceed representable coordinate/derivative range')
    return out
