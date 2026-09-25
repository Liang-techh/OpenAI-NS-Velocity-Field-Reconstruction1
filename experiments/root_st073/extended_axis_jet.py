"""Implicit bivariate axis jets beyond the old |eta|<=0.5 expansion limit.

This is an experimental alternative to the Lagrange m-series in full_radial.
It solves the implicit source-coordinate equation degree by degree, avoiding
the geometric convergence restriction of that particular expansion.
"""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from radial_continuation import FullRadialField
from full_radial import LD, Parameters, multiply


def _power_next(q, p, i, j, beta):
    """Coefficient of q**beta using q * d(p) = beta * p * d(q)."""
    qa = q[:i+1, :j+1]
    pa = p[:i+1, :j+1]
    qr = q[i::-1, j::-1]
    pr = p[i::-1, j::-1]
    if i:
        weight = (i-np.arange(i+1, dtype=LD))[:, None]
        degree = i
    else:
        weight = (j-np.arange(j+1, dtype=LD))[None, :]
        degree = j
    return (beta*np.sum(pa*qr*weight)-np.sum(qa*pr*weight))/degree


def q_power_implicit_jet(beta, eta, shape, h=.005):
    """Taylor coefficients of Q**beta at (eta,1-eta**2).

    Q-(eta+zeta)^2 Q**(2h)=1-eta**2+theta, Q(0,0)=1.
    """
    e, gamma = LD(eta), LD(2*h)
    if abs(e) >= 1:
        raise ValueError('Only |eta|<1 is supported')
    if 1-gamma*e*e <= 0:
        raise ValueError('Degenerate implicit derivative')
    q = np.zeros(shape, dtype=LD)
    pg = np.zeros(shape, dtype=LD)
    q[0, 0] = pg[0, 0] = 1
    for degree in range(1, sum(shape)-1):
        for i in range(max(0, degree-shape[1]+1), min(shape[0]-1, degree)+1):
            j = degree-i
            known_pg = _power_next(q, pg, i, j, gamma)
            rhs = e*e*known_pg
            if i >= 1:
                rhs += 2*e*pg[i-1, j]
            if i >= 2:
                rhs += pg[i-2, j]
            if i == 0 and j == 1:
                rhs += 1
            q[i, j] = rhs/(1-gamma*e*e)
            pg[i, j] = gamma*q[i, j]+known_pg
    p = np.zeros(shape, dtype=LD)
    p[0, 0] = 1
    for degree in range(1, sum(shape)-1):
        for i in range(max(0, degree-shape[1]+1), min(shape[0]-1, degree)+1):
            j = degree-i
            p[i, j] = _power_next(q, p, i, j, LD(beta))
    return p


@dataclass(frozen=True)
class ExtendedParameters(Parameters):
    eta_max: float = .85

    def check(self):
        if not (1 <= self.order <= 18 and self.jet_extra >= 3
                and self.axis_terms >= 60):
            raise ValueError('Unsupported numerical settings')
        if not (0 < self.h < .01 and self.nu > 0 and self.axis_swirl > 0
                and 0 < self.eta_max < 1 and 0 < self.X_max <= .25):
            raise ValueError('Invalid extended physical settings')


class ExtendedFullRadialField(FullRadialField):
    @classmethod
    def load(cls, path):
        import json
        from pathlib import Path
        return cls(ExtendedParameters(**json.loads(Path(path).read_text())['parameters']))

    @lru_cache(maxsize=3000)
    def axis_jets(self, e):
        p = self.p
        shape = (2*p.order+p.jet_extra+2, p.order+p.jet_extra+1)
        powers = {b: q_power_implicit_jet(b, e, shape, p.h)
                  for b in [-1-p.h, -1., -self.A, -2*self.A, -2.]}
        z = np.zeros(shape, dtype=LD)
        z[0, 0] = LD(e)
        z[1, 0] = 1
        B = LD(p.axis_swirl)*powers[-1-p.h]
        C = LD(p.axial_slope)*multiply(z, powers[-1.], shape) \
            + LD(p.axial_bias)*powers[-self.A]
        P = -LD(p.pressure)*powers[-2*self.A] \
            + LD(.5*p.pressure)*multiply(multiply(z, z, shape),
                                          powers[-2.], shape)
        return B, C, P
