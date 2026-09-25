"""Measure all five outgoing moment changes after the new M-step restores."""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def moments(field, eta, tau, order=64):
    # Split at every join and at both endpoints of the compensating step.
    edges = (0., field.compact.joined.join_X, 1.5, 1.75, 3., 3.5)
    g, w = leggauss(order)
    xs = np.concatenate([(a+b)/2+(b-a)*g/2
                         for a, b in zip(edges[:-1], edges[1:])])
    ws = np.concatenate([(b-a)*w/2
                         for a, b in zip(edges[:-1], edges[1:])])
    U, E = profile(field, xs, eta, tau)
    H = np.sqrt(2*xs)*E
    return np.array([ws@U, ws@H, ws@(U*H),
                     ws@(U**2-E**2/2), ws@(E**2/(2*xs))])


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    field = RadialMomentStep(base, -.5)
    rows = []
    for eta in (.2, .3):
        original = moments(base, eta, tau)
        changed = moments(field, eta, tau)
        row = dict(eta=eta, base=original.tolist(),
                   changed=changed.tolist(),
                   delta=(changed-original).tolist())
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, x_out=3.5, quadrature_per_piece=64,
                  rows=rows,
                  scope='Five normalized radial moment differences at X=3.5, after exact physical field restoration at X>=3. The M step can control the old X=1.5 exit yet leaves nonlinear J and S moment defects that a later correction must remove. Finite axial/time slices, no global matching or PDE acceptance.',
                  accepted=False)
    (ROOT/'moment_shear_outgoing_moments.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
