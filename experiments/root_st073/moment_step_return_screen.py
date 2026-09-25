"""Compare earlier streamfunction returns after the local physical-cone point."""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from extended_physical_cone_map import physical_cone
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


RETURNS = ((1.05, 1.3), (1.05, 1.5), (1.1, 1.5),
           (1.25, 1.6), (1.5, 2.0), (1.75, 3.0))


def moments(field, eta, tau, start, end, order=64):
    edges = sorted(set((0., 3.5, field.compact.joined.join_X,
                        1.5, start, end)))
    g, w = leggauss(order)
    X = np.concatenate([(a+b)/2+(b-a)*g/2
                        for a, b in zip(edges[:-1], edges[1:])])
    weights = np.concatenate([(b-a)*w/2
                              for a, b in zip(edges[:-1], edges[1:])])
    U, E = profile(field, X, eta, tau)
    return moment_vector(U, E, X, weights)


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    reference = RadialMomentStep(base, -.5)
    at_cone = [physical_cone(reference, 1., eta, tau, order=64)
               for eta in (.2, .3)]
    rows = []
    for start, end in RETURNS:
        field = RadialMomentStep(base, -.5, start, end)
        moments_delta = []
        for eta in (.2, .3):
            delta = (moments(field, eta, tau, start, end)
                     -moments(base, eta, tau, start, end))
            moments_delta.append(dict(eta=eta, delta=delta.tolist()))
        row = dict(restore_start=start, restore_end=end,
                   cone=at_cone, outgoing=moments_delta)
        rows.append(row)
        print(json.dumps(dict(restore_start=start, restore_end=end,
                              strict_pass=[c['strict_pass'] for c in at_cone],
                              outgoing=moments_delta)), flush=True)
    report = dict(tau=tau, shear=2., moment=-.5, rows=rows,
                  scope='Finite cone points and normalized outgoing moments for alternative solenoidal return geometries. A return beginning after X=1 leaves the field at X=1 unchanged. No full cone or momentum acceptance.',
                  accepted=False)
    (ROOT/'moment_step_return_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
