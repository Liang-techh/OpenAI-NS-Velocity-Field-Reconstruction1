"""Test wider axial taper for the fixed-moment physical lift."""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, residual
from radial_continuation import ROOT


def run():
    slice_name = 'delayed005_wide04_multitime_tangent.json'
    xs = (1.015, 1.05, 2.975)
    etas = (.2, .3, .31, .35, .4, .45)
    times = (.5*2**(-5.5), .5*2**(-5.25))
    rows = []
    for fall_end in (.4, .45, .5):
        field = CoupledMomentPhysicalLift(slice_filename=slice_name)
        field.axial_fall_end = fall_end
        screens = []
        for tau in times:
            points, X, eta = nodes(field, xs, etas, tau)
            R, divergence = residual(field, points, tau)
            norms = np.linalg.norm(R, axis=1)
            screen = dict(tau=tau, max=float(np.max(norms)),
                          max_divergence=float(np.max(np.abs(divergence))),
                          nodes=[dict(X=float(x), eta=float(e),
                                      residual=float(n))
                                 for x, e, n in zip(X, eta, norms)])
            screens.append(screen)
        rows.append(dict(axial_fall_end=fall_end, screens=screens))
        print(json.dumps(dict(axial_fall_end=fall_end,
                              maxima=[s['max'] for s in screens])),
              flush=True)
    report = dict(slice_filename=slice_name, X=xs, eta=etas, rows=rows,
                  scope='Three radial nodes across six axial samples and '
                        'two times. The reference eta=.2,.3 slice '
                        'profiles remain unchanged. Local full momentum '
                        'only, no continuous cone or volume certificate.',
                  accepted=False)
    (ROOT/'delayed_axial_support_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
