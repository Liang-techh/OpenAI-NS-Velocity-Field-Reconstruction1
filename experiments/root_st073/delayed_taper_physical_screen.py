"""Same-node physical residual comparison of original and delayed U starts."""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    tau = .5*2**(-5.5)
    xs = (1.005, 1.01, 1.025, 1.1, 2., 2.975, 2.99)
    etas = (.2, .25, .3)
    old = CoupledMomentPhysicalLift(
        slice_filename='wide_taper_curvature_optimize.json')
    new = CoupledMomentPhysicalLift(
        slice_filename='delayed_taper_curvature_optimize.json')
    X, eta = np.meshgrid(xs, etas, indexing='ij')
    points = new.compact.joined.inner.from_similarity(
        X.ravel(), eta.ravel(), tau)
    rows = []
    for label, field in (('original', old), ('delayed', new)):
        residual, divergence = independent_fd(
            field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
        norms = np.linalg.norm(residual, axis=1)
        top = int(np.argmax(norms))
        row = dict(name=label, node_count=len(points),
                   max_full_momentum=float(np.max(norms)),
                   max_divergence=float(np.max(np.abs(divergence))),
                   hotspot=dict(X=float(X.ravel()[top]),
                                eta=float(eta.ravel()[top]),
                                residual=residual[top].tolist()),
                   nodes=[dict(X=float(x), eta=float(e),
                               residual_norm=float(n))
                          for x, e, n in zip(X.ravel(), eta.ravel(), norms)])
        rows.append(row)
        print(json.dumps({key: row[key] for key in
                          ('name', 'node_count', 'max_full_momentum',
                           'max_divergence', 'hotspot')}), flush=True)
    report = dict(tau=tau, sample_X=xs, sample_eta=etas, rows=rows,
                  scope='Same 21 physical finite-difference nodes at one '
                        'late time. Two fitted slice moments only; no '
                        'volume L2, continuous cone, or PDE admission.',
                  accepted=False)
    (ROOT/'delayed_taper_physical_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
