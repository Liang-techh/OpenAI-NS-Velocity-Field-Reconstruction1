"""Nearby-time same-node comparison for the delayed wide-taper lift."""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    tau = .5*2**(-5.25)
    xs = (1.005, 1.01, 1.0125, 1.015, 1.0175,
          1.02, 1.025, 2.6, 2.975)
    X, eta = np.meshgrid(xs, (.2, .3), indexing='ij')
    fields = (
        ('original', 'wide_taper_curvature_optimize.json'),
        ('delayed_wide04', 'delayed005_wide04_curvature_optimize.json'))
    rows = []
    for name, filename in fields:
        field = CoupledMomentPhysicalLift(slice_filename=filename)
        points = field.compact.joined.inner.from_similarity(
            X.ravel(), eta.ravel(), tau)
        residual, divergence = independent_fd(
            field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
        norms = np.linalg.norm(residual, axis=1)
        top = int(np.argmax(norms))
        row = dict(name=name, max_full_momentum=float(np.max(norms)),
                   max_divergence=float(np.max(np.abs(divergence))),
                   hotspot=dict(X=float(X.ravel()[top]),
                                eta=float(eta.ravel()[top]),
                                residual=residual[top].tolist()),
                   nodes=[dict(X=float(x), eta=float(e),
                               residual_norm=float(n))
                          for x, e, n in zip(X.ravel(), eta.ravel(), norms)])
        rows.append(row)
        print(json.dumps({key: row[key] for key in
                          ('name', 'max_full_momentum', 'hotspot')}),
              flush=True)
    report = dict(tau=tau, k=5.25, sample_X=xs, sample_eta=[.2, .3],
                  rows=rows,
                  scope='Eighteen same physical similarity nodes at one '
                        'nearby time. No uniform-time or volume L2 '
                        'admission.', accepted=False)
    (ROOT/'delayed_wide04_time_holdout.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
