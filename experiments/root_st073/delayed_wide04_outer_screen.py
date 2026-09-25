"""Screen the wider delayed lift across the middle and return regions."""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    field = CoupledMomentPhysicalLift(
        slice_filename='delayed005_wide04_curvature_optimize.json')
    tau = .5*2**(-5.5)
    xs = (1.05, 1.25, 1.5, 2., 2.5, 2.6, 2.7, 2.8,
          2.9, 2.95, 2.975, 2.99)
    X, eta = np.meshgrid(xs, (.2, .3), indexing='ij')
    points = field.compact.joined.inner.from_similarity(
        X.ravel(), eta.ravel(), tau)
    residual, divergence = independent_fd(
        field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
    norms = np.linalg.norm(residual, axis=1)
    top = int(np.argmax(norms))
    report = dict(tau=tau, sample_X=xs, sample_eta=[.2, .3],
                  max_full_momentum=float(np.max(norms)),
                  max_divergence=float(np.max(np.abs(divergence))),
                  hotspot=dict(X=float(X.ravel()[top]),
                               eta=float(eta.ravel()[top]),
                               residual=residual[top].tolist()),
                  nodes=[dict(X=float(x), eta=float(e),
                              residual_norm=float(n))
                         for x, e, n in zip(X.ravel(), eta.ravel(), norms)],
                  scope='Twenty-four physical nodes covering middle and '
                        'return regions at one time. No continuous '
                        'maximum or volume-L2 admission.', accepted=False)
    (ROOT/'delayed_wide04_outer_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in
                      ('max_full_momentum', 'max_divergence', 'hotspot')}),
          flush=True)


if __name__ == '__main__':
    run()
