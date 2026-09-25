"""Original curvature lift on the delayed-candidate dense comparison grid."""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    field = CoupledMomentPhysicalLift(
        slice_filename='wide_taper_curvature_optimize.json')
    tau = .5*2**(-5.5)
    xs = (1.005, 1.01, 1.0125, 1.015, 1.0175,
          1.02, 1.025, 1.1, 2., 2.975)
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
                  scope='Original curvature lift on the same 20 physical '
                        'nodes and FD step as delayed005 screens. '
                        'No continuous maximum or volume-L2 admission.',
                  accepted=False)
    (ROOT/'original_dense_onset_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in
                      ('max_full_momentum', 'max_divergence', 'hotspot')}),
          flush=True)


if __name__ == '__main__':
    run()
