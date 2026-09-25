"""Spatial finite-difference refinement of the delayed wide-taper peak."""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    field = CoupledMomentPhysicalLift(
        slice_filename='delayed005_wide04_curvature_optimize.json')
    tau = .5*2**(-5.5)
    xs = np.array([1.0125, 1.015])
    points = field.compact.joined.inner.from_similarity(
        xs, np.full(len(xs), .3), tau)
    rows = []
    for factor in (.001, .0005, .00025, .000125):
        residual, divergence = independent_fd(
            field, points, tau, factor*np.sqrt(field.nu*tau),
            .00025*tau)
        row = dict(spatial_step_factor=factor,
                   nodes=[dict(X=float(x),
                               residual=R.tolist(),
                               norm=float(np.linalg.norm(R)),
                               divergence=float(d))
                          for x, R, d in zip(xs, residual, divergence)])
        rows.append(row)
        print(json.dumps(dict(factor=factor,
                              norms=[item['norm'] for item in row['nodes']],
                              divergence=[item['divergence']
                                          for item in row['nodes']])),
              flush=True)
    report = dict(tau=tau, eta=.3, rows=rows,
                  scope='Two physical nodes around the sampled peak; '
                        'spatial FD refinement only. No global maximum '
                        'or volume-L2 admission.', accepted=False)
    (ROOT/'delayed_wide04_fd_refinement.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
