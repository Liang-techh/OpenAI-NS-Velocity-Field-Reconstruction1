"""Resolve the narrow residual spike just outside the delayed U start."""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    tau = .5*2**(-5.5)
    xs = np.array([1.0125, 1.015, 1.0175, 1.02])
    eta = np.full(len(xs), .3)
    rows = []
    for label, filename in (
        ('original', 'wide_taper_curvature_optimize.json'),
        ('delayed', 'delayed_taper_curvature_optimize.json')):
        field = CoupledMomentPhysicalLift(slice_filename=filename)
        points = field.compact.joined.inner.from_similarity(xs, eta, tau)
        for factor in (.001, .0005, .00025):
            residual, divergence = independent_fd(
                field, points, tau, factor*np.sqrt(field.nu*tau),
                .00025*tau)
            row = dict(name=label, spatial_step_factor=factor,
                       nodes=[dict(X=float(x),
                                   norm=float(np.linalg.norm(R)),
                                   residual=R.tolist(),
                                   divergence=float(d))
                              for x, R, d in zip(xs, residual, divergence)])
            rows.append(row)
            print(json.dumps(dict(name=label, factor=factor,
                                  norms=[node['norm'] for node in row['nodes']],
                                  max_divergence=max(abs(node['divergence'])
                                                     for node in row['nodes']))),
                  flush=True)
    report = dict(tau=tau, eta=.3, rows=rows,
                  scope='Four narrow physical X nodes with spatial '
                        'step refinement. No continuous maximum or '
                        'volume-L2 admission.', accepted=False)
    (ROOT/'delayed_taper_hotspot_refinement.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
