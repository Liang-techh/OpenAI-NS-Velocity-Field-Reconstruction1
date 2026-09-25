"""Nearby-time same-node comparison for the delayed wide-taper lift."""
import argparse
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from radial_continuation import ROOT


def run(slice_filename='delayed005_wide04_curvature_optimize.json',
        output_name='delayed_wide04_time_holdout.json'):
    tau = .5*2**(-5.25)
    xs = (1.005, 1.01, 1.0125, 1.015, 1.0175,
          1.02, 1.025, 2.6, 2.975)
    X, eta = np.meshgrid(xs, (.2, .3), indexing='ij')
    fields = (
        ('original', 'wide_taper_curvature_optimize.json'),
        ('delayed_wide04', slice_filename))
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
    (ROOT/output_name).write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--slice-filename',
                        default='delayed005_wide04_curvature_optimize.json')
    parser.add_argument('--output-name',
                        default='delayed_wide04_time_holdout.json')
    arguments = parser.parse_args()
    run(arguments.slice_filename, arguments.output_name)
