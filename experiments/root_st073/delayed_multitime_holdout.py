"""Independent space/time screen for a delayed physical lift candidate."""

import argparse
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, residual
from radial_continuation import ROOT


def run(candidate_name='delayed005_wide04_multitime_tangent.json',
        output_name='delayed_multitime_holdout.json'):
    names = ('delayed005_wide04_reoptimized_E.json',
             candidate_name)
    xs = (1.01125, 1.01625, 1.0225, 1.07, 2.975)
    etas = (.29, .31)
    times = (.5*2**(-5.4), .5*2**(-5.15))
    fields = [CoupledMomentPhysicalLift(slice_filename=name)
              for name in names]
    rows = []
    for tau in times:
        points, X, eta = nodes(fields[0], xs, etas, tau)
        norms = []
        for field in fields:
            R, _ = residual(field, points, tau)
            norms.append(np.linalg.norm(R, axis=1))
        group = []
        for lo, hi, label in ((0, 6, 'entrance'),
                              (6, 8, 'middle'),
                              (8, 10, 'far_return')):
            before, after = (float(np.max(n[lo:hi])) for n in norms)
            group.append(dict(region=label, before_max=before,
                              after_max=after, ratio=after/before))
        row = dict(tau=tau, groups=group,
                   before_max=float(np.max(norms[0])),
                   after_max=float(np.max(norms[1])),
                   nodes=[dict(X=float(x), eta=float(e),
                               before=float(b), after=float(a))
                          for x, e, b, a in zip(X, eta, *norms)])
        rows.append(row)
        print(json.dumps({k: row[k] for k in
                          ('tau', 'before_max', 'after_max', 'groups')}),
              flush=True)
    report = dict(slice_filenames=names, X=xs, eta=etas,
                  rows=rows, scope='Two unseen times and disjoint spatial '
                        'nodes in entrance/middle/far return; local full '
                        'Cartesian finite-difference momentum only.',
                  accepted=False)
    (ROOT/output_name).write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--candidate-name',
                        default='delayed005_wide04_multitime_tangent.json')
    parser.add_argument('--output-name',
                        default='delayed_multitime_holdout.json')
    args = parser.parse_args()
    run(args.candidate_name, args.output_name)
