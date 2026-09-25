"""Physical-space audit of the remote-patched mean-field increment."""
import json

import numpy as np

from joined_field import ROOT, independent_fd
from moment_matched_joined_field import MomentMatchedJoinedField
from paper_moment_bridge import slice_data


def run():
    field = MomentMatchedJoinedField()
    base = field.base
    rows = []
    for tau in (.0084, .032):
        for eta in (0., .3):
            X = np.array([5., 10., 15.])
            points = field.inner.from_similarity(X, eta, tau)
            hs = .0005*np.sqrt(field.nu*tau)
            ht = .0001*tau
            baseline_residual, baseline_div = independent_fd(
                base, points, tau, hs, ht)
            corrected_residual, corrected_div = independent_fd(
                field, points, tau, hs, ht)
            u0, _ = base.fields(points, tau)
            u1, _ = field.fields(points, tau)
            row = {'tau': tau, 'eta': eta, 'X': X.tolist(),
                   'baseline_residual_norms': np.linalg.norm(
                       baseline_residual, axis=1).tolist(),
                   'corrected_residual_norms': np.linalg.norm(
                       corrected_residual, axis=1).tolist(),
                   'baseline_divergence_max': float(np.max(np.abs(baseline_div))),
                   'corrected_divergence_max': float(np.max(np.abs(corrected_div))),
                   'correction_speed_max': float(np.max(np.linalg.norm(u1-u0, axis=1)))}
            rows.append(row)
            print(json.dumps(row), flush=True)
    moment_holdout = []
    for tau in (.012, .064):
        for eta in (-.35, -.15, .15, .35):
            row = slice_data(field, eta, tau, 48,
                             field.patch_start, field.patch_end)
            moment_holdout.append({
                'tau': tau, 'eta': eta,
                'angular_moment': row['baseline_angular_moment'],
                'kinetic_moment': row['kinetic_moment_baseline']})
    report = {'rows': rows, 'moment_holdout': moment_holdout,
              'moment_holdout_max_angular': max(abs(row['angular_moment'])
                                                for row in moment_holdout),
              'moment_holdout_max_kinetic': max(abs(row['kinetic_moment'])
                                                for row in moment_holdout),
              'scope': 'Physical full-momentum and divergence finite differences at selected remote-patch points. The correction realizes slice moment coefficients through a smooth eta polynomial and a streamfunction. Large residual is expected before stress realization. Joined heat exterior is not finite-energy in z, and no critical-time forcing or global field is claimed.',
              'accepted': False}
    (ROOT/'moment_matched_joined_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'moment_holdout_max_angular': report['moment_holdout_max_angular'],
                      'moment_holdout_max_kinetic': report['moment_holdout_max_kinetic']}),
          flush=True)


if __name__ == '__main__':
    run()
