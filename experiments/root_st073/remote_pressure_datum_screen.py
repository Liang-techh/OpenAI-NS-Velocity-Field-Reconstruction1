"""Restore centrifugal radial pressure and expose the datum tradeoff."""
import json

import numpy as np

from joined_field import ROOT, independent_fd
from moment_matched_joined_field import MomentMatchedJoinedField
from paper_moment_bridge import XI


def centrifugal_check(field, eta, tau):
    point = field.inner.from_similarity(np.array([10.]), eta, tau)
    radius = float(point[0, 0])
    step = 1e-5*radius
    plus, minus = point.copy(), point.copy()
    plus[0, 0] += step
    minus[0, 0] -= step
    def dp(p):
        return float((field.fields(p, tau)[1]-field.base.fields(p, tau)[1])[0])
    radial_derivative = (dp(plus)-dp(minus))/(2*step)
    corrected, _ = field.fields(point, tau)
    baseline, _ = field.base.fields(point, tau)
    centrifugal_change = float((corrected[0, 1]**2-baseline[0, 1]**2)/radius)
    return {'pressure_radial_derivative': radial_derivative,
            'centrifugal_change': centrifugal_change,
            'difference': radial_derivative-centrifugal_change}


def run():
    tau = .0084
    fields = (('unchanged', MomentMatchedJoinedField()),
              ('inner_datum', MomentMatchedJoinedField(
                  pressure_correction=True, pressure_datum='inner')),
              ('outer_datum', MomentMatchedJoinedField(
                  pressure_correction=True, pressure_datum='outer')))
    X = np.array([XI/2, 6., 10., 14., 18.])
    rows = []
    for eta in (0., .2, .35):
        for name, field in fields:
            points = field.inner.from_similarity(X, eta, tau)
            hs = .0005*np.sqrt(field.nu*tau)
            ht = .0001*tau
            residual, divergence = independent_fd(field, points, tau, hs, ht)
            _, pressure = field.fields(points, tau)
            _, baseline_pressure = field.base.fields(points, tau)
            row = {'eta': eta, 'variant': name,
                   'X': X.tolist(), 'residual': residual.tolist(),
                   'residual_norm': np.linalg.norm(residual, axis=1).tolist(),
                   'pressure_increment': (pressure-baseline_pressure).tolist(),
                   'max_abs_divergence': float(np.max(np.abs(divergence)))}
            rows.append(row)
            print(json.dumps({'eta': eta, 'variant': name,
                              'inner_residual': row['residual_norm'][0],
                              'patch_max_residual': max(row['residual_norm'][1:4]),
                              'outer_residual': row['residual_norm'][4]}),
                  flush=True)
    checks = {name: centrifugal_check(field, .2, tau)
              for name, field in fields[1:]}
    report = {'tau': tau, 'rows': rows,
              'centrifugal_checks': checks,
              'scope': 'Radial centrifugal pressure increment integrated from the remote swirl patch with either preserved inner or preserved outer radial datum. Full momentum and divergence finite differences at five radial locations. Pressure correction alone does not solve axial/poloidal residuals or match all five paper moments.',
              'accepted': False}
    (ROOT/'remote_pressure_datum_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
