"""Screen whether current momentum residual trends toward a bounded force.

The paper permits a smooth compact force across the singular time. A force
defined tautologically as the candidate residual on a truncated slab does
not establish that extension. This samples matched relative collar points.
"""
import json
import sys

import numpy as np

from joined_field import ROOT
from joint_collar_fit import kinematics
from local_poloidal_basis_screen import compact_base, load_robust_candidate
from radial_peak_cone import operator


def run():
    field = load_robust_candidate(.1)
    extrapolate = '--extrapolate' in sys.argv
    if extrapolate:
        compact = compact_base(field)
        compact.experimental_time_extension = True
        compact.base.experimental_time_extension = True
        times = [.000525, .00105, .0021, .0042, .0084,
                 .012, .024, .048]
    else:
        times = [.0084, .01, .012, .016, .024, .032, .064, .128]
    radial_fractions = [.3, .5]
    collar_fractions = [.35, .65]
    rows = []
    for tau in times:
        _, radius, zflat, zsupport = field.support(tau)
        points = np.array([
            [rf*radius, 0., sign*(zflat+s*(zsupport-zflat))]
            for sign in (-1, 1)
            for rf in radial_fractions
            for s in collar_fractions])
        if extrapolate:
            hs = .0005*np.sqrt(field.nu*tau)
            ht = .0001*tau
            velocity, gradient, part = kinematics(
                field, points, tau, hs, ht, time_min=0.)
            residual = part+np.einsum('nij,nj->ni', gradient, velocity)
        else:
            velocity, _, residual = operator(field, points, tau)
        norms = np.linalg.norm(residual, axis=1)
        velocities = np.linalg.norm(velocity, axis=1)
        worst = int(np.argmax(norms))
        row = {'tau': tau, 'radial_support': float(radius),
               'zflat': float(zflat), 'axial_support': float(zsupport),
               'sample_max_residual': float(norms[worst]),
               'sample_max_speed': float(np.max(velocities)),
               'component_max_abs': np.max(np.abs(residual), axis=0).tolist(),
               'worst_point': points[worst].tolist(),
               'worst_residual': residual[worst].tolist(),
               'point_norms': norms.tolist()}
        rows.append(row)
        print(json.dumps({k: row[k] for k in
                          ('tau', 'sample_max_residual', 'sample_max_speed',
                           'worst_point')}), flush=True)
    small = rows[:5]
    exponent = float(np.polyfit(np.log([row['tau'] for row in small]),
                                np.log([row['sample_max_residual'] for row in small]), 1)[0])
    report = {'times': times, 'experimental_extrapolation': extrapolate,
              'radial_fractions': radial_fractions,
              'collar_fractions': collar_fractions, 'rows': rows,
              'small_tau_log_log_slope': exponent,
              'scope': ('Eight relative collar points at each of eight positive backward times. Extrapolation explicitly bypasses the ST073 registered finite-slab guards and uses central time differences; formulas and fitted coefficients are unvalidated there. A fitted trend neither proves divergence nor smooth extendibility of a whole-space force at tau=0.' if extrapolate else 'Eight relative collar points at each of eight registered positive backward times only. A fitted trend neither proves divergence nor smooth extendibility of a whole-space force at tau=0. Current candidate is registered only for tau>=0.5/64.'),
              'accepted': False}
    path = ROOT/'compact_potential'/(
        'forcing_extension_extrapolation.json' if extrapolate else
        'forcing_extension_screen.json')
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'small_tau_log_log_slope': exponent}), flush=True)


if __name__ == '__main__':
    run()
