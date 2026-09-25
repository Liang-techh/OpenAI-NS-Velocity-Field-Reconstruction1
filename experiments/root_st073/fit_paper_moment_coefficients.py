"""Fit smooth eta-dependent remote-patch coefficients for both radial moments."""
import json

import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval

from joined_field import JoinedField, ROOT
from paper_moment_bridge import slice_data


PATCH_START = 4.
PATCH_END = 16.


def predicted_moments(row, aE, aU):
    angular = (row['baseline_angular_moment']
               + aE*row['swirl_bump_response'])
    kinetic = (row['kinetic_moment_baseline']
               + aE*row['kinetic_swirl_linear_response']
               + aE*aE*row['kinetic_swirl_quadratic_response']
               + aU*row['meridional_linear_response']
               + aU*aU*row['meridional_quadratic_response'])
    return float(angular), float(kinetic)


def run():
    field = JoinedField()
    training = [slice_data(field, float(eta), tau, 32,
                           PATCH_START, PATCH_END)
                for tau in (.5/64, .032, .128)
                for eta in np.linspace(-.4, .4, 9)]
    eta = np.array([row['eta'] for row in training])
    cE = chebfit(eta/.5,
                [row['swirl_amplitude'] for row in training], 4)
    cU = chebfit(eta/.5,
                [row['meridional_amplitude'] for row in training], 4)
    holdout = []
    for tau in (.012, .064):
        for e in (-.35, -.15, .15, .35):
            row = slice_data(field, e, tau, 48,
                             PATCH_START, PATCH_END)
            aE, aU = float(chebval(e/.5, cE)), float(chebval(e/.5, cU))
            angular, kinetic = predicted_moments(row, aE, aU)
            holdout.append({'tau': tau, 'eta': e,
                            'fitted_swirl_amplitude': aE,
                            'exact_slice_swirl_amplitude': row['swirl_amplitude'],
                            'fitted_meridional_amplitude': aU,
                            'exact_slice_meridional_amplitude': row['meridional_amplitude'],
                            'angular_moment_defect': angular,
                            'kinetic_moment_defect': kinetic})
    report = {'patch_start': PATCH_START, 'patch_end': PATCH_END,
              'eta_chebyshev_scale': .5,
              'swirl_chebyshev_coefficients': cE.tolist(),
              'meridional_chebyshev_coefficients': cU.tolist(),
              'training_times': [.5/64, .032, .128],
              'training_eta': np.linspace(-.4, .4, 9).tolist(),
              'holdout': holdout,
              'holdout_max_angular_defect': max(abs(row['angular_moment_defect'])
                                                 for row in holdout),
              'holdout_max_kinetic_defect': max(abs(row['kinetic_moment_defect'])
                                                 for row in holdout),
              'scope': 'Smooth eta-polynomial approximation of exact slice-level moment coefficients in remote X patch. Time-independent fit across three registered times. Two moment holdouts only; no physical momentum, stress cone, axial localization or global acceptance.',
              'accepted': False}
    (ROOT/'paper_moment_coefficients.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'swirl_coefficients': report['swirl_chebyshev_coefficients'],
                      'meridional_coefficients': report['meridional_chebyshev_coefficients'],
                      'holdout_max_angular_defect': report['holdout_max_angular_defect'],
                      'holdout_max_kinetic_defect': report['holdout_max_kinetic_defect']}),
          flush=True)


if __name__ == '__main__':
    run()
