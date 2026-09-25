"""Screen how widening U support changes optimistic J/S capacity.

Uses the fixed-slice positive E profiles. U is allowed to vary as an
arbitrary L2 function on each proposed interval while preserving M and
J; smoothness, solenoidal lifting, and full momentum are omitted.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import current_mean
from delayed_remote_positive_e import quadrature
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT


SUPPORTS = ((1.04, 2.98), (1.04, 3.0),
            (1.04, 3.25), (1.04, 3.5),
            (1.03, 3.5), (1.02, 3.5),
            (1.0, 3.5), (.5, 3.5))


def capacity(U, E, target, X, weights, support):
    inside = (X >= support[0]) & (X <= support[1])
    outside = ~inside
    H = np.sqrt(2*X)*E
    wi, hi = weights[inside], H[inside]
    gram = np.array([[np.sum(wi), wi@hi],
                     [wi@hi, wi@(hi**2)]])
    required = np.array([target[0]-weights[outside]@U[outside],
                         target[2]-weights[outside]@(
                             U[outside]*H[outside])])
    minimum = (float(weights[outside]@(U[outside]**2))
               +float(required@np.linalg.solve(gram, required))
               -float(weights@(E**2))/2)
    return dict(S_min=minimum, S_slack=float(target[3]-minimum),
                gram_condition=float(np.linalg.cond(gram)))


def run():
    source = json.loads((ROOT/'delayed_remote_positive_e.json').read_text())
    radial_intervals = [tuple(interval)
                        for interval in source['radial_intervals']]
    source_eta_interval = tuple(source['mode_eta_interval'])
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target_field = make_field(16, 2.)
    X, weights = quadrature()
    radial = np.array([bump(X, *interval)[0]
                       for interval in radial_intervals])
    tau = source['tau']
    rows = []
    for prior in source['rows']:
        eta = prior['eta']
        U, E = profile(mean, X, eta, tau)
        target_U, target_E = profile(target_field, X, eta, tau)
        target = moment_vector(target_U, target_E, X, weights)
        axial = bump(np.array([eta]), *source_eta_interval)[0][0]
        corrected_E = E+axial*np.asarray(prior['coefficients'])@radial
        rows.append(dict(eta=eta,
                         supports=[dict(interval=support,
                                        **capacity(U, corrected_E, target,
                                                   X, weights, support))
                                   for support in SUPPORTS]))
    report = dict(source='delayed_remote_positive_e.json', rows=rows,
                  scope='Optimistic fixed-E, fixed-slice lower bounds '
                        'for arbitrary U on proposed intervals. '
                        'S slack is necessary, not sufficient, for '
                        'smooth exact-solenoidal J/S repair or NS '
                        'momentum acceptance.', accepted=False)
    (ROOT/'delayed_remote_support_capacity.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps([(row['eta'],
                       [round(s['S_slack'], 7) for s in row['supports']])
                      for row in rows]), flush=True)


if __name__ == '__main__':
    run()
