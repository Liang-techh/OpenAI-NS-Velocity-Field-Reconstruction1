"""Test smooth finite-dimensional U capacity after positive E repair.

This solves fixed-slice M/J constraints by a weighted KKT projection
and measures remaining S slack. Only a nonnegative slack can be spent
on a null direction to achieve S. No physical streamfunction is built.
"""

import json

import numpy as np
from numpy.polynomial.legendre import legval

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import current_mean
from delayed_remote_positive_e import quadrature
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import septic_step


def u_basis(X, start, end, width, degree):
    up = septic_step((X-start)/width)[0]
    down = septic_step((X-(end-width))/width)[0]
    window = up*(1-down)
    argument = np.clip(2*(X-start)/(end-start)-1, -1, 1)
    return np.array([window*legval(argument, [0.]*order+[1.])
                     for order in range(degree+1)])


def capacity(U, E, target, X, weights, B):
    # Weighted QR avoids squaring the polynomial-basis condition number.
    orth, factor = np.linalg.qr((B*np.sqrt(weights)).T)
    Q = orth.T/np.sqrt(weights)
    H = np.sqrt(2*X)*E
    constraints = np.vstack((Q@weights, Q@(weights*H)))
    current = moment_vector(U, E, X, weights)
    required = np.array([target[0]-current[0],
                         target[2]-current[2]])
    unconstrained = -(Q@(weights*U))
    dual = constraints@constraints.T
    coefficients = (unconstrained+constraints.T@np.linalg.solve(
        dual, required-constraints@unconstrained))
    corrected_U = U+coefficients@Q
    achieved = moment_vector(corrected_U, E, X, weights)
    return dict(S_min=float(achieved[3]),
                target_minus_S_min=float(target[3]-achieved[3]),
                M_J_defect=(achieved-target)[[0, 2]].tolist(),
                polynomial_condition=float(np.linalg.cond(factor)),
                constraint_condition=float(np.linalg.cond(dual)),
                max_abs_U=float(np.max(np.abs(corrected_U))))


def run():
    source = json.loads((ROOT/'delayed_remote_swirl_floor_capacity.json').read_text())
    row = next(entry['result'] for entry in source['rows']
               if entry['relative_E_floor'] == .05)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target_field = make_field(16, 2.)
    X, weights = quadrature()
    eta = .3
    tau = .5*2**(-5.5)
    U, E = profile(mean, X, eta, tau)
    U0, E0 = profile(target_field, X, eta, tau)
    target = moment_vector(U0, E0, X, weights)
    first = json.loads((ROOT/'delayed_remote_positive_e.json').read_text())
    radial = np.array([bump(X, *interval)[0]
                       for interval in first['radial_intervals']])
    axial = bump(np.array([eta]), *first['mode_eta_interval'])[0][0]
    E_new = E+axial*np.asarray(row['coefficients'])@radial
    results = []
    for width in (.002, .005, .01, .02, .04):
        for degree in (11, 19, 31, 47, 63):
            B = u_basis(X, 1.03, 3.5, width, degree)
            results.append(dict(width=width, degree=degree,
                                **capacity(U, E_new, target, X, weights, B)))
    report = dict(source='delayed_remote_swirl_floor_capacity.json',
                  eta=eta, relative_E_floor=.05,
                  E_min_relative_to_target=float(np.min(E_new/E0)),
                  results=results,
                  scope='Fixed-slice tapered-polynomial U basis '
                        'capacity under exact linear M/J. A positive '
                        'S gap is needed before quadratic completion. '
                        'No eta/time lift or complete momentum check.',
                  accepted=False)
    (ROOT/'delayed_remote_u_basis_capacity.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps([(r['width'], r['degree'],
                       r['target_minus_S_min'],
                       r['polynomial_condition']) for r in results]),
          flush=True)


if __name__ == '__main__':
    run()
