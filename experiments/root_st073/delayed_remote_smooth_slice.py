"""Complete the smoother eta=.3 fixed-slice M/J/S algebra.

The U correction begins at X=1.01 and therefore perturbs the inner
cone. This reports its exact slice moments and cone-node U changes;
it is not a physical eta/time lift or a pressure-cone admission.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import current_mean
from delayed_remote_positive_e import quadrature
from delayed_remote_u_basis_capacity import u_basis
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT


def run():
    floor_source = json.loads((ROOT/'delayed_remote_swirl_floor_capacity.json').read_text())
    e_row = next(entry['result'] for entry in floor_source['rows']
                 if entry['relative_E_floor'] == .05)
    e_source = json.loads((ROOT/'delayed_remote_positive_e.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target_field = make_field(16, 2.)
    X, weights = quadrature()
    eta = .3
    tau = .5*2**(-5.5)
    U, E = profile(mean, X, eta, tau)
    U0, E0 = profile(target_field, X, eta, tau)
    target = moment_vector(U0, E0, X, weights)
    radial = np.array([bump(X, *interval)[0]
                       for interval in e_source['radial_intervals']])
    axial = bump(np.array([eta]), *e_source['mode_eta_interval'])[0][0]
    E_new = E+axial*np.asarray(e_row['coefficients'])@radial
    start, end, width, degree = 1.01, 3.5, .02, 11
    B = u_basis(X, start, end, width, degree)
    orth, factor = np.linalg.qr((B*np.sqrt(weights)).T)
    Q = orth.T/np.sqrt(weights)
    H = np.sqrt(2*X)*E_new
    constraints = np.vstack((Q@weights, Q@(weights*H)))
    current = moment_vector(U, E_new, X, weights)
    required = np.array([target[0]-current[0],
                         target[2]-current[2]])
    unconstrained = -(Q@(weights*U))
    dual = constraints@constraints.T
    minimum = (unconstrained+constraints.T@np.linalg.solve(
        dual, required-constraints@unconstrained))
    at_minimum = moment_vector(U+minimum@Q, E_new, X, weights)
    gap = float(target[3]-at_minimum[3])
    if gap < 0:
        raise RuntimeError(f'finite basis S capacity negative: {gap}')
    options = []
    for index in range(len(Q)):
        trial = np.eye(len(Q))[index]
        null = trial-constraints.T@np.linalg.solve(
            dual, constraints@trial)
        norm = float(null@null)
        if norm < 1e-12:
            continue
        for sign in (-1., 1.):
            coeff_q = minimum+sign*np.sqrt(gap/norm)*null
            new_U = U+coeff_q@Q
            options.append((float(np.max(np.abs(new_U))), coeff_q,
                            new_U, index, sign))
    peak, coeff_q, U_new, index, sign = min(options,
        key=lambda row: row[0])
    coeff_b = np.linalg.solve(factor, coeff_q)
    achieved = moment_vector(U_new, E_new, X, weights)
    cone_X = np.array([1.008, 1.016, 1.02, 1.03])
    cone_delta_U = coeff_b@u_basis(cone_X, start, end, width, degree)
    report = dict(source='delayed_remote_swirl_floor_capacity.json',
                  eta=eta, start=start, end=end, width=width,
                  degree=degree, coefficients=coeff_b.tolist(),
                  selected_null_index=index, selected_null_sign=sign,
                  S_slack_before_completion=gap,
                  moment_defect=(achieved-target).tolist(),
                  max_abs_corrected_U=peak,
                  cone_X=cone_X.tolist(),
                  cone_delta_U=cone_delta_U.tolist(),
                  max_abs_cone_delta_U=float(np.max(np.abs(cone_delta_U))),
                  min_relative_E=float(np.min(E_new/E0)),
                  scope='Exact fixed-eta five-moment algebra with a '
                        'smooth radial U basis entering the current '
                        'stress cone. No eta/time-dependent physical '
                        'lift, cone re-audit, momentum or volume L2 '
                        'acceptance.', accepted=False)
    (ROOT/'delayed_remote_smooth_slice.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(moment_defect=report['moment_defect'],
                          gap=gap, cone_delta_U=report['cone_delta_U'],
                          max_abs_coeff=float(np.max(np.abs(coeff_b))))),
          flush=True)


if __name__ == '__main__':
    run()
