"""KKT S capacity with selected inner-cone U values kept unchanged.

Value constraints are only sampled safeguards. Full stress geometry
also depends on gradients and radial primitives, so positive slack
here is not pressure-cone admissibility.
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


def kkt_capacity(U, E, target, X, weights, B,
                 cone_X, start, end, width, degree,
                 pin_gradient=False):
    H = np.sqrt(2*X)*E
    current = moment_vector(U, E, X, weights)
    gram = (B*weights)@B.T
    gradient = B@(weights*U)
    constraints = [B@weights, B@(weights*H)]
    required = [target[0]-current[0], target[2]-current[2]]
    cone_basis = u_basis(cone_X, start, end, width, degree)
    for column in cone_basis.T:
        if np.linalg.norm(column) > 1e-12:
            constraints.append(column)
            required.append(0.)
    if pin_gradient:
        eps = 1e-5
        at_plus = u_basis(cone_X+eps, start, end, width, degree)
        at_minus = u_basis(cone_X-eps, start, end, width, degree)
        for column in ((at_plus-at_minus)/(2*eps)).T:
            if np.linalg.norm(column) > 1e-12:
                constraints.append(column)
                required.append(0.)
    C = np.array(constraints)
    d = np.array(required)
    kkt = np.block([[gram, C.T],
                    [C, np.zeros((len(C), len(C)))]])
    solution = np.linalg.solve(kkt, np.r_[-gradient, d])
    coefficient = solution[:len(B)]
    achieved = moment_vector(U+coefficient@B, E, X, weights)
    return dict(S_min=float(achieved[3]),
                S_slack=float(target[3]-achieved[3]),
                linear_constraint_error=float(np.max(np.abs(
                    C@coefficient-d))),
                cone_delta_U=(coefficient@cone_basis).tolist(),
                kkt_condition=float(np.linalg.cond(kkt)))


def run():
    source = json.loads((ROOT/'delayed_remote_swirl_floor_capacity.json').read_text())
    e_row = next(entry['result'] for entry in source['rows']
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
    cone_X = np.array([1.016, 1.02, 1.03])
    rows = []
    for start in (1.01, 1.02):
        for width in (.005, .01, .02, .04):
            for degree in (11, 19, 31):
                B = u_basis(X, start, 3.5, width, degree)
                for pin_gradient in (False, True):
                    try:
                        result = kkt_capacity(U, E_new, target, X,
                            weights, B, cone_X, start, 3.5,
                            width, degree, pin_gradient)
                    except np.linalg.LinAlgError:
                        result = dict(error='singular KKT')
                    rows.append(dict(start=start, width=width,
                                     degree=degree,
                                     pin_gradient=pin_gradient,
                                     **result))
    report = dict(source='delayed_remote_swirl_floor_capacity.json',
                  eta=eta, cone_X=cone_X.tolist(), rows=rows,
                  scope='Fixed-slice quadratic S capacity with sampled '
                        'cone U values and optionally U_X pinned. '
                        'No full cone stress primitive or PDE residual '
                        'is certified.', accepted=False)
    (ROOT/'delayed_remote_cone_null_capacity.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps([dict(start=start,pin_gradient=pin,
                           best=max((r for r in rows
                                     if r['start']==start and
                                     r['pin_gradient']==pin and
                                     'S_slack' in r),
                                    key=lambda r:r['S_slack']))
                      for start in (1.01,1.02)
                      for pin in (False,True)]), flush=True)


if __name__ == '__main__':
    run()
