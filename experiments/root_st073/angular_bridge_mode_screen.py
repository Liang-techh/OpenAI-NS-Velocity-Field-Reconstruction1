"""Four-mode screen of the pressure-invariant bridge momentum component.

This searches the existing two swirl and two poloidal velocity modes under
a necessary fifth-moment capacity constraint. It is a finite-sample search,
not a global residual or continuous-cone construction.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import differential_evolution

from bridge_swirl_moment_balance import cp_tail, make_field
from extended_compact_join import load_extended_heated_candidate
from extended_relaxed_cone_screen import cone_point, profile
from joined_field import independent_fd
from joint_collar_fit import kinematics
from radial_continuation import ROOT


def response_model(points, tau):
    fields = [make_field(0., 0., 0., 0.),
              make_field(1., 0., 0., 0.),
              make_field(0., 1., 0., 0.),
              make_field(0., 0., 1., 0.),
              make_field(0., 0., 0., 1.)]
    hs = .001*np.sqrt(fields[0].nu*tau)
    data = [kinematics(field, points, tau, hs, .00025*tau)
            for field in fields]
    u0, j0, p0 = data[0]
    du = np.stack([row[0]-u0 for row in data[1:]])
    dj = np.stack([row[1]-j0 for row in data[1:]])
    dp = np.stack([row[2]-p0 for row in data[1:]])

    def residual(c):
        c = np.asarray(c, float)
        u = u0+np.einsum('k,knj->nj', c, du)
        jac = j0+np.einsum('k,knij->nij', c, dj)
        part = p0+np.einsum('k,knj->nj', c, dp)
        return part+np.einsum('nij,nj->ni', jac, u)

    return residual


def cp_response(base, etas, tau):
    radial = base.compact.joined
    xi, xo = radial.join_X, radial.join_X*radial.outer_ratio**2
    nodes, weights = leggauss(64)
    xs = (xi+xo)/2+(xo-xi)*nodes/2
    weights = (xo-xi)*weights/(4*xs)
    y = (np.sqrt(xs/xi)-1)/(radial.outer_ratio-1)
    symmetric = 64*y**3*(1-y)**3
    outer = symmetric/.421875*(y/.75)**8
    profiles = [profile(base, xs, eta, tau)[1] for eta in etas]
    capacity = np.array([cp_tail(base, xo, eta, tau) for eta in etas])

    def delta(c):
        correction = c[0]*symmetric+c[1]*outer
        return np.array([weights@((E+correction)**2-E**2)
                         for E in profiles])

    return delta, capacity


def run():
    tau = .5*2**(-5.5)
    base = load_extended_heated_candidate()
    train_x, train_eta = np.meshgrid(
        [.5859375, .75, 1., 1.25], [.2, .3], indexing='ij')
    train = base.compact.joined.inner.from_similarity(
        train_x.ravel(), train_eta.ravel(), tau)
    residual = response_model(train, tau)
    delta_cp, capacity = cp_response(base, [.2, .3, .65], tau)
    old = np.array([-.27311627313724157, 1., -2.75, 4.5])

    def objective(c):
        excess = np.maximum(delta_cp(c)-capacity, 0.)
        return (float(np.max(np.abs(residual(c)[:, 1])))
                +1e5*float(np.max(excess)))

    opt = differential_evolution(
        objective, bounds=[(-1., .1), (.75, 1.5), (-4., -1.), (2., 7.)],
        maxiter=60, popsize=8, seed=7304, polish=True)
    rows = []
    for label, c in [('prior', old), ('optimized', opt.x)]:
        field = make_field(*c)
        train_R, _ = independent_fd(
            field, train, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
        hold_x, hold_eta = np.meshgrid(
            [.68, .9, 1.1, 1.34], [.25, .35], indexing='ij')
        hold = field.compact.joined.inner.from_similarity(
            hold_x.ravel(), hold_eta.ravel(), tau)
        hold_R, _ = independent_fd(
            field, hold, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
        cones = [cone_point(field, 1., eta, tau, order=12)
                 for eta in (.2, .3)]
        rows.append(dict(label=label, coefficients=c.tolist(),
                         delta_cp=delta_cp(c).tolist(),
                         train_angular_max=float(np.max(np.abs(train_R[:, 1]))),
                         holdout_angular_max=float(np.max(np.abs(hold_R[:, 1]))),
                         train_full_max=float(np.max(np.linalg.norm(train_R, axis=1))),
                         holdout_full_max=float(np.max(np.linalg.norm(hold_R, axis=1))),
                         surrogate_fd_discrepancy=float(np.max(np.linalg.norm(
                             train_R-residual(c), axis=1))),
                         cones=[dict(eta=p['eta'], Pc=p.get('Pc'),
                                     margin=(p['upper']-p['vs']
                                             if p.get('upper') is not None
                                             else None),
                                     relaxed_pass=p.get('relaxed_pass'))
                                for p in cones]))
        print(json.dumps(rows[-1]), flush=True)
    report = dict(tau=tau, downstream_cp_capacity=capacity.tolist(),
                  optimizer_message=str(opt.message), rows=rows,
                  scope='Four-mode angular-momentum screen at eight train and eight holdout physical points, with necessary Cp capacity and two snapshot cone checks. No continuous/global acceptance.',
                  accepted=False)
    (ROOT/'angular_bridge_mode_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
