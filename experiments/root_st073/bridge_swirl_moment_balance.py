"""Necessary fifth-moment screen and an internal swirl-balance experiment.

The downstream swirl contribution to Cp is nonnegative. Thus a bridge
correction that adds more Cp than the entire baseline downstream tail cannot
be repaired by another swirl correction outside the bridge while keeping the
same exterior and outgoing fifth moment.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from bridge_poloidal_mode import BridgePoloidalMode
from extended_compact_join import load_extended_heated_candidate
from extended_relaxed_cone_screen import cone_point, moment_data, profile
from heat_exterior import tail_moments
from joined_field import independent_fd
from radial_continuation import ROOT


def cp_tail(field, x_out, eta, tau, order=64):
    q = tau/(1-eta**2)
    x_heat = float(field.attachment_radius(tau)**2/(2*field.nu*q))
    nodes, weights = leggauss(order)
    xs = (x_out+x_heat)/2+(x_heat-x_out)*nodes/2
    _, E = profile(field, xs, eta, tau)
    near = float((x_heat-x_out)/2*weights@(E**2/(2*xs)))
    # tail_moments already uses the normalized E=c X^(-A) H profile.
    far = float(tail_moments(
        x_heat, eta, c=field.heat.amplitude,
        h=field.heat.h)[0])
    return near+far


def balanced_symmetric_amplitudes(base, outer_amplitude, eta, tau):
    radial = base.compact.joined
    xi = radial.join_X
    xo = radial.outer_ratio**2*xi
    nodes, weights = leggauss(64)
    xs = (xi+xo)/2+(xo-xi)*nodes/2
    weights = (xo-xi)*weights/2
    _, e0 = profile(base, xs, eta, tau)
    y = (np.sqrt(xs/xi)-1)/(radial.outer_ratio-1)
    symmetric = 64*y**3*(1-y)**3
    outer = symmetric/.421875*(y/.75)**8
    fixed = e0+outer_amplitude*outer
    weight = weights/(2*xs)
    quadratic = np.array([
        weight@(symmetric**2),
        2*weight@(fixed*symmetric),
        weight@(fixed**2-e0**2),
    ])
    roots = np.roots(quadratic)
    return [float(x.real) for x in roots if abs(x.imag) < 1e-9], quadratic


def make_field(symmetric, outer=1., minimum=-2.75, moment=4.5):
    base = load_extended_heated_candidate(
        swirl_bubble_amplitude=symmetric,
        outer_swirl_bubble_amplitude=outer)
    return BridgePoloidalMode(
        BridgePoloidalMode(base, minimum, shape='minimum_curvature'),
        moment, shape='moment')


def run():
    tau = .5*2**(-5.5)
    base = load_extended_heated_candidate()
    prior = make_field(0.)
    x_out = base.compact.joined.outer_ratio**2*base.compact.joined.join_X
    rows = []
    for eta in (.2, .3, .65):
        m0, _, _ = moment_data(base, x_out, eta, tau, 32)
        m1, _, _ = moment_data(prior, x_out, eta, tau, 32)
        tail = cp_tail(base, x_out, eta, tau)
        roots, coefficients = balanced_symmetric_amplitudes(
            base, 1., eta, tau)
        rows.append(dict(eta=eta, delta_moments=(m1-m0).tolist(),
                         downstream_cp_capacity=tail,
                         unavoidable_cp_defect=max(0., m1[4]-m0[4]-tail),
                         symmetric_cp_balance_roots=roots,
                         cp_quadratic=coefficients.tolist()))
    # One constant amplitude brackets the two endpoint slice roots.
    roots = [min(row['symmetric_cp_balance_roots'], key=abs)
             for row in rows]
    symmetric = (roots[0]+roots[-1])/2
    balanced = make_field(symmetric)
    X, eta = np.meshgrid([.5859375, .75, 1., 1.25], [.2, .3],
                         indexing='ij')
    X, eta = X.ravel(), eta.ravel()
    points = base.compact.joined.inner.from_similarity(X, eta, tau)
    residual, divergence = independent_fd(
        balanced, points, tau, .001*np.sqrt(base.nu*tau), .00025*tau)
    cones = [cone_point(balanced, 1., e, tau, order=16)
             for e in (.2, .3)]
    for row in rows:
        m, _, _ = moment_data(balanced, x_out, row['eta'], tau, 32)
        m0, _, _ = moment_data(base, x_out, row['eta'], tau, 32)
        row['balanced_delta_moments'] = (m-m0).tolist()
    report = dict(tau=tau, x_out=x_out, rows=rows,
                  balanced_symmetric_amplitude=symmetric,
                  balanced_momentum_max=float(np.max(np.linalg.norm(
                      residual, axis=1))),
                  balanced_divergence_max=float(np.max(np.abs(divergence))),
                  balanced_cones=[dict(eta=c['eta'], Pc=c.get('Pc'),
                                       margin=(c['upper']-c['vs']
                                               if c.get('upper') is not None
                                               else None),
                                       relaxed_pass=c.get('relaxed_pass'))
                                  for c in cones],
                  scope='Necessary downstream fifth-moment screen and one internal swirl-balance candidate. Finite sample only; no global momentum or cone acceptance.',
                  accepted=False)
    (ROOT/'bridge_swirl_moment_balance.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run()
