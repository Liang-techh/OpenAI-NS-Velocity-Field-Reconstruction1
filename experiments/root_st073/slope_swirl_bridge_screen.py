"""Probe a value-zero, shear-changing swirl mode in the ST073 bridge.

The mode vanishes at X=1 and at both radial bridge boundaries. It follows
the paper's distinction between local shear control and cumulative moments,
but is only a low-frequency numerical diagnostic, not the paper's shear loop.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from bridge_poloidal_mode import BridgePoloidalMode
from bridge_swirl_moment_balance import cp_tail
from extended_compact_join import load_extended_heated_candidate
from extended_relaxed_cone_screen import cone_point, profile
from joined_field import independent_fd
from joint_collar_fit import kinematics
from radial_continuation import ROOT


SYMMETRIC = -.27311627313724157


def make_field(slope):
    field = load_extended_heated_candidate(
        swirl_bubble_amplitude=SYMMETRIC,
        outer_swirl_bubble_amplitude=1.,
        slope_swirl_bubble_amplitude=slope)
    field = BridgePoloidalMode(field, -2.75, shape='minimum_curvature')
    return BridgePoloidalMode(field, 4.5, shape='moment')


def complete_response(points, tau):
    fields = [make_field(0.), make_field(1.)]
    data = [kinematics(field, points, tau,
                       .001*np.sqrt(field.nu*tau), .00025*tau)
            for field in fields]
    u0, j0, p0 = data[0]
    du, dj, dp = [a-b for a, b in zip(data[1], data[0])]

    def residual(amplitude):
        u = u0+amplitude*du
        jac = j0+amplitude*dj
        part = p0+amplitude*dp
        return part+np.einsum('nij,nj->ni', jac, u)

    return residual


def cp_response(etas, tau):
    base = load_extended_heated_candidate()
    radial = base.compact.joined
    xi, xo = radial.join_X, radial.join_X*radial.outer_ratio**2
    nodes, weights = leggauss(64)
    xs = (xi+xo)/2+(xo-xi)*nodes/2
    weights = (xo-xi)*weights/(4*xs)
    y = (np.sqrt(xs/xi)-1)/(radial.outer_ratio-1)
    y0 = (np.sqrt(1/xi)-1)/(radial.outer_ratio-1)
    symmetric = 64*y**3*(1-y)**3
    outer = symmetric/.421875*(y/.75)**8
    slope_shape = symmetric*(y-y0)
    e0 = [profile(base, xs, eta, tau)[1] for eta in etas]
    capacity = [cp_tail(base, xo, eta, tau) for eta in etas]

    def delta(amplitude):
        dE = SYMMETRIC*symmetric+outer+amplitude*slope_shape
        return [float(weights@((E+dE)**2-E**2)) for E in e0]

    return delta, capacity


def run():
    tau = .5*2**(-5.5)
    field = make_field(0.)
    x, eta = np.meshgrid([.5859375, .75, 1., 1.25], [.2, .3],
                         indexing='ij')
    points = field.compact.joined.inner.from_similarity(
        x.ravel(), eta.ravel(), tau)
    residual = complete_response(points, tau)
    etas = [.2, .3, .65]
    delta_cp, capacity = cp_response(etas, tau)
    rows = []
    for amplitude in (-4., -2., -1., 0., 1., 2., 4.):
        candidate = make_field(amplitude)
        cones = [cone_point(candidate, 1., e, tau, order=8)
                 for e in (.2, .3)]
        R = residual(amplitude)
        row = dict(amplitude=amplitude,
                   cp_delta=delta_cp(amplitude),
                   cp_capacity_pass=all(a <= b for a, b in
                                        zip(delta_cp(amplitude), capacity)),
                   sampled_full_max=float(np.max(np.linalg.norm(R, axis=1))),
                   sampled_angular_max=float(np.max(np.abs(R[:, 1]))),
                   cones=[dict(eta=e['eta'], Pc=e.get('Pc'),
                               margin=(e['upper']-e['vs']
                                       if e.get('upper') is not None else None),
                               relaxed_pass=e.get('relaxed_pass'))
                          for e in cones])
        rows.append(row)
        print(json.dumps(row), flush=True)
    feasible = [row for row in rows if row['cp_capacity_pass']
                and all(c['relaxed_pass'] for c in row['cones'])]
    selected = min(feasible, key=lambda row: row['sampled_full_max'])
    candidate = make_field(selected['amplitude'])
    direct, divergence = independent_fd(
        candidate, points, tau, .001*np.sqrt(candidate.nu*tau),
        .00025*tau)
    selected['direct_full_max'] = float(np.max(np.linalg.norm(direct, axis=1)))
    selected['surrogate_fd_discrepancy'] = float(np.max(np.linalg.norm(
        direct-residual(selected['amplitude']), axis=1)))
    selected['divergence_max'] = float(np.max(np.abs(divergence)))
    report = dict(tau=tau, cp_capacity=capacity, rows=rows,
                  selected_amplitude=selected['amplitude'],
                  scope='Seven-point slope-mode screen at eight physical momentum points and two snapshot cone points. Cp capacity is a one-sided necessary screen for downstream deletion, not a fifth-moment match. No continuous cone or global PDE acceptance.',
                  accepted=False)
    (ROOT/'slope_swirl_bridge_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
