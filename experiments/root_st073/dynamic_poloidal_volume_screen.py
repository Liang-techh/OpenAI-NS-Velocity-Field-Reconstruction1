"""Whole-support quadrature screen for the staged ST073 correction."""
import json

import numpy as np

from axial_swirl_multitime_screen import load_dense_candidate
from dynamic_poloidal_multitime_screen import load_dynamic_candidate, full_residual
from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from local_poloidal_basis_screen import load_robust_candidate


def increment(modes, amplitudes, points, tau, hs, ht):
    velocity = np.zeros_like(points)
    gradient = np.zeros((len(points), 3, 3))
    part = np.zeros_like(points)
    for mode, amplitude in zip(modes, amplitudes):
        if amplitude:
            du, dJ, dpart = kinematics(mode, points, tau, hs, ht)
            velocity += amplitude*du
            gradient += amplitude*dJ
            part += amplitude*dpart
    return velocity, gradient, part


def add(state, delta):
    return tuple(a+b for a, b in zip(state, delta))


def scaled(delta, factor):
    return tuple(factor*a for a in delta)


def metrics(state, points, weights, support):
    u, J, part = state
    residual = full_residual(u, J, part)
    norms = np.linalg.norm(residual, axis=1)
    worst = int(np.argmax(norms))
    rflat, _, zflat, _ = support
    radius = np.hypot(points[:, 0], points[:, 1])
    axial_collar = np.abs(points[:, 2]) > zflat
    radial_collar = radius > rflat
    zones = {'interior': ~(axial_collar | radial_collar),
             'axial_collar': axial_collar & ~radial_collar,
             'radial_collar': radial_collar & ~axial_collar,
             'corner': axial_collar & radial_collar}
    return {'max': float(norms[worst]),
            'physical_volume_l2': float(np.sqrt(weights@norms**2)),
            'worst_point': points[worst].tolist(),
            'worst_residual': residual[worst].tolist(),
            'zone_max': {name: float(np.max(norms[mask])) if np.any(mask) else None
                         for name, mask in zones.items()},
            'max_speed': float(np.max(np.linalg.norm(u, axis=1)))}


def run():
    base = load_robust_candidate(.1)
    swirl = load_dense_candidate(base)
    candidate = load_dynamic_candidate()
    rows = []
    for order in (6, 8):
        tau = .0084
        points, weights = nodes(base, tau, order)
        hs = .0005*np.sqrt(base.nu*tau)
        ht = .0001*tau
        state0 = kinematics(base, points, tau, hs, ht)
        state1 = add(state0, increment(swirl.modes, swirl.amplitudes,
                                       points, tau, hs, ht))
        poloidal_delta = increment(candidate.modes, candidate.amplitudes,
                                   points, tau, hs, ht)
        state2 = add(state1, poloidal_delta)
        scale_screen = []
        for factor in (0., .1, .25, .5, .75, 1.):
            result = metrics(add(state1, scaled(poloidal_delta, factor)),
                             points, weights, base.support(tau))
            scale_screen.append({'factor': factor, 'max': result['max'],
                                 'physical_volume_l2': result['physical_volume_l2']})
        row = {'tau': tau, 'order': order, 'point_count': len(points),
               'robust_base': metrics(state0, points, weights, base.support(tau)),
               'dense_swirl': metrics(state1, points, weights, base.support(tau)),
               'dynamic_poloidal': metrics(state2, points, weights, base.support(tau)),
               'poloidal_strength_screen': scale_screen}
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = {'rows': rows,
              'scope': 'Whole compact-support cylindrical quadrature at one registered time, Gauss orders 6 and 8. Exact additive finite-difference kinematics and full nonlinear Cartesian momentum. These sparse orders are diagnostics, not certified continuum max or converged physical-volume L2.',
              'accepted': False}
    (ROOT/'compact_potential'/'dynamic_poloidal_volume_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
