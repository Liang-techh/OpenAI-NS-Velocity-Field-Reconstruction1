"""Fit localized exact-solenoidal collar streamfunctions to momentum peaks."""
import json

import numpy as np
from numpy.polynomial.legendre import Legendre
from scipy.optimize import least_squares, minimize

from compact_potential import CompactPotentialField
from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from radial_pressure_volume_constrained import load_candidate as load_base


class LocalPoloidalMode:
    """Axisymmetric u=curl of a compact azimuthal vector potential."""

    def __init__(self, base, radial_degree, axial_degree, parity,
                 reference_tau=.5/64, temporal_power=2.):
        self.base = base
        self.nu = base.nu
        self.radial_degree = radial_degree
        self.axial_degree = axial_degree
        self.parity = parity
        self.reference_tau = reference_tau
        self.temporal_power = temporal_power

    def fields(self, points, tau):
        points = np.asarray(points)
        velocity = np.zeros_like(points)
        _, radius, zflat, zsupport = self.base.support(tau)
        r = np.hypot(points[:, 0], points[:, 1])
        x = (r/radius)**2
        s = (np.abs(points[:, 2])-zflat)/(zsupport-zflat)
        inside = (x < 1) & (s > 0) & (s < 1)
        if not np.any(inside):
            return velocity, np.zeros(len(points))
        xx, ss = x[inside], s[inside]
        radial = Legendre.basis(self.radial_degree)
        axial = Legendre.basis(self.axial_degree)
        rv = radial(2*xx-1)
        f = (1-xx)**5*rv
        fx = -5*(1-xx)**4*rv+2*(1-xx)**5*radial.deriv()(2*xx-1)
        av = axial(2*ss-1)
        bubble = 1024*ss**5*(1-ss)**5
        bubble_s = 5120*ss**4*(1-ss)**4*(1-2*ss)
        g = bubble*av
        gs = bubble_s*av+2*bubble*axial.deriv()(2*ss-1)
        sign = np.sign(points[inside, 2])
        radial_parity = sign if self.parity == 'even' else np.ones(len(sign))
        axial_parity = np.ones(len(sign)) if self.parity == 'even' else sign
        factor = (self.reference_tau/tau)**self.temporal_power
        velocity[inside, 0] = (-15*points[inside, 0]*f*gs*radial_parity
                               /(zsupport-zflat)*factor)
        velocity[inside, 1] = (-15*points[inside, 1]*f*gs*radial_parity
                               /(zsupport-zflat)*factor)
        velocity[inside, 2] = (30*(f+xx*fx)*g*axial_parity*factor)
        return velocity, np.zeros(len(points))


class LocalPoloidalCandidate:
    def __init__(self, base, modes, amplitudes):
        self.base = base
        self.nu = base.nu
        self.modes = modes
        self.amplitudes = np.asarray(amplitudes)

    def support(self, tau):
        return self.base.support(tau)

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        for amplitude, mode in zip(self.amplitudes, self.modes):
            if amplitude:
                velocity += amplitude*mode.fields(points, tau)[0]
        return velocity, pressure


def compact_base(field):
    node = field
    while not isinstance(node, CompactPotentialField):
        node = node.base
    return node


def load_robust_candidate(strength=1.):
    report = json.loads((ROOT/'compact_potential'/'local_poloidal_basis_screen.json').read_text())
    if not report['robust_optimizer_success']:
        raise ValueError('No robust poloidal fit is registered')
    base = load_base()
    compact = compact_base(base)
    modes = [LocalPoloidalMode(compact, item['radial_degree'],
                               item['axial_degree'], item['parity'])
             for item in report['modes']]
    amplitudes = strength*np.asarray(report['scales'])*np.asarray(
        report['robust_relative_amplitudes'])
    return LocalPoloidalCandidate(base, modes, amplitudes)


def run():
    field = load_base()
    modes = [LocalPoloidalMode(compact_base(field), i, j, parity)
             for parity in ('even', 'odd') for i in (0, 1)
             for j in (0, 1)]
    mode_ids = [{'parity': mode.parity,
                 'radial_degree': mode.radial_degree,
                 'axial_degree': mode.axial_degree} for mode in modes]
    tau = .0084
    blocks = []
    max_base_speed = 0.
    max_mode_speed = np.zeros(len(modes))
    for order in (6, 8, 10):
        points, weights = nodes(field, tau, order)
        hs = .0005*np.sqrt(field.nu*tau)
        ht = .0001*tau
        u0, J0, part0 = kinematics(field, points, tau, hs, ht)
        pieces = [kinematics(mode, points, tau, hs, ht) for mode in modes]
        du = np.stack([item[0] for item in pieces], axis=-1)
        dJ = np.stack([item[1] for item in pieces], axis=-1)
        dpart = np.stack([item[2] for item in pieces], axis=-1)
        R0 = part0+np.einsum('nia,na->ni', J0, u0)
        baseline_norm = np.linalg.norm(R0, axis=1)
        baseline_l2 = float(np.sqrt(weights@baseline_norm**2))
        max_base_speed = max(max_base_speed,
                             float(np.max(np.linalg.norm(u0, axis=1))))
        max_mode_speed = np.maximum(max_mode_speed,
            np.max(np.linalg.norm(du, axis=1), axis=0))
        blocks.append({'order': order, 'points': points, 'weights': weights,
                       'u0': u0, 'J0': J0, 'part0': part0,
                       'du': du, 'dJ': dJ, 'dpart': dpart,
                       'baseline_l2': baseline_l2,
                       'baseline_max': float(np.max(baseline_norm))})
        print(json.dumps({'loaded_order': order,
                          'baseline_l2': baseline_l2,
                          'baseline_max': blocks[-1]['baseline_max']}),
              flush=True)
    scales = .25*max_base_speed/np.maximum(max_mode_speed, 1e-12)
    def evaluate(block, x):
        a = scales*x
        u = block['u0']+np.einsum('nik,k->ni', block['du'], a)
        J = block['J0']+np.einsum('niak,k->nia', block['dJ'], a)
        part = block['part0']+np.einsum('nik,k->ni', block['dpart'], a)
        return part+np.einsum('nia,na->ni', J, u)
    def objective(x):
        parts = []
        for block in blocks:
            residual = evaluate(block, x)
            parts.append((np.sqrt(block['weights'])[:, None]*residual
                          /block['baseline_l2']).ravel())
        parts.append(.01*x)
        return np.concatenate(parts)
    fit = least_squares(objective, np.zeros(len(modes)),
                        bounds=(-2*np.ones(len(modes)),
                                 2*np.ones(len(modes))),
                        max_nfev=500, xtol=1e-10, ftol=1e-10, gtol=1e-10)
    def robust_objective(x):
        error = objective(x)
        return .5*(error@error)
    constraints = []
    for block in blocks:
        constraints.append({
            'type': 'ineq',
            'fun': lambda x, item=block:
                1.-(item['weights']@np.sum(evaluate(item, x)**2, axis=1))
                /item['baseline_l2']**2})
        constraints.append({
            'type': 'ineq',
            'fun': lambda x, item=block:
                1.-np.sum(evaluate(item, x)**2, axis=1)
                /item['baseline_max']**2})
    robust_fit = minimize(robust_objective, np.zeros(len(modes)),
                          method='SLSQP', constraints=constraints,
                          bounds=[(-2, 2)]*len(modes),
                          options={'maxiter': 1000, 'ftol': 1e-12})
    rows = []
    for strength in (0., .25, .5, .75, 1.):
        x = strength*fit.x
        metrics = []
        for block in blocks:
            residual = evaluate(block, x)
            norm = np.linalg.norm(residual, axis=1)
            metrics.append({'order': block['order'],
                            'physical_volume_l2': float(np.sqrt(block['weights']@norm**2)),
                            'max': float(np.max(norm))})
        rows.append({'strength': strength,
                     'amplitudes': (scales*x).tolist(),
                     'metrics': metrics,
                     'worst_l2_ratio': max(metric['physical_volume_l2']/block['baseline_l2']
                                           for metric, block in zip(metrics, blocks)),
                     'worst_max_ratio': max(metric['max']/block['baseline_max']
                                            for metric, block in zip(metrics, blocks))})
    report = {'tau': tau, 'modes': mode_ids, 'scales': scales.tolist(),
              'optimizer_success': bool(fit.success),
              'optimizer_message': fit.message,
              'relative_amplitudes': fit.x.tolist(),
              'robust_optimizer_success': bool(robust_fit.success),
              'robust_optimizer_message': robust_fit.message,
              'robust_relative_amplitudes': robust_fit.x.tolist(),
              'rows': rows,
              'scope': 'Eight compact exact-solenoidal poloidal streamfunctions, full nonlinear physical momentum at one time on Gauss6/8/10. No cone or time holdout validation yet; numerical derivatives are finite-difference approximations.',
              'accepted': False}
    path = ROOT/'compact_potential'/'local_poloidal_basis_screen.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    robust_metrics = []
    for block in blocks:
        norm = np.linalg.norm(evaluate(block, robust_fit.x), axis=1)
        robust_metrics.append({
            'order': block['order'],
            'l2_ratio': float(np.sqrt(block['weights']@norm**2)
                              /block['baseline_l2']),
            'max_ratio': float(np.max(norm)/block['baseline_max'])})
    report['robust_metrics'] = robust_metrics
    report['robust_speed_coordinates_norm'] = float(np.linalg.norm(robust_fit.x))
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'optimizer_success': report['optimizer_success'],
                      'robust_optimizer_success': report['robust_optimizer_success'],
                      'robust_metrics': robust_metrics,
                      'robust_speed_coordinates_norm': report['robust_speed_coordinates_norm'],
                      'rows': [{'strength': row['strength'],
                                'worst_l2_ratio': row['worst_l2_ratio'],
                                'worst_max_ratio': row['worst_max_ratio'],
                                'metrics': row['metrics']} for row in rows]}),
          flush=True)


if __name__ == '__main__':
    run()
