"""Fit compact axial-collar swirl transport across registered times.

Axisymmetric pure swirl is solenoidal. Its azimuthal momentum equation is
linear when the meridional velocity is fixed, but the full radial momentum
still changes through the centrifugal term and is checked after fitting.
"""
import json
import sys

import numpy as np
from numpy.polynomial.legendre import Legendre

from joined_field import ROOT
from joint_collar_fit import kinematics
from local_poloidal_basis_screen import load_robust_candidate


class AxialSwirlMode:
    def __init__(self, base, radial_degree, axial_degree,
                 reference_tau=.0084, temporal_power=.5):
        self.base = base
        self.nu = base.nu
        self.radial_degree = radial_degree
        self.axial_degree = axial_degree
        self.reference_tau = reference_tau
        self.temporal_power = temporal_power

    def fields(self, points, tau):
        points = np.asarray(points, float)
        velocity = np.zeros_like(points)
        _, radius, zflat, zsupport = self.base.support(tau)
        r = np.hypot(points[:, 0], points[:, 1])
        y = r/radius
        s = (np.abs(points[:, 2])-zflat)/(zsupport-zflat)
        inside = (y < 1) & (s > 0) & (s < 1) & (r > 0)
        if np.any(inside):
            yy, ss = y[inside], s[inside]
            radial = (yy*(1-yy**2)**5
                      * Legendre.basis(self.radial_degree)(2*yy**2-1))
            axial = (1024*ss**5*(1-ss)**5
                     * Legendre.basis(self.axial_degree)(2*ss-1))
            theta = ((self.reference_tau/tau)**self.temporal_power
                     * radial*axial)
            velocity[inside, 0] = -theta*points[inside, 1]/r[inside]
            velocity[inside, 1] = theta*points[inside, 0]/r[inside]
        return velocity, np.zeros(len(points))


class SwirlCandidate:
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


def load_dense_candidate(base=None):
    """Load the multi-time dense diagnostic candidate (not accepted)."""
    report = json.loads((ROOT/'compact_potential'/'axial_swirl_dense.json').read_text())
    if not report['dense']:
        raise ValueError('Expected dense axial-swirl report')
    base = base if base is not None else load_robust_candidate(.1)
    modes = [AxialSwirlMode(base, item['radial_degree'], item['axial_degree'],
                            report['reference_tau'], report['temporal_power'])
             for item in report['mode_ids']]
    return SwirlCandidate(base, modes, report['amplitudes'])


def momentum(field, points, tau):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, grad, part = kinematics(field, points, tau, hs, ht)
    return u, grad, part + np.einsum('nij,nj->ni', grad, u)


def run():
    low_order = '--low-order' in sys.argv
    medium_order = '--medium-order' in sys.argv
    dense = '--dense' in sys.argv
    if sum((low_order, medium_order, dense)) > 1:
        raise ValueError('Choose one basis order')
    source = json.loads((ROOT/'compact_potential'/'forcing_extension_screen.json').read_text())
    base = load_robust_candidate(.1)
    modes = [AxialSwirlMode(base, i, j)
             for i in ((0,) if low_order or medium_order else (0, 1))
             for j in ((0,) if low_order else
                       (0, 1) if medium_order else (0, 1, 2))]
    fit_blocks = []
    cases = []
    for item in source['rows'][:5]:
        tau = item['tau']
        _, radius, zflat, zsupport = base.support(tau)
        radial_fractions = (.25, .4, .55) if dense else (.3, .5)
        axial_fractions = (.2, .35, .5, .65, .8) if dense else (.35, .65)
        points = np.array([[rf*radius, 0., sign*(zflat+s*(zsupport-zflat))]
                           for sign in (-1, 1) for rf in radial_fractions
                           for s in axial_fractions])
        u, J, residual = momentum(base, points, tau)
        columns = []
        for mode in modes:
            du, dJ, dpart = kinematics(
                mode, points, tau, .0005*np.sqrt(base.nu*tau), .0001*tau)
            columns.append((du, dJ, dpart))
        matrix = np.stack([
            dpart[:, 1] + np.einsum('nij,nj->ni', dJ, u)[:, 1]
            + np.einsum('nij,nj->ni', J, du)[:, 1]
            for du, dJ, dpart in columns], axis=1)
        scale = np.max(np.abs(residual[:, 1]))
        fit_blocks.append((matrix/scale, -residual[:, 1]/scale))
        cases.append({'tau': tau, 'points': points,
                      'baseline': residual})
    matrix = np.concatenate([item[0] for item in fit_blocks])
    target = np.concatenate([item[1] for item in fit_blocks])
    amplitudes, *_ = np.linalg.lstsq(matrix, target, rcond=1e-10)
    candidate = SwirlCandidate(base, modes, amplitudes)
    rows = []
    for case in cases:
        _, _, corrected = momentum(candidate, case['points'], case['tau'])
        before, after = case['baseline'], corrected
        rows.append({'tau': case['tau'],
                     'baseline_sample_max': float(np.max(np.linalg.norm(before, axis=1))),
                     'corrected_sample_max': float(np.max(np.linalg.norm(after, axis=1))),
                     'baseline_angular_max': float(np.max(np.abs(before[:, 1]))),
                     'corrected_angular_max': float(np.max(np.abs(after[:, 1]))),
                     'baseline_radial_max': float(np.max(np.abs(before[:, 0]))),
                     'corrected_radial_max': float(np.max(np.abs(after[:, 0])))})
        print(json.dumps(rows[-1]), flush=True)
    report = {'mode_ids': [{'radial_degree': mode.radial_degree,
                            'axial_degree': mode.axial_degree}
                           for mode in modes],
              'reference_tau': .0084, 'temporal_power': .5,
              'low_order': low_order, 'medium_order': medium_order,
              'dense': dense,
              'amplitudes': amplitudes.tolist(), 'rows': rows,
              'scope': 'Moving axial-collar points at each of five registered times (30 per time for dense, eight otherwise). Least-squares fit of azimuthal momentum with full-operator check at the same points; separate holdout required before interpreting gains.',
              'accepted': False}
    filename = ('axial_swirl_low_order.json' if low_order else
                'axial_swirl_medium_order.json' if medium_order else
                'axial_swirl_dense.json' if dense else
                'axial_swirl_multitime_screen.json')
    (ROOT/'compact_potential'/filename).write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
