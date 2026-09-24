"""Fit annular solenoidal modes against multiple grids and peak residuals."""
import json

import numpy as np
from scipy.optimize import minimize

from annular_poloidal_multigrid_fit import WINDOWS, REFERENCE_TAU, precompute, recombine
from annular_poloidal_volume_fit import AnnularPoloidalMode
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


MODE_SPECS = tuple((lo, hi, degree)
                   for lo, hi in WINDOWS for degree in (0, 1, 2))


class RobustAnnularCandidate:
    def __init__(self, base, amplitudes):
        self.base = base
        self.nu = base.nu
        self.amplitudes = np.asarray(amplitudes)
        self.modes = [AnnularPoloidalMode(base.base, lo, hi, REFERENCE_TAU,
                                          axial_degree=degree)
                      for lo, hi, degree in MODE_SPECS]

    def support(self, tau):
        return self.base.support(tau)

    def fields(self, points, tau):
        u, p = self.base.fields(points, tau)
        for amplitude, mode in zip(self.amplitudes, self.modes):
            if amplitude:
                u += amplitude*mode.fields(points, tau)[0]
        return u, p


def load_candidate():
    report = json.loads((ROOT/'compact_potential'/'annular_robust_fit.json').read_text())
    return RobustAnnularCandidate(current_field(), report['amplitudes'])


def patch_data(field, modes):
    radii = np.array([.005, .0065, .008, .0095, .011])
    heights = np.array([.00285, .0031, .00335, .0036, .00385])
    points = np.array([[r, 0., sign*z]
                       for sign in (-1, 1) for z in heights for r in radii])
    weights = (2*np.pi*points[:, 0]*(radii[1]-radii[0])
               *(heights[1]-heights[0]))
    tau = REFERENCE_TAU
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u0, J0, part0 = kinematics(field, points, tau, hs, ht)
    pieces = [kinematics(mode, points, tau, hs, ht) for mode in modes]
    du = np.stack([part[0] for part in pieces], axis=-1)
    dJ = np.stack([part[1] for part in pieces], axis=-1)
    dpart = np.stack([part[2] for part in pieces], axis=-1)
    baseline = part0+np.einsum('nij,nj->ni', J0, u0)
    base_norm = np.linalg.norm(baseline, axis=1)
    data = {'k': 6., 'order': 'off_gauss_patch',
            'weights': weights, 'u0': u0, 'J0': J0, 'part0': part0,
            'du': du, 'dJ': dJ, 'dpart': dpart,
            'baseline_l2': float(np.sqrt(weights @ base_norm**2)),
            'baseline_max': float(np.max(base_norm))}
    print(json.dumps({'loaded_patch': len(points),
                      'baseline_l2': data['baseline_l2']}), flush=True)
    return data


def row(data, a, held_out=False):
    R, _ = recombine(data, a)
    norms = np.linalg.norm(R, axis=1)
    return {'k': data['k'], 'order': data['order'],
            'held_out': held_out,
            'baseline_max': data['baseline_max'],
            'baseline_l2': data['baseline_l2'],
            'candidate_max': float(np.max(norms)),
            'candidate_l2': float(np.sqrt(data['weights'] @ norms**2))}


def run():
    field = current_field()
    modes = [AnnularPoloidalMode(field.base, lo, hi, REFERENCE_TAU,
                                 axial_degree=degree)
             for lo, hi, degree in MODE_SPECS]
    train = [precompute(field, modes, k, order)
             for k, order in ((6., 6), (6., 8), (6., 10), (5.5, 6))]
    train.append(patch_data(field, modes))
    peak_weight = .3
    def objective(a):
        value = 0.
        gradient = np.zeros(len(modes))
        for data in train:
            R, tangent = recombine(data, a)
            weights = data['weights']
            l2factor = 1/(len(train)*data['baseline_l2']**2)
            value += .5*l2factor*np.sum(weights[:, None]*R**2)
            gradient += l2factor*np.einsum('n,ni,nik->k', weights, R, tangent)
            norm2 = np.sum(R**2, axis=1)
            peakfactor = peak_weight/(len(train)*len(R)*data['baseline_max']**4)
            value += .25*peakfactor*np.sum(norm2**2)
            gradient += peakfactor*np.einsum('n,ni,nik->k', norm2, R, tangent)
        return value, gradient
    previous = json.loads((ROOT/'compact_potential'/'annular_poloidal_multigrid_fit.json').read_text())
    initial = np.zeros(len(modes))
    # Previous eight-mode solution supplies a useful warm start; new degree-2 modes start at zero.
    for j, amplitude in enumerate(previous['amplitudes']):
        initial[3*(j//2)+(j%2)] = amplitude
    candidates = [minimize(objective, start, jac=True, method='L-BFGS-B',
                           bounds=[(-.1, .1)]*len(modes),
                           options={'maxiter': 500, 'ftol': 1e-13})
                  for start in (np.zeros(len(modes)), initial)]
    solution = min(candidates, key=lambda fit: fit.fun)
    rows = [row(data, solution.x) for data in train]
    holdouts = [precompute(field, modes, k, order)
                for k, order in ((6., 7), (5.75, 7))]
    rows.extend(row(data, solution.x, True) for data in holdouts)
    report = {'mode_specs': MODE_SPECS, 'reference_tau': REFERENCE_TAU,
              'temporal_power': 2., 'peak_weight': peak_weight,
              'amplitudes': solution.x.tolist(),
              'optimizer_success': bool(solution.success),
              'optimizer_message': solution.message,
              'objective': float(solution.fun), 'rows': rows,
              'scope': 'Twelve compact solenoidal modes fit to Gauss6/8/10 k6, Gauss6 k5.5, and off-Gauss hotspot patch, with quartic sampled-peak penalty. Gauss7 k6 and k5.75 held out. No continuum, cone, or PDE certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_robust_fit.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'amplitudes': report['amplitudes'],
                      'objective': report['objective'], 'rows': rows}), flush=True)


if __name__ == '__main__':
    run()
