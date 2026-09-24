"""Fit compact solenoidal collar modes on several volume grids and times."""
import json

import numpy as np
from scipy.optimize import minimize

from annular_poloidal_volume_fit import AnnularPoloidalMode
from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from radial_peak_cone import current_field


WINDOWS = ((.12, .36), (.28, .55), (.12, .65), (.45, .80))
MODE_SPECS = tuple((lo, hi, degree)
                   for lo, hi in WINDOWS for degree in (0, 1))
REFERENCE_TAU = .5/64


class MultigridAnnularCandidate:
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
    report = json.loads((ROOT/'compact_potential'/'annular_poloidal_multigrid_fit.json').read_text())
    return MultigridAnnularCandidate(current_field(), report['amplitudes'])


def precompute(field, modes, k, order):
    tau = .5*2**(-k)
    points, weights = nodes(field, tau, order)
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u0, J0, part0 = kinematics(field, points, tau, hs, ht)
    pieces = [kinematics(mode, points, tau, hs, ht) for mode in modes]
    du = np.stack([part[0] for part in pieces], axis=-1)
    dJ = np.stack([part[1] for part in pieces], axis=-1)
    dpart = np.stack([part[2] for part in pieces], axis=-1)
    baseline = part0+np.einsum('nij,nj->ni', J0, u0)
    baseline_l2 = float(np.linalg.norm(np.sqrt(weights)[:, None]*baseline))
    print(json.dumps({'loaded_k': k, 'order': order,
                      'baseline_l2': baseline_l2}), flush=True)
    return {'k': k, 'order': order, 'weights': weights,
            'u0': u0, 'J0': J0, 'part0': part0,
            'du': du, 'dJ': dJ, 'dpart': dpart,
            'baseline_l2': baseline_l2,
            'baseline_max': float(np.max(np.linalg.norm(baseline, axis=1)))}


def recombine(data, a):
    u = data['u0']+np.einsum('nik,k->ni', data['du'], a)
    J = data['J0']+np.einsum('nijk,k->nij', data['dJ'], a)
    part = data['part0']+np.einsum('nik,k->ni', data['dpart'], a)
    R = part+np.einsum('nij,nj->ni', J, u)
    tangent = (data['dpart']+np.einsum('nijk,nj->nik', data['dJ'], u)
               +np.einsum('nij,njk->nik', J, data['du']))
    return R, tangent


def run():
    field = current_field()
    modes = [AnnularPoloidalMode(field.base, lo, hi, REFERENCE_TAU,
                                 axial_degree=degree)
             for lo, hi, degree in MODE_SPECS]
    train = [precompute(field, modes, k, order)
             for k, order in ((6., 8), (6., 10), (5.5, 6))]
    def objective(a):
        value = 0.
        gradient = np.zeros(len(modes))
        for data in train:
            R, tangent = recombine(data, a)
            factor = 1/(len(train)*data['baseline_l2']**2)
            w = data['weights']
            value += .5*factor*np.sum(w[:, None]*R**2)
            gradient += factor*np.einsum('n,ni,nik->k', w, R, tangent)
        return value, gradient
    solution = minimize(objective, np.zeros(len(modes)), jac=True,
                        method='L-BFGS-B', bounds=[(-.1, .1)]*len(modes),
                        options={'maxiter': 400, 'ftol': 1e-13})
    rows = []
    for data in train:
        R, _ = recombine(data, solution.x)
        norms = np.linalg.norm(R, axis=1)
        rows.append({'k': data['k'], 'order': data['order'],
                     'baseline_max': data['baseline_max'],
                     'baseline_l2': data['baseline_l2'],
                     'candidate_max': float(np.max(norms)),
                     'candidate_l2': float(np.sqrt(data['weights'] @ norms**2))})
    holdout = precompute(field, modes, 6., 6)
    R, _ = recombine(holdout, solution.x)
    norms = np.linalg.norm(R, axis=1)
    rows.append({'k': 6., 'order': 6,
                 'held_out': True,
                 'baseline_max': holdout['baseline_max'],
                 'baseline_l2': holdout['baseline_l2'],
                 'candidate_max': float(np.max(norms)),
                 'candidate_l2': float(np.sqrt(holdout['weights'] @ norms**2))})
    report = {'windows': WINDOWS, 'mode_specs': MODE_SPECS,
              'reference_tau': REFERENCE_TAU,
              'temporal_power': 2., 'amplitudes': solution.x.tolist(),
              'optimizer_success': bool(solution.success),
              'optimizer_message': solution.message, 'rows': rows,
              'scope': 'Four compact solenoidal modes fit jointly on Gauss8/Gauss10 at k6 and Gauss6 at k5.5, with k6 Gauss6 held out. Finite grids only, no local cone or continuum certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_poloidal_multigrid_fit.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'amplitudes': report['amplitudes'], 'rows': rows}), flush=True)


if __name__ == '__main__':
    run()
