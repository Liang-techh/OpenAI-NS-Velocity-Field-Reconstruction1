"""Fit time-scaled solenoidal meridional collar modes to full momentum."""
import json

import numpy as np
from scipy.optimize import least_squares

from axial_swirl_multitime_screen import load_dense_candidate
from joined_field import ROOT
from joint_collar_fit import kinematics
from local_poloidal_basis_screen import (LocalPoloidalCandidate,
                                         LocalPoloidalMode, compact_base)


def make_modes(base):
    compact = compact_base(base)
    return [LocalPoloidalMode(compact, i, j, 'odd',
                              reference_tau=.0084, temporal_power=.5)
            for i in (0, 1) for j in (0, 1, 2)]


def load_dynamic_candidate():
    report = json.loads((ROOT/'compact_potential'/'dynamic_poloidal_multitime.json').read_text())
    base = load_dense_candidate()
    return LocalPoloidalCandidate(base, make_modes(base), report['amplitudes'])


def full_residual(u, J, part):
    return part + np.einsum('nij,nj->ni', J, u)


def run():
    base = load_dense_candidate()
    modes = make_modes(base)
    blocks = []
    max_base_speed = 0.
    max_mode_speed = np.zeros(len(modes))
    for tau in (.0084, .012, .024):
        _, radius, zflat, zsupport = base.support(tau)
        points = np.array([[rf*radius, 0., sign*(zflat+s*(zsupport-zflat))]
                           for sign in (-1, 1) for rf in (.25, .4, .55)
                           for s in (.2, .35, .5, .65, .8)])
        hs = .0005*np.sqrt(base.nu*tau)
        ht = .0001*tau
        u0, J0, part0 = kinematics(base, points, tau, hs, ht)
        pieces = [kinematics(mode, points, tau, hs, ht) for mode in modes]
        du = np.stack([item[0] for item in pieces], axis=-1)
        dJ = np.stack([item[1] for item in pieces], axis=-1)
        dpart = np.stack([item[2] for item in pieces], axis=-1)
        baseline = full_residual(u0, J0, part0)
        baseline_norms = np.linalg.norm(baseline, axis=1)
        scale = float(np.max(baseline_norms))
        max_base_speed = max(max_base_speed,
                             float(np.max(np.linalg.norm(u0, axis=1))))
        max_mode_speed = np.maximum(max_mode_speed,
            np.max(np.linalg.norm(du, axis=1), axis=0))
        blocks.append({'tau': tau, 'u0': u0, 'J0': J0, 'part0': part0,
                       'du': du, 'dJ': dJ, 'dpart': dpart,
                       'scale': scale,
                       'baseline_max': scale,
                       'baseline_radial_max': float(np.max(np.abs(baseline[:, 0]))),
                       'baseline_angular_max': float(np.max(np.abs(baseline[:, 1]))),
                       'baseline_axial_max': float(np.max(np.abs(baseline[:, 2])))})
        print(json.dumps({'loaded_tau': tau, 'points': len(points),
                          'baseline_max': scale}), flush=True)
    amplitude_scales = .5*max_base_speed/np.maximum(max_mode_speed, 1e-12)

    def evaluate(block, x):
        a = amplitude_scales*x
        u = block['u0'] + np.einsum('nik,k->ni', block['du'], a)
        J = block['J0'] + np.einsum('niak,k->nia', block['dJ'], a)
        part = block['part0'] + np.einsum('nik,k->ni', block['dpart'], a)
        return full_residual(u, J, part)

    def objective(x):
        residuals = [(evaluate(block, x)/block['scale']).ravel()
                     for block in blocks]
        return np.concatenate(residuals + [.01*x])

    fit = least_squares(objective, np.zeros(len(modes)), bounds=(-1., 1.),
                        max_nfev=200, xtol=1e-9, ftol=1e-9, gtol=1e-9)
    amplitudes = amplitude_scales*fit.x
    rows = []
    for block in blocks:
        residual = evaluate(block, fit.x)
        norms = np.linalg.norm(residual, axis=1)
        rows.append({'tau': block['tau'], 'point_count': len(norms),
                     'baseline_max': block['baseline_max'],
                     'corrected_max': float(np.max(norms)),
                     'baseline_radial_max': block['baseline_radial_max'],
                     'corrected_radial_max': float(np.max(np.abs(residual[:, 0]))),
                     'baseline_angular_max': block['baseline_angular_max'],
                     'corrected_angular_max': float(np.max(np.abs(residual[:, 1]))),
                     'baseline_axial_max': block['baseline_axial_max'],
                     'corrected_axial_max': float(np.max(np.abs(residual[:, 2])))})
        print(json.dumps(rows[-1]), flush=True)
    report = {'times': [block['tau'] for block in blocks],
              'mode_ids': [{'radial_degree': mode.radial_degree,
                            'axial_degree': mode.axial_degree,
                            'parity': mode.parity} for mode in modes],
              'reference_tau': .0084, 'temporal_power': .5,
              'amplitude_scales': amplitude_scales.tolist(),
              'amplitudes': amplitudes.tolist(),
              'optimizer_success': bool(fit.success),
              'optimizer_message': fit.message,
              'normalized_coordinates': fit.x.tolist(),
              'rows': rows,
              'scope': 'Full nonlinear Cartesian momentum at 30 moving collar points at each of three registered times. Fit is local and requires independent points, broad spatial checks, critical-time forcing and exterior matching before acceptance.',
              'accepted': False}
    (ROOT/'compact_potential'/'dynamic_poloidal_multitime.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
