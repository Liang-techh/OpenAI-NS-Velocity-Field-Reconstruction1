"""Refit the solenoidal collar against two full-support quadratures."""
import json

import numpy as np
from scipy.optimize import least_squares

from axial_swirl_multitime_screen import load_dense_candidate
from dynamic_poloidal_multitime_screen import make_modes, full_residual
from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from local_poloidal_basis_screen import LocalPoloidalCandidate


def load_volume_refit_candidate(base=None):
    report = json.loads((ROOT/'compact_potential'/'dynamic_poloidal_volume_refit.json').read_text())
    base = base if base is not None else load_dense_candidate()
    return LocalPoloidalCandidate(base, make_modes(base), report['amplitudes'])


def build_block(base, modes, tau, order):
    points, weights = nodes(base, tau, order)
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    u0, J0, part0 = kinematics(base, points, tau, hs, ht)
    pieces = [kinematics(mode, points, tau, hs, ht) for mode in modes]
    du = np.stack([item[0] for item in pieces], axis=-1)
    dJ = np.stack([item[1] for item in pieces], axis=-1)
    dpart = np.stack([item[2] for item in pieces], axis=-1)
    residual = full_residual(u0, J0, part0)
    norm = np.linalg.norm(residual, axis=1)
    return {'order': order, 'points': points, 'weights': weights,
            'u0': u0, 'J0': J0, 'part0': part0,
            'du': du, 'dJ': dJ, 'dpart': dpart,
            'baseline_max': float(np.max(norm)),
            'baseline_l2': float(np.sqrt(weights@norm**2))}


def evaluate(block, amplitudes):
    u = block['u0'] + np.einsum('nik,k->ni', block['du'], amplitudes)
    J = block['J0'] + np.einsum('niak,k->nia', block['dJ'], amplitudes)
    part = block['part0'] + np.einsum('nik,k->ni', block['dpart'], amplitudes)
    return full_residual(u, J, part)


def metrics(block, amplitudes):
    residual = evaluate(block, amplitudes)
    norms = np.linalg.norm(residual, axis=1)
    worst = int(np.argmax(norms))
    return {'max': float(norms[worst]),
            'physical_volume_l2': float(np.sqrt(block['weights']@norms**2)),
            'worst_point': block['points'][worst].tolist(),
            'worst_residual': residual[worst].tolist()}


def run():
    base = load_dense_candidate()
    modes = make_modes(base)
    tau = .0084
    blocks = []
    for order in (6, 8):
        block = build_block(base, modes, tau, order)
        blocks.append(block)
        print(json.dumps({'loaded_order': order,
                          'baseline_max': block['baseline_max'],
                          'baseline_l2': block['baseline_l2']}), flush=True)
    max_base_speed = max(float(np.max(np.linalg.norm(block['u0'], axis=1)))
                         for block in blocks)
    max_mode_speed = np.maximum.reduce([
        np.max(np.linalg.norm(block['du'], axis=1), axis=0)
        for block in blocks])
    scales = .5*max_base_speed/np.maximum(max_mode_speed, 1e-12)

    def objective(x):
        parts = []
        for block in blocks:
            residual = evaluate(block, scales*x)
            norm2 = np.sum(residual**2, axis=1)
            count = len(norm2)
            parts.append((np.sqrt(block['weights'])[:, None]*residual
                          /block['baseline_l2']).ravel())
            parts.append((residual/block['baseline_max']/np.sqrt(count)).ravel())
            parts.append(.5*norm2/block['baseline_max']**2/np.sqrt(count))
        parts.append(.01*x)
        return np.concatenate(parts)

    fit = least_squares(objective, np.zeros(len(modes)), bounds=(-1., 1.),
                        max_nfev=300, xtol=1e-9, ftol=1e-9, gtol=1e-9)
    amplitudes = scales*fit.x
    rows = []
    for block in blocks:
        row = {'order': block['order'],
               'baseline': metrics(block, np.zeros(len(modes))),
               'refit': metrics(block, amplitudes)}
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = {'tau': tau, 'train_orders': [block['order'] for block in blocks],
              'mode_ids': [{'radial_degree': mode.radial_degree,
                            'axial_degree': mode.axial_degree,
                            'parity': mode.parity} for mode in modes],
              'reference_tau': .0084, 'temporal_power': .5,
              'amplitudes': amplitudes.tolist(),
              'normalized_coordinates': fit.x.tolist(),
              'optimizer_success': bool(fit.success),
              'optimizer_message': fit.message,
              'rows': rows,
              'scope': 'Full nonlinear physical momentum refit on whole-support Gauss6 and Gauss8 at one registered time. Objective combines volume L2, point residual and fourth-power peak penalty. Gauss10 and off-grid holdouts are not fitted; this is not an acceptance certificate.',
              'accepted': False}
    (ROOT/'compact_potential'/'dynamic_poloidal_volume_refit.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
