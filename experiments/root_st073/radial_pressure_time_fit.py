"""Fit constant or affine-time compact pressure increments on a sampled box."""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import linprog

from joined_field import ROOT
from radial_peak_cone import operator
from radial_pressure_patch_screen import (load_dense_candidate,
                                          pressure_axial_gradient)


def solve(A, b):
    scales = 1/np.maximum(np.max(np.abs(A), axis=0), 1e-12)
    fit = linprog(np.zeros(A.shape[1]), A_ub=A*scales, b_ub=b,
                  bounds=[(None, None)]*A.shape[1], method='highs')
    if not fit.success:
        return {'feasible': False, 'message': fit.message}
    amplitudes = scales*fit.x
    return {'feasible': True, 'message': fit.message,
            'amplitudes': amplitudes.tolist(),
            'maximum_constraint_violation': float(np.max(A@amplitudes-b))}


def run():
    field = load_dense_candidate()
    source = json.loads((ROOT/'compact_potential'/'radial_pressure_time_screen.json').read_text())
    degrees = [tuple(pair) for pair in json.loads(
        (ROOT/'compact_potential'/'radial_pressure_patch_dense.json').read_text())['mode_degrees']]
    rows = source['rows']
    A, b, diagnostics = [], [], []
    qnodes, qweights = leggauss(24)
    centers = np.array([[row['r'], 0., row['z']] for row in rows])
    for tau in source['times']:
        at_time = [i for i, row in enumerate(rows) if row['tau'] == tau]
        u, J, _ = operator(field, centers[at_time], tau)
        _, radius_support, zflat, zsupport = field.support(tau)
        for i, ui, Ji in zip(at_time, u, J):
            row = rows[i]
            r, z = row['r'], row['z']
            F = ui[1]/r
            shear = np.array([Ji[1, 0]-F, Ji[2, 0]])
            N = shear/np.linalg.norm(shear)
            K = np.array([-N[1], N[0]])
            lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
            if lam2 <= 0:
                raise ValueError(f'Nonpositive lambda squared at {(tau, r, z)}')
            c = np.sqrt(lam2)/(2*F*N[0])
            qr = r*(qnodes+1)/2
            qw = r*qweights/2
            targets = np.array([
                [0., -np.dot(qw*qr, pressure_axial_gradient(
                    qr, z, radius_support, zflat, zsupport, a, d))/r]
                for a, d in degrees]).T
            target0 = np.asarray(row['target'])
            dt = (tau-.0084)/.00015
            for sign in (-1, 1):
                direction = .8*N+sign*c*K
                increment = direction@targets
                A.append(np.r_[increment, dt*increment])
                b.append(-1.-direction@target0)
            diagnostics.append({'tau': tau, 'r': r, 'z': z,
                                'baseline_pass': row['strict_local_pass'],
                                'baseline_ratio': row.get('cone_ratio')})
    A, b = np.asarray(A), np.asarray(b)
    constant = solve(A[:, :len(degrees)], b)
    affine = solve(A, b)
    report = {'times': source['times'], 'radii': source['radii'],
              'heights': source['heights'], 'mode_degrees': degrees,
              'constant_fit': constant, 'affine_time_fit': affine,
              'baseline_pass_count': source['passing_count'],
              'total_count': len(rows), 'rows': diagnostics,
              'scope': 'Linear pressure correction at three sampled times on fixed 7x3 spatial grid. Affine-time amplitudes are not yet checked at intermediate times or for full momentum/pressure-Poisson consistency.',
              'accepted': False}
    path = ROOT/'compact_potential'/'radial_pressure_time_fit.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'constant_feasible': constant['feasible'],
                      'affine_time_feasible': affine['feasible'],
                      'affine_max_violation': affine.get('maximum_constraint_violation')}),
          flush=True)


if __name__ == '__main__':
    run()
