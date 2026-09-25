"""Search exact-solenoidal axial shapes for local paper-style cone feasibility.

This is a frozen physical-field diagnostic. The paper's normalized stress
cone, supported waves, moment repair and smooth forcing are not implemented.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import differential_evolution

from curl_wave_cone_parameter_screen import cone_row
from joined_field import ROOT, coordinates
from joint_collar_fit import kinematics
from local_poloidal_basis_screen import LocalPoloidalMode
from radial_peak_cone import current_field


def run():
    current = current_field()
    compact = current.current.joint.base
    tau = .5/64
    radius = .0056890761915166545
    heights = (.0025, .00275, .003, .00325,
               .003432627453438968)
    modes = [LocalPoloidalMode(compact, i, j, 'odd',
                               reference_tau=tau, temporal_power=2.)
             for i in (0, 1) for j in (0, 1, 2)]
    nodes, weights = leggauss(10)
    blocks = []
    for z in heights:
        q = float(coordinates(0., z/np.sqrt(current.nu), tau,
                              compact.base.inner.h)['q'])
        inner_radius = np.sqrt(current.nu)*np.sqrt(2*q*3/64)
        quad_r, quad_w = [], []
        for lo, hi in zip((0., inner_radius), (inner_radius, radius)):
            quad_r.extend((lo+hi)/2+(hi-lo)/2*nodes)
            quad_w.extend((hi-lo)/2*weights)
        quad_r, quad_w = np.asarray(quad_r), np.asarray(quad_w)
        rr = np.r_[quad_r, radius]
        points = np.column_stack((rr, np.zeros(len(rr)),
                                  np.full(len(rr), z)))
        blocks.append((z, quad_r, quad_w, points))
    points = np.vstack([block[3] for block in blocks])
    hs = .0005*np.sqrt(current.nu*tau)
    ht = .0001*tau
    u0, J0, part0 = kinematics(current, points, tau, hs, ht)
    pieces = [kinematics(mode, points, tau, hs, ht) for mode in modes]
    du = np.stack([item[0] for item in pieces], axis=-1)
    dJ = np.stack([item[1] for item in pieces], axis=-1)
    dpart = np.stack([item[2] for item in pieces], axis=-1)
    max_base_speed = float(np.max(np.linalg.norm(u0, axis=1)))
    max_mode_speed = np.max(np.linalg.norm(du, axis=1), axis=0)
    scales = .5*max_base_speed/np.maximum(max_mode_speed, 1e-12)

    def evaluate(x):
        amplitudes = scales*np.asarray(x)
        u = u0+np.einsum('nik,k->ni', du, amplitudes)
        J = J0+np.einsum('nijk,k->nij', dJ, amplitudes)
        part = part0+np.einsum('nik,k->ni', dpart, amplitudes)
        residual = part+np.einsum('nij,nj->ni', J, u)
        rows, offset = [], 0
        for z, quad_r, quad_w, pts in blocks:
            sl = slice(offset, offset+len(pts))
            rows.append(cone_row(radius, z, u[sl][-1], J[sl][-1],
                                 residual[sl], quad_r, quad_w))
            offset += len(pts)
        return rows, float(np.max(np.linalg.norm(residual, axis=1))), \
            float(np.max(np.linalg.norm(u, axis=1)))

    initial_rows, initial_residual, _ = evaluate(np.zeros(len(modes)))
    best_by_pass = {'pass_count': sum(row['strict_local_pass']
                                      for row in initial_rows),
                    'normalized_coordinates': np.zeros(len(modes)).tolist(),
                    'max_sample_residual': initial_residual,
                    'rows': initial_rows}

    def objective(x):
        rows, max_residual, max_speed = evaluate(x)
        pass_count = sum(row['strict_local_pass'] for row in rows)
        if (pass_count > best_by_pass['pass_count'] or
                (pass_count == best_by_pass['pass_count'] and
                 max_residual < best_by_pass['max_sample_residual'])):
            best_by_pass.update(pass_count=pass_count,
                                normalized_coordinates=np.asarray(x).tolist(),
                                max_sample_residual=max_residual,
                                rows=rows)
        penalty = 0.
        for row in rows:
            if row['lambda_squared'] <= 0:
                penalty += 100.
            penalty += 100.*max(0., row['target_dot_N']/100.+.05)**2
            ratio = row.get('cone_ratio')
            penalty += 100.*max(0., (ratio if ratio is not None else 10.)-.9)**2
        penalty += .001*(max_residual/initial_residual)**2
        penalty += .001*(max_speed/max_base_speed)**2
        return penalty

    fit = differential_evolution(objective, [(-4., 4.)]*len(modes),
                                 seed=1, maxiter=120, popsize=8,
                                 polish=True, tol=1e-7)
    rows, max_residual, max_speed = evaluate(fit.x)
    report = {'tau': tau, 'radius': radius, 'heights': heights,
              'mode_ids': [{'radial_degree': mode.radial_degree,
                            'axial_degree': mode.axial_degree,
                            'parity': 'odd', 'temporal_power': 2.}
                           for mode in modes],
              'amplitude_scales': scales.tolist(),
              'normalized_coordinates': fit.x.tolist(),
              'amplitudes': (scales*fit.x).tolist(),
              'optimizer_success': bool(fit.success),
              'optimizer_message': fit.message,
              'initial_rows': initial_rows,
              'optimized_rows': rows,
              'initial_max_sample_residual': initial_residual,
              'optimized_max_sample_residual': max_residual,
              'optimized_max_sample_speed': max_speed,
              'all_five_local_cones_pass': all(row['strict_local_pass'] for row in rows),
              'best_by_pass_count': best_by_pass,
              'scope': 'Fixed radius, five positive heights, one registered time; physical-field analog of the paper local cone. Six exact-solenoidal radial/axial shapes optimized with soft feasibility penalties. Not the normalized paper cone, no supported nonaxisymmetric wave, full-volume momentum gate or force extension.',
              'accepted': False}
    (ROOT/'compact_potential'/'cone_axial_shape_search.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'initial_pass_count': sum(r['strict_local_pass'] for r in initial_rows),
                      'optimized_pass_count': sum(r['strict_local_pass'] for r in rows),
                      'best_observed_pass_count': best_by_pass['pass_count'],
                      'amplitudes': report['amplitudes'],
                      'max_sample_residual': max_residual,
                      'rows': rows}), flush=True)


if __name__ == '__main__':
    run()
