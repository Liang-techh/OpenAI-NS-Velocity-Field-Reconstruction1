"""Fit compact radial/axial pressure shapes to a sampled 2D cone patch.

The velocity is the registered annular candidate. Pressure changes its full
physical momentum residual by its gradient, but does not change the velocity
or the shear cone. This is a finite fixed-time design screen only.
"""
import json
import sys

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.optimize import linprog

from annular_pressure_scale_screen import load_candidate
from curl_wave_cone_region import evaluate_height
from joined_field import ROOT
from joint_collar_fit import nodes
from radial_peak_cone import operator


def pressure_axial_gradient(r, z, radius, zflat, zsupport, radial_degree,
                            axial_degree):
    """Unit compact pressure basis derivative in positive z collar."""
    r = np.asarray(r)
    y = r/radius
    s = (z-zflat)/(zsupport-zflat)
    if not (0 < s < 1):
        return np.zeros_like(r)
    bubble = 256*s**4*(1-s)**4
    bubble_prime = 1024*s**3*(1-s)**3*(1-2*s)
    radial = (1-y*y)**5*Legendre.basis(radial_degree)(2*y*y-1)
    axial = Legendre.basis(axial_degree)
    axial_derivative = (bubble_prime*axial(2*s-1)
                        + 2*bubble*axial.deriv()(2*s-1))/(zsupport-zflat)
    return radial*axial_derivative


class RadialPressureCandidate:
    def __init__(self, base, degrees, amplitudes):
        self.base = base
        self.nu = base.nu
        self.degrees = degrees
        self.amplitudes = np.asarray(amplitudes)

    def support(self, tau):
        return self.base.support(tau)

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        points = np.asarray(points)
        _, radius, zflat, zsupport = self.support(tau)
        y = np.hypot(points[:, 0], points[:, 1])/radius
        s = (np.abs(points[:, 2])-zflat)/(zsupport-zflat)
        inside = (y < 1) & (s > 0) & (s < 1)
        if np.any(inside):
            yy, ss = y[inside], s[inside]
            bubble = 256*ss**4*(1-ss)**4
            for (i, j), amplitude in zip(self.degrees, self.amplitudes):
                if amplitude:
                    pressure[inside] += (amplitude*(1-yy*yy)**5
                                         * Legendre.basis(i)(2*yy*yy-1)
                                         * bubble*Legendre.basis(j)(2*ss-1))
        return velocity, pressure


def load_dense_candidate():
    report = json.loads((ROOT/'compact_potential'/'radial_pressure_patch_dense.json').read_text())
    if not report['feasible'] or report.get('holdout_pass_count') != len(report.get('holdout_rows', [])):
        raise ValueError('No passing dense pressure candidate is registered')
    return RadialPressureCandidate(load_candidate(),
                                   report['mode_degrees'], report['amplitudes'])


def run():
    field = load_candidate()
    report0 = json.loads((ROOT/'compact_potential'/'annular_cone_area_map.json').read_text())
    tau = report0['tau']
    dense = '--dense' in sys.argv
    radii = ([.0075, .00825, .009, .010, .011, .012, .013]
             if dense else [.0075, .009, .011, .013])
    heights = ([.00355, .003625, .0037]
               if dense else [.003, .0032, .0034, .00355, .0037])
    _, radius_support, zflat, zsupport = field.support(tau)
    mode_degrees = [(i, j) for i in range(4) for j in range(4)]
    rows = ([row for z in heights for row in evaluate_height(
                field, np.asarray(radii), z, tau, radius_support, zsupport,
                quadrature_order=8)] if dense else
            [row for row in report0['rows']
             if row['r'] in radii and row['z'] in heights])
    centers = np.array([[row['r'], 0., row['z']] for row in rows])
    velocity, gradient, _ = operator(field, centers, tau)
    qnodes, qweights = leggauss(24)
    A, b, detail = [], [], []
    for row, u, J in zip(rows, velocity, gradient):
        r, z = row['r'], row['z']
        F = u[1]/r
        shear = np.array([J[1, 0]-F, J[2, 0]])
        N = shear/np.linalg.norm(shear)
        K = np.array([-N[1], N[0]])
        lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
        if lam2 <= 0:
            raise ValueError(f'Nonpositive lambda squared at {(r, z)}')
        c = np.sqrt(lam2)/(2*F*N[0])
        qr = r*(qnodes+1)/2
        qw = r*qweights/2
        pressure_targets = np.array([
            [0., -np.dot(qw*qr, pressure_axial_gradient(
                qr, z, radius_support, zflat, zsupport, i, j))/r]
            for i, j in mode_degrees]).T
        target0 = np.asarray(row['target'])
        for sign in (-1, 1):
            direction = .8*N+sign*c*K
            A.append(direction @ pressure_targets)
            b.append(-1.-direction @ target0)
        detail.append({'r': r, 'z': z, 'lambda_squared': float(lam2),
                       'baseline_cone_ratio': row.get('cone_ratio'),
                       'baseline_strict_pass': row['strict_local_pass'],
                       'N': N.tolist(), 'K': K.tolist(), 'c': float(c),
                       'target0': target0.tolist(),
                       'pressure_targets': pressure_targets.tolist()})
    A, b = np.asarray(A), np.asarray(b)
    scales = 1/np.maximum(np.max(np.abs(A), axis=0), 1e-12)
    fit = linprog(np.zeros(len(mode_degrees)), A_ub=A*scales,
                  b_ub=b, bounds=[(None, None)]*len(mode_degrees),
                  method='highs')
    out = {'tau': tau, 'radii': radii, 'heights': heights,
           'mode_degrees': mode_degrees, 'feasible': bool(fit.success),
           'solver_message': fit.message, 'rows': detail,
           'scope': 'Pressure-only finite 4x5 point cone screen on one annular candidate at one time. No held-out, continuum, time-persistence, pressure-Poisson, or global momentum certificate.',
           'accepted': False}
    if fit.success:
        amplitudes = scales*fit.x
        out['amplitudes'] = amplitudes.tolist()
        out['maximum_constraint_violation'] = float(np.max(A@amplitudes-b))
        out['corrected_rows'] = []
        for row in detail:
            target = np.asarray(row['target0'])+np.asarray(row['pressure_targets'])@amplitudes
            N, K = np.asarray(row['N']), np.asarray(row['K'])
            out['corrected_rows'].append({
                'r': row['r'], 'z': row['z'],
                'target_dot_N': float(target@N),
                'cone_ratio': float(abs(row['c']*(target@K)/(target@N)))})
        if dense:
            candidate = RadialPressureCandidate(field, mode_degrees, amplitudes)
            holdout = []
            for z in (.0035875, .0036625):
                holdout.extend(evaluate_height(
                    candidate, np.array([.007875, .008625, .0095,
                                         .0105, .0115, .0125]),
                    z, tau, radius_support, zsupport,
                    quadrature_order=8))
            out['holdout_rows'] = holdout
            out['holdout_pass_count'] = sum(
                row['strict_local_pass'] and row.get('cone_ratio', np.inf) < .8
                for row in holdout)
            volume = []
            for order in (6, 8):
                points, weights = nodes(field, tau, order)
                before = operator(field, points, tau)[2]
                after = operator(candidate, points, tau)[2]
                volume.append({
                    'order': order,
                    'before_l2': float(np.sqrt(weights@np.sum(before**2, axis=1))),
                    'after_l2': float(np.sqrt(weights@np.sum(after**2, axis=1))),
                    'before_max': float(np.max(np.linalg.norm(before, axis=1))),
                    'after_max': float(np.max(np.linalg.norm(after, axis=1)))})
            out['volume_rows'] = volume
    else:
        feasible_rectangles = []
        for ri in range(len(radii)-1):
            for rj in range(ri+1, len(radii)):
                for zi in range(len(heights)-1):
                    for zj in range(zi+1, len(heights)):
                        selected = [n for n, row in enumerate(detail)
                                    if radii[ri] <= row['r'] <= radii[rj]
                                    and heights[zi] <= row['z'] <= heights[zj]]
                        constraints = np.array([k for n in selected
                                                for k in (2*n, 2*n+1)])
                        subfit = linprog(np.zeros(len(mode_degrees)),
                                         A_ub=A[constraints]*scales,
                                         b_ub=b[constraints],
                                         bounds=[(None, None)]*len(mode_degrees),
                                         method='highs')
                        if subfit.success:
                            feasible_rectangles.append({
                                'r_interval': [radii[ri], radii[rj]],
                                'z_interval': [heights[zi], heights[zj]],
                                'node_count': len(selected),
                                'area': (radii[rj]-radii[ri])*(heights[zj]-heights[zi])})
        feasible_rectangles.sort(key=lambda row: row['area'], reverse=True)
        out['feasible_subrectangles'] = feasible_rectangles[:20]
        if feasible_rectangles:
            best = feasible_rectangles[0]
            selected = [n for n, row in enumerate(detail)
                        if best['r_interval'][0] <= row['r'] <= best['r_interval'][1]
                        and best['z_interval'][0] <= row['z'] <= best['z_interval'][1]]
            constraints = np.array([k for n in selected for k in (2*n, 2*n+1)])
            subfit = linprog(np.zeros(len(mode_degrees)),
                             A_ub=A[constraints]*scales,
                             b_ub=b[constraints],
                             bounds=[(None, None)]*len(mode_degrees),
                             method='highs')
            amplitudes = scales*subfit.x
            best['amplitudes'] = amplitudes.tolist()
            best['maximum_constraint_violation'] = float(np.max(A[constraints]@amplitudes-b[constraints]))
            best['corrected_full_grid'] = []
            for row in detail:
                target = np.asarray(row['target0'])+np.asarray(row['pressure_targets'])@amplitudes
                N, K = np.asarray(row['N']), np.asarray(row['K'])
                projection = float(target@N)
                best['corrected_full_grid'].append({
                    'r': row['r'], 'z': row['z'],
                    'target_dot_N': projection,
                    'cone_ratio': float(abs(row['c']*(target@K)/projection))
                    if projection else None,
                    'passes': bool(projection < 0 and abs(row['c']*(target@K)/projection) < .8)})
            candidate = RadialPressureCandidate(field, mode_degrees, amplitudes)
            holdout = []
            for z in (.00355, .003625, .0037):
                holdout.extend(evaluate_height(candidate,
                                               np.array([.00825, .010, .012]),
                                               z, tau, radius_support, zsupport,
                                               quadrature_order=8))
            best['holdout_rows'] = holdout
            best['holdout_pass_count'] = sum(
                row['strict_local_pass'] and row.get('cone_ratio', np.inf) < .8
                for row in holdout)
    path = ROOT/'compact_potential'/(
        'radial_pressure_patch_dense.json' if dense else
        'radial_pressure_patch_screen.json')
    path.write_bytes((json.dumps(out, indent=2)+'\n').encode())
    print(json.dumps({'feasible': out['feasible'],
                      'solver_message': out['solver_message'],
                      'maximum_constraint_violation': out.get('maximum_constraint_violation'),
                      'best_subrectangle': {
                          k: v for k, v in out.get('feasible_subrectangles', [{}])[0].items()
                          if k in ('r_interval', 'z_interval', 'node_count',
                                   'area', 'holdout_pass_count')},
                      'holdout_pass_count': out.get('holdout_pass_count'),
                      'volume_rows': out.get('volume_rows'),
                      'max_cone_ratio': max((x['cone_ratio'] for x in out.get('corrected_rows', [])), default=None)}),
          flush=True)


if __name__ == '__main__':
    run()
