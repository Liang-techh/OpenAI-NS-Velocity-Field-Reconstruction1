"""Minimize sampled physical momentum subject to the spacetime cone margins."""
import json

import numpy as np
from numpy.polynomial.legendre import Legendre
from scipy.optimize import minimize

from joined_field import ROOT
from joint_collar_fit import nodes
from radial_peak_cone import operator
from radial_pressure_patch_screen import (RadialPressureCandidate,
                                          load_dense_candidate)


def pressure_gradient_columns(points, tau, field, degrees):
    """Analytic Cartesian gradients of compact radial/axial pressure modes."""
    points = np.asarray(points)
    columns = np.zeros((len(points), 3, len(degrees)))
    _, radius, zflat, zsupport = field.support(tau)
    rr = np.hypot(points[:, 0], points[:, 1])
    y = rr/radius
    s = (np.abs(points[:, 2])-zflat)/(zsupport-zflat)
    inside = (y < 1) & (s > 0) & (s < 1)
    if not np.any(inside):
        return columns
    yy, ss = y[inside], s[inside]
    xx = 2*yy*yy-1
    zz = 2*ss-1
    bubble = 256*ss**4*(1-ss)**4
    bubble_prime = 1024*ss**3*(1-ss)**3*(1-2*ss)
    cart_r = np.zeros((len(yy), 2))
    nonzero = rr[inside] > 0
    cart_r[nonzero] = points[inside, :2][nonzero]/rr[inside][nonzero, None]
    for col, (i, j) in enumerate(degrees):
        radial = Legendre.basis(i)
        axial = Legendre.basis(j)
        rv = radial(xx)
        dr = (-10*yy*(1-yy*yy)**4*rv
              +4*yy*(1-yy*yy)**5*radial.deriv()(xx))/radius
        av = axial(zz)
        az = (bubble_prime*av+2*bubble*axial.deriv()(zz))/(zsupport-zflat)
        dpdr = dr*bubble*av
        dpdz = (1-yy*yy)**5*rv*az*np.sign(points[inside, 2])
        columns[inside, 0, col] = dpdr*cart_r[:, 0]
        columns[inside, 1, col] = dpdr*cart_r[:, 1]
        columns[inside, 2, col] = dpdz
    return columns


def run():
    base = load_dense_candidate()
    source = json.loads((ROOT/'compact_potential'/'radial_pressure_time_fit.json').read_text())
    degrees = source['mode_degrees']
    A = np.asarray(source['constant_cone_matrix'])
    b = np.asarray(source['cone_rhs'])
    initial = np.asarray(source['constant_fit']['amplitudes'])
    grids = []
    weighted_columns, weighted_residuals = [], []
    for tau, order in ((.00825, 6), (.0084, 6), (.00855, 6),
                       (.0084, 8), (.0084, 10)):
        points, weights = nodes(base, tau, order)
        residual = operator(base, points, tau)[2]
        columns = pressure_gradient_columns(points, tau, base, degrees)
        baseline_l2 = float(np.sqrt(weights@np.sum(residual**2, axis=1)))
        root_weights = np.sqrt(weights)/baseline_l2
        weighted_residuals.append((root_weights[:, None]*residual).ravel())
        weighted_columns.append((root_weights[:, None, None]*columns).reshape(-1, len(degrees)))
        grids.append({'tau': tau, 'order': order, 'baseline_l2': baseline_l2,
                      'weights': weights, 'residual': residual,
                      'columns': columns})
        print(json.dumps({'loaded_tau': tau, 'order': order,
                          'baseline_l2': baseline_l2}), flush=True)
    rhs = np.concatenate(weighted_residuals)
    mat = np.vstack(weighted_columns)
    ridge = 1e-9
    def objective(x):
        error = rhs+mat@x
        return .5*(error@error+ridge*(x@x))
    def jac(x):
        return mat.T@(rhs+mat@x)+ridge*x
    constraints = [{'type': 'ineq', 'fun': lambda x: b-A@x,
                    'jac': lambda x: -A}]
    for block_mat, block_rhs in zip(weighted_columns, weighted_residuals):
        constraints.append({
            'type': 'ineq',
            'fun': lambda x, G=block_mat, r=block_rhs:
                1.-np.dot(r+G@x, r+G@x),
            'jac': lambda x, G=block_mat, r=block_rhs:
                -2*G.T@(r+G@x)})
    fit = minimize(objective, initial, jac=jac,
                   constraints=constraints, method='SLSQP',
                   options={'maxiter': 3000, 'ftol': 1e-12})
    amplitudes = fit.x
    def grid_metrics(x):
        rows = []
        for item in grids:
            residual = item['residual']+np.einsum('nik,k->ni', item['columns'], x)
            norm = np.linalg.norm(residual, axis=1)
            rows.append({'tau': item['tau'], 'order': item['order'],
                         'baseline_l2': item['baseline_l2'],
                         'corrected_l2': float(np.sqrt(item['weights']@norm**2)),
                         'corrected_max': float(np.max(norm))})
        return rows
    report = {'degrees': degrees, 'initial_amplitudes': initial.tolist(),
              'optimized_amplitudes': amplitudes.tolist(),
              'optimizer_success': bool(fit.success),
              'optimizer_message': fit.message,
              'maximum_constraint_violation': float(np.max(A@amplitudes-b)),
              'maximum_volume_ratio': max(row['corrected_l2']/row['baseline_l2']
                                          for row in grid_metrics(amplitudes)),
              'initial_objective': float(objective(initial)),
              'optimized_objective': float(objective(amplitudes)),
              'initial_rows': grid_metrics(initial),
              'optimized_rows': grid_metrics(amplitudes),
              'scope': 'Convex pressure-gradient least-squares objective on Gauss6 at three times and Gauss8/10 at center time, constrained by 63-node local cone inequalities and non-worsening L2 on each of the five grids. Velocity unchanged. No continuum, pressure-Poisson, or converged global momentum certificate.',
              'accepted': False}
    path = ROOT/'compact_potential'/'radial_pressure_volume_constrained.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({key: report[key] for key in
                      ('optimizer_success', 'maximum_constraint_violation',
                       'maximum_volume_ratio',
                       'initial_objective', 'optimized_objective',
                       'optimized_rows')}), flush=True)


def load_candidate():
    report = json.loads((ROOT/'compact_potential'/'radial_pressure_volume_constrained.json').read_text())
    if (not report['optimizer_success'] or
            report['maximum_constraint_violation'] > 1e-6 or
            report['maximum_volume_ratio'] > 1+1e-6):
        raise ValueError('No feasible volume-constrained candidate registered')
    return RadialPressureCandidate(load_dense_candidate(), report['degrees'],
                                   report['optimized_amplitudes'])


if __name__ == '__main__':
    run()
