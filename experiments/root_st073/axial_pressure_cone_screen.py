"""Screen independent compact axial-pressure shapes against the local cone.

Pressure leaves the solenoidal velocity unchanged and enters the complete
momentum residual linearly. This is a one-time, one-radius design experiment.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss, Legendre
from scipy.optimize import linprog, minimize

from curl_wave_cone_parameter_screen import cone_row
from joined_field import ROOT, coordinates
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


class AxialPressureMode:
    def __init__(self, base, degree):
        self.base = base
        self.nu = base.nu
        self.degree = degree
        self.poly = Legendre.basis(degree)

    def fields(self, points, tau):
        points = np.asarray(points)
        times = np.broadcast_to(tau, (len(points),))
        p = np.zeros(len(points))
        for i, (point, t) in enumerate(zip(points, times)):
            radius = np.hypot(point[0], point[1])
            _, rsupp, zflat, zsupp = self.base.support(t)
            zabs = abs(point[2])
            if radius >= rsupp or zabs <= zflat or zabs >= zsupp:
                continue
            y = radius / rsupp
            s = (zabs-zflat)/(zsupp-zflat)
            p[i] = ((1-y*y)**5 * 256*s**4*(1-s)**4
                    * self.poly(2*s-1))
        return np.zeros_like(points), p


class AxialPressureCandidate:
    """Finite-slab field with the screened pressure and unchanged velocity."""

    def __init__(self, base, amplitudes):
        self.base = base
        self.nu = base.nu
        self.amplitudes = np.asarray(amplitudes)
        self.modes = [AxialPressureMode(base.base, degree)
                      for degree in range(len(self.amplitudes))]

    def support(self, tau):
        return self.base.support(tau)

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        for amplitude, mode in zip(self.amplitudes, self.modes):
            if amplitude:
                pressure += amplitude*mode.fields(points, tau)[1]
        return velocity, pressure


def load_candidate():
    report = json.loads((ROOT/'compact_potential'/'axial_pressure_cone_screen.json').read_text())
    if not report['linear_feasible'] or not report['optimizer_success']:
        raise ValueError('No optimized axial-pressure candidate is registered')
    return AxialPressureCandidate(current_field(), report['amplitudes'])


def make_blocks(field, tau, radius, heights):
    nodes, weights = leggauss(10)
    blocks = []
    for z in heights:
        q = float(coordinates(0., z/np.sqrt(field.nu), tau,
                              field.base.base.inner.h)['q'])
        inner_radius = np.sqrt(field.nu)*np.sqrt(2*q*3/64)
        edges = sorted(set([0., inner_radius, radius]))
        qr, qw = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            qr.extend((lo+hi)/2+(hi-lo)/2*nodes)
            qw.extend((hi-lo)/2*weights)
        qr, qw = np.asarray(qr), np.asarray(qw)
        points = np.column_stack((np.r_[qr, radius],
                                  np.zeros(len(qr)+1),
                                  np.full(len(qr)+1, z)))
        blocks.append((z, qr, qw, points))
    return blocks


def evaluate_rows(field, modes, amplitudes, tau, radius, heights, hs, ht):
    blocks = make_blocks(field, tau, radius, heights)
    points = np.vstack([block[3] for block in blocks])
    u, J, part = kinematics(field, points, tau, hs, ht)
    residual = part+np.einsum('nij,nj->ni', J, u)
    for amplitude, mode in zip(amplitudes, modes):
        if amplitude:
            residual += amplitude*kinematics(mode, points, tau, hs, ht)[2]
    rows, offset = [], 0
    for z, qr, qw, block_points in blocks:
        sl = slice(offset, offset+len(block_points))
        offset += len(block_points)
        rows.append(cone_row(radius, z, u[sl][-1], J[sl][-1],
                             residual[sl], qr, qw))
    return rows


def run():
    field = current_field()
    tau = .5/64
    radius = .0056890761915166545
    heights = tuple(np.linspace(.0025, .003432627453438968, 13))
    blocks = make_blocks(field, tau, radius, heights)
    points = np.vstack([block[3] for block in blocks])
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, J, part = kinematics(field, points, tau, hs, ht)
    residual0 = part + np.einsum('nij,nj->ni', J, u)
    modes = [AxialPressureMode(field.base, degree) for degree in range(4)]
    columns = np.stack([kinematics(mode, points, tau, hs, ht)[2]
                        for mode in modes], axis=-1)
    target0, target_columns, cone_data = [], [], []
    offset = 0
    for z, qr, qw, block_points in blocks:
        sl = slice(offset, offset+len(block_points))
        offset += len(block_points)
        def target(residual):
            return np.array([-np.dot(qw*qr**2, residual[:-1, 1])/radius**2,
                             -np.dot(qw*qr, residual[:-1, 2])/radius])
        t0 = target(residual0[sl])
        tc = np.stack([target(columns[sl, :, j])
                       for j in range(len(modes))], axis=-1)
        F = u[sl][-1, 1]/radius
        shear = np.array([J[sl][-1, 1, 0]-F, J[sl][-1, 2, 0]])
        N = shear/np.linalg.norm(shear)
        K = np.array([-N[1], N[0]])
        lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
        c = np.sqrt(lam2)/(2*F*N[0])
        target0.append(t0)
        target_columns.append(tc)
        cone_data.append((N, K, c))
    target_cone_ratio = .8
    A, b = [], []
    for t0, tc, (N, K, c) in zip(target0, target_columns, cone_data):
        for sign in (-1, 1):
            direction = target_cone_ratio*N+sign*c*K
            A.append(direction @ tc)
            b.append(-1.-direction @ t0)
    A, b = np.asarray(A), np.asarray(b)
    scale = 1/np.maximum(np.max(np.abs(A), axis=0), 1e-12)
    feasibility = linprog(np.zeros(len(modes)), A_ub=A*scale,
                          b_ub=b, bounds=[(None, None)]*len(modes),
                          method='highs')
    report = {'tau': tau, 'radius': radius, 'heights': heights,
              'target_cone_ratio': target_cone_ratio,
              'mode_degrees': list(range(len(modes))),
              'linear_feasible': bool(feasibility.success),
              'linear_solver_message': feasibility.message,
              'scope': 'Independent compact pressure correction only; velocity unchanged. Local cone constrained at five axial points; weighted residual sampled on radial quadrature at one time/radius. No global NS or pressure-Poisson certificate.',
              'accepted': False}
    if feasibility.success:
        base_norm = np.linalg.norm(residual0)
        mat = (columns*scale[None, None, :]).reshape(-1, len(modes))/base_norm
        rhs = residual0.reshape(-1)/base_norm
        def objective(v):
            e = rhs+mat @ v
            return .5*(e @ e + 1e-8*(v @ v))
        def jac(v):
            return mat.T @ (rhs+mat @ v)+1e-8*v
        solution = minimize(objective, feasibility.x, jac=jac,
                            constraints=[{'type': 'ineq',
                                          'fun': lambda v: b-A @ (scale*v),
                                          'jac': lambda v: -A*scale}],
                            method='SLSQP',
                            options={'maxiter': 1000, 'ftol': 1e-12})
        v = solution.x if np.max(A @ (scale*solution.x)-b) <= 1e-5 else feasibility.x
        amplitudes = scale*v
        residual = residual0+np.einsum('nik,k->ni', columns, amplitudes)
        rows0, rows1 = [], []
        offset = 0
        for z, qr, qw, block_points in blocks:
            sl = slice(offset, offset+len(block_points))
            offset += len(block_points)
            rows0.append(cone_row(radius, z, u[sl][-1], J[sl][-1],
                                  residual0[sl], qr, qw))
            rows1.append(cone_row(radius, z, u[sl][-1], J[sl][-1],
                                  residual[sl], qr, qw))
        report.update({'amplitudes': amplitudes.tolist(),
                       'optimizer_success': bool(solution.success),
                       'optimizer_message': solution.message,
                       'baseline': {'pass_count': sum(x['strict_local_pass'] for x in rows0),
                                    'max_center_residual': max(x['center_residual_norm'] for x in rows0),
                                    'sample_residual_rms': float(np.sqrt(np.mean(residual0**2)))},
                       'candidate': {'pass_count': sum(x['strict_local_pass'] for x in rows1),
                                     'max_center_residual': max(x['center_residual_norm'] for x in rows1),
                                     'sample_residual_rms': float(np.sqrt(np.mean(residual**2)))},
                       'rows': rows1})
        holdout_heights = tuple((left+right)/2
                                for left, right in zip(heights[:-1], heights[1:]))
        holdout_rows = evaluate_rows(field, modes, amplitudes, tau, radius,
                                     holdout_heights, hs, ht)
        report['holdout'] = {'heights': holdout_heights,
                             'pass_count': sum(row['strict_local_pass'] for row in holdout_rows),
                             'max_cone_ratio': max(row['cone_ratio'] for row in holdout_rows),
                             'rows': holdout_rows}
    out = ROOT/'compact_potential'/'axial_pressure_cone_screen.json'
    out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report.get(key) for key in
                      ('linear_feasible', 'amplitudes', 'baseline',
                       'candidate', 'optimizer_success')})
          , flush=True)


if __name__ == '__main__':
    run()
