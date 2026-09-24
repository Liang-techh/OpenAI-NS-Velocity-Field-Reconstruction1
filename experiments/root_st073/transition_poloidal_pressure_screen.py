"""Joint finite-sample screen of a solenoidal collar mode and pressure shapes."""
import json

import numpy as np
from scipy.optimize import linprog, minimize

from axial_pressure_cone_screen import AxialPressureMode, make_blocks
from curl_wave_cone_parameter_screen import cone_row
from joined_field import ROOT
from joint_collar_fit import TransitionPoloidalMode, kinematics
from radial_peak_cone import current_field


class TransitionPressureCandidate:
    def __init__(self, base, poloidal_amplitude, pressure_amplitudes,
                 reference_tau=.5/64):
        self.base = base
        self.nu = base.nu
        self.poloidal_amplitude = poloidal_amplitude
        self.pressure_amplitudes = np.asarray(pressure_amplitudes)
        compact_base = base.base
        self.poloidal = TransitionPoloidalMode(compact_base, 2., reference_tau)
        self.pressure_modes = [AxialPressureMode(compact_base, i)
                               for i in range(len(self.pressure_amplitudes))]

    def support(self, tau):
        return self.base.support(tau)

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        if self.poloidal_amplitude:
            velocity += self.poloidal_amplitude*self.poloidal.fields(points, tau)[0]
        for amplitude, mode in zip(self.pressure_amplitudes, self.pressure_modes):
            if amplitude:
                pressure += amplitude*mode.fields(points, tau)[1]
        return velocity, pressure


def load_candidate():
    report = json.loads((ROOT/'compact_potential'/'transition_poloidal_pressure_screen.json').read_text())
    best = report['best']
    return TransitionPressureCandidate(current_field(), best['poloidal_amplitude'],
                                       best['pressure_amplitudes'])


def target(residual, radius, qr, qw):
    return np.array([-np.dot(qw*qr**2, residual[:-1, 1])/radius**2,
                     -np.dot(qw*qr, residual[:-1, 2])/radius])


def solve_pressure(u, J, residual, pressure_columns, blocks, radius, ratio=.8,
                   objective_residual=None, objective_columns=None,
                   objective_weights=None):
    A, b = [], []
    offset = 0
    for z, qr, qw, block_points in blocks:
        sl = slice(offset, offset+len(block_points))
        offset += len(block_points)
        t0 = target(residual[sl], radius, qr, qw)
        tc = np.stack([target(pressure_columns[sl, :, j], radius, qr, qw)
                       for j in range(pressure_columns.shape[-1])], axis=-1)
        F = u[sl][-1, 1]/radius
        shear = np.array([J[sl][-1, 1, 0]-F, J[sl][-1, 2, 0]])
        N = shear/np.linalg.norm(shear)
        K = np.array([-N[1], N[0]])
        lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
        if lam2 <= 0:
            return None
        c = np.sqrt(lam2)/(2*F*N[0])
        for sign in (-1, 1):
            direction = ratio*N+sign*c*K
            A.append(direction @ tc)
            b.append(-1.-direction @ t0)
    A, b = np.asarray(A), np.asarray(b)
    scale = 1/np.maximum(np.max(np.abs(A), axis=0), 1e-12)
    feasibility = linprog(np.zeros(pressure_columns.shape[-1]),
                          A_ub=A*scale, b_ub=b,
                          bounds=[(None, None)]*pressure_columns.shape[-1],
                          method='highs')
    if not feasibility.success:
        return None
    if objective_residual is None:
        objective_residual = residual
        objective_columns = pressure_columns
        objective_weights = np.ones(len(residual))
    root_weights = np.sqrt(objective_weights)
    weighted_residual = root_weights[:, None]*objective_residual
    weighted_columns = root_weights[:, None, None]*objective_columns
    norm = np.linalg.norm(weighted_residual)
    mat = (weighted_columns*scale[None, None, :]).reshape(-1, len(scale))/norm
    rhs = weighted_residual.reshape(-1)/norm
    def objective(v):
        e = rhs+mat @ v
        return .5*(e @ e+1e-8*(v @ v))
    def jac(v):
        return mat.T @ (rhs+mat @ v)+1e-8*v
    solution = minimize(objective, feasibility.x, jac=jac,
                        constraints=[{'type': 'ineq',
                                      'fun': lambda v: b-A @ (scale*v),
                                      'jac': lambda v: -A*scale}],
                        method='SLSQP', options={'maxiter': 1000, 'ftol': 1e-12})
    v = solution.x if np.max(A @ (scale*solution.x)-b) <= 1e-5 else feasibility.x
    amplitudes = scale*v
    corrected = residual+np.einsum('nik,k->ni', pressure_columns, amplitudes)
    corrected_objective = (objective_residual
                           + np.einsum('nik,k->ni', objective_columns, amplitudes))
    rows, offset = [], 0
    for z, qr, qw, block_points in blocks:
        sl = slice(offset, offset+len(block_points))
        offset += len(block_points)
        rows.append(cone_row(radius, z, u[sl][-1], J[sl][-1],
                             corrected[sl], qr, qw))
    return {'pressure_amplitudes': amplitudes.tolist(),
            'pass_count': sum(row['strict_local_pass'] for row in rows),
            'max_cone_ratio': max(row['cone_ratio'] for row in rows),
            'sample_residual_rms': float(np.sqrt(np.mean(corrected**2))),
            'objective_weighted_l2': float(np.linalg.norm(
                root_weights[:, None]*corrected_objective)),
            'max_center_residual': max(row['center_residual_norm'] for row in rows),
            'max_poloidal_speed': None,
            'rows': rows}


def run():
    field = current_field()
    tau = .5/64
    radius = .0056890761915166545
    heights = tuple(np.linspace(.0025, .003432627453438968, 13))
    blocks = make_blocks(field, tau, radius, heights)
    points = np.vstack([block[3] for block in blocks])
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u0, J0, part0 = kinematics(field, points, tau, hs, ht)
    poloidal = TransitionPoloidalMode(field.base, 2., tau)
    du, dJ, dpart = kinematics(poloidal, points, tau, hs, ht)
    pressure_modes = [AxialPressureMode(field.base, i) for i in range(4)]
    pressure_columns = np.stack([kinematics(mode, points, tau, hs, ht)[2]
                                 for mode in pressure_modes], axis=-1)
    amplitudes_to_screen = (0., .002, .004, .006, .008, .010,
                            .012, .014, .016, .018, .020)
    results = []
    for amplitude in amplitudes_to_screen:
        u = u0+amplitude*du
        J = J0+amplitude*dJ
        part = part0+amplitude*dpart
        residual = part+np.einsum('nij,nj->ni', J, u)
        result = solve_pressure(u, J, residual, pressure_columns,
                                blocks, radius)
        if result is not None:
            result['poloidal_amplitude'] = amplitude
            result['max_poloidal_speed'] = float(np.max(np.linalg.norm(amplitude*du, axis=1)))
            results.append(result)
    results.sort(key=lambda row: (row['sample_residual_rms'],
                                  row['max_center_residual']))
    report = {'tau': tau, 'radius': radius, 'heights': heights,
              'screened_poloidal_amplitudes': amplitudes_to_screen,
              'target_cone_ratio': .8,
              'results': results,
              'best': results[0] if results else None,
              'scope': 'One time and one radius. Complete nonlinear momentum recombination for the poloidal velocity; pressure changes residual linearly. No global NS or pressure-Poisson certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'transition_poloidal_pressure_screen.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'feasible': len(results),
                      'best': {k: results[0][k] for k in
                               ('poloidal_amplitude', 'pressure_amplitudes',
                                'pass_count', 'max_cone_ratio',
                                'sample_residual_rms', 'max_center_residual')}
                      if results else None}), flush=True)


if __name__ == '__main__':
    run()
