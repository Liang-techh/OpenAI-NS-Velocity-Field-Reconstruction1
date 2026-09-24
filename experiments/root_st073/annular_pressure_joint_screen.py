"""Fit pressure cone after scaling robust solenoidal annular correction."""
import json

import numpy as np

from annular_cone_scale_screen import CorrectionOnly
from annular_robust_fit import load_candidate as load_robust
from axial_pressure_cone_screen import AxialPressureMode, make_blocks
from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from radial_peak_cone import current_field
from transition_poloidal_pressure_screen import solve_pressure


class ScaledAnnularField:
    def __init__(self, base, correction, scale):
        self.base = base
        self.correction = correction
        self.scale = scale
        self.nu = base.nu

    def support(self, tau):
        return self.base.support(tau)

    def fields(self, points, tau):
        u, p = self.base.fields(points, tau)
        du, _ = self.correction.fields(points, tau)
        return u+self.scale*du, p


class AnnularPressureCandidate(ScaledAnnularField):
    def __init__(self, base, correction, scale, pressure_amplitudes):
        super().__init__(base, correction, scale)
        self.pressure_amplitudes = np.asarray(pressure_amplitudes)
        self.pressure_modes = [AxialPressureMode(base.base, degree)
                               for degree in range(len(self.pressure_amplitudes))]

    def fields(self, points, tau):
        u, p = super().fields(points, tau)
        for amplitude, mode in zip(self.pressure_amplitudes, self.pressure_modes):
            if amplitude:
                p += amplitude*mode.fields(points, tau)[1]
        return u, p


def load_candidate():
    report = json.loads((ROOT/'compact_potential'/'annular_pressure_joint_screen.json').read_text())
    if not report['linear_feasible']:
        raise ValueError('No registered joint pressure candidate')
    robust = load_robust()
    base = robust.base
    return AnnularPressureCandidate(base, CorrectionOnly(robust),
                                    report['velocity_scale'],
                                    report['pressure_amplitudes'])


def kinematic_residual(field, points, tau, hs, ht):
    u, J, part = kinematics(field, points, tau, hs, ht)
    return u, J, part+np.einsum('nij,nj->ni', J, u)


def run():
    base = current_field()
    robust = load_robust()
    scale = .4
    field = ScaledAnnularField(base, CorrectionOnly(robust), scale)
    tau = .5/64
    radius = .0056890761915166545
    heights = tuple(np.linspace(.0025, .003432627453438968, 13))
    blocks = make_blocks(base, tau, radius, heights)
    points = np.vstack([block[3] for block in blocks])
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    line_u, line_J, line_R = kinematic_residual(field, points, tau, hs, ht)
    pressure_modes = [AxialPressureMode(base.base, i) for i in range(4)]
    line_columns = np.stack([kinematics(mode, points, tau, hs, ht)[2]
                             for mode in pressure_modes], axis=-1)
    volume = []
    for order in (6, 8, 10):
        grid, weights = nodes(base, tau, order)
        _, _, residual = kinematic_residual(field, grid, tau, hs, ht)
        columns = np.stack([kinematics(mode, grid, tau, hs, ht)[2]
                            for mode in pressure_modes], axis=-1)
        baseline_l2 = float(np.sqrt(weights @ np.linalg.norm(residual, axis=1)**2))
        volume.append({'order': order, 'weights': weights,
                       'residual': residual, 'columns': columns,
                       'scaled_velocity_l2': baseline_l2})
        print(json.dumps({'loaded_order': order,
                          'scaled_velocity_l2': baseline_l2}), flush=True)
    objective_residual = np.vstack([data['residual'] for data in volume])
    objective_columns = np.concatenate([data['columns'] for data in volume], axis=0)
    objective_weights = np.concatenate([data['weights']/data['scaled_velocity_l2']**2
                                        for data in volume])
    fit = solve_pressure(line_u, line_J, line_R, line_columns,
                         blocks, radius, ratio=.8,
                         objective_residual=objective_residual,
                         objective_columns=objective_columns,
                         objective_weights=objective_weights)
    rows = []
    if fit is not None:
        for data in volume:
            residual = data['residual']+np.einsum(
                'nik,k->ni', data['columns'], fit['pressure_amplitudes'])
            norms = np.linalg.norm(residual, axis=1)
            rows.append({'order': data['order'],
                         'scaled_velocity_l2': data['scaled_velocity_l2'],
                         'joint_l2': float(np.sqrt(data['weights'] @ norms**2)),
                         'joint_max': float(np.max(norms))})
    report = {'tau': tau, 'velocity_scale': scale, 'radius': radius,
              'heights': heights, 'linear_feasible': fit is not None,
              'pressure_amplitudes': fit['pressure_amplitudes'] if fit else None,
              'local_cone_pass_count': fit['pass_count'] if fit else 0,
              'local_max_cone_ratio': fit['max_cone_ratio'] if fit else None,
              'volume_rows': rows,
              'scope': 'Scaled solenoidal annular velocity plus four compact pressure modes. Local cone constrained at 13 axial points; pressure fit to Gauss6/8/10 relative volume residual at one time. No radial-time wave or continuum certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_pressure_joint_screen.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({key: report[key] for key in
                      ('linear_feasible', 'pressure_amplitudes',
                       'local_cone_pass_count', 'local_max_cone_ratio',
                       'volume_rows')}), flush=True)


if __name__ == '__main__':
    run()
