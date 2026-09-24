"""Late-time incremental full-momentum fit at the remaining radial peak."""
import json

import numpy as np

from compact_potential import CompactPotentialField
from joined_field import ROOT
from joint_collar_fit import CollarMode
from joint_collar_scale import ScaledJointCollarField
from radial_swirl_fit import RadialSwirlRepairedField, linear_screen, metrics


class IncrementalRadialRepairField:
    """Callable finite-slab field with conservative late-time increment."""

    def __init__(self, current, amplitudes, reference_tau):
        self.current = current
        self.base = current.base
        self.nu = current.nu
        self.amplitudes = np.asarray(amplitudes, float)
        self.modes = [CollarMode(self.base, i, temporal_power=2.,
                                 reference_tau=reference_tau) for i in range(6)]

    def support(self, tau):
        return self.base.support(tau)

    def fields(self, points, tau):
        velocity, pressure = self.current.fields(points, tau)
        for amplitude, mode in zip(self.amplitudes, self.modes):
            if amplitude:
                du, dp = mode.fields(points, tau)
                velocity += amplitude*du
                pressure += amplitude*dp
        return velocity, pressure


def run():
    first = json.loads((ROOT/'compact_potential'/'joint_collar_fit.json').read_text())
    second = json.loads((ROOT/'compact_potential'/'radial_swirl_fit.json').read_text())
    base = CompactPotentialField()
    reference_tau = first['tau']
    joint = ScaledJointCollarField(.1*np.array(first['linear_fit_amplitudes']),
                                   reference_tau, base=base)
    current = RadialSwirlRepairedField(joint,
                                      second['linear_fit_amplitudes'], reference_tau)
    modes = [CollarMode(base, i, temporal_power=2., reference_tau=reference_tau)
             for i in range(6)]
    train = linear_screen(current, modes, .5/64, 8)
    _, weights, residual, columns, _, _ = train
    matrix = (np.sqrt(weights)[:, None, None]*columns).reshape(-1, 6)
    target = (-np.sqrt(weights)[:, None]*residual).reshape(-1)
    norms = np.linalg.norm(matrix, axis=0)
    scaled, *_ = np.linalg.lstsq(matrix/norms, target, rcond=1e-8)
    direction = scaled/norms
    rows = []
    for k, order in ((6., 8), (6., 6), (5.5, 6), (4., 6), (.4, 6)):
        tau = .5*2**(-k)
        data = train if (k == 6. and order == 8) else linear_screen(current, modes, tau, order)
        for strength in (0., .03, .1, .3, 1.):
            row = {'k': k, 'order': order, 'strength': strength,
                   'metrics': metrics(data, strength*direction)}
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = {'reference_tau': reference_tau, 'training_k': 6.,
              'temporal_power': 2., 'mode_names': CollarMode.names,
              'increment_amplitudes': direction.tolist(), 'rows': rows,
              'scope': 'Incremental six-mode solenoidal velocity/pressure fit on top of axial-collar plus radial-swirl candidate. Full nonlinear physical momentum; Gauss8 training and Gauss6 time/grid holdouts. Not a global norm certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'radial_repair_refit.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
