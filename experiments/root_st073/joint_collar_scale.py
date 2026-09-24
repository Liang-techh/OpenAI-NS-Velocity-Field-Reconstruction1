"""Screen a growing time weight for the fitted joint collar correction."""
import json

import numpy as np

from compact_potential import CompactPotentialField
from joined_field import ROOT
from joint_collar_fit import CollarMode, metrics, screen


class ScaledJointCollarField:
    """Callable experimental field retaining the compact baseline support."""

    def __init__(self, amplitudes, reference_tau, temporal_power=1.5, base=None):
        self.base = base or CompactPotentialField()
        self.nu = self.base.nu
        self.amplitudes = np.asarray(amplitudes, float)
        self.modes = [CollarMode(self.base, i, temporal_power, reference_tau)
                      for i in range(6)]

    def support(self, tau):
        return self.base.support(tau)

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        for amplitude, mode in zip(self.amplitudes, self.modes):
            if amplitude:
                du, dp = mode.fields(points, tau)
                velocity += amplitude*du
                pressure += amplitude*dp
        return velocity, pressure


def run():
    fit = json.loads((ROOT/'compact_potential'/'joint_collar_fit.json').read_text())
    direction = np.array(fit['linear_fit_amplitudes'])
    reference_tau = fit['tau']
    strength = .1
    temporal_power = 1.5
    base = CompactPotentialField()
    rows = []
    for k in (.4, 4., 5.5, 6.):
        tau = .5*2**(-k)
        data = screen(base, tau, 6, temporal_power, reference_tau)
        row = {
            'k': k, 'tau': tau,
            'amplitude_factor': float((reference_tau/tau)**temporal_power),
            'baseline': metrics(data, np.zeros(6)),
            'candidate': metrics(data, strength*direction),
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = {
        'reference_tau': reference_tau, 'strength': strength,
        'temporal_power': temporal_power, 'rows': rows,
        'scope': 'The same six fitted coefficients multiplied by (reference_tau/tau)^1.5. Time derivative of the multiplier is included by fourth-order physical-time finite differences. Independent Gauss6 grid at four registered times; sampled screen only.',
        'accepted': False,
    }
    out = ROOT/'compact_potential'/'joint_collar_scale.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
