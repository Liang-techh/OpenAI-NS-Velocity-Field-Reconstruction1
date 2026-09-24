"""Direct midpoint residual of the first evolved exact-curl patch interval.

The velocity interpolates potential coefficients linearly between the first
two collocation nodes. This tests the actual continuous-time field between
nodes, rather than only the instantaneous projected endpoint residual.
"""
import json

import numpy as np

from curl_wave_patch_evolution import FrozenPotentialField, polynomial_data
from curl_wave_prototype import LocalizedCurlWave
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


class FirstIntervalField:
    def __init__(self, base, wave, amplitude, stage):
        self.base, self.wave, self.amplitude = base, wave, amplitude
        self.nu = base.nu
        self.tau0 = stage['tau']
        self.harmonic_slopes = [np.array([complex(re, im) for re, im in row])
                                for row in stage['harmonic_slope_coefficients']]
        self.mean_slope = np.array(stage['mean_slope_coefficients'])

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        if not np.all(ts == ts[0]):
            raise ValueError('Evaluate one time per spatial batch')
        t = float(ts[0])
        interval = t-self.tau0
        frozen = FrozenPotentialField(self.base, self.wave, self.amplitude,
                                      [interval*s[:27] for s in self.harmonic_slopes],
                                      interval*self.mean_slope[:18])
        velocity, pressure = frozen.fields(pts, t)
        for i, (x, y, z) in enumerate(pts):
            r = np.hypot(x, y)
            if r == 0 or abs(r-self.wave.radius) >= self.wave.radial_halfwidth or abs(z-self.wave.zcenter) >= self.wave.axial_halfwidth:
                continue
            theta = np.arctan2(y, x)
            for mode, slope in zip(self.wave.waves, self.harmonic_slopes):
                m = mode['m']
                kr, _, kz = mode['normal']
                p = polynomial_data(self.wave, r, z, slope[27:], 1)[0][0]
                carrier = np.exp(1j*(m*theta+kr*(r-self.wave.radius)
                                     +kz*(z-self.wave.zcenter)
                                     +mode['omega']*(t-self.wave.tau0)))
                pressure[i] += (carrier*p).real
            pressure[i] += polynomial_data(self.wave, r, z,
                                           self.mean_slope[18:], 1)[0][0]
        return velocity, pressure


def residual_metrics(field, points, tau):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, grad, part = kinematics(field, points, tau, hs, ht)
    residual = part+np.einsum('nij,nj->ni', grad, u)
    norm = np.linalg.norm(residual, axis=1)
    return {'max_momentum': float(norm.max()),
            'rms_momentum': float(np.sqrt(np.mean(norm**2))),
            'max_fd_divergence': float(np.max(np.abs(np.trace(grad,
                                                             axis1=1,
                                                             axis2=2)))),
            'max_speed': float(np.max(np.linalg.norm(u, axis=1)))}


def run():
    evolution = json.loads((ROOT/'compact_potential'/'curl_wave_patch_evolution.json').read_text())
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    base = current_field()
    wave = LocalizedCurlWave(source)
    field = FirstIntervalField(base, wave, evolution['amplitude'],
                               evolution['stages'][0])
    tau = wave.tau0+evolution['time_step']/2
    angles = np.arange(8)*2*np.pi/8
    axis = np.array([-.45, -.15, .15, .45])
    points = np.vstack([np.column_stack((r*np.cos(angles), r*np.sin(angles),
                                         np.full(len(angles), z)))
                        for xi in axis for eta in axis
                        for r, z in [(wave.radius+wave.radial_halfwidth*xi,
                                      wave.zcenter+wave.axial_halfwidth*eta)]])
    report = {'tau': tau, 'time_step': evolution['time_step'],
              'heldout_spatial_nodes': len(axis)**2, 'angles_per_node': len(angles),
              'baseline': residual_metrics(base, points, tau),
              'evolved_first_interval': residual_metrics(field, points, tau),
              'scope': 'Direct finite-difference full momentum at the midpoint of the first piecewise-linear potential-coefficient interval. This is a spatially held-out sample, not an inter-node/global bound or a matched pulse.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_patch_trajectory.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
