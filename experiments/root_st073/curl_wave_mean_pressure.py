"""Local compact mean-pressure proxy for the curl wave's quadratic radial term.

This cancels the mean radial derivative at one chosen point only. It is not the
paper's global pressure primitive or a momentum solution.
"""
import json

import numpy as np

from compact_control import cutoff
from curl_wave_prototype import LocalizedCurlWave, bump, cylindrical_residual
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


class PressureProxyField:
    def __init__(self, base, wave, amplitude, radial_mean_coefficient):
        self.base = base
        self.wave = wave
        self.amplitude = amplitude
        self.coefficient = radial_mean_coefficient
        self.nu = base.nu

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity, pressure = self.base.fields(pts, ts)
        delta, _ = self.wave.fields(pts, ts)
        velocity += self.amplitude*delta
        pressure = pressure.copy()
        for i, (point, t) in enumerate(zip(pts, ts)):
            r = np.hypot(point[0], point[1])
            z = point[2]
            br, _ = bump(r, self.wave.radius, self.wave.radial_halfwidth)
            bz, _ = bump(z, self.wave.zcenter, self.wave.axial_halfwidth)
            if br and bz:
                time_cut = float(cutoff((t-.012)/.008)[0])
                pressure[i] -= (self.amplitude**2*self.coefficient
                                *(r-self.wave.radius)*(br*bz*time_cut)**2)
        return velocity, pressure


def run():
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    amplitudes = json.loads((ROOT/'compact_potential'/'curl_wave_amplitude.json').read_text())
    coefficient = amplitudes['mean_quadratic_cylindrical'][0]
    base = current_field()
    wave = LocalizedCurlWave(source)
    tau = source['tau']
    radius, _, z = source['point']
    angles = np.arange(16)*2*np.pi/16
    points = np.column_stack((radius*np.cos(angles), radius*np.sin(angles),
                              np.full(len(angles), z)))
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    rows = []
    for amplitude in (.1, 1.):
        field = PressureProxyField(base, wave, amplitude, coefficient)
        u, grad, part = kinematics(field, points, tau, hs, ht)
        residual = part + np.einsum('nij,nj->ni', grad, u)
        cylindrical = cylindrical_residual(residual, points)
        norms = np.linalg.norm(cylindrical, axis=1)
        rows.append({'amplitude': amplitude,
                     'covariance_fraction': amplitude**2,
                     'mean_residual_cylindrical': cylindrical.mean(axis=0).tolist(),
                     'max_momentum': float(norms.max()),
                     'rms_momentum': float(np.sqrt(np.mean(norms**2)))})
    report = {'tau': tau, 'point': source['point'],
              'local_radial_pressure_gradient': -coefficient,
              'rows': rows,
              'scope': 'Compact pressure bubble whose radial derivative cancels wave quadratic mean radial residual at one point. No exact global pressure reconstruction or amplitude equation. Full nonlinear physical momentum at 16 angles.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_mean_pressure.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
