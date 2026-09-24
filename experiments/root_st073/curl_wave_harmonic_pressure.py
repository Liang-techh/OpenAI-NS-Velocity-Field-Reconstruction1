"""Projected local harmonic pressure for the exact-curl periodic wave.

At the selected point it removes the phase-normal component of the linear
nonzero Fourier residual. Transverse amplitude dynamics remain unsolved.
"""
import json

import numpy as np

from compact_control import cutoff
from curl_wave_mean_pressure import PressureProxyField
from curl_wave_prototype import LocalizedCurlWave, bump, cylindrical_residual
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


class HarmonicPressureField(PressureProxyField):
    def __init__(self, base, wave, amplitude, mean_coefficient, harmonic_coefficients):
        super().__init__(base, wave, amplitude, mean_coefficient)
        self.harmonic_coefficients = harmonic_coefficients

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity, pressure = super().fields(pts, ts)
        for i, (point, t) in enumerate(zip(pts, ts)):
            x, ycart, z = point
            r = np.hypot(x, ycart)
            if r == 0:
                continue
            br, _ = bump(r, self.wave.radius, self.wave.radial_halfwidth)
            bz, _ = bump(z, self.wave.zcenter, self.wave.axial_halfwidth)
            if not br or not bz:
                continue
            envelope = br*bz*float(cutoff((t-.012)/.008)[0])
            theta = np.arctan2(ycart, x)
            for coeff, mode in zip(self.harmonic_coefficients, self.wave.waves):
                kr, _, kz = mode['normal']
                phase = (mode['m']*theta + kr*(r-self.wave.radius)
                         + kz*(z-self.wave.zcenter)
                         + mode['omega']*(t-self.wave.tau0))
                pressure[i] += self.amplitude*envelope*np.real(coeff*np.exp(1j*phase))
        return velocity, pressure


def full_operator(field, points, tau):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, grad, part = kinematics(field, points, tau, hs, ht)
    return part + np.einsum('nij,nj->ni', grad, u)


def run():
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    previous = json.loads((ROOT/'compact_potential'/'curl_wave_amplitude.json').read_text())
    base = current_field()
    wave = LocalizedCurlWave(source)
    tau = source['tau']
    radius, _, z = source['point']
    angles = np.arange(16)*2*np.pi/16
    points = np.column_stack((radius*np.cos(angles), radius*np.sin(angles),
                              np.full(len(angles), z)))
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    u, grad, _ = kinematics(base, points, tau, hs, ht)
    w, wgrad, wpart = kinematics(wave, points, tau, hs, ht)
    linear = (wpart + np.einsum('nij,nj->ni', grad, w)
              + np.einsum('nij,nj->ni', wgrad, u))
    linear_cyl = cylindrical_residual(linear, points)
    coefficients = []
    projections = []
    for mode in wave.waves:
        m = mode['m']
        harmonic = 2*np.mean(linear_cyl*np.exp(-1j*m*angles)[:, None], axis=0)
        normal = mode['normal']
        coefficient = 1j*np.dot(normal, harmonic)/np.dot(normal, normal)
        longitudinal = normal*np.dot(normal, harmonic)/np.dot(normal, normal)
        coefficients.append(coefficient)
        projections.append({'m': m,
                            'linear_harmonic_magnitude': float(np.linalg.norm(harmonic)),
                            'longitudinal_fraction': float(np.linalg.norm(longitudinal)/np.linalg.norm(harmonic)),
                            'pressure_coefficient': [float(coefficient.real), float(coefficient.imag)]})
    rows = []
    for amplitude, mean_coefficient in ((1., 0.), (1., previous['mean_quadratic_cylindrical'][0]),
                                        (.1, previous['mean_quadratic_cylindrical'][0])):
        field = HarmonicPressureField(base, wave, amplitude,
                                     mean_coefficient, coefficients)
        residual = cylindrical_residual(full_operator(field, points, tau), points)
        norms = np.linalg.norm(residual, axis=1)
        rows.append({'amplitude': amplitude, 'with_mean_pressure': bool(mean_coefficient),
                     'mean_residual_cylindrical': residual.mean(axis=0).tolist(),
                     'max_momentum': float(norms.max()),
                     'rms_momentum': float(np.sqrt(np.mean(norms**2)))})
    report = {'tau': tau, 'point': source['point'],
              'harmonic_projections': projections, 'rows': rows,
              'scope': 'Local pressure harmonic chosen from phase-normal projection of the full linear wave residual. Compact pressure follows each integer angular phase and envelope. Full physical nonlinear momentum sampled at 16 angles. Transverse amplitude equation, global pressure primitive, and support moment conditions remain unsolved.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_harmonic_pressure.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
