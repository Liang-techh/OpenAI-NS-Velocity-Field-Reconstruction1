"""One-step transverse amplitude evolution for the compact periodic wave.

The time-slope vector potential cancels the selected-point linear Fourier
residual at tau0. It is a local Taylor step, not the paper's pulse ODE solve.
"""
import json

import numpy as np

from compact_control import cutoff
from curl_wave_harmonic_pressure import HarmonicPressureField, full_operator
from curl_wave_prototype import LocalizedCurlWave, bump, cylindrical_residual
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


def curl_matrix(normal, radius):
    matrix = np.column_stack([1j*np.cross(normal, np.eye(3)[:, i])
                              for i in range(3)])
    matrix[2, 1] += 1/radius
    return matrix


class TaylorAmplitudeField:
    def __init__(self, pressure_field, corrections):
        self.pressure_field = pressure_field
        self.wave = pressure_field.wave
        self.nu = pressure_field.nu
        self.corrections = corrections

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity, pressure = self.pressure_field.fields(pts, ts)
        for i, (point, t) in enumerate(zip(pts, ts)):
            s = t-self.wave.tau0
            if s < 0 or s >= .0001:
                continue
            x, ycart, z = point
            r = np.hypot(x, ycart)
            if r == 0:
                continue
            br, br_r = bump(r, self.wave.radius, self.wave.radial_halfwidth)
            bz, bz_z = bump(z, self.wave.zcenter, self.wave.axial_halfwidth)
            if not br or not bz:
                continue
            taper = float(cutoff((s-.00005)/.00005)[0])
            envelope = br*bz*taper
            er = br_r*bz*taper
            ez = br*bz_z*taper
            theta = np.arctan2(ycart, x)
            radial = angular = axial = 0.
            for mode, coeff in zip(self.wave.waves, self.corrections):
                cr, ct, cz = coeff
                m = mode['m']
                kr, _, kz = mode['normal']
                phase = (m*theta + kr*(r-self.wave.radius)
                         + kz*(z-self.wave.zcenter)
                         + mode['omega']*(t-self.wave.tau0))
                carrier = np.exp(1j*phase)
                wr = (1j*(m*cz/r-kz*ct)*envelope -ct*ez)*carrier
                wt = (1j*(kz*cr-kr*cz)*envelope +(cr*ez-cz*er))*carrier
                wz = (1j*(kr*ct-m*cr/r)*envelope +ct*(er+envelope/r))*carrier
                radial += wr.real
                angular += wt.real
                axial += wz.real
            ca, sa = x/r, ycart/r
            scale = self.pressure_field.amplitude*s
            velocity[i] += scale*np.array([radial*ca-angular*sa,
                                           radial*sa+angular*ca, axial])
        return velocity, pressure


def run():
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    amplitude_report = json.loads((ROOT/'compact_potential'/'curl_wave_amplitude.json').read_text())
    harmonic_report = json.loads((ROOT/'compact_potential'/'curl_wave_harmonic_pressure.json').read_text())
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
    pressure_coeffs, correction_coeffs, mode_rows = [], [], []
    for mode in wave.waves:
        m = mode['m']
        normal = mode['normal']
        harmonic = 2*np.mean(linear_cyl*np.exp(-1j*m*angles)[:, None], axis=0)
        pressure = 1j*np.dot(normal, harmonic)/np.dot(normal, normal)
        transverse = harmonic+1j*normal*pressure
        matrix = curl_matrix(normal, radius)
        correction, *_ = np.linalg.lstsq(matrix, transverse, rcond=None)
        pressure_coeffs.append(pressure)
        correction_coeffs.append(correction)
        mode_rows.append({'m': m,
                          'transverse_target_magnitude': float(np.linalg.norm(transverse)),
                          'curl_replay_relative_error': float(np.linalg.norm(matrix@correction-transverse)/np.linalg.norm(transverse)),
                          'potential_time_slope': [[float(v.real), float(v.imag)] for v in correction]})
    mean_coeff = amplitude_report['mean_quadratic_cylindrical'][0]
    rows = []
    direct = {}
    for amplitude in (.1, 1.):
        pressure_field = HarmonicPressureField(base, wave, amplitude,
                                               mean_coeff, pressure_coeffs)
        taylor_field = TaylorAmplitudeField(pressure_field, correction_coeffs)
        for name, field in (('pressure_only', pressure_field),
                            ('taylor_amplitude', taylor_field)):
            residual = cylindrical_residual(full_operator(field, points, tau), points)
            direct[(amplitude, name)] = residual
            norms = np.linalg.norm(residual, axis=1)
            rows.append({'amplitude': amplitude, 'field': name,
                         'mean_cylindrical': residual.mean(axis=0).tolist(),
                         'max_momentum': float(norms.max()),
                         'rms_momentum': float(np.sqrt(np.mean(norms**2)))})
    baseline = cylindrical_residual(full_operator(base, points, tau), points)
    negative_pressure = HarmonicPressureField(base, wave, -1., mean_coeff,
                                              pressure_coeffs)
    negative = TaylorAmplitudeField(negative_pressure, correction_coeffs)
    negative_residual = cylindrical_residual(full_operator(negative, points, tau), points)
    positive_residual = direct[(1., 'taylor_amplitude')]
    linear_coefficient = (positive_residual-negative_residual)/2
    quadratic_coefficient = (positive_residual+negative_residual)/2-baseline
    def polynomial_metrics(amplitude):
        result = baseline+amplitude*linear_coefficient+amplitude**2*quadratic_coefficient
        norms = np.linalg.norm(result, axis=1)
        return {'amplitude': float(amplitude),
                'covariance_fraction': float(amplitude**2),
                'max_momentum': float(norms.max()),
                'rms_momentum': float(np.sqrt(np.mean(norms**2)))}
    grid = np.linspace(-.3, .3, 1201)
    best = min((polynomial_metrics(x) for x in grid),
               key=lambda item: item['rms_momentum'])
    best_max = min((polynomial_metrics(x) for x in grid),
                   key=lambda item: item['max_momentum'])
    best_pressure = HarmonicPressureField(base, wave, best['amplitude'],
                                          mean_coeff, pressure_coeffs)
    best_field = TaylorAmplitudeField(best_pressure, correction_coeffs)
    best_direct = cylindrical_residual(full_operator(best_field, points, tau), points)
    best_reconstruction = (baseline+best['amplitude']*linear_coefficient
                           +best['amplitude']**2*quadratic_coefficient)
    report = {'tau': tau, 'point': source['point'],
              'mode_rows': mode_rows, 'rows': rows,
              'baseline': polynomial_metrics(0.),
              'best_local_rms': best,
              'best_local_max': best_max,
              'best_polynomial_direct_max_difference': float(np.max(np.abs(best_direct-best_reconstruction))),
              'scope': 'Exact-curl time-slope vector potential cancels the selected-point transverse linear Fourier residual at the registered endpoint; compact taper limits it to a 1e-4 tau interval. Pressure harmonics and local mean pressure proxy are retained. This is a local Taylor step, not a solved supported pulse amplitude PDE; full momentum sampled at 16 angles.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_taylor_amplitude.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
