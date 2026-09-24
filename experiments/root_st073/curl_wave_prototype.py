"""Exact-curl two-harmonic compact wave prototype at the current peak.

Integer angular modes make the field periodic. This freezes local Kelvin
polarizations and has no amplitude/pressure PDE solve, so full momentum is the
decisive rejection or continuation test.
"""
import json

import numpy as np
from scipy.optimize import nnls

from compact_control import cutoff
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


def bump(value, center, halfwidth):
    s = (value-center+halfwidth)/(2*halfwidth)
    if s <= 0 or s >= 1:
        return 0., 0.
    return (1024*s**5*(1-s)**5,
            5120*s**4*(1-s)**4*(1-2*s)/(2*halfwidth))


class LocalizedCurlWave:
    def __init__(self, source):
        self.radius, _, self.zcenter = source['point']
        self.tau0 = source['tau']
        self.nu = .01
        self.radial_halfwidth = .0025
        self.axial_halfwidth = .00075
        self.target = np.array(source['local_tangential_stress_primitive'])
        self.waves = []
        columns = []
        for angular_mode, item in zip((1, 2), source['covariance_nnls']['selected']):
            pulse = source['kelvin_pulses'][item['index']]
            nold = np.array(pulse['peak_wavevector'])
            normal = np.array([nold[0], angular_mode/self.radius, nold[2]])
            amplitude = np.array(pulse['peak_amplitude'])
            amplitude -= normal*np.dot(normal, amplitude)/np.dot(normal, normal)
            columns.append(.5*amplitude[0]*amplitude[1:])
            potential = np.cross(normal, amplitude)/np.dot(normal, normal)
            omega = float(np.dot(source['velocity'], normal))
            self.waves.append({'m': angular_mode, 'normal': normal,
                               'amplitude': amplitude, 'potential': potential,
                               'omega': omega})
        matrix = np.array(columns).T
        weights, error = nnls(matrix, self.target)
        self.weights = weights
        self.covariance_error = float(error/np.linalg.norm(self.target))

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity = np.zeros_like(pts)
        pressure = np.zeros(len(pts))
        for i, (point, t) in enumerate(zip(pts, ts)):
            x, ycart, z = point
            r = np.hypot(x, ycart)
            if r == 0 or t < .5/64 or t > .5:
                continue
            br, br_r = bump(r, self.radius, self.radial_halfwidth)
            bz, bz_z = bump(z, self.zcenter, self.axial_halfwidth)
            if br == 0 or bz == 0:
                continue
            time_cut = float(cutoff((t-.012)/.008)[0])
            envelope = br*bz*time_cut
            er = br_r*bz*time_cut
            ez = br*bz_z*time_cut
            theta = np.arctan2(ycart, x)
            wr = wt = wz = 0.
            for weight, wave in zip(self.weights, self.waves):
                m = wave['m']
                kr, _, kz = wave['normal']
                cr, ct, cz = wave['potential']
                phase = (m*theta + kr*(r-self.radius) + kz*(z-self.zcenter)
                         + wave['omega']*(t-self.tau0))
                sn, cs = np.sin(phase), np.cos(phase)
                factor = np.sqrt(weight)
                wr += factor*((-m*cz/r+kz*ct)*envelope*sn -ct*ez*cs)
                wt += factor*((-kz*cr+kr*cz)*envelope*sn
                              +(cr*ez-cz*er)*cs)
                wz += factor*((-kr*ct+m*cr/r)*envelope*sn
                              +ct*(er+envelope/r)*cs)
            ca, sa = x/r, ycart/r
            velocity[i] = [wr*ca-wt*sa, wr*sa+wt*ca, wz]
        return velocity, pressure


class WavePerturbedField:
    def __init__(self, base, wave):
        self.base = base
        self.wave = wave
        self.nu = base.nu

    def fields(self, points, tau):
        u, p = self.base.fields(points, tau)
        w, _ = self.wave.fields(points, tau)
        return u+w, p


def cylindrical_residual(residual, points):
    theta = np.arctan2(points[:, 1], points[:, 0])
    ca, sa = np.cos(theta), np.sin(theta)
    return np.column_stack((ca*residual[:, 0]+sa*residual[:, 1],
                            -sa*residual[:, 0]+ca*residual[:, 1],
                            residual[:, 2]))


def run():
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    base = current_field()
    wave = LocalizedCurlWave(source)
    perturbed = WavePerturbedField(base, wave)
    angles = np.arange(16)*2*np.pi/16
    r, _, z = source['point']
    points = np.column_stack((r*np.cos(angles), r*np.sin(angles),
                              np.full(len(angles), z)))
    w = wave.fields(points, source['tau'])[0]
    wc = np.column_stack((np.cos(angles)*w[:, 0]+np.sin(angles)*w[:, 1],
                          -np.sin(angles)*w[:, 0]+np.cos(angles)*w[:, 1],
                          w[:, 2]))
    covariance = np.mean(wc[:, 0, None]*wc[:, 1:], axis=0)
    hs = .0005*np.sqrt(base.nu*source['tau'])
    ht = .0001*source['tau']
    def operator(field):
        u, grad, part = kinematics(field, points, source['tau'], hs, ht)
        residual = part + np.einsum('nij,nj->ni', grad, u)
        return cylindrical_residual(residual, points), np.trace(grad, axis1=1, axis2=2)
    old, old_div = operator(base)
    new, new_div = operator(perturbed)
    report = {
        'point': source['point'], 'tau': source['tau'],
        'target_stress': wave.target.tolist(),
        'periodicized_wavevectors': [x['normal'].tolist() for x in wave.waves],
        'transverse_amplitudes': [x['amplitude'].tolist() for x in wave.waves],
        'positive_weights': wave.weights.tolist(),
        'ideal_covariance_relative_error': wave.covariance_error,
        'sampled_exact_curl_covariance': covariance.tolist(),
        'sampled_covariance_relative_error': float(np.linalg.norm(covariance-wave.target)/np.linalg.norm(wave.target)),
        'baseline_mean_residual_cylindrical': old.mean(axis=0).tolist(),
        'wave_mean_residual_cylindrical': new.mean(axis=0).tolist(),
        'baseline_max_residual': float(np.max(np.linalg.norm(old, axis=1))),
        'wave_max_residual': float(np.max(np.linalg.norm(new, axis=1))),
        'wave_max_velocity': float(np.max(np.linalg.norm(w, axis=1))),
        'baseline_max_fd_divergence': float(np.max(np.abs(old_div))),
        'wave_max_fd_divergence': float(np.max(np.abs(new_div))),
        'support': {'radial_halfwidth': wave.radial_halfwidth,
                    'axial_halfwidth': wave.axial_halfwidth},
        'localization': {
            'carrier_magnitudes': [float(np.linalg.norm(x['normal'])) for x in wave.waves],
            'inverse_axial_halfwidth': 1/wave.axial_halfwidth,
            'inverse_radial_halfwidth': 1/wave.radial_halfwidth,
            'axial_cutoff_diffusion_scale': base.nu/wave.axial_halfwidth**2,
            'physical_analogue_growth_rate': float(np.sqrt(source['cone']['lambda_squared'])),
            'note': 'Dimensional scale comparison only; the paper uses normalized chart estimates and an evolving amplitude equation.'},
        'scope': 'Two integer angular harmonics from projected frozen Kelvin peak vectors. Exact analytic curl of a C4 compact vector potential; 16-angle local sample at one late point. No amplitude/pressure equation or stress moment repair. Full physical momentum including nonlinear wave terms is evaluated.',
        'accepted': False,
    }
    out = ROOT/'compact_potential'/'curl_wave_prototype.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
