"""Axisymmetric compact mean update for the patch-collocated curl wave.

This fits the mean endpoint time derivative using two vector-potential
components and pressure. It screens the complete nonlinear momentum residual
at held-out nodes, but does not solve the paper's time-dependent mean system.
"""
import json

import numpy as np

from compact_control import cutoff
from curl_wave_patch_collocation import (PatchTaylorField, basis, fit_mode,
                                         sample)
from curl_wave_prototype import LocalizedCurlWave, bump, cylindrical_residual
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


def mean_basis(wave, r, z):
    """Columns for -curl(A_theta e_theta+A_z e_z)+grad(p)."""
    dr, dz = wave.radial_halfwidth, wave.axial_halfwidth
    xi, eta = (r-wave.radius)/dr, (z-wave.zcenter)/dz
    br, brr = bump(r, wave.radius, dr)
    bz, bzz = bump(z, wave.zcenter, dz)
    envelope, er, ez = br*bz, brr*bz, br*bzz
    matrix = np.zeros((3, 27))
    for a in range(3):
        for b in range(3):
            index = 3*a+b
            poly = xi**a*eta**b
            q = envelope*poly
            qr = er*poly+(a/dr*envelope*xi**(a-1)*eta**b if a else 0.)
            qz = ez*poly+(b/dz*envelope*xi**a*eta**(b-1) if b else 0.)
            matrix[:, index] = [qz, 0., -qr-q/r]
            matrix[:, 9+index] = [0., qr, 0.]
            matrix[:, 18+index] = [qr, 0., qz]
    return matrix


def fit_mean(rows, wave, amplitude, regularization=1e-4):
    matrix = np.vstack([mean_basis(wave, row['r'], row['z']) for row in rows])
    target = -np.concatenate([np.mean(row['baseline']+amplitude*row['linear']
                                       +amplitude**2*row['quadratic'], axis=0)
                              for row in rows])
    scale = np.maximum(np.linalg.norm(matrix, axis=0), 1e-30)
    normalized = matrix/scale
    ridge = regularization*np.linalg.norm(normalized, ord=2)
    augmented = np.vstack((normalized, ridge*np.eye(matrix.shape[1])))
    coefficients = np.linalg.lstsq(augmented,
                                   np.r_[target, np.zeros(matrix.shape[1])],
                                   rcond=None)[0]/scale
    return coefficients


class MeanPatchField(PatchTaylorField):
    """Callable exact-curl mean update at the registered endpoint."""
    def __init__(self, base, wave, coefficients, amplitude, mean_coefficients,
                 temporal_cutoff=True):
        super().__init__(base, wave, coefficients, amplitude, temporal_cutoff)
        self.mean_coefficients = mean_coefficients

    def fields(self, points, tau):
        velocity, pressure = super().fields(points, tau)
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        for i, ((x, y, z), t) in enumerate(zip(pts, ts)):
            r = np.hypot(x, y)
            if r == 0 or abs(r-self.wave.radius) >= self.wave.radial_halfwidth or abs(z-self.wave.zcenter) >= self.wave.axial_halfwidth:
                continue
            dt = t-self.wave.tau0
            if dt < 0 or (self.temporal_cutoff and dt >= .0001):
                continue
            taper = (float(cutoff((dt-.00005)/.00005)[0])
                     if self.temporal_cutoff else 1.)
            xi = (r-self.wave.radius)/self.wave.radial_halfwidth
            eta = (z-self.wave.zcenter)/self.wave.axial_halfwidth
            br, brr = bump(r, self.wave.radius, self.wave.radial_halfwidth)
            bz, bzz = bump(z, self.wave.zcenter, self.wave.axial_halfwidth)
            envelope, er, ez = br*bz, brr*bz, br*bzz
            values = np.zeros(3)
            radial = np.zeros(3)
            axial = np.zeros(3)
            for a in range(3):
                for b in range(3):
                    index = 3*a+b
                    poly = xi**a*eta**b
                    q = envelope*poly
                    qr = er*poly+(a/self.wave.radial_halfwidth*envelope*xi**(a-1)*eta**b if a else 0.)
                    qz = ez*poly+(b/self.wave.axial_halfwidth*envelope*xi**a*eta**(b-1) if b else 0.)
                    for comp in range(3):
                        c = self.mean_coefficients[9*comp+index]
                        values[comp] += c*q
                        radial[comp] += c*qr
                        axial[comp] += c*qz
            atheta, az, p = values
            curl = np.array([-axial[0], -radial[1], radial[0]+atheta/r])
            ca, sa = x/r, y/r
            wc = dt*taper*curl
            velocity[i] += [ca*wc[0]-sa*wc[1], sa*wc[0]+ca*wc[1], wc[2]]
            pressure[i] += taper*p
        return velocity, pressure


def metrics(rows, wave, harmonics, mean_coefficients, amplitude):
    nangle = rows[0]['linear'].shape[0]
    angles = np.arange(nangle)*2*np.pi/nangle
    baseline, uncorrected, harmonic_only, corrected = [], [], [], []
    mean_before, mean_after = [], []
    for row in rows:
        harmonic_response = np.zeros_like(row['linear'])
        for mode, coeff in zip(wave.waves, harmonics):
            complex_response = basis(wave, mode, row['r'], row['z'])@coeff
            harmonic_response += (complex_response[None, :]
                                  *np.exp(1j*mode['m']*angles)[:, None]).real
        mean_response = mean_basis(wave, row['r'], row['z'])@mean_coefficients
        old = row['baseline']+amplitude*row['linear']+amplitude**2*row['quadratic']
        new = old+amplitude*harmonic_response+mean_response
        baseline.extend(np.linalg.norm(row['baseline'], axis=1))
        uncorrected.extend(np.linalg.norm(old, axis=1))
        harmonic_only.extend(np.linalg.norm(old+amplitude*harmonic_response, axis=1))
        corrected.extend(np.linalg.norm(new, axis=1))
        mean_before.append(np.linalg.norm(np.mean(old, axis=0)))
        mean_after.append(np.linalg.norm(np.mean(new, axis=0)))
    def stats(items):
        items = np.asarray(items)
        return {'max': float(items.max()), 'rms': float(np.sqrt(np.mean(items**2)))}
    return {'baseline': stats(baseline), 'wave_only': stats(uncorrected),
            'harmonics_only': stats(harmonic_only), 'with_mean': stats(corrected),
            'mean_before': stats(mean_before), 'mean_after': stats(mean_after)}


def run():
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    wave = LocalizedCurlWave(source)
    base = current_field()
    amplitude = .1
    angles = np.arange(8)*2*np.pi/8
    train_axis = np.linspace(-.6, .6, 5)
    test_axis = np.array([-.45, -.15, .15, .45])
    train = sample(base, wave, [(x, y) for x in train_axis for y in train_axis], angles)
    test = sample(base, wave, [(x, y) for x in test_axis for y in test_axis], angles)
    harmonics = [fit_mode(train, wave, i) for i in range(len(wave.waves))]
    amplitude_scan = []
    for trial_amplitude in (0., .005, .01, .02, .05, .1):
        trial_mean = fit_mean(train, wave, trial_amplitude)
        amplitude_scan.append({'amplitude': trial_amplitude,
                               'train': metrics(train, wave, harmonics, trial_mean,
                                                trial_amplitude),
                               'heldout': metrics(test, wave, harmonics, trial_mean,
                                                  trial_amplitude)})
    mean_coefficients = fit_mean(train, wave, amplitude)
    held = test[10]
    angle = angles[1]
    point = np.array([[held['r']*np.cos(angle), held['r']*np.sin(angle), held['z']]])
    hs = .0005*np.sqrt(base.nu*wave.tau0)
    ht = .0001*wave.tau0
    field = MeanPatchField(base, wave, harmonics, amplitude, mean_coefficients)
    u, j, part = kinematics(field, point, wave.tau0, hs, ht)
    direct = cylindrical_residual(part+np.einsum('nij,nj->ni', j, u), point)[0]
    correction = mean_basis(wave, held['r'], held['z'])@mean_coefficients
    for mode, coeff in zip(wave.waves, harmonics):
        correction += amplitude*((basis(wave, mode, held['r'], held['z'])@coeff)
                                 *np.exp(1j*mode['m']*angle)).real
    predicted = (held['baseline'][1]+amplitude*held['linear'][1]
                 +amplitude**2*held['quadratic'][1]+correction)
    probe_amplitude = .005
    probe_mean = fit_mean(train, wave, probe_amplitude)
    probe_field = MeanPatchField(base, wave, harmonics, probe_amplitude, probe_mean)
    continued_field = MeanPatchField(base, wave, harmonics, probe_amplitude,
                                     probe_mean, temporal_cutoff=False)
    mean_only_field = MeanPatchField(base, wave, harmonics, 0.,
                                     fit_mean(train, wave, 0.),
                                     temporal_cutoff=False)
    harmonic_only_field = PatchTaylorField(base, wave, harmonics,
                                            probe_amplitude,
                                            temporal_cutoff=False)
    time_probe = []
    for offset in (0., .00001, .00005, .00009, .00011):
        t = wave.tau0+offset
        probe_u, probe_j, probe_part = kinematics(probe_field, point, t, hs, ht)
        probe_residual = probe_part+np.einsum('nij,nj->ni', probe_j, probe_u)
        continued_u, continued_j, continued_part = kinematics(continued_field,
                                                              point, t, hs, ht)
        continued_residual = (continued_part
                              +np.einsum('nij,nj->ni', continued_j, continued_u))
        mean_u, mean_j, mean_part = kinematics(mean_only_field, point, t, hs, ht)
        mean_residual = mean_part+np.einsum('nij,nj->ni', mean_j, mean_u)
        harmonic_u, harmonic_j, harmonic_part = kinematics(harmonic_only_field,
                                                            point, t, hs, ht)
        harmonic_residual = (harmonic_part
                             +np.einsum('nij,nj->ni', harmonic_j, harmonic_u))
        base_u, base_j, base_part = kinematics(base, point, t, hs, ht)
        base_residual = base_part+np.einsum('nij,nj->ni', base_j, base_u)
        time_probe.append({'offset': offset,
                           'baseline_norm': float(np.linalg.norm(base_residual)),
                           'corrected_norm': float(np.linalg.norm(probe_residual)),
                           'continued_norm': float(np.linalg.norm(continued_residual)),
                           'mean_only_norm': float(np.linalg.norm(mean_residual)),
                           'harmonic_only_norm': float(np.linalg.norm(harmonic_residual)),
                           'continued_speed': float(np.linalg.norm(continued_u)),
                           'baseline_speed': float(np.linalg.norm(base_u))})
    report = {'tau': wave.tau0, 'amplitude': amplitude,
              'train_nodes': len(train), 'heldout_nodes': len(test),
              'mean_coefficient_norm': float(np.linalg.norm(mean_coefficients)),
              'train': metrics(train, wave, harmonics, mean_coefficients, amplitude),
              'heldout': metrics(test, wave, harmonics, mean_coefficients, amplitude),
              'amplitude_scan': amplitude_scan,
              'time_probe_amplitude': probe_amplitude,
              'time_probe': time_probe,
              'direct_endpoint_replay': {'predicted': predicted.tolist(),
                                         'direct': direct.tolist(),
                                         'max_absolute_difference': float(np.max(np.abs(predicted-direct)))},
              'scope': 'Endpoint patch collocation of exact-curl axisymmetric mean and two nonaxisymmetric harmonics. No time-integrated amplitude equation, pressure Poisson solve, or global momentum bound.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_patch_mean.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
