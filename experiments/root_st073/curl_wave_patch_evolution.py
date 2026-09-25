"""Two-stage time-marching collocation for a compact curl-potential patch.

At each time, evaluate the full nonlinear NS residual of the frozen current
potential, solve for its time derivative and pressure, then advance the
potential coefficients. This is an exploratory semi-discrete ODE, not a
certified solution or a completed pulse with temporal endpoint matching.
"""
import argparse
import json

import numpy as np

from curl_wave_patch_collocation import basis
from curl_wave_patch_mean import mean_basis
from curl_wave_prototype import LocalizedCurlWave, bump, cylindrical_residual
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


def polynomial_data(wave, r, z, coefficient, components):
    dr, dz = wave.radial_halfwidth, wave.axial_halfwidth
    xi, eta = (r-wave.radius)/dr, (z-wave.zcenter)/dz
    br, brr = bump(r, wave.radius, dr)
    bz, bzz = bump(z, wave.zcenter, dz)
    envelope, er, ez = br*bz, brr*bz, br*bzz
    values = np.zeros(components, dtype=coefficient.dtype)
    radial = np.zeros_like(values)
    axial = np.zeros_like(values)
    for a in range(3):
        for b in range(3):
            index = 3*a+b
            poly = xi**a*eta**b
            q = envelope*poly
            qr = er*poly+(a/dr*envelope*xi**(a-1)*eta**b if a else 0.)
            qz = ez*poly+(b/dz*envelope*xi**a*eta**(b-1) if b else 0.)
            for comp in range(components):
                c = coefficient[9*comp+index]
                values[comp] += c*q
                radial[comp] += c*qr
                axial[comp] += c*qz
    return values, radial, axial


class FrozenPotentialField:
    def __init__(self, base, wave, amplitude, harmonic_potentials,
                 mean_potential):
        self.base = base
        self.wave = wave
        self.amplitude = amplitude
        self.harmonic_potentials = harmonic_potentials
        self.mean_potential = mean_potential
        self.nu = base.nu

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity, pressure = self.base.fields(pts, ts)
        wave_velocity, _ = self.wave.fields(pts, ts)
        velocity += self.amplitude*wave_velocity
        for i, ((x, y, z), t) in enumerate(zip(pts, ts)):
            r = np.hypot(x, y)
            if r == 0 or abs(r-self.wave.radius) >= self.wave.radial_halfwidth or abs(z-self.wave.zcenter) >= self.wave.axial_halfwidth:
                continue
            theta = np.arctan2(y, x)
            cyl = np.zeros(3)
            for mode, coefficient in zip(self.wave.waves,
                                         self.harmonic_potentials):
                m = mode['m']
                kr, _, kz = mode['normal']
                values, radial, axial = polynomial_data(self.wave, r, z,
                                                          coefficient, 3)
                ar, at, az = values
                carrier = np.exp(1j*(m*theta+kr*(r-self.wave.radius)
                                     +kz*(z-self.wave.zcenter)
                                     +mode['omega']*(t-self.wave.tau0)))
                curl = np.array([1j*m*az/r-axial[1]-1j*kz*at,
                                 axial[0]+1j*kz*ar-radial[2]-1j*kr*az,
                                 radial[1]+at/r+1j*kr*at-1j*m*ar/r])
                cyl += (carrier*curl).real
            values, radial, axial = polynomial_data(self.wave, r, z,
                                                      self.mean_potential, 2)
            atheta = values[0]
            cyl += [-axial[0], -radial[1], radial[0]+atheta/r]
            ca, sa = x/r, y/r
            velocity[i] += [ca*cyl[0]-sa*cyl[1], sa*cyl[0]+ca*cyl[1], cyl[2]]
        return velocity, pressure


def sample_residual(field, wave, tau, grid, angles, time_step=None,
                    time_min=.5/64):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = min(.0001*tau, time_step/8) if time_step is not None else .0001*tau
    nodes = []
    points_blocks = []
    for xi, eta in grid:
        r = wave.radius+wave.radial_halfwidth*xi
        z = wave.zcenter+wave.axial_halfwidth*eta
        nodes.append((r, z))
        points_blocks.append(np.column_stack((r*np.cos(angles),
                                              r*np.sin(angles),
                                              np.full(len(angles), z))))
    points = np.vstack(points_blocks)
    u, grad, part = kinematics(field, points, tau, hs, ht,
                               time_min=time_min)
    residual = cylindrical_residual(part+np.einsum('nij,nj->ni', grad, u),
                                    points).reshape(len(nodes), len(angles), 3)
    return [{'r': r, 'z': z, 'residual': residual[i]}
            for i, (r, z) in enumerate(nodes)]


def ridge_solve(matrix, target, regularization=1e-4):
    scale = np.maximum(np.linalg.norm(matrix, axis=0), 1e-30)
    normalized = matrix/scale
    ridge = regularization*np.linalg.norm(normalized, ord=2)
    augmented = np.vstack((normalized, ridge*np.eye(matrix.shape[1])))
    rhs = np.r_[target, np.zeros(matrix.shape[1], dtype=target.dtype)]
    return np.linalg.lstsq(augmented, rhs, rcond=None)[0]/scale


def fit_slope(rows, wave, tau, angles, regularization=1e-4):
    harmonic_fits = []
    for mode in wave.waves:
        m = mode['m']
        phase = np.exp(1j*mode['omega']*(tau-wave.tau0))
        matrix = np.vstack([basis(wave, mode, row['r'], row['z'])*phase
                            for row in rows])
        target = -np.concatenate([2*np.mean(row['residual']
                                            *np.exp(-1j*m*angles)[:, None],
                                            axis=0) for row in rows])
        harmonic_fits.append(ridge_solve(matrix, target, regularization))
    mean_matrix = np.vstack([mean_basis(wave, row['r'], row['z'])
                             for row in rows])
    mean_target = -np.concatenate([row['residual'].mean(axis=0)
                                   for row in rows])
    mean_fit = ridge_solve(mean_matrix, mean_target, regularization)
    return harmonic_fits, mean_fit


def metrics(rows, wave, tau, angles, harmonic_fits, mean_fit):
    before, after = [], []
    for row in rows:
        delta = np.broadcast_to(mean_basis(wave, row['r'], row['z'])@mean_fit,
                                row['residual'].shape).copy()
        for mode, coefficients in zip(wave.waves, harmonic_fits):
            phase = np.exp(1j*mode['omega']*(tau-wave.tau0))
            harmonic = basis(wave, mode, row['r'], row['z'])@coefficients*phase
            delta += (harmonic[None, :]
                      *np.exp(1j*mode['m']*angles)[:, None]).real
        before.extend(np.linalg.norm(row['residual'], axis=1))
        after.extend(np.linalg.norm(row['residual']+delta, axis=1))
    def stats(values):
        values = np.asarray(values)
        return {'max': float(values.max()),
                'rms': float(np.sqrt(np.mean(values**2)))}
    return {'frozen': stats(before), 'projected_derivative_and_pressure': stats(after)}


def run(dt=1e-5, output_name='curl_wave_patch_evolution.json',
        stage_count=2):
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    base = current_field()
    wave = LocalizedCurlWave(source)
    amplitude = .005
    angles = np.arange(8)*2*np.pi/8
    train_axis = np.linspace(-.6, .6, 5)
    test_axis = np.array([-.45, -.15, .15, .45])
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    test_grid = [(x, y) for x in test_axis for y in test_axis]
    harmonic_potentials = [np.zeros(27, complex) for _ in wave.waves]
    mean_potential = np.zeros(18)
    stages = []
    for stage in range(stage_count):
        tau = wave.tau0+stage*dt
        frozen = FrozenPotentialField(base, wave, amplitude,
                                      harmonic_potentials, mean_potential)
        train = sample_residual(frozen, wave, tau, train_grid, angles)
        heldout = sample_residual(frozen, wave, tau, test_grid, angles)
        harmonic_fit, mean_fit = fit_slope(train, wave, tau, angles)
        stages.append({'stage': stage, 'tau': tau,
                       'train': metrics(train, wave, tau, angles,
                                        harmonic_fit, mean_fit),
                       'heldout': metrics(heldout, wave, tau, angles,
                                          harmonic_fit, mean_fit),
                       'harmonic_slope_norms': [float(np.linalg.norm(x[:27]))
                                                for x in harmonic_fit],
                       'mean_slope_norm': float(np.linalg.norm(mean_fit[:18])),
                       'pressure_norms': [float(np.linalg.norm(x[27:]))
                                          for x in harmonic_fit]
                                          +[float(np.linalg.norm(mean_fit[18:]))],
                       'harmonic_slope_coefficients': [[[float(v.real),
                                                          float(v.imag)] for v in fit]
                                                        for fit in harmonic_fit],
                       'mean_slope_coefficients': mean_fit.tolist()})
        harmonic_potentials = [current+dt*fit[:27]
                               for current, fit in zip(harmonic_potentials,
                                                       harmonic_fit)]
        mean_potential = mean_potential+dt*mean_fit[:18]
    report = {'amplitude': amplitude, 'time_step': dt,
              'stage_count': stage_count,
              'train_nodes': len(train_grid), 'heldout_nodes': len(test_grid),
              'stages': stages,
              'scope': 'Explicit-Euler coefficient steps from the full frozen-field residual, with instantaneous projected pressure. The separate trajectory screen directly evaluates interval interiors; no global residual bound or pulse endpoint matching.',
              'accepted': False}
    out = ROOT/'compact_potential'/output_name
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--time-step', type=float, default=1e-5)
    parser.add_argument('--output-name', default='curl_wave_patch_evolution.json')
    parser.add_argument('--stages', type=int, default=2)
    args = parser.parse_args()
    if args.time_step <= 0:
        parser.error('--time-step must be positive')
    if args.stages < 2:
        parser.error('--stages must be at least 2')
    run(args.time_step, args.output_name, args.stages)
