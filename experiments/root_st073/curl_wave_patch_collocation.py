"""Patch-wide first Taylor step for the compact exact-curl wave.

Fits the time derivative of a compact vector potential and harmonic pressure
over an r-z patch. This is a local collocation diagnostic, not an evolved
amplitude equation or a Navier--Stokes solution.
"""
import json

import numpy as np

from curl_wave_prototype import LocalizedCurlWave, bump, cylindrical_residual
from compact_control import cutoff
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


def basis(wave, mode, r, z):
    """Complex Fourier residual columns for -curl(A_tau)+grad(p)."""
    dr, dz = wave.radial_halfwidth, wave.axial_halfwidth
    xi, eta = (r-wave.radius)/dr, (z-wave.zcenter)/dz
    br, brr = bump(r, wave.radius, dr)
    bz, bzz = bump(z, wave.zcenter, dz)
    envelope, er, ez = br*bz, brr*bz, br*bzz
    m = mode['m']
    kr, _, kz = mode['normal']
    phase = np.exp(1j*(kr*(r-wave.radius)+kz*(z-wave.zcenter)))
    matrix = np.zeros((3, 36), complex)
    for a in range(3):
        for b in range(3):
            index = 3*a+b
            polynomial = xi**a*eta**b
            q = envelope*polynomial
            qr = er*polynomial+(a/dr*envelope*xi**(a-1)*eta**b if a else 0.)
            qz = ez*polynomial+(b/dz*envelope*xi**a*eta**(b-1) if b else 0.)
            # Curl in cylindrical coordinates, including the connection term.
            matrix[1, index] = -(qz+1j*kz*q)
            matrix[2, index] = 1j*m*q/r
            matrix[0, 9+index] = qz+1j*kz*q
            matrix[2, 9+index] = -(qr+q/r+1j*kr*q)
            matrix[0, 18+index] = -1j*m*q/r
            matrix[1, 18+index] = qr+1j*kr*q
            matrix[:, 27+index] = [qr+1j*kr*q, 1j*m*q/r, qz+1j*kz*q]
    # The first 27 columns are minus curl, since u_tau enters as -u_tau.
    matrix[:, :27] *= phase
    matrix[:, 27:] *= phase
    return matrix


def sample(base, wave, grid, angles):
    tau = wave.tau0
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    rows = []
    for xi, eta in grid:
        r = wave.radius+wave.radial_halfwidth*xi
        z = wave.zcenter+wave.axial_halfwidth*eta
        central = np.array([[r, 0., z]])
        u0, j0, part0 = kinematics(base, central, tau, hs, ht)
        points = np.column_stack((r*np.cos(angles), r*np.sin(angles),
                                  np.full(len(angles), z)))
        w, jw, partw = kinematics(wave, points, tau, hs, ht)
        ca, sa = np.cos(angles), np.sin(angles)
        rotation = np.zeros((len(angles), 3, 3))
        rotation[:, 0, 0] = ca
        rotation[:, 0, 1] = -sa
        rotation[:, 1, 0] = sa
        rotation[:, 1, 1] = ca
        rotation[:, 2, 2] = 1.
        u = np.einsum('nij,j->ni', rotation, u0[0])
        j = np.einsum('nik,kl,njl->nij', rotation, j0[0], rotation)
        part = np.einsum('nij,j->ni', rotation, part0[0])
        baseline = part+np.einsum('nij,nj->ni', j, u)
        linear = (partw+np.einsum('nij,nj->ni', j, w)
                  +np.einsum('nij,nj->ni', jw, u))
        quadratic = np.einsum('nij,nj->ni', jw, w)
        linear_cyl = cylindrical_residual(linear, points)
        harmonics = [2*np.mean(linear_cyl*np.exp(-1j*mode['m']*angles)[:, None],
                               axis=0) for mode in wave.waves]
        rows.append({'xi': xi, 'eta': eta, 'r': r, 'z': z,
                     'baseline': cylindrical_residual(baseline, points),
                     'linear': linear_cyl,
                     'quadratic': cylindrical_residual(quadratic, points),
                     'harmonics': harmonics})
    return rows


def fit_mode(rows, wave, mode_index, regularization=1e-4):
    blocks = [basis(wave, wave.waves[mode_index], row['r'], row['z'])
              for row in rows]
    matrix = np.vstack(blocks)
    target = -np.concatenate([row['harmonics'][mode_index] for row in rows])
    scale = np.maximum(np.linalg.norm(matrix, axis=0), 1e-30)
    normalized = matrix/scale
    ridge = regularization*np.linalg.norm(normalized, ord=2)
    augmented = np.vstack((normalized, ridge*np.eye(matrix.shape[1])))
    rhs = np.r_[target, np.zeros(matrix.shape[1], complex)]
    coefficients = np.linalg.lstsq(augmented, rhs, rcond=None)[0]/scale
    return coefficients


class PatchTaylorField:
    """Real callable velocity/pressure realizing the collocated endpoint jet."""
    def __init__(self, base, wave, coefficients, amplitude):
        self.base, self.wave = base, wave
        self.coefficients = coefficients
        self.amplitude = amplitude
        self.nu = base.nu

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity, pressure = self.base.fields(pts, ts)
        w, _ = self.wave.fields(pts, ts)
        velocity += self.amplitude*w
        for i, ((x, y, z), t) in enumerate(zip(pts, ts)):
            r = np.hypot(x, y)
            if r == 0 or abs(r-self.wave.radius) >= self.wave.radial_halfwidth or abs(z-self.wave.zcenter) >= self.wave.axial_halfwidth:
                continue
            dt = t-self.wave.tau0
            if dt < 0 or dt >= .0001:
                continue
            taper = float(cutoff((dt-.00005)/.00005)[0])
            theta = np.arctan2(y, x)
            xi = (r-self.wave.radius)/self.wave.radial_halfwidth
            eta = (z-self.wave.zcenter)/self.wave.axial_halfwidth
            br, brr = bump(r, self.wave.radius, self.wave.radial_halfwidth)
            bz, bzz = bump(z, self.wave.zcenter, self.wave.axial_halfwidth)
            e, er, ez = br*bz, brr*bz, br*bzz
            for mode, coeff in zip(self.wave.waves, self.coefficients):
                m = mode['m']
                kr, _, kz = mode['normal']
                carrier = np.exp(1j*(m*theta+kr*(r-self.wave.radius)
                                     +kz*(z-self.wave.zcenter)
                                     +mode['omega']*dt))
                values = np.zeros(4, complex)
                radial = np.zeros(4, complex)
                axial = np.zeros(4, complex)
                for a in range(3):
                    for b in range(3):
                        index = 3*a+b
                        poly = xi**a*eta**b
                        q = e*poly
                        qr = er*poly+(a/self.wave.radial_halfwidth*e*xi**(a-1)*eta**b if a else 0.)
                        qz = ez*poly+(b/self.wave.axial_halfwidth*e*xi**a*eta**(b-1) if b else 0.)
                        for comp in range(4):
                            c = coeff[9*comp+index]
                            values[comp] += c*q
                            radial[comp] += c*qr
                            axial[comp] += c*qz
                ar, at, az, p = values
                curl = np.array([1j*m*az/r-axial[1]-1j*kz*at,
                                 axial[0]+1j*kz*ar-radial[2]-1j*kr*az,
                                 radial[1]+at/r+1j*kr*at-1j*m*ar/r])
                wc = (dt*taper*carrier*curl).real
                ca, sa = x/r, y/r
                velocity[i] += self.amplitude*np.array([ca*wc[0]-sa*wc[1],
                                                        sa*wc[0]+ca*wc[1], wc[2]])
                pressure[i] += self.amplitude*(taper*carrier*p).real
        return velocity, pressure


def metrics(rows, wave, coefficients, amplitude):
    baseline, before, after, harmonic_before, harmonic_after = [], [], [], [], []
    angles = np.arange(rows[0]['linear'].shape[0])*2*np.pi/rows[0]['linear'].shape[0]
    for row in rows:
        correction = np.zeros_like(row['linear'])
        for mode, coeff, harmonic in zip(wave.waves, coefficients, row['harmonics']):
            complex_correction = basis(wave, mode, row['r'], row['z'])@coeff
            correction += (complex_correction[None, :]
                           *np.exp(1j*mode['m']*angles)[:, None]).real
            harmonic_before.extend(np.abs(harmonic))
            harmonic_after.extend(np.abs(harmonic+complex_correction))
        before.extend(np.linalg.norm(row['baseline']+amplitude*row['linear']
                                     +amplitude**2*row['quadratic'], axis=1))
        baseline.extend(np.linalg.norm(row['baseline'], axis=1))
        after.extend(np.linalg.norm(row['baseline']+amplitude*(row['linear']+correction)
                                    +amplitude**2*row['quadratic'], axis=1))
    def stats(values):
        values = np.asarray(values)
        return {'max': float(values.max()), 'rms': float(np.sqrt(np.mean(values**2)))}
    return {'amplitude': amplitude, 'baseline': stats(baseline),
            'before': stats(before), 'after': stats(after),
            'harmonic_before': stats(harmonic_before),
            'harmonic_after': stats(harmonic_after)}


def run():
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    wave = LocalizedCurlWave(source)
    base = current_field()
    angles = np.arange(8)*2*np.pi/8
    train_axis = np.linspace(-.6, .6, 5)
    test_axis = np.array([-.45, -.15, .15, .45])
    train = sample(base, wave, [(x, y) for x in train_axis for y in train_axis], angles)
    test = sample(base, wave, [(x, y) for x in test_axis for y in test_axis], angles)
    coefficients = [fit_mode(train, wave, i) for i in range(len(wave.waves))]
    held = test[10]
    check_angle = angles[1]
    check_point = np.array([[held['r']*np.cos(check_angle),
                             held['r']*np.sin(check_angle), held['z']]])
    check_amplitude = .1
    field = PatchTaylorField(base, wave, coefficients, check_amplitude)
    hs = .0005*np.sqrt(base.nu*wave.tau0)
    ht = .0001*wave.tau0
    du, dj, dp = kinematics(field, check_point, wave.tau0, hs, ht)
    direct = cylindrical_residual(dp+np.einsum('nij,nj->ni', dj, du),
                                  check_point)[0]
    correction = np.zeros(3)
    for mode, coeff in zip(wave.waves, coefficients):
        correction += ((basis(wave, mode, held['r'], held['z'])@coeff)
                       *np.exp(1j*mode['m']*check_angle)).real
    prediction = (held['baseline'][1]
                  +check_amplitude*(held['linear'][1]+correction)
                  +check_amplitude**2*held['quadratic'][1])
    report = {'tau': wave.tau0, 'center': [wave.radius, wave.zcenter],
              'train_nodes': len(train), 'heldout_nodes': len(test),
              'coefficient_norms': [float(np.linalg.norm(x)) for x in coefficients],
              'train': [metrics(train, wave, coefficients, alpha) for alpha in (.1, 1.)],
              'heldout': [metrics(test, wave, coefficients, alpha) for alpha in (.1, 1.)],
              'direct_endpoint_replay': {'predicted': prediction.tolist(),
                                         'direct': direct.tolist(),
                                         'max_absolute_difference': float(np.max(np.abs(prediction-direct)))},
              'scope': 'Patch collocation of one endpoint Taylor slope and harmonic pressure. Exact-curl basis preserves divergence when integrated in time, but no amplitude PDE or mean pressure/flux correction is solved.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_patch_collocation.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
