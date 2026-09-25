"""Compare frozen and locally transported periodic phases for ST073 waves.

For backward time tau, D=-d_tau+u.grad. At an axisymmetric moving center,
the phase m*theta+kr(tau)*(r-rc(tau))+kz(tau)*(z-zc(tau))+phi0(tau)
obeys the eikonal equation through first spatial order when its center and
covector satisfy the ODEs below. This is a local physical-coordinate analogue
of the paper's phase transport, not the paper's normalized chart phase.
"""
import json

import numpy as np
from scipy.integrate import solve_ivp

from curl_wave_prototype import LocalizedCurlWave
from joined_field import ROOT
from local_poloidal_basis_screen import load_robust_candidate


def coefficients(field, r, z, tau):
    h = .0005*np.sqrt(field.nu*tau)
    points = np.array([[r, 0., z], [r+h, 0., z], [r-h, 0., z],
                       [r, 0., z+h], [r, 0., z-h]])
    velocity, _ = field.fields(points, tau)
    ur, ut, uz = velocity[0]
    ur_r = (velocity[1, 0]-velocity[2, 0])/(2*h)
    ur_z = (velocity[3, 0]-velocity[4, 0])/(2*h)
    uz_r = (velocity[1, 2]-velocity[2, 2])/(2*h)
    uz_z = (velocity[3, 2]-velocity[4, 2])/(2*h)
    F = ut/r
    F_r = (velocity[1, 1]/(r+h)-velocity[2, 1]/(r-h))/(2*h)
    F_z = (velocity[3, 1]-velocity[4, 1])/(2*h*r)
    return np.array([ur, uz, F, ur_r, ur_z, uz_r, uz_z, F_r, F_z])


def derivative(field, mode, tau, state):
    r, z, kr, kz, _ = state
    ur, uz, F, ur_r, ur_z, uz_r, uz_z, F_r, F_z = coefficients(
        field, r, z, tau)
    return np.array([-ur, -uz,
                     ur_r*kr+uz_r*kz+mode*F_r,
                     ur_z*kr+uz_z*kz+mode*F_z,
                     mode*F])


def flow(field, mode, tau0, state0, tau_low, tau_high):
    options = {'rtol': 1e-9, 'atol': 1e-10, 'max_step': 1.5e-5,
               'dense_output': True}
    left = solve_ivp(lambda tau, state: derivative(field, mode, tau, state),
                     (tau0, tau_low), state0, **options)
    right = solve_ivp(lambda tau, state: derivative(field, mode, tau, state),
                      (tau0, tau_high), state0, **options)
    if not (left.success and right.success):
        raise RuntimeError(f'Characteristic integration failed: {left.message}; {right.message}')
    return lambda tau: (left.sol(tau) if tau < tau0 else right.sol(tau))


def evaluate(field, mode, initial, trajectory, times, radial_offsets,
             axial_offsets):
    r0, z0, tau0, kr0, kz0, omega0 = initial
    rows = []
    for tau in times:
        state = trajectory(tau)
        rc, zc, kr, kz, _ = state
        rates = derivative(field, mode, tau, state)
        radii = rc+radial_offsets
        heights = zc+axial_offsets
        points = np.array([[r, 0., z] for r in radii for z in heights])
        u, _ = field.fields(points, tau)
        for point, velocity in zip(points, u):
            r, _, z = point
            ur, ut, uz = velocity
            F = ut/r
            frozen = -omega0+ur*kr0+F*mode+uz*kz0
            dphi = (rates[2]*(r-rc)+rates[3]*(z-zc)
                    -kr*rates[0]-kz*rates[1]+rates[4])
            transported = -dphi+ur*kr+F*mode+uz*kz
            rows.append({'tau': float(tau), 'r': float(r), 'z': float(z),
                         'frozen_defect': float(frozen),
                         'transported_defect': float(transported)})
    frozen = np.array([row['frozen_defect'] for row in rows])
    transported = np.array([row['transported_defect'] for row in rows])
    return {'rows': rows,
            'frozen_max_abs': float(np.max(np.abs(frozen))),
            'transported_max_abs': float(np.max(np.abs(transported))),
            'frozen_rms': float(np.sqrt(np.mean(frozen**2))),
            'transported_rms': float(np.sqrt(np.mean(transported**2))),
            'trajectory': [{'tau': float(tau),
                            'state': trajectory(tau).tolist()}
                           for tau in times]}


def run():
    field = load_robust_candidate(.1)
    source = json.loads((ROOT/'compact_potential'/'local_poloidal_10pct_source.json').read_text())
    wave = LocalizedCurlWave(source, radial_halfwidth=.00275,
                             axial_halfwidth=.000075,
                             time_halfwidth=.00015)
    tau0 = source['tau']
    times = [tau0-.00012, tau0-.00006, tau0,
             tau0+.00006, tau0+.00012]
    radial_offsets = np.linspace(-.0011, .0011, 5)
    axial_offsets = np.linspace(-.00003, .00003, 5)
    modes = []
    for item in wave.waves:
        m = item['m']
        kr0, _, kz0 = item['normal']
        state0 = np.array([wave.radius, wave.zcenter, kr0, kz0, 0.])
        trajectory = flow(field, m, tau0, state0,
                          min(times), max(times))
        metrics = evaluate(field, m,
                           (wave.radius, wave.zcenter, tau0,
                            kr0, kz0, item['omega']), trajectory,
                           times, radial_offsets, axial_offsets)
        modes.append({'m': m, 'initial_normal': item['normal'].tolist(),
                      'initial_omega': item['omega'], **metrics})
        print(json.dumps({'m': m,
                          'frozen_max_abs': metrics['frozen_max_abs'],
                          'transported_max_abs': metrics['transported_max_abs'],
                          'frozen_rms': metrics['frozen_rms'],
                          'transported_rms': metrics['transported_rms']}),
              flush=True)
    report = {'source': 'local_poloidal_10pct_source.json',
              'times': times,
              'radial_offsets_from_transport_center': radial_offsets.tolist(),
              'axial_offsets_from_transport_center': axial_offsets.tolist(),
              'modes': modes,
              'scope': 'Local cylindrical characteristic phase with integer angular mode, compared against frozen phase on a 5x5x5 moving sample. No amplitude transport, exact-curl field with this phase, full momentum test, or paper normalized estimate.',
              'accepted': False}
    path = ROOT/'compact_potential'/'transported_phase_screen.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
