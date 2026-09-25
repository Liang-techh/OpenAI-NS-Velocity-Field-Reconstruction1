"""Physical-coordinate Kelvin amplitude along the transported ST073 phase.

This projects the local linearized amplitude equation onto k-perpendicular
vectors, including shear and nu*|k|^2. It is a centerline analogue of the
paper's normalized pulse equation, not a supported amplitude/pressure solve.
"""
import json

import numpy as np
from scipy.integrate import solve_ivp

from curl_wave_prototype import LocalizedCurlWave
from joined_field import ROOT
from local_poloidal_basis_screen import load_robust_candidate
from transported_phase_screen import coefficients


def derivative(field, mode, tau, state):
    r, z, kr, kz, _, ar, at, az = state
    ur, uz, F, ur_r, ur_z, uz_r, uz_z, F_r, F_z = coefficients(
        field, r, z, tau)
    J = np.array([[ur_r, -F, ur_z],
                  [F+r*F_r, ur/r, r*F_z],
                  [uz_r, 0., uz_z]])
    amplitude = np.array([ar, at, az])
    normal = np.array([kr, mode/r, kz])
    Ja = J@amplitude
    rotation = F*np.array([-at, ar, 0.])
    pressure_projection = -2*normal*np.dot(normal, Ja)/np.dot(normal, normal)
    amp_tau = rotation+Ja+field.nu*np.dot(normal, normal)*amplitude+pressure_projection
    return np.r_[-ur, -uz,
                 ur_r*kr+uz_r*kz+mode*F_r,
                 ur_z*kr+uz_z*kz+mode*F_z,
                 mode*F, amp_tau]


def solve(field, mode, tau0, initial, times):
    options = {'rtol': 1e-9, 'atol': 1e-10, 'max_step': 1.5e-5,
               'dense_output': True}
    left = solve_ivp(lambda t, y: derivative(field, mode, t, y),
                     (tau0, min(times)), initial, **options)
    right = solve_ivp(lambda t, y: derivative(field, mode, t, y),
                      (tau0, max(times)), initial, **options)
    if not (left.success and right.success):
        raise RuntimeError(f'Amplitude integration failed: {left.message}; {right.message}')
    return lambda t: left.sol(t) if t < tau0 else right.sol(t)


def run():
    field = load_robust_candidate(.1)
    source = json.loads((ROOT/'compact_potential'/'local_poloidal_10pct_source.json').read_text())
    wave = LocalizedCurlWave(source, .00275, .000075,
                             time_halfwidth=.00015)
    tau0 = source['tau']
    times = [tau0-.00012, tau0-.00006, tau0,
             tau0+.00006, tau0+.00012]
    modes = []
    for item in wave.waves:
        kr, _, kz = item['normal']
        initial = np.r_[wave.radius, wave.zcenter, kr, kz, 0.,
                        item['amplitude']]
        trajectory = solve(field, item['m'], tau0, initial, times)
        initial_norm = np.linalg.norm(item['amplitude'])
        rows = []
        for tau in times:
            state = trajectory(tau)
            r, z, kr, kz, _ = state[:5]
            amplitude = state[5:]
            normal = np.array([kr, item['m']/r, kz])
            ur, uz, F, ur_r, ur_z, uz_r, uz_z, F_r, F_z = coefficients(
                field, r, z, tau)
            J = np.array([[ur_r, -F, ur_z],
                          [F+r*F_r, ur/r, r*F_z],
                          [uz_r, 0., uz_z]])
            shear_rate = float(np.dot(amplitude, J@amplitude)
                               /np.dot(amplitude, amplitude))
            backward_viscous_rate = float(field.nu*np.dot(normal, normal))
            rows.append({'tau': tau,
                         'center': [float(r), float(z)],
                         'normal': normal.tolist(),
                         'amplitude': amplitude.tolist(),
                         'amplitude_norm_ratio': float(np.linalg.norm(amplitude)/initial_norm),
                         'relative_transversality_error': float(abs(np.dot(normal, amplitude))
                                                              /(np.linalg.norm(normal)*np.linalg.norm(amplitude))),
                         'backward_shear_rate': shear_rate,
                         'backward_viscous_rate': backward_viscous_rate})
        modes.append({'m': item['m'], 'rows': rows})
        print(json.dumps({'m': item['m'],
                          'amplitude_norm_ratios': [x['amplitude_norm_ratio'] for x in rows],
                          'max_transversality_error': max(x['relative_transversality_error'] for x in rows),
                          'viscous_rates': [x['backward_viscous_rate'] for x in rows]}),
              flush=True)
    report = {'source': 'local_poloidal_10pct_source.json',
              'times': times, 'modes': modes,
              'scope': 'Unforced physical-coordinate projected Kelvin amplitude along a single moving center per harmonic, with viscosity and local velocity gradient. No spatially supported amplitude, pulse cutoff, pressure field, exact-curl momentum acceptance, or normalized-chart theorem.',
              'accepted': False}
    path = ROOT/'compact_potential'/'kelvin_amplitude_transport.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
