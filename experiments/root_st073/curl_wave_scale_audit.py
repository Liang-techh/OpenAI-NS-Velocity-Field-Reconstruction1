"""Physical scale audit for the current compact exact-curl wave.

This is a diagnostic of one prototype, not a verification of the normalized
small-parameter hypotheses in the OpenAI Navier--Stokes paper.
"""
import json
import argparse

import numpy as np

from curl_wave_prototype import LocalizedCurlWave, bump
from joined_field import ROOT


def pieces(wave, mode, r, z):
    br, brr = bump(r, wave.radius, wave.radial_halfwidth)
    bz, bzz = bump(z, wave.zcenter, wave.axial_halfwidth)
    envelope = br*bz
    er, ez = brr*bz, br*bzz
    cr, ct, cz = mode['potential']
    kr, _, kz = mode['normal']
    normal = np.array([kr, mode['m']/r, kz])
    leading = 1j*np.cross(normal, mode['potential'])*envelope
    remainder = np.array([-ct*ez, cr*ez-cz*er,
                          ct*(er+envelope/r)])
    return leading, remainder


def run(source_name='radial_peak_cone.json',
        output_name='curl_wave_scale_audit.json',
        radial_halfwidth=.0025, axial_halfwidth=.00075,
        dt=2.5e-6):
    source = json.loads((ROOT/'compact_potential'/source_name).read_text())
    wave = LocalizedCurlWave(source, radial_halfwidth, axial_halfwidth)
    axis = np.linspace(-.6, .6, 5)
    rows = []
    for mode in wave.waves:
        samples = []
        leading_sq = remainder_sq = 0.
        for xi in axis:
            for eta in axis:
                r = wave.radius+wave.radial_halfwidth*xi
                z = wave.zcenter+wave.axial_halfwidth*eta
                leading, remainder = pieces(wave, mode, r, z)
                lnorm = np.linalg.norm(leading)
                rnorm = np.linalg.norm(remainder)
                leading_sq += lnorm**2
                remainder_sq += rnorm**2
                samples.append({'xi': float(xi), 'eta': float(eta),
                                'remainder_over_carrier': float(rnorm/lnorm)})
        center_leading, center_remainder = pieces(wave, mode,
                                                  wave.radius, wave.zcenter)
        normal_norm = np.linalg.norm(mode['normal'])
        damping = wave.nu*normal_norm**2*dt
        ratio = max(item['remainder_over_carrier'] for item in samples)
        minimum_multiplier = int(np.ceil(ratio))
        maximum_multiplier = int(np.floor(1/np.sqrt(damping)))
        rows.append({'m': mode['m'], 'normal_norm': float(normal_norm),
                     'radial_phase_radians_per_halfwidth': float(abs(mode['normal'][0])*wave.radial_halfwidth),
                     'axial_phase_radians_per_halfwidth': float(abs(mode['normal'][2])*wave.axial_halfwidth),
                     'center_remainder_over_carrier': float(np.linalg.norm(center_remainder)/np.linalg.norm(center_leading)),
                     'grid_max_remainder_over_carrier': float(ratio),
                     'grid_rms_remainder_over_carrier': float(np.sqrt(remainder_sq/leading_sq)),
                     'viscous_damping_exponent_per_time_step': float(damping),
                     'integer_multiplier_for_remainder_at_most_carrier': minimum_multiplier,
                     'integer_multiplier_for_damping_at_most_one': maximum_multiplier,
                     'both_lenient_conditions_feasible': minimum_multiplier <= maximum_multiplier,
                     'largest_time_step_at_minimum_multiplier': float(1/(wave.nu*(minimum_multiplier*normal_norm)**2)),
                     'multiplier_screen': [{'multiplier': multiplier,
                                            'grid_max_remainder_over_carrier': float(ratio/multiplier),
                                            'damping_exponent_per_time_step': float(damping*multiplier**2)}
                                           for multiplier in (1, 2, 4, 8, 16, 32, 64)],
                     'samples': samples})
    report = {'tau': wave.tau0, 'nu': wave.nu, 'time_step': dt,
              'source_name': source_name,
              'field_id': source.get('field_id', 'current'),
              'selected_pulse_indices': wave.pulse_indices,
              'radius': wave.radius, 'zcenter': wave.zcenter,
              'radial_halfwidth': wave.radial_halfwidth,
              'axial_halfwidth': wave.axial_halfwidth,
              'grid_normalized_axis': axis.tolist(),
              'modes': rows,
              'interpretation': 'Scale all integer angular modes and radial/axial wavevectors by one positive integer while scaling vector potentials inversely. Then the curl envelope remainder falls as 1/L but the physical viscous damping exponent grows as L^2. Thresholds remainder/carrier <= 1 and damping exponent <= 1 are deliberately lenient diagnostics, not the paper normalized hypotheses.',
              'accepted': False}
    out = ROOT/'compact_potential'/output_name
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'source_name': source_name, 'time_step': dt,
                      'modes': [{key: row[key] for key in
                                 ('m', 'grid_max_remainder_over_carrier',
                                  'integer_multiplier_for_remainder_at_most_carrier',
                                  'integer_multiplier_for_damping_at_most_one',
                                  'both_lenient_conditions_feasible')}
                                for row in rows]}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-name', default='radial_peak_cone.json')
    parser.add_argument('--output-name', default='curl_wave_scale_audit.json')
    parser.add_argument('--radial-halfwidth', type=float, default=.0025)
    parser.add_argument('--axial-halfwidth', type=float, default=.00075)
    parser.add_argument('--time-step', type=float, default=2.5e-6)
    args = parser.parse_args()
    run(args.source_name, args.output_name, args.radial_halfwidth,
        args.axial_halfwidth, args.time_step)
