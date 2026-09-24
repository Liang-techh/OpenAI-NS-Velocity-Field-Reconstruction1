"""Neighboring-point full momentum screen of the local Taylor wave correction."""
import json

import numpy as np

from curl_wave_harmonic_pressure import HarmonicPressureField, full_operator
from curl_wave_prototype import LocalizedCurlWave
from curl_wave_taylor_amplitude import TaylorAmplitudeField
from joined_field import ROOT
from radial_peak_cone import current_field


def run():
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    harmonic = json.loads((ROOT/'compact_potential'/'curl_wave_harmonic_pressure.json').read_text())
    taylor = json.loads((ROOT/'compact_potential'/'curl_wave_taylor_amplitude.json').read_text())
    amplitude_report = json.loads((ROOT/'compact_potential'/'curl_wave_amplitude.json').read_text())
    base = current_field()
    wave = LocalizedCurlWave(source)
    pressure_coeffs = [complex(*row['pressure_coefficient'])
                       for row in harmonic['harmonic_projections']]
    corrections = [np.array([complex(*pair) for pair in row['potential_time_slope']])
                   for row in taylor['mode_rows']]
    field = TaylorAmplitudeField(
        HarmonicPressureField(base, wave, .1,
                              amplitude_report['mean_quadratic_cylindrical'][0],
                              pressure_coeffs), corrections)
    radius, _, z = source['point']
    positions = [('center', radius, z), ('radial_in', radius-.0005, z),
                 ('radial_out', radius+.0005, z), ('axial_down', radius, z-.0002),
                 ('axial_up', radius, z+.0002)]
    angles = np.arange(8)*2*np.pi/8
    points = np.array([[r*np.cos(theta), r*np.sin(theta), zz]
                       for _, r, zz in positions for theta in angles])
    before = np.linalg.norm(full_operator(base, points, source['tau']), axis=1)
    after = np.linalg.norm(full_operator(field, points, source['tau']), axis=1)
    rows = []
    for i, (name, r, zz) in enumerate(positions):
        sl = slice(i*len(angles), (i+1)*len(angles))
        rows.append({'position': name, 'r': r, 'z': zz,
                     'baseline_max': float(before[sl].max()),
                     'candidate_max': float(after[sl].max()),
                     'baseline_rms': float(np.sqrt(np.mean(before[sl]**2))),
                     'candidate_rms': float(np.sqrt(np.mean(after[sl]**2)))})
    report = {'tau': source['tau'], 'amplitude': .1, 'rows': rows,
              'scope': 'Eight angular points at center and four nearby meridional positions, with complete physical momentum and endpoint one-sided time differences. Diagnostic only.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_taylor_neighbor.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
