"""Amplitude and divergence diagnostics for the transported exact-curl wave."""
import json

import numpy as np

from curl_wave_prototype import LocalizedCurlWave
from joined_field import ROOT
from local_poloidal_basis_screen import load_robust_candidate
from transported_curl_wave import TransportedCurlWave


def divergence_fd(wave, points, tau, h):
    divergence = np.zeros(len(points))
    for axis in range(3):
        step = h*np.eye(3)[axis]
        up = wave.fields(points+step, tau)[0]
        um = wave.fields(points-step, tau)[0]
        divergence += (up[:, axis]-um[:, axis])/(2*h)
    return float(np.max(np.abs(divergence)))


def run():
    source = json.loads((ROOT/'compact_potential'/'local_poloidal_10pct_source.json').read_text())
    base = load_robust_candidate(.1)
    frozen = LocalizedCurlWave(source, .00275, .000075,
                               time_halfwidth=.00015)
    transported = TransportedCurlWave(source, base)
    radius, _, z = source['point']
    angles = np.arange(16)*2*np.pi/16
    points = np.column_stack((radius*np.cos(angles), radius*np.sin(angles),
                              np.full(len(angles), z)))
    rows = []
    for tau in (source['tau'], source['tau']+.00006):
        row = {'tau': tau}
        for name, wave in (('frozen', frozen), ('transported', transported)):
            velocity = wave.fields(points, tau)[0]
            row[name] = {'max_wave_speed': float(np.max(np.linalg.norm(velocity, axis=1))),
                         'max_fd_divergence': [divergence_fd(wave, points, tau, h)
                                               for h in (4e-6, 2e-6, 1e-6)]}
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = {'rows': rows,
              'scope': '16-angle ring speed and second-order finite-difference divergence at three spacings. Exact incompressibility follows from the spatial curl formula, not this sample.',
              'accepted': False}
    path = ROOT/'compact_potential'/'transported_wave_diagnostics.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
