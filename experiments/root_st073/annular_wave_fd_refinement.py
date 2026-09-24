"""Step-size check for complete momentum of the interior exact-curl pulse."""
import json

import numpy as np

from annular_pressure_scale_screen import load_candidate
from curl_wave_prototype import LocalizedCurlWave, WavePerturbedField
from joined_field import ROOT
from joint_collar_fit import kinematics


def metrics(field, points, tau, hs, ht):
    u, J, part = kinematics(field, points, tau, hs, ht)
    residual = part+np.einsum('nij,nj->ni', J, u)
    norms = np.linalg.norm(residual, axis=1)
    return {'max': float(np.max(norms)),
            'rms': float(np.sqrt(np.mean(norms**2))),
            'max_fd_divergence': float(np.max(np.abs(np.trace(J, axis1=1, axis2=2))))}


def run():
    source = json.loads((ROOT/'compact_potential'/'annular_pressure_scale_cone_t0084_z0034.json').read_text())
    base = load_candidate()
    wave = LocalizedCurlWave(source, .0007, .00005,
                             time_halfwidth=.0003)
    field = WavePerturbedField(base, wave)
    radius, _, z = source['point']
    theta = np.arange(16)*2*np.pi/16
    points = np.column_stack((radius*np.cos(theta), radius*np.sin(theta),
                              np.full(len(theta), z)))
    rows = []
    for fraction in (1., .5, .25):
        hs = fraction*.0005*np.sqrt(base.nu*source['tau'])
        ht = fraction*.0001*source['tau']
        row = {'step_fraction': fraction, 'hs': hs, 'ht': ht,
               'baseline': metrics(base, points, source['tau'], hs, ht),
               'wave': metrics(field, points, source['tau'], hs, ht)}
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = {'source_name': 'annular_pressure_scale_cone_t0084_z0034.json',
              'rows': rows,
              'scope': 'Fourth-order Cartesian/time finite-difference refinement at 16 angles on one ring; not a global error bound or analytic derivative certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_wave_fd_refinement.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
