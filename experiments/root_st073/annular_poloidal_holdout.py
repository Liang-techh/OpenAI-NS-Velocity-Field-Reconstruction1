"""Independent spatial/time quadrature checks for the annular correction."""
import json

import numpy as np

from annular_poloidal_volume_fit import load_candidate
from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from radial_peak_cone import current_field


def metrics(field, tau, order):
    points, weights = nodes(field, tau, order)
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, J, part = kinematics(field, points, tau, hs, ht)
    R = part+np.einsum('nij,nj->ni', J, u)
    norms = np.linalg.norm(R, axis=1)
    return {'max': float(np.max(norms)),
            'physical_volume_l2': float(np.sqrt(weights @ norms**2))}


def run():
    current = current_field()
    candidate = load_candidate()
    rows = []
    for k, order in ((6., 6), (6., 10), (5.5, 6)):
        tau = .5*2**(-k)
        row = {'k': k, 'tau': tau, 'order': order,
               'current': metrics(current, tau, order),
               'candidate': metrics(candidate, tau, order)}
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = {'rows': rows,
              'scope': 'Held-out Gauss6/Gauss10 spatial grids and earlier time for Gauss8-trained annular mode. No continuum or full-time certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_poloidal_holdout.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
