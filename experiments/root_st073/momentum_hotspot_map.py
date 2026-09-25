"""Locate dominant full-momentum residual components of the current candidate."""
import json

import numpy as np

from joined_field import ROOT
from joint_collar_fit import nodes
from radial_peak_cone import operator
from radial_pressure_volume_constrained import load_candidate


def run():
    field = load_candidate()
    tau = .0084
    rows = []
    for order in (8, 10):
        points, weights = nodes(field, tau, order)
        velocity, gradient, residual = operator(field, points, tau)
        norm = np.linalg.norm(residual, axis=1)
        ranked = np.argsort(norm)[::-1]
        component_l2 = np.sqrt(np.sum(weights[:, None]*residual**2, axis=0))
        row = {'order': order, 'tau': tau,
               'max': float(np.max(norm)),
               'physical_volume_l2': float(np.sqrt(weights@norm**2)),
               'component_l2': component_l2.tolist(),
               'component_l2_fraction_squared': (component_l2**2/sum(component_l2**2)).tolist(),
               'worst_nodes': [
                   {'point': points[i].tolist(),
                    'residual': residual[i].tolist(),
                    'norm': float(norm[i]),
                    'velocity': velocity[i].tolist(),
                    'weight': float(weights[i])}
                   for i in ranked[:10]]}
        rows.append(row)
        print(json.dumps({k: row[k] for k in ('order', 'max',
                                             'physical_volume_l2',
                                             'component_l2_fraction_squared',
                                             'worst_nodes')}), flush=True)
    report = {'rows': rows,
              'scope': 'Gauss8/10 Cartesian physical momentum at one time; diagnostic of hotspots, not converged extrema or domain certificate.',
              'accepted': False}
    path = ROOT/'compact_potential'/'momentum_hotspot_map.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
