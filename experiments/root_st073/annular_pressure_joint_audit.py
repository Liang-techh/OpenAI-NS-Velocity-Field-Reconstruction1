"""Midpoint cone and independent-time volume checks for joint candidate."""
import json

import numpy as np

from annular_pressure_joint_screen import load_candidate
from annular_robust_cone_audit import rows as cone_rows
from axial_pressure_cone_screen import make_blocks
from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from radial_peak_cone import current_field


def metrics(field, grid_base, tau, order):
    points, weights = nodes(grid_base, tau, order)
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, J, part = kinematics(field, points, tau, hs, ht)
    residual = part+np.einsum('nij,nj->ni', J, u)
    norm = np.linalg.norm(residual, axis=1)
    return {'max': float(np.max(norm)),
            'physical_volume_l2': float(np.sqrt(weights @ norm**2))}


def run():
    base = current_field()
    candidate = load_candidate()
    fit = json.loads((ROOT/'compact_potential'/'annular_pressure_joint_screen.json').read_text())
    heights = fit['heights']
    midpoints = tuple((a+b)/2 for a, b in zip(heights[:-1], heights[1:]))
    blocks = make_blocks(base, fit['tau'], fit['radius'], midpoints)
    checked = cone_rows(candidate, blocks, fit['tau'], fit['radius'])
    volume = []
    for k, order in ((6., 7), (5.75, 7)):
        tau = .5*2**(-k)
        item = {'k': k, 'order': order,
                'baseline': metrics(base, base, tau, order),
                'candidate': metrics(candidate, base, tau, order)}
        volume.append(item)
        print(json.dumps(item), flush=True)
    report = {'midpoint_pass_count': sum(row['strict_local_pass'] for row in checked),
              'midpoint_max_cone_ratio': max(row['cone_ratio'] for row in checked),
              'midpoint_rows': checked, 'volume_holdouts': volume,
              'scope': 'Twelve unfitted axial midpoints and two independent Gauss7 spacetime samples. No radial cone patch, continuum, or exact-curl wave certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_pressure_joint_audit.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'midpoint_pass_count': report['midpoint_pass_count'],
                      'midpoint_max_cone_ratio': report['midpoint_max_cone_ratio']}),
          flush=True)


if __name__ == '__main__':
    run()
