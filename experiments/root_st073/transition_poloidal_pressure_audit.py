"""Independent axial midpoint and coarse full-support audit of joint candidate."""
import json

import numpy as np

from axial_pressure_cone_screen import make_blocks
from curl_wave_cone_parameter_screen import cone_row
from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from radial_peak_cone import current_field
from transition_poloidal_pressure_screen import load_candidate


def residual(field, points, tau):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, J, part = kinematics(field, points, tau, hs, ht)
    return u, J, part+np.einsum('nij,nj->ni', J, u)


def full_support_metrics(field, tau, order=6):
    points, weights = nodes(field.base, tau, order)
    _, _, R = residual(field, points, tau)
    norm = np.linalg.norm(R, axis=1)
    return {'max': float(np.max(norm)),
            'physical_volume_l2': float(np.sqrt(weights @ norm**2))}


def midpoint_rows(field, heights, tau, radius):
    blocks = make_blocks(field.base, tau, radius, heights)
    points = np.vstack([block[3] for block in blocks])
    u, J, R = residual(field, points, tau)
    rows, offset = [], 0
    for z, qr, qw, block_points in blocks:
        sl = slice(offset, offset+len(block_points))
        offset += len(block_points)
        rows.append(cone_row(radius, z, u[sl][-1], J[sl][-1],
                             R[sl], qr, qw))
    return rows


def run():
    candidate = load_candidate()
    current = current_field()
    report_path = ROOT/'compact_potential'/'transition_poloidal_pressure_screen.json'
    fit = json.loads(report_path.read_text())
    heights = fit['heights']
    midpoints = tuple((a+b)/2 for a, b in zip(heights[:-1], heights[1:]))
    rows = midpoint_rows(candidate, midpoints, fit['tau'], fit['radius'])
    metrics = {}
    for k in (6., 5.5):
        tau = .5*2**(-k)
        metrics[str(k)] = {'tau': tau,
                           'current': full_support_metrics(current, tau),
                           'candidate': full_support_metrics(candidate, tau)}
    report = {'midpoint_heights': midpoints,
              'midpoint_pass_count': sum(row['strict_local_pass'] for row in rows),
              'midpoint_max_cone_ratio': max(row['cone_ratio'] for row in rows),
              'midpoint_rows': rows,
              'full_support_gauss6': metrics,
              'scope': 'Twelve unfitted axial midpoints plus coarse Gauss6 physical-volume screen at two times. No global supremum or spacetime L2 certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'transition_poloidal_pressure_audit.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({key: report[key] for key in
                      ('midpoint_pass_count', 'midpoint_max_cone_ratio',
                       'full_support_gauss6')}), flush=True)


if __name__ == '__main__':
    run()
