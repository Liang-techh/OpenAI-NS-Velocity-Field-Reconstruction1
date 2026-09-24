"""Check wave-cone compatibility of the multigrid annular correction."""
import json

import numpy as np

from annular_robust_fit import load_candidate
from axial_pressure_cone_screen import make_blocks
from curl_wave_cone_parameter_screen import cone_row
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


def rows(field, blocks, tau, radius):
    points = np.vstack([block[3] for block in blocks])
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, J, part = kinematics(field, points, tau, hs, ht)
    residual = part+np.einsum('nij,nj->ni', J, u)
    answer, offset = [], 0
    for z, qr, qw, block_points in blocks:
        sl = slice(offset, offset+len(block_points))
        offset += len(block_points)
        answer.append(cone_row(radius, z, u[sl][-1], J[sl][-1],
                               residual[sl], qr, qw))
    return answer


def run():
    current = current_field()
    candidate = load_candidate()
    tau = .5/64
    radius = .0056890761915166545
    heights = tuple(np.linspace(.0025, .003432627453438968, 13))
    blocks = make_blocks(current, tau, radius, heights)
    baseline = rows(current, blocks, tau, radius)
    corrected = rows(candidate, blocks, tau, radius)
    report = {'tau': tau, 'radius': radius, 'heights': heights,
              'baseline_pass_count': sum(x['strict_local_pass'] for x in baseline),
              'candidate_pass_count': sum(x['strict_local_pass'] for x in corrected),
              'baseline_rows': baseline, 'candidate_rows': corrected,
              'scope': 'Thirteen axial points at one radius/time. Local physical-field cone analogue only; no uniform wave or paper normalized-cone proof.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_robust_cone_audit.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'baseline_pass_count': report['baseline_pass_count'],
                      'candidate_pass_count': report['candidate_pass_count'],
                      'candidate_rows': corrected}), flush=True)


if __name__ == '__main__':
    run()
