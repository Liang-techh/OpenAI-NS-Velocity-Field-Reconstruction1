"""Screen existing compact correction groups for a wider local cone.

Group kinematics are precomputed once. The full residual for any three group
scales is then evaluated exactly as a quadratic polynomial in those scales.
This is a fixed-radius, finite-sample design screen, not a global NS fit.
"""
import itertools
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from joined_field import ROOT, coordinates
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


class ModeGroup:
    def __init__(self, modes, amplitudes, nu):
        self.modes = modes
        self.amplitudes = np.asarray(amplitudes)
        self.nu = nu

    def fields(self, points, tau):
        points = np.asarray(points)
        velocity = np.zeros_like(points)
        pressure = np.zeros(len(points))
        for amplitude, mode in zip(self.amplitudes, self.modes):
            if amplitude:
                du, dp = mode.fields(points, tau)
                velocity += amplitude*du
                pressure += amplitude*dp
        return velocity, pressure


def cone_row(radius, z, u, J, residual, quad_r, quad_w):
    target = np.array([-np.dot(quad_w*quad_r**2, residual[:-1, 1])/radius**2,
                       -np.dot(quad_w*quad_r, residual[:-1, 2])/radius])
    F = u[1]/radius
    shear = np.array([J[1, 0]-F, J[2, 0]])
    N = shear/np.linalg.norm(shear)
    K = np.array([-N[1], N[0]])
    lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
    row = {'z': z, 'target': target.tolist(),
           'target_dot_N': float(target@N),
           'lambda_squared': float(lam2),
           'center_residual_norm': float(np.linalg.norm(residual[-1]))}
    if lam2 > 0 and abs(2*F*N[0]) > 1e-14 and abs(target@N) > 1e-14:
        c = np.sqrt(lam2)/(2*F*N[0])
        row['cone_ratio'] = float(abs(c*(target@K)/(target@N)))
        row['strict_local_pass'] = bool(target@N < 0 and row['cone_ratio'] < 1)
    else:
        row['strict_local_pass'] = False
    return row


def run():
    current = current_field()
    swirl = current.current
    joint = swirl.joint
    base = joint.base
    groups = [ModeGroup(joint.modes, joint.amplitudes, current.nu),
              ModeGroup(swirl.modes, swirl.amplitudes, current.nu),
              ModeGroup(current.modes, current.amplitudes, current.nu)]
    tau = .5/64
    radius = .0056890761915166545
    heights = (.0025, .00275, .003, .00325,
               .003432627453438968)
    nodes, weights = leggauss(10)
    points_blocks, quad_data = [], []
    for z in heights:
        q = float(coordinates(0., z/np.sqrt(current.nu), tau,
                              base.base.inner.h)['q'])
        inner_radius = np.sqrt(current.nu)*np.sqrt(2*q*3/64)
        edges = sorted(set([0., inner_radius, radius]))
        quad_r, quad_w = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            quad_r.extend((lo+hi)/2+(hi-lo)/2*nodes)
            quad_w.extend((hi-lo)/2*weights)
        quad_r, quad_w = np.asarray(quad_r), np.asarray(quad_w)
        block_r = np.r_[quad_r, radius]
        points_blocks.append(np.column_stack((block_r, np.zeros(len(block_r)),
                                              np.full(len(block_r), z))))
        quad_data.append((quad_r, quad_w, len(block_r)))
    points = np.vstack(points_blocks)
    hs = .0005*np.sqrt(current.nu*tau)
    ht = .0001*tau
    fields = [base]+groups
    kins = [kinematics(field, points, tau, hs, ht) for field in fields]
    base_u, base_J, base_part = kins[0]
    group_u = np.stack([data[0] for data in kins[1:]], axis=-1)
    group_J = np.stack([data[1] for data in kins[1:]], axis=-1)
    group_part = np.stack([data[2] for data in kins[1:]], axis=-1)
    def evaluate(scales):
        scales = np.asarray(scales)
        u = base_u+np.einsum('nik,k->ni', group_u, scales)
        J = base_J+np.einsum('nijk,k->nij', group_J, scales)
        part = base_part+np.einsum('nik,k->ni', group_part, scales)
        residual = part+np.einsum('nij,nj->ni', J, u)
        rows, offset = [], 0
        for z, (quad_r, quad_w, block_length) in zip(heights, quad_data):
            sl = slice(offset, offset+block_length)
            rows.append(cone_row(radius, z, u[sl][-1], J[sl][-1],
                                 residual[sl], quad_r, quad_w))
            offset += block_length
        return {'scales': list(map(float, scales)),
                'pass_count': sum(row['strict_local_pass'] for row in rows),
                'max_center_residual': max(row['center_residual_norm'] for row in rows),
                'rows': rows}
    current_result = evaluate((1., 1., 1.))
    factors = {'scaled_joint_collar': (0., .5, .75, 1., 1.25, 1.5, 2.),
               'radial_swirl': (1.,),
               'incremental_collar': (0., .5, 1., 1.25, 1.5, 1.75,
                                      2., 2.5, 3.)}
    candidates = [evaluate(scales) for scales in itertools.product(
        factors['scaled_joint_collar'], factors['radial_swirl'],
        factors['incremental_collar'])]
    candidates.sort(key=lambda row: (-row['pass_count'],
                                     row['max_center_residual']))
    report = {'tau': tau, 'radius': radius, 'heights': heights,
              'group_names': ['scaled_joint_collar', 'radial_swirl',
                              'incremental_collar'],
              'factors': factors, 'current': current_result,
              'group_max_speed_on_sample': [float(np.max(np.linalg.norm(group_u[:, :, i], axis=1)))
                                             for i in range(3)],
              'best_by_pass_then_center_residual': candidates[:20],
              'all_summary': [{'scales': row['scales'],
                               'pass_count': row['pass_count'],
                               'max_center_residual': row['max_center_residual']}
                              for row in candidates],
              'scope': 'Existing group amplitudes multiplied by three independent scalars. Residual is exact quadratic recombination of precomputed group kinematics at one time/radius. Cone tests and residual norms are sampled only; parameter improvement is not a global NS or support certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_cone_parameter_screen.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'current': current_result,
                      'best': candidates[:5]}), flush=True)


if __name__ == '__main__':
    run()
