"""Batched fixed-time (r,z) screen for local wave stress-cone locations.

This uses the physical residual primitive of the current ST073 baseline. It
does not certify the paper's normalized leading cone or construct a pulse.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from compact_potential import CompactPotentialField
from joined_field import ROOT, coordinates
from radial_peak_cone import current_field, inner_similarity_exponent, operator


def evaluate_height(field, radii, z, tau, r_support, z_support,
                    quadrature_order=8):
    q = float(coordinates(0., z/np.sqrt(field.nu), tau,
                          inner_similarity_exponent(field))['q'])
    inner_radius = np.sqrt(field.nu)*np.sqrt(2*q*3/64)
    edges = sorted(set([0., inner_radius]+list(radii)))
    nodes, weights = leggauss(quadrature_order)
    segment_radii, segment_weights, segment_ends = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi > max(radii):
            break
        rr = (lo+hi)/2+(hi-lo)/2*nodes
        segment_radii.extend(rr)
        segment_weights.extend((hi-lo)/2*weights)
        segment_ends.append(hi)
    segment_radii = np.asarray(segment_radii)
    segment_weights = np.asarray(segment_weights)
    all_radii = np.r_[segment_radii, radii]
    points = np.column_stack((all_radii, np.zeros(len(all_radii)),
                              np.full(len(all_radii), z)))
    velocity, gradient, residual = operator(field, points, tau)
    primitive = {}
    angular_integral = axial_integral = 0.
    for i, end in enumerate(segment_ends):
        sl = slice(i*quadrature_order, (i+1)*quadrature_order)
        rr = segment_radii[sl]
        ww = segment_weights[sl]
        angular_integral += float(np.dot(ww*rr**2, residual[sl, 1]))
        axial_integral += float(np.dot(ww*rr, residual[sl, 2]))
        primitive[end] = np.array([-angular_integral/end**2,
                                   -axial_integral/end])
    offset = len(segment_radii)
    rows = []
    for i, radius in enumerate(radii):
        u = velocity[offset+i]
        J = gradient[offset+i]
        target = primitive[radius]
        F = u[1]/radius
        shear = np.array([J[1, 0]-F, J[2, 0]])
        shear_norm = np.linalg.norm(shear)
        N = shear/shear_norm
        K = np.array([-N[1], N[0]])
        lam2 = -2*F*N[0]*(2*F*N[0]+shear_norm)
        row = {'r': float(radius), 'z': float(z),
               'radial_symmetric_support_halfwidth_ceiling': float(min(radius, r_support-radius)),
               'axial_symmetric_support_halfwidth_ceiling': float(z_support-abs(z)),
               'residual_norm': float(np.linalg.norm(residual[offset+i])),
               'target': target.tolist(),
               'lambda_squared': float(lam2),
               'target_dot_N': float(target@N),
               'target_dot_K': float(target@K)}
        if lam2 > 0 and abs(2*F*N[0]) > 1e-14 and abs(target@N) > 1e-14:
            c = np.sqrt(lam2)/(2*F*N[0])
            row['cone_ratio'] = float(abs(c*(target@K)/(target@N)))
            row['strict_local_pass'] = bool(target@N < 0 and row['cone_ratio'] < 1)
        else:
            row['strict_local_pass'] = False
        rows.append(row)
    return rows


def run():
    field = current_field()
    tau = .5/64
    original_r = .0056890761915166545
    original_z = .003432627453438968
    radii = np.array([.003, .004, original_r, .0075, .0095, .012, .016])
    heights = np.array([0., .0015, .0025, .00275, .003, .00325,
                        original_z, .0036])
    _, r_support, _, z_support = CompactPotentialField().support(tau)
    rows = []
    for z in heights:
        rows.extend(evaluate_height(field, radii, z, tau,
                                    r_support, z_support))
    passing = [row for row in rows if row['strict_local_pass']]
    passing.sort(key=lambda row: (row['axial_symmetric_support_halfwidth_ceiling'],
                                   row['radial_symmetric_support_halfwidth_ceiling']),
                 reverse=True)
    reference = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    center = next(row for row in rows
                  if row['r'] == original_r and row['z'] == original_z)
    reference_target = np.array(reference['local_tangential_stress_primitive'])
    report = {'tau': tau, 'radii': radii.tolist(), 'heights': heights.tolist(),
              'r_support': float(r_support), 'z_support': float(z_support),
              'quadrature_order_per_radial_segment': 8,
              'reference_target_relative_difference': float(np.linalg.norm(np.array(center['target'])-reference_target)/np.linalg.norm(reference_target)),
              'passing_count': len(passing), 'total_count': len(rows),
              'passing_by_axial_room': passing,
              'rows': rows,
              'scope': 'Coarse physical local cone screen at one time. Support ceilings are geometric only; cone passing at a point does not establish a supported positive covariance through a patch, exact pressure/mean correction, or a global NS solution.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_cone_region.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'passing_count': len(passing),
                      'total_count': len(rows),
                      'reference_target_relative_difference': report['reference_target_relative_difference'],
                      'best': passing[:8]}), flush=True)


if __name__ == '__main__':
    run()
