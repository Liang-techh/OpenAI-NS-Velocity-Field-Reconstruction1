"""Physical stress primitives of the remote moment-matched mean field."""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from joined_field import ROOT
from joint_collar_fit import kinematics
from moment_matched_joined_field import MomentMatchedJoinedField
from paper_moment_bridge import XI, XB


def cone_data(radius, velocity, gradient, target):
    F = velocity[1]/radius
    shear = np.array([gradient[1, 0]-F, gradient[2, 0]])
    shear_norm = float(np.linalg.norm(shear))
    if shear_norm == 0:
        return {'F': float(F), 'shear': shear.tolist(),
                'lambda_squared': None, 'target_dot_N': None,
                'cone_ratio': None, 'strict_local_pass': False}
    N = shear/shear_norm
    K = np.array([-N[1], N[0]])
    lam2 = -2*F*N[0]*(2*F*N[0]+shear_norm)
    dotN = float(target@N)
    out = {'F': float(F), 'shear': shear.tolist(),
           'lambda_squared': float(lam2),
           'target_dot_N': dotN, 'target_dot_K': float(target@K)}
    if lam2 > 0 and abs(2*F*N[0]) > 1e-14 and abs(dotN) > 1e-14:
        c = np.sqrt(lam2)/(2*F*N[0])
        out['cone_ratio'] = float(abs(c*(target@K)/dotN))
        out['strict_local_pass'] = bool(dotN < 0 and out['cone_ratio'] < 1)
    else:
        out['cone_ratio'] = None
        out['strict_local_pass'] = False
    return out


def points_for_eta(field, eta, tau, order=8):
    q = tau/(1-eta**2)
    radius = lambda X: np.sqrt(2*field.nu*q*X)
    z = np.sqrt(field.nu)*q**field.inner.D*eta
    g, w = leggauss(order)
    edges = (0., XI, XB, 4., 6., 10., 14.)
    parts = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        a, b = radius(lo), radius(hi)
        radii = (a+b)/2+(b-a)/2*g
        weights = (b-a)/2*w
        parts.append((hi, radii, weights))
    all_r = np.r_[*[part[1] for part in parts],
                  radius(6.), radius(10.), radius(14.)]
    points = np.column_stack((all_r, np.zeros(len(all_r)),
                              np.full(len(all_r), z)))
    return parts, points


def summarize(field, parts, points, tau):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    velocity, gradient, linear = kinematics(field, points, tau, hs, ht)
    residual = linear+np.einsum('nij,nj->ni', gradient, velocity)
    rows = []
    integration_count = len(points)-3
    for j, targetX in enumerate((6., 10., 14.)):
        target_radius = points[integration_count+j, 0]
        theta_integral, axial_integral = 0., 0.
        offset = 0
        for upperX, radii, weights in parts:
            count = len(radii)
            if upperX <= targetX:
                theta_integral += float(np.dot(weights*radii**2,
                                               residual[offset:offset+count, 1]))
                axial_integral += float(np.dot(weights*radii,
                                               residual[offset:offset+count, 2]))
            offset += count
        target = np.array([-theta_integral/target_radius**2,
                           -axial_integral/target_radius])
        row = {'X': targetX, 'radius': float(target_radius),
               'stress_target': target.tolist(),
               'center_residual_norm': float(np.linalg.norm(
                   residual[integration_count+j]))}
        row.update(cone_data(target_radius,
                             velocity[integration_count+j],
                             gradient[integration_count+j], target))
        rows.append(row)
    return rows, float(np.max(np.linalg.norm(residual, axis=1)))


def run():
    patched = MomentMatchedJoinedField()
    tau = .0084
    rows = []
    for eta in (0., .2, .35):
        parts, points = points_for_eta(patched, eta, tau)
        baseline, base_max = summarize(patched.base, parts, points, tau)
        corrected, patched_max = summarize(patched, parts, points, tau)
        row = {'eta': eta, 'tau': tau, 'baseline': baseline,
               'moment_patched': corrected,
               'baseline_sample_max_residual': base_max,
               'patched_sample_max_residual': patched_max}
        rows.append(row)
        print(json.dumps({'eta': eta,
                          'baseline_pass_count': sum(x['strict_local_pass'] for x in baseline),
                          'patched_pass_count': sum(x['strict_local_pass'] for x in corrected),
                          'patched_targets': [x['stress_target'] for x in corrected]}),
              flush=True)
    report = {'rows': rows, 'radius_X': [6., 10., 14.],
              'scope': 'Physical full-residual radial stress primitives and paper-inspired local cone analog in the remote patch at one registered time. Piecewise eight-node radial quadrature, three eta slices. This is not the paper normalized leading stress cone, no wave or force realization, and no continuum certificate.',
              'accepted': False}
    (ROOT/'remote_moment_stress_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
