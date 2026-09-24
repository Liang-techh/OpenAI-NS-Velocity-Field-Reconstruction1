"""Local stress-cone and Kelvin covariance screen at the current radial peak.

The cone is evaluated from the physical corrected field, unlike the paper's
normalized leading background. Passing it would not establish a supported wave.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import nnls

from affine_pulse import affine_pulse
from compact_potential import CompactPotentialField
from joined_field import ROOT, coordinates
from joint_collar_fit import kinematics
from joint_collar_scale import ScaledJointCollarField
from radial_repair_refit import IncrementalRadialRepairField
from radial_swirl_fit import RadialSwirlRepairedField


def current_field():
    first = json.loads((ROOT/'compact_potential'/'joint_collar_fit.json').read_text())
    second = json.loads((ROOT/'compact_potential'/'radial_swirl_fit.json').read_text())
    third = json.loads((ROOT/'compact_potential'/'radial_repair_refit.json').read_text())
    base = CompactPotentialField()
    joint = ScaledJointCollarField(.1*np.array(first['linear_fit_amplitudes']),
                                   first['tau'], base=base)
    swirl = RadialSwirlRepairedField(joint,
                                    second['linear_fit_amplitudes'], first['tau'])
    return IncrementalRadialRepairField(swirl,
                                       .1*np.array(third['increment_amplitudes']),
                                       first['tau'])


def operator(field, points, tau):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    velocity, gradient, part = kinematics(field, points, tau, hs, ht)
    return velocity, gradient, part + np.einsum('nij,nj->ni', gradient, velocity)


def stress_primitive(field, radius, z, tau, order=12):
    q = float(coordinates(0., z/np.sqrt(field.nu), tau,
                          field.base.base.inner.h)['q'])
    ri = np.sqrt(field.nu)*np.sqrt(2*q*3/64)
    edges = sorted(set(np.clip([0., ri, radius], 0., radius)))
    g, w = leggauss(order)
    theta_integral = 0.
    axial_integral = 0.
    for lo, hi in zip(edges[:-1], edges[1:]):
        rr = (hi+lo)/2 + (hi-lo)/2*g
        weights = (hi-lo)/2*w
        points = np.column_stack((rr, np.zeros(order), np.full(order, z)))
        _, _, residual = operator(field, points, tau)
        theta_integral += float(np.dot(weights*rr**2, residual[:, 1]))
        axial_integral += float(np.dot(weights*rr, residual[:, 2]))
    return np.array([-theta_integral/radius**2, -axial_integral/radius])


def run():
    field = current_field()
    tau = .5/64
    radius = .0056890761915166545
    z = .003432627453438968
    point = np.array([[radius, 0., z]])
    velocity, gradient, residual = operator(field, point, tau)
    u = velocity[0]
    J = gradient[0]
    target = stress_primitive(field, radius, z, tau)
    F = u[1]/radius
    shear = np.array([J[1, 0]-F, J[2, 0]])
    shear_norm = np.linalg.norm(shear)
    N = shear/shear_norm
    K = np.array([-N[1], N[0]])
    lambda_squared = -2*F*N[0]*(2*F*N[0]+shear_norm)
    cone = {'F': float(F), 'shear': shear.tolist(),
            'lambda_squared': float(lambda_squared),
            'target_dot_N': float(target@N),
            'target_dot_K': float(target@K)}
    if lambda_squared > 0 and abs(2*F*N[0]) > 1e-14 and abs(target@N) > 1e-14:
        c = np.sqrt(lambda_squared)/(2*F*N[0])
        cone['c'] = float(c)
        cone['ratio'] = float(abs(c*(target@K)/(target@N)))
        cone['strict_pass'] = bool((target@N) < 0 and cone['ratio'] < 1)
    else:
        cone['strict_pass'] = False
    duration = .1*tau
    rows = []
    for mode in (1, 4, 16):
        for axial in (-1, 0, 1):
            for radial in (-1, 0, 1):
                tangential = np.array([mode/radius, axial*mode/radius])
                normal = np.r_[radial*np.linalg.norm(tangential), tangential]
                pulse = affine_pulse(J, field.nu, normal, duration)
                rows.append({'mode': mode, 'axial_ratio': axial,
                             'radial_ratio': radial, **pulse})
    matrix = np.array([r['mean_radial_tangential_covariance'] for r in rows]).T
    weights, error = nnls(matrix, target)
    selected = np.flatnonzero(weights > 1e-8)
    report = {'tau': tau, 'point': point[0].tolist(),
              'velocity': u.tolist(), 'gradient': J.tolist(),
              'residual': residual[0].tolist(),
              'local_tangential_stress_primitive': target.tolist(),
              'cone': cone,
              'kelvin_pulses': rows,
              'covariance_nnls': {
                  'relative_error': float(error/np.linalg.norm(target)),
                  'selected': [{'index': int(i), 'weight': float(weights[i])}
                               for i in selected],
                  'achieved': (matrix@weights).tolist(),
              },
              'scope': 'Physical local analogue of paper Section 7 cone, not its normalized leading theorem. Primitive integrates current full residual only from axis to selected radius; nonzero total radial moments prevent global compact stress. Frozen Cartesian Kelvin rays omit spatial support, exact curl, phase periodicity, amplitude PDE and nonlinear full-field momentum.',
              'pde_validated': False, 'global_field_ready': False}
    out = ROOT/'compact_potential'/'radial_peak_cone.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({key: report[key] for key in
                      ('point', 'residual', 'local_tangential_stress_primitive',
                       'cone', 'covariance_nnls')}), flush=True)


if __name__ == '__main__':
    run()
