"""Coarse axial relocation screen for the physical local stress cone.

The cone here is the same local analogue as radial_peak_cone.py, not the
normalized leading cone in the paper. These points are design candidates,
not supported pulses or globally admissible stresses.
"""
import json

import numpy as np

from compact_potential import CompactPotentialField
from joined_field import ROOT
from radial_peak_cone import current_field, operator, stress_primitive


def run():
    field = current_field()
    tau = .5/64
    radius = .0056890761915166545
    original_z = .003432627453438968
    _, _, _, zsupport = CompactPotentialField().support(tau)
    rows = []
    for z in (0., .0015, .0025, .00275, .003, .00325, original_z):
        point = np.array([[radius, 0., z]])
        u, gradient, residual = operator(field, point, tau)
        target = stress_primitive(field, radius, z, tau, order=10)
        F = u[0, 1]/radius
        shear = np.array([gradient[0, 1, 0]-F, gradient[0, 2, 0]])
        N = shear/np.linalg.norm(shear)
        K = np.array([-N[1], N[0]])
        lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
        row = {'z': z, 'axial_symmetric_support_halfwidth_ceiling': float(zsupport-abs(z)),
               'velocity': u[0].tolist(),
               'residual_norm': float(np.linalg.norm(residual[0])),
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
    report = {'tau': tau, 'radius': radius, 'base_axial_support_halfwidth': float(zsupport),
              'rows': rows,
              'scope': 'Fixed-radius coarse axial screen of the physical local stress-cone analogue. Uses the full corrected baseline residual primitive, not the paper normalized leading stress. Local cone passing does not imply a realizable supported wave or global Navier--Stokes solution.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_location_screen.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
