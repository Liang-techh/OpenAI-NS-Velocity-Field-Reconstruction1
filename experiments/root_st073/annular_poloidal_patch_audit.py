"""Off-quadrature radial-axial patch check of multigrid annular candidate."""
import json

import numpy as np

from annular_poloidal_multigrid_fit import load_candidate
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


def evaluate(field, points, tau):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, J, part = kinematics(field, points, tau, hs, ht)
    return part+np.einsum('nij,nj->ni', J, u)


def run():
    tau = .5/64
    radii = np.array([.005, .0065, .008, .0095, .011])
    heights = np.array([.00285, .0031, .00335, .0036, .00385])
    points = np.array([[r, 0., sign*z]
                       for sign in (-1, 1) for z in heights for r in radii])
    current_R = evaluate(current_field(), points, tau)
    candidate_R = evaluate(load_candidate(), points, tau)
    dr, dz = radii[1]-radii[0], heights[1]-heights[0]
    weights = 2*np.pi*points[:, 0]*dr*dz
    def metrics(R):
        n = np.linalg.norm(R, axis=1)
        return {'max': float(np.max(n)),
                'physical_patch_l2': float(np.sqrt(weights @ n**2)),
                'rms': float(np.sqrt(np.mean(n**2))),
                'max_radial': float(np.max(np.abs(R[:, 0])))}
    report = {'tau': tau, 'radii': radii.tolist(),
              'heights': heights.tolist(),
              'baseline': metrics(current_R),
              'candidate': metrics(candidate_R),
              'pointwise_improved_count': int(sum(np.linalg.norm(candidate_R, axis=1)
                                                  < np.linalg.norm(current_R, axis=1))),
              'total_points': len(points),
              'scope': 'Off-Gauss tensor patch with rectangle physical-volume weights. Sampled diagnostic, not an integral/supremum certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_poloidal_patch_audit.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
