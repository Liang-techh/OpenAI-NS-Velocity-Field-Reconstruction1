"""Independent full-operator point check of incremental radial repair."""
import json

import numpy as np

from compact_potential import CompactPotentialField
from joined_field import ROOT
from joint_collar_fit import kinematics
from joint_collar_scale import ScaledJointCollarField
from radial_repair_refit import IncrementalRadialRepairField
from radial_swirl_fit import RadialSwirlRepairedField


def full_residual(field, points, tau):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    velocity, grad, part = kinematics(field, points, tau, hs, ht)
    return part + np.einsum('nij,nj->ni', grad, velocity), np.trace(grad, axis1=1, axis2=2)


def run():
    first = json.loads((ROOT/'compact_potential'/'joint_collar_fit.json').read_text())
    second = json.loads((ROOT/'compact_potential'/'radial_swirl_fit.json').read_text())
    third = json.loads((ROOT/'compact_potential'/'radial_repair_refit.json').read_text())
    base = CompactPotentialField()
    joint = ScaledJointCollarField(.1*np.array(first['linear_fit_amplitudes']),
                                   first['tau'], base=base)
    current = RadialSwirlRepairedField(joint,
                                      second['linear_fit_amplitudes'], first['tau'])
    candidate = IncrementalRadialRepairField(current,
                                            .1*np.array(third['increment_amplitudes']),
                                            first['tau'])
    tau = .5/64
    points = np.array([[r, 0., z] for z in (.003432627453438968, -.003432627453438968)
                       for r in (.0024380588399712513, .0056890761915166545,
                                 .009790979681828347)])
    before, div0 = full_residual(current, points, tau)
    after, div1 = full_residual(candidate, points, tau)
    report = {
        'tau': tau, 'points': points.tolist(),
        'before_norms': np.linalg.norm(before, axis=1).tolist(),
        'after_norms': np.linalg.norm(after, axis=1).tolist(),
        'before_radial': before[:, 0].tolist(),
        'after_radial': after[:, 0].tolist(),
        'before_divergence': div0.tolist(),
        'after_divergence': div1.tolist(),
        'scope': 'Independent full Cartesian finite differences at six late axial-collar points. Fourth-order one-sided time derivative at registered endpoint. Not a global gate.',
        'accepted': False,
    }
    out = ROOT/'compact_potential'/'radial_repair_direct.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
