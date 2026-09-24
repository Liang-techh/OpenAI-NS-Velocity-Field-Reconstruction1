"""Independent physical finite differences in the radial cutoff collar."""
import json

import numpy as np

from compact_potential import CompactPotentialField
from joined_field import ROOT, independent_fd
from joint_collar_scale import ScaledJointCollarField
from radial_swirl_fit import RadialSwirlRepairedField


def run():
    fit = json.loads((ROOT/'compact_potential'/'joint_collar_fit.json').read_text())
    swirl = json.loads((ROOT/'compact_potential'/'radial_swirl_fit.json').read_text())
    base = CompactPotentialField()
    joint = ScaledJointCollarField(.1*np.array(fit['linear_fit_amplitudes']),
                                   fit['tau'], base=base)
    candidate = RadialSwirlRepairedField(joint,
                                        swirl['linear_fit_amplitudes'], fit['tau'])
    tau = fit['tau']
    rflat, rsupp, zflat, _ = base.support(tau)
    points = np.array([[rflat+y*(rsupp-rflat), 0., z]
                       for z in (0., .7*zflat) for y in (.25, .5, .75)])
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    before, div0 = independent_fd(joint, points, tau, hs, ht)
    after, div1 = independent_fd(candidate, points, tau, hs, ht)
    report = {
        'tau': tau, 'points': points.tolist(),
        'before_momentum_norms': np.linalg.norm(before, axis=1).tolist(),
        'after_momentum_norms': np.linalg.norm(after, axis=1).tolist(),
        'before_angular_residual': before[:, 1].tolist(),
        'after_angular_residual': after[:, 1].tolist(),
        'before_divergence': div0.tolist(),
        'after_divergence': div1.tolist(),
        'scope': 'Independent fourth-order Cartesian full momentum and divergence at six radial-collar points. Not a global max or L2 gate.',
        'accepted': False,
    }
    out = ROOT/'compact_potential'/'radial_swirl_direct.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
