"""Independent full finite-difference spot check of callable collar candidate."""
import json

import numpy as np

from compact_potential import CompactPotentialField
from joined_field import ROOT, independent_fd
from joint_collar_scale import ScaledJointCollarField


def run():
    fit = json.loads((ROOT/'compact_potential'/'joint_collar_fit.json').read_text())
    amplitudes = .1*np.array(fit['linear_fit_amplitudes'])
    base = CompactPotentialField()
    field = ScaledJointCollarField(amplitudes, fit['tau'], base=base)
    tau = fit['tau']
    _, radius, zflat, zsupp = base.support(tau)
    z = (zflat+zsupp)/2
    points = np.array([[fraction*radius, 0., z] for fraction in (.1, .5, .9)])
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    baseline, div0 = independent_fd(base, points, tau, hs, ht)
    candidate, div1 = independent_fd(field, points, tau, hs, ht)
    report = {
        'tau': tau, 'points': points.tolist(),
        'baseline_momentum_norms': np.linalg.norm(baseline, axis=1).tolist(),
        'candidate_momentum_norms': np.linalg.norm(candidate, axis=1).tolist(),
        'baseline_divergence': div0.tolist(),
        'candidate_divergence': div1.tolist(),
        'scope': 'Independent fourth-order Cartesian momentum and divergence spot check of the callable candidate at three collar points. Not a full-domain gate.',
        'accepted': False,
    }
    out = ROOT/'compact_potential'/'joint_collar_direct.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
