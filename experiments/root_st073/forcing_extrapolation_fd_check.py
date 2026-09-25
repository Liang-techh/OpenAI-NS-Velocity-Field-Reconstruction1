"""FD sensitivity at the smallest extrapolated time and sampled hotspot."""
import json

import numpy as np

from joined_field import ROOT
from joint_collar_fit import kinematics
from local_poloidal_basis_screen import compact_base, load_robust_candidate


def run():
    field = load_robust_candidate(.1)
    compact = compact_base(field)
    compact.experimental_time_extension = True
    compact.base.experimental_time_extension = True
    source = json.loads((ROOT/'compact_potential'/'forcing_extension_extrapolation.json').read_text())
    row = source['rows'][0]
    tau = row['tau']
    points = np.array([row['worst_point']])
    rows = []
    for factor in (1., .5, .25):
        hs = factor*.0005*np.sqrt(field.nu*tau)
        ht = factor*.0001*tau
        u, J, part = kinematics(field, points, tau, hs, ht,
                                time_min=0.)
        residual = part+np.einsum('nij,nj->ni', J, u)
        rows.append({'factor': factor,
                     'residual': residual[0].tolist(),
                     'norm': float(np.linalg.norm(residual[0]))})
    report = {'tau': tau, 'point': points[0].tolist(),
              'rows': rows,
              'scope': 'One extrapolated point, three finite-difference spacings. Stability does not validate out-of-domain model formulas or an asymptotic force limit.',
              'accepted': False}
    path = ROOT/'compact_potential'/'forcing_extrapolation_fd_check.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
