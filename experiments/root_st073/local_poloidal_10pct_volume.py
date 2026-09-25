"""Direct physical-momentum check of the cone-compatible 10% poloidal fit."""
import json

import numpy as np

from joined_field import ROOT
from joint_collar_fit import nodes
from local_poloidal_basis_screen import load_robust_candidate
from radial_peak_cone import operator


def run():
    field = load_robust_candidate(.1)
    source = json.loads((ROOT/'compact_potential'/'local_poloidal_basis_screen.json').read_text())
    baseline = {row['order']: row for row in source['rows'][0]['metrics']}
    tau = .0084
    rows = []
    for order in (6, 8, 10):
        points, weights = nodes(field, tau, order)
        residual = operator(field, points, tau)[2]
        norms = np.linalg.norm(residual, axis=1)
        row = {'tau': tau, 'order': order,
               'max': float(np.max(norms)),
               'physical_volume_l2': float(np.sqrt(weights@norms**2)),
               'baseline': baseline[order]}
        row['max_ratio'] = row['max']/baseline[order]['max']
        row['l2_ratio'] = row['physical_volume_l2']/baseline[order]['physical_volume_l2']
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = {'strength': .1, 'rows': rows,
              'scope': 'Gauss6/8/10 full physical finite-difference momentum at one time. Not a converged spatial or spacetime bound.',
              'accepted': False}
    path = ROOT/'compact_potential'/'local_poloidal_10pct_volume.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
