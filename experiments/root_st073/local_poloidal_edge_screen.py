"""Find a nonzero collar correction strength preserving the cone edge."""
import json

import numpy as np

from curl_wave_cone_region import evaluate_height
from joined_field import ROOT
from local_poloidal_basis_screen import load_robust_candidate


def run():
    rows = []
    for strength in (0., .1, .2, .3, .5, .75, 1.):
        field = load_robust_candidate(strength)
        for tau in (.0084, .008475, .00855):
            _, r_support, _, z_support = field.support(tau)
            row = evaluate_height(field, np.array([.012]), .0037, tau,
                                  r_support, z_support, quadrature_order=8)[0]
            rows.append({'strength': strength, 'tau': tau,
                         'cone_ratio': row.get('cone_ratio'),
                         'strict_local_pass': row['strict_local_pass']})
            print(json.dumps(rows[-1]), flush=True)
    report = {'rows': rows,
              'scope': 'One previously failed spatial edge at three times. Passing here does not verify the entire spacetime grid or global momentum.',
              'accepted': False}
    path = ROOT/'compact_potential'/'local_poloidal_edge_screen.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
