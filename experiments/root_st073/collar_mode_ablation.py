"""Compare baseline and local-poloidal correction at registered collar hotspots."""
import json

import numpy as np

from joined_field import ROOT
from joint_collar_fit import kinematics
from local_poloidal_basis_screen import load_robust_candidate


def run():
    source = json.loads((ROOT/'compact_potential'/'forcing_extension_screen.json').read_text())
    rows = []
    for sample in source['rows'][:5]:
        tau = sample['tau']
        points = np.array([sample['worst_point']])
        cases = {}
        for strength in (0., .1):
            field = load_robust_candidate(strength)
            hs = .0005*np.sqrt(field.nu*tau)
            ht = .0001*tau
            _, _, _, terms = kinematics(field, points, tau, hs, ht,
                                        return_terms=True)
            residual = sum(terms.values())[0]
            cases[str(strength)] = {
                'residual': residual.tolist(),
                'residual_norm': float(np.linalg.norm(residual)),
                'viscosity': terms['viscosity'][0].tolist(),
            }
        rows.append({'tau': tau, 'point': points[0].tolist(),
                     'cases': cases,
                     'residual_norm_ratio_corrected_to_base':
                         cases['0.1']['residual_norm']/cases['0.0']['residual_norm']})
        print(json.dumps({'tau': tau, 'base': cases['0.0']['residual_norm'],
                          'corrected': cases['0.1']['residual_norm']}), flush=True)
    report = {'rows': rows, 'scope': 'Same selected collar points for baseline and robust local-poloidal strength 0.1. The points were selected using corrected-field residual, so this is an ablation, not a fair maximum comparison.',
              'accepted': False}
    (ROOT/'compact_potential'/'collar_mode_ablation.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
