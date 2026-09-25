"""Decompose full momentum at the registered moving collar hotspots."""
import json

import numpy as np

from joined_field import ROOT
from joint_collar_fit import kinematics
from local_poloidal_basis_screen import load_robust_candidate


def run():
    field = load_robust_candidate(.1)
    source = json.loads((ROOT/'compact_potential'/'forcing_extension_screen.json').read_text())
    rows = []
    for sample in source['rows']:
        tau = sample['tau']
        points = np.array([sample['worst_point']])
        hs = .0005*np.sqrt(field.nu*tau)
        ht = .0001*tau
        _, _, _, terms, viscosity_axes = kinematics(
            field, points, tau, hs, ht, return_terms=True,
            return_viscosity_axes=True)
        vectors = {key: value[0].tolist() for key, value in terms.items()}
        residual = sum(terms.values())[0]
        row = {'tau': tau, 'point': points[0].tolist(),
               'terms': vectors, 'residual': residual.tolist(),
               'viscosity_axes_xyz': viscosity_axes[:, 0, :].tolist(),
               'residual_norm': float(np.linalg.norm(residual)),
               'pressure_to_nonpressure_x': float(
                   terms['pressure'][0, 0]/sum(
                       terms[key][0, 0] for key in
                       ('time', 'convection', 'viscosity')))}
        rows.append(row)
        print(json.dumps({'tau': tau, 'residual': row['residual'],
                          'terms_x': {key: value[0] for key, value in vectors.items()}}),
              flush=True)
    report = {'rows': rows, 'scope': 'Selected worst points from the earlier registered eight-point collar screen. Finite-difference decomposition of the physical momentum operator; not a global maximum or proof of asymptotic behavior.',
              'accepted': False}
    path = ROOT/'compact_potential'/'collar_residual_terms.json'
    path.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
