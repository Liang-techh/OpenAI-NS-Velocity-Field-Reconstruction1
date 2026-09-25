"""Attribute the delayed lift's entrance momentum to E and U changes."""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    tau = .5*2**(-5.5)
    xs = np.array([1.01, 1.0125, 1.015, 1.0175, 1.02])
    points = CoupledMomentPhysicalLift(
        slice_filename='delayed005_wide04_curvature_optimize.json'
    ).compact.joined.inner.from_similarity(xs, np.full(len(xs), .3), tau)
    rows = []
    for label, use_E, use_U in (
        ('base', False, False), ('E_only', True, False),
        ('U_only', False, True), ('coupled', True, True)):
        field = CoupledMomentPhysicalLift(
            slice_filename='delayed005_wide04_curvature_optimize.json')
        if not use_E:
            field.e_rows[:] = 0.
        if not use_U:
            field.u_rows[:] = 0.
        residual, divergence = independent_fd(
            field, points, tau, .0005*np.sqrt(field.nu*tau), .00025*tau)
        row = dict(name=label,
                   nodes=[dict(X=float(x), residual=R.tolist(),
                               norm=float(np.linalg.norm(R)),
                               divergence=float(d))
                          for x, R, d in zip(xs, residual, divergence)])
        rows.append(row)
        print(json.dumps(dict(name=label,
                              norms=[item['norm'] for item in row['nodes']],
                              hotspot_residual=row['nodes'][2]['residual'])),
              flush=True)
    report = dict(tau=tau, eta=.3, rows=rows,
                  scope='Five physical entrance nodes with E/U corrections '
                        'separately toggled. Off-manifold variants are '
                        'diagnostic only and do not match five moments.',
                  accepted=False)
    (ROOT/'delayed_wide04_component_decomposition.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
