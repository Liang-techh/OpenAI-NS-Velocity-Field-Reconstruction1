"""Measure an independent M-moment degree of freedom and its PDE cost."""
import json

import numpy as np

from extended_physical_cone_map import physical_cone
from extended_relaxed_cone_screen import moment_data
from high_frequency_shear_screen import make_field
from joined_field import independent_fd
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 0.)
    one = RadialMomentStep(base, 1.)
    x_out = one.rise_end
    moment_rows = []
    for eta in (.2, .3):
        m0, _, _ = moment_data(base, x_out, eta, tau, 32)
        m1, _, _ = moment_data(one, x_out, eta, tau, 32)
        predicted = (1-eta**2)**(1-(.5+one.heat.h))
        moment_rows.append(dict(eta=eta, x_out=x_out,
                                delta_moments=(m1-m0).tolist(),
                                analytic_delta_M=predicted))
    cone_rows = []
    for amplitude in (-2., -1., 0., 1., 2.):
        field = RadialMomentStep(base, amplitude)
        row = dict(amplitude=amplitude,
                   **physical_cone(field, 1., .2, tau, order=32))
        cone_rows.append(row)
        print(json.dumps(row), flush=True)
    selected = min(cone_rows, key=lambda row: row['target_dot_N'])
    chosen = RadialMomentStep(base, selected['amplitude'])
    X, eta = np.meshgrid([.75, 1., 1.25, 2.25, 2.65], [.2, .3],
                         indexing='ij')
    points = base.compact.joined.inner.from_similarity(
        X.ravel(), eta.ravel(), tau)
    results = []
    for label, field in [('base', base), ('selected', chosen)]:
        R, div = independent_fd(
            field, points, tau, .001*np.sqrt(field.nu*tau),
            .00025*tau)
        results.append(dict(field=label,
                            sampled_full_max=float(np.max(np.linalg.norm(
                                R, axis=1))),
                            sampled_angular_max=float(np.max(np.abs(R[:, 1]))),
                            divergence_max=float(np.max(np.abs(div)))))
    report = dict(tau=tau, moment_rows=moment_rows,
                  physical_cone_rows=cone_rows,
                  selected_amplitude=selected['amplitude'],
                  residual_rows=results,
                  scope='Two-stage solenoidal streamfunction step changes M at the old radial exit and restores its physical streamfunction by X=3. Finite samples only; other moments, continuous cone, full momentum and volume-L2 gates remain open.',
                  accepted=False)
    (ROOT/'radial_moment_step_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(moment_rows=moment_rows,
                          selected_amplitude=selected['amplitude'],
                          residual_rows=results)), flush=True)


if __name__ == '__main__':
    run()
