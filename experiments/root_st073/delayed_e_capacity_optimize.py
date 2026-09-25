"""Reoptimize swirl E for the delayed, wide, degree-31 U repair basis."""
import json

from delayed_taper_capacity_screen import grid
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep
from taper_width_capacity import optimize_slice


def run():
    tau = .5*2**(-5.5)
    width, degree, start_X = .4, 31, 1.005
    X, weights = grid(order=64)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    old = json.loads((ROOT/'delayed005_wide04_curvature_optimize.json').read_text())
    rows = []
    for prior in old['rows']:
        row = optimize_slice(base, changed, X, weights, prior['eta'],
                             tau, width, prior['e_coefficients'],
                             degree=degree, start_X=start_X)
        rows.append(row)
        print(json.dumps(dict(eta=row['eta'],
                              new_S_slack=row['finite_basis_S_slack'],
                              min_relative_E=row['min_relative_E'],
                              feasible=row['feasible'],
                              optimizer_success=row['optimizer_success'])),
              flush=True)
    report = dict(tau=tau, width=width, degree=degree,
                  start_X=start_X, quadrature_per_piece=64,
                  rows=rows,
                  scope='Fixed-slice E redistribution maximizing S slack '
                        'for delayed U support, with I/Cp equalities and '
                        'positive E floor. No selected physical lift or '
                        'PDE admission.', accepted=False)
    (ROOT/'delayed_e_capacity_optimize.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
