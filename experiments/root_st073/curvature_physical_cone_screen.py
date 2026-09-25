"""Physical stress-cone continuity check on the curvature-optimized lift.

This is a necessary local wave-feasibility screen, not a constructed wave or
a continuous cone certificate. The target uses complete physical residuals.
"""
import json

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from extended_physical_cone_map import physical_cone
from radial_continuation import ROOT


def run():
    field = CoupledMomentPhysicalLift(
        slice_filename='wide_taper_curvature_optimize.json')
    tau = .5*2**(-5.5)
    rows = []
    for eta in (.2, .3):
        for X in (1., 1.005, 1.01, 1.025):
            row = physical_cone(field, X, eta, tau, order=12)
            rows.append(row)
            print(json.dumps({key: row[key] for key in
                              ('X', 'eta', 'lambda_squared', 'target_dot_N',
                               'ratio', 'strict_pass', 'residual_norm')}),
                  flush=True)
    refinement = []
    for eta in (.2, .3):
        for X in (1.005, 1.01):
            row = physical_cone(field, X, eta, tau, order=24)
            refinement.append(row)
            print(json.dumps(dict(refined_order=24, X=X, eta=eta,
                                  target_dot_N=row['target_dot_N'],
                                  ratio=row['ratio'],
                                  strict_pass=row['strict_pass'])), flush=True)
    report = dict(tau=tau, quadrature_order=12, rows=rows,
                  refinement_order=24, refinement=refinement,
                  scope='Eight sampled physical stress-cone nodes on the '
                        'curvature-optimized physical lift. The target is '
                        'an axis-to-radius primitive of complete physical '
                        'momentum residual; no continuous cone, wave, '
                        'mean closure, or PDE acceptance.',
                  accepted=False)
    (ROOT/'curvature_physical_cone_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
