"""Check physical stress cone on the delayed-start five-moment lift."""
import argparse
import json

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from extended_physical_cone_map import physical_cone
from radial_continuation import ROOT


def run(slice_filename='delayed_taper_curvature_optimize.json',
        output_name='delayed_taper_cone_screen.json'):
    field = CoupledMomentPhysicalLift(
        slice_filename=slice_filename)
    tau = .5*2**(-5.5)
    rows = []
    for eta in (.2, .3):
        for X in (1., 1.005, 1.01, 1.015, 1.025):
            order = 24 if X in (1.005, 1.01, 1.015) else 12
            row = physical_cone(field, X, eta, tau, order=order)
            rows.append(row)
            print(json.dumps({key: row[key] for key in
                              ('X', 'eta', 'order', 'lambda_squared',
                               'target_dot_N', 'ratio', 'strict_pass',
                               'residual_norm')}), flush=True)
    report = dict(tau=tau, start_X=field.start_X, rows=rows,
                  scope='Ten sampled physical full-residual cone nodes '
                        'near delayed U start. Order 24 on sensitive '
                        'nodes, with quadrature split at U/E support '
                        'edges. No continuous cone or wave admission.',
                  accepted=False)
    (ROOT/output_name).write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--slice-filename',
                        default='delayed_taper_curvature_optimize.json')
    parser.add_argument('--output-name',
                        default='delayed_taper_cone_screen.json')
    arguments = parser.parse_args()
    run(arguments.slice_filename, arguments.output_name)
