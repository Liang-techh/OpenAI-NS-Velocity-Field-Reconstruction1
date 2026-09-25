"""Find the smallest moment-step magnitude retaining the physical cone."""
import json

from extended_physical_cone_map import physical_cone
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def run():
    tau = .5*2**(-5.5)
    rows = []
    for amplitude in (-.5, -.45, -.4, -.35, -.3):
        field = RadialMomentStep(make_field(16, 2.), amplitude)
        for eta in (.2, .3):
            row = dict(amplitude=amplitude,
                       **physical_cone(field, 1., eta, tau, order=24))
            rows.append(row)
            print(json.dumps({key: row[key] for key in
                              ('amplitude', 'eta', 'strict_pass',
                               'target_dot_N', 'ratio',
                               'lambda_squared')}), flush=True)
    report = dict(tau=tau, rows=rows,
                  scope='Ten physical cone samples at X=1 with Gauss24 '
                        'stress integration, before five-moment repair. '
                        'No continuous cone or full PDE admission.',
                  accepted=False)
    (ROOT/'radial_step_cone_amplitude.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
