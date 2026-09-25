"""High-order cone check at the newly shaped U onset."""

import json

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from extended_physical_cone_map import physical_cone
from radial_continuation import ROOT


def run():
    slice_name = 'delayed005_rise146_degree31.json'
    tau = .5*2**(-5.5)
    field = CoupledMomentPhysicalLift(slice_filename=slice_name)
    rows = [physical_cone(field, 1.005, .2, tau, order=64)]
    rows += [physical_cone(field, 1.005, .3, tau, order=order)
             for order in (64, 96, 128)]
    report = dict(slice_filename=slice_name, tau=tau, rows=rows,
                  scope='Two fixed-slice onset points with piecewise '
                        'Gauss64 and eta=.3 order escalation to 128 for '
                        'the physical residual primitive. Passing '
                        'here does not imply a continuous cone or wave '
                        'admission.', accepted=False)
    (ROOT/'delayed_rise146_cone_refine.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(rows), flush=True)


if __name__ == '__main__':
    run()
