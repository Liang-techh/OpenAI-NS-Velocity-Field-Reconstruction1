"""Refine full-momentum finite differences across the narrow U onset."""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, stats
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import current_mean
import delayed_remote_u_physical as physical
from joined_field import independent_fd
from radial_continuation import ROOT


def run(smooth=False):
    if smooth:
        physical.START, physical.WIDTH, physical.DEGREE = 1.01, .02, 11
    source_name = ('delayed_remote_u_physical_smooth.json' if smooth
                   else 'delayed_remote_u_physical.json')
    source = json.loads((ROOT/source_name).read_text())
    e_source = json.loads((ROOT/'delayed_remote_positive_e.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    field = physical.RemoteUPatchedField(mean, base,
        [tuple(pair) for pair in e_source['radial_intervals']],
        source['slice']['e_coefficients'],
        source['slice']['u_coefficients'])
    tau = .5*2**(-5.5)
    xs = ((1.015,1.02) if smooth else (1.031,1.04))
    points, X, eta = nodes(base, xs, (.3,), tau)
    rows = []
    for scale in (1., .5, .25, .125):
        hs = scale*.001*np.sqrt(base.nu*tau)
        ht = scale*.00025*tau
        R, div = independent_fd(field, points, tau, hs, ht)
        rows.append(dict(scale=scale, hs=hs, ht=ht,
                         momentum=stats(R),
                         norms=np.linalg.norm(R,axis=1).tolist(),
                         divergence=div.tolist()))
    report = dict(source=source_name,
                  X=X.tolist(), eta=eta.tolist(), rows=rows,
                  scope='Finite-difference refinement at two narrow '
                        'onset nodes. It tests discretization '
                        'sensitivity but does not establish a '
                        'continuous residual bound.', accepted=False)
    output = ('delayed_remote_u_fd_refinement_smooth.json' if smooth
              else 'delayed_remote_u_fd_refinement.json')
    (ROOT/output).write_bytes(
        (json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps([(row['scale'],row['norms'],row['divergence'])
                      for row in rows]),flush=True)


if __name__ == '__main__':
    import sys
    run(smooth='--smooth' in sys.argv[1:])
