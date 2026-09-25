"""Check spatial-stencil convergence for the rejected two-mode field."""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode
from delayed_similarity_curl_screen import CurlPatchedLift, SimilarityCurlMode
from delayed_swirl_cone_response import BASE_NAME, ETA_INTERVAL
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    source = json.loads((ROOT/'delayed_coupled_cone_response.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    modes = [SimilaritySwirlMode(
                 base, source['swirl_interval'], eta_interval=ETA_INTERVAL),
             SimilarityCurlMode(
                 base, source['poloidal_interval'],
                 eta_interval=ETA_INTERVAL)]
    chosen = source['selected']
    field = CurlPatchedLift(base, modes,
                            (chosen['swirl'], chosen['poloidal']))
    tau = .5*2**(-5.4)
    points = base.compact.joined.inner.from_similarity(
        np.array([1.025]), np.array([.22]), tau)
    spatial = .001*np.sqrt(base.nu*tau)
    temporal = .00025*tau
    rows = []
    for factor in (1., .5, .25, .125):
        R, divergence = independent_fd(field, points, tau,
                                        factor*spatial, temporal)
        rows.append(dict(h_factor=factor,
                         divergence=float(divergence[0]),
                         residual_norm=float(np.linalg.norm(R[0]))))
    report = dict(source='delayed_coupled_cone_response.json',
                  X=1.025, eta=.22, tau=tau, rows=rows,
                  scope='Fourth-order spatial finite-difference '
                        'convergence at one holdout node only.',
                  accepted=False)
    (ROOT/'delayed_coupled_divergence_refine.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(rows), flush=True)


if __name__ == '__main__':
    run()
