"""Independent space/time momentum check of the coupled cone screen."""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_similarity_curl_screen import CurlPatchedLift, SimilarityCurlMode
from delayed_swirl_cone_response import BASE_NAME, ETA_INTERVAL
from radial_continuation import ROOT


def run():
    source = json.loads((ROOT/'delayed_coupled_cone_response.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    modes = [SimilaritySwirlMode(
                 base, source['swirl_interval'], eta_interval=ETA_INTERVAL),
             SimilarityCurlMode(
                 base, source['poloidal_interval'],
                 eta_interval=ETA_INTERVAL)]
    selected = source['selected']
    patched = CurlPatchedLift(base, modes,
                              (selected['swirl'], selected['poloidal']))
    tau = .5*2**(-5.4)
    points, X, eta = nodes(base, (1.01, 1.015, 1.025),
                           (.22, .28, .32), tau)
    before, _ = residual(base, points, tau)
    after, divergence = residual(patched, points, tau)
    report = dict(source='delayed_coupled_cone_response.json',
                  tau=tau, X=X.tolist(), eta=eta.tolist(),
                  before=stats(before), after=stats(after),
                  node_before=np.linalg.norm(before, axis=1).tolist(),
                  node_after=np.linalg.norm(after, axis=1).tolist(),
                  max_abs_fd_divergence=float(np.max(np.abs(divergence))),
                  scope='Disjoint space/time full Cartesian momentum '
                        'and finite-difference divergence of the selected '
                        'two-mode field. No cone interval or PDE gate.',
                  accepted=False)
    (ROOT/'delayed_coupled_cone_holdout.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(before=report['before'], after=report['after'],
                          max_abs_fd_divergence=(
                              report['max_abs_fd_divergence']))), flush=True)


if __name__ == '__main__':
    run()
