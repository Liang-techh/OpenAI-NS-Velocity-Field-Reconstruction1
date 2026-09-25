"""Disjoint Cartesian momentum check of the selected eight-mode mean field.

These nodes and this time are absent from both the cone fit and its
nearby-time objective. This is a sampled diagnostic, not a PDE acceptance.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_multimode_cone_fit import BASE_NAME, make_modes
from delayed_similarity_curl_screen import CurlPatchedLift
from radial_continuation import ROOT


def run(tradeoff=False):
    source_name = ('delayed_multimode_tradeoff_fit.json' if tradeoff
                   else 'delayed_multimode_cone_fit.json')
    source = json.loads((ROOT/source_name).read_text())
    coefficients = (source['searches'][0]['coefficients'] if tradeoff
                    else source['coefficients'])
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    modes, _ = make_modes(base)
    patched = CurlPatchedLift(base, modes, coefficients)
    tau = .5*2**(-5.35)
    points, X, eta = nodes(base, (1.011, 1.017, 1.023),
                           (.215, .265, .315), tau)
    before, _ = residual(base, points, tau)
    after, divergence = residual(patched, points, tau)
    node_before = np.linalg.norm(before, axis=1)
    node_after = np.linalg.norm(after, axis=1)
    report = dict(source=source_name,
                  selected_weight=.05 if tradeoff else None,
                  tau=tau, X=X.tolist(), eta=eta.tolist(),
                  before=stats(before), after=stats(after),
                  node_before=node_before.tolist(),
                  node_after=node_after.tolist(),
                  max_abs_fd_divergence=float(np.max(np.abs(divergence))),
                  scope='Nine previously unused space/time points, full '
                        'Cartesian momentum and fourth-order finite-'
                        'difference divergence. No continuous cone, '
                        'exact moments, volume L2, wave, or PDE gate.',
                  accepted=False)
    output = ('delayed_multimode_tradeoff_disjoint_holdout.json' if tradeoff
              else 'delayed_multimode_disjoint_holdout.json')
    (ROOT/output).write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(before=report['before'], after=report['after'],
                          worst_node=dict(X=float(X[np.argmax(node_after)]),
                                          eta=float(eta[np.argmax(node_after)]),
                                          after=float(np.max(node_after))),
                          max_abs_fd_divergence=(
                              report['max_abs_fd_divergence']))), flush=True)


if __name__ == '__main__':
    import sys
    run(tradeoff='--tradeoff' in sys.argv[1:])
