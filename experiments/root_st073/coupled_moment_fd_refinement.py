"""Resolve FD divergence and momentum at the lifted return hotspot."""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    tau = .5*2**(-5.5)
    field = CoupledMomentPhysicalLift()
    point = field.compact.joined.inner.from_similarity(
        np.array([2.99]), np.array([.2]), tau)
    rows = []
    for factor in (.001, .0005, .00025, .000125):
        residual, divergence = independent_fd(
            field, point, tau, factor*np.sqrt(field.nu*tau),
            .00025*tau)
        row = dict(spatial_factor=factor, residual=residual[0].tolist(),
                   residual_norm=float(np.linalg.norm(residual[0])),
                   divergence=float(divergence[0]))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, X=2.99, eta=.2, rows=rows,
                  scope='Fourth-order FD refinement at the largest sampled lifted return hotspot. Divergence convergence is a local numerical check of the streamfunction lift. The large complete momentum persists under refinement; no full-domain or volume-L2 acceptance.',
                  accepted=False)
    (ROOT/'coupled_moment_fd_refinement.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
