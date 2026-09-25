"""Check the radial bridge of the selected extended-core candidate."""
import json
import numpy as np

from extended_compact_join import load_extended_heated_candidate
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    field = load_extended_heated_candidate()
    tau = .5*2**(-5.5)
    X = np.full(3, .5859375)
    eta = np.array([.3, .65, .8])
    points = field.compact.joined.inner.from_similarity(X, eta, tau)
    R, divergence = independent_fd(
        field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
    rows = [dict(eta=float(eta[i]), X=float(X[i]),
                 residual_norm=float(np.linalg.norm(R[i])),
                 divergence=float(divergence[i])) for i in range(len(eta))]
    report = dict(tau=tau, rows=rows,
                  scope='Three sampled radial-bridge points; no full-volume acceptance.',
                  accepted=False)
    (ROOT/'extended_radial_bridge_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run()
