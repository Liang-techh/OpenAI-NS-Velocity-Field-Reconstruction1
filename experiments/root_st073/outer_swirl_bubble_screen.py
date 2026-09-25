"""Test an outer-biased C2 swirl bubble against cone and full momentum."""
import json
import numpy as np

from extended_compact_join import load_extended_heated_candidate
from extended_relaxed_cone_screen import cone_point
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    tau = .5*2**(-5.5)
    X = np.array([1., 1.25, 1., 1.25])
    eta = np.array([.2, .2, .3, .3])
    rows = []
    for amplitude in (0., .76):
        field = load_extended_heated_candidate(
            outer_swirl_bubble_amplitude=amplitude)
        points = field.compact.joined.inner.from_similarity(X, eta, tau)
        R, divergence = independent_fd(
            field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
        cone = [cone_point(field, 1., e, tau, order=12) for e in (.2, .3)]
        row = dict(amplitude=amplitude,
                   complete_residual_norms=np.linalg.norm(R, axis=1).tolist(),
                   divergence_max=float(np.max(np.abs(divergence))),
                   cone_X1=[dict(eta=p['eta'], a=p.get('a'),
                                 vs=p.get('vs'), Pc=p.get('Pc'),
                                 upper=p.get('upper'),
                                 relaxed_pass=p.get('relaxed_pass'))
                            for p in cone])
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, X=X.tolist(), eta=eta.tolist(), rows=rows,
                  scope='Four momentum points and two cone points per variant; outer-biased swirl-only mode preserves radial C2 but does not solve bridge momentum or relaxed cone.',
                  accepted=False)
    (ROOT/'outer_swirl_bubble_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
