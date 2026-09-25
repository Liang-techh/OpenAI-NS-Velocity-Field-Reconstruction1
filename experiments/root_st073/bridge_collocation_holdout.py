"""Directly verify the best sampled three-mode bridge compromise."""
import json
import numpy as np

from bridge_poloidal_mode import BridgePoloidalMode
from extended_compact_join import load_extended_heated_candidate
from extended_relaxed_cone_screen import cone_point
from joined_field import independent_fd
from radial_continuation import ROOT


def load_candidate():
    swirl = load_extended_heated_candidate(
        outer_swirl_bubble_amplitude=1.)
    minimum = BridgePoloidalMode(
        swirl, -2.75, shape='minimum_curvature')
    return BridgePoloidalMode(minimum, 4.5, shape='moment')


def run():
    field = load_candidate()
    tau = .5*2**(-5.5)
    X, eta = np.meshgrid([.5859375, .75, 1., 1.25], [.2, .3],
                         indexing='ij')
    X, eta = X.ravel(), eta.ravel()
    points = field.compact.joined.inner.from_similarity(X, eta, tau)
    R, divergence = independent_fd(
        field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
    cones = []
    for k in (3., 5.5, 6.):
        t = .5*2**(-k)
        for e in (.2, .3):
            xs = (.75, .9, 1., 1.1, 1.25) if k == 5.5 else (1.,)
            for x in xs:
                d = cone_point(field, x, e, t, order=16)
                cones.append(dict(k=k, X=x, eta=e,
                                  Pc=d.get('Pc'), vs=d.get('vs'),
                                  upper=d.get('upper'),
                                  margin=(d['upper']-d['vs']
                                          if d.get('upper') is not None else None),
                                  relaxed_pass=d.get('relaxed_pass')))
    report = dict(coefficients=dict(outer_swirl=1.,
                                    minimum_curvature=-2.75,
                                    moment=4.5),
                  tau=tau, X=X.tolist(), eta=eta.tolist(),
                  direct_momentum_norms=np.linalg.norm(R, axis=1).tolist(),
                  direct_momentum_max=float(np.max(np.linalg.norm(R, axis=1))),
                  divergence_max=float(np.max(np.abs(divergence))),
                  cone_rows=cones,
                  scope='Direct eight-point physical momentum and limited space/time cone holdout. The cone fails outside a narrow radius band and momentum remains far above 1e-3; no global acceptance.',
                  accepted=False)
    (ROOT/'bridge_collocation_holdout.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(direct_momentum_max=report['direct_momentum_max'],
                          cone_pass_count=sum(bool(d['relaxed_pass'])
                                              for d in cones),
                          cone_point_count=len(cones))))


if __name__ == '__main__':
    run()
