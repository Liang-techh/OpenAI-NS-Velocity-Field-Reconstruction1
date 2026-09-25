"""Compare local relaxed-cone gain with the actual momentum cost."""
import json
import numpy as np

from bridge_poloidal_mode import BridgePoloidalMode
from extended_compact_join import load_extended_heated_candidate
from extended_relaxed_cone_screen import cone_point
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    base = load_extended_heated_candidate()
    coupled = BridgePoloidalMode(
        load_extended_heated_candidate(outer_swirl_bubble_amplitude=1.1),
        amplitude=-2.5, shape='curvature')
    tau = .5*2**(-5.5)
    X, eta = np.meshgrid([.5859375, .75, 1., 1.25], [.2, .3],
                         indexing='ij')
    X, eta = X.ravel(), eta.ravel()
    points = base.compact.joined.inner.from_similarity(X, eta, tau)
    rows = []
    for name, field in (('baseline', base), ('coupled', coupled)):
        R, div = independent_fd(
            field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
        rows.append(dict(name=name,
                         momentum_norms=np.linalg.norm(R, axis=1).tolist(),
                         momentum_max=float(np.max(np.linalg.norm(R, axis=1))),
                         finite_difference_divergence_max=float(np.max(np.abs(div)))))
    cone_rows = []
    for k in (3., 5.5, 6.):
        t = .5*2**(-k)
        for e in (.2, .3):
            xs = (.5859375, .75, 1., 1.25, 1.5) if k == 5.5 else (1.,)
            for x in xs:
                d = cone_point(coupled, x, e, t, order=12)
                cone_rows.append(dict(k=k, eta=e, X=x, a=d.get('a'),
                                      vs=d.get('vs'), Pc=d.get('Pc'),
                                      upper=d.get('upper'),
                                      margin=(d['upper']-d['vs']
                                              if d.get('upper') is not None else None),
                                      relaxed_pass=d.get('relaxed_pass')))
    report = dict(tau=tau, X=X.tolist(), eta=eta.tolist(), rows=rows,
                  swirl_outer_bubble_amplitude=1.1,
                  poloidal_curvature_amplitude=-2.5,
                  cone_rows=cone_rows,
                  scope='Sampled local cone passes are outweighed by much larger complete momentum residual and failures at neighboring radii. Diagnostic candidate rejected; no global acceptance.',
                  accepted=False)
    (ROOT/'coupled_bridge_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(momentum_max_by_variant={r['name']: r['momentum_max'] for r in rows},
                          cone_pass_count=sum(bool(r['relaxed_pass']) for r in cone_rows),
                          cone_point_count=len(cone_rows))), flush=True)


if __name__ == '__main__':
    run()
