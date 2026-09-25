"""Test minimum-third-derivative poloidal shape against momentum and cone."""
import json
import numpy as np
from numpy.polynomial import Polynomial
from numpy.polynomial.legendre import leggauss

from bridge_poloidal_mode import BridgePoloidalMode, minimum_curvature_coefficients
from extended_compact_join import load_extended_heated_candidate
from extended_relaxed_cone_screen import cone_point
from joined_field import independent_fd
from radial_continuation import ROOT


def curvature_diagnostics(join_X, ratio):
    y0, coefficients = minimum_curvature_coefficients(join_X, ratio)
    y = Polynomial([0., 1.])
    common = y**4*(1-y)**4*(y-y0)**2
    minimum = sum((coefficients[k]*common*(2*y-1)**k
                   for k in range(len(coefficients))), Polynomial([0.]))
    old = 8*1024*y**5*(1-y)**5*(y-y0)**2
    g, w = leggauss(64)
    nodes, weights = (g+1)/2, w/2
    return dict(y0=float(y0), coefficients=coefficients.tolist(),
                old_third_derivative_l2=float(np.sqrt(
                    weights@(old.deriv(3)(nodes)**2))),
                minimum_third_derivative_l2=float(np.sqrt(
                    weights@(minimum.deriv(3)(nodes)**2))),
                old_second_at_reference=float(old.deriv(2)(y0)),
                minimum_second_at_reference=float(minimum.deriv(2)(y0)))


def run():
    base = load_extended_heated_candidate()
    swirl = load_extended_heated_candidate(outer_swirl_bubble_amplitude=1.1)
    field = BridgePoloidalMode(swirl, -2.5, shape='minimum_curvature')
    tau = .5*2**(-5.5)
    X, eta = np.meshgrid([.5859375, .75, 1., 1.25], [.2, .3],
                         indexing='ij')
    X, eta = X.ravel(), eta.ravel()
    points = base.compact.joined.inner.from_similarity(X, eta, tau)
    rows = []
    for name, candidate in (('baseline', base), ('minimum_curvature', field)):
        R, div = independent_fd(
            candidate, points, tau, .001*np.sqrt(candidate.nu*tau),
            .00025*tau)
        norms = np.linalg.norm(R, axis=1)
        rows.append(dict(name=name, complete_residual_norms=norms.tolist(),
                         sampled_max=float(np.max(norms)),
                         divergence_max=float(np.max(np.abs(div)))))
    cone_rows = []
    for weight in (0., .3, .6, .9, 1.):
        candidate = BridgePoloidalMode(
            swirl, -2.5, shape='mixed', curvature_weight=weight)
        points_cone = [cone_point(candidate, 1., e, tau, order=12)
                       for e in (.2, .3)]
        cone_rows.append(dict(curvature_weight=weight,
                              points=[dict(eta=p['eta'], Pc=p.get('Pc'),
                                           vs=p.get('vs'), upper=p.get('upper'),
                                           margin=(p['upper']-p['vs']
                                                   if p.get('upper') is not None else None),
                                           relaxed_pass=p.get('relaxed_pass'))
                                      for p in points_cone]))
    radial = base.compact.joined
    report = dict(tau=tau, X=X.tolist(), eta=eta.tolist(),
                  shape=curvature_diagnostics(radial.join_X, radial.outer_ratio),
                  rows=rows, cone_rows=cone_rows,
                  scope='The variational shape minimizes integral |phi_yyy|^2 in the chosen finite-dimensional polynomial family at fixed local curvature. It reduces sampled momentum cost but fails the relaxed cone; no global acceptance.',
                  accepted=False)
    (ROOT/'minimum_curvature_bridge_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(shape_third_derivative_l2=[
        report['shape']['old_third_derivative_l2'],
        report['shape']['minimum_third_derivative_l2']],
        sampled_max_by_variant={row['name']: row['sampled_max'] for row in rows},
        minimum_curvature_cone=cone_rows[0]['points'])), flush=True)


if __name__ == '__main__':
    run()
