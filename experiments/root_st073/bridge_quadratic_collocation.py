"""Quadratic complete-momentum surrogate for three bridge correction modes.

Velocity and all linear PDE terms are affine in the mode amplitudes;
convection is exactly quadratic. The surrogate accelerates candidate
search but only on its sampled physical points.
"""
import json
import numpy as np

from bridge_poloidal_mode import BridgePoloidalMode
from extended_compact_join import load_extended_heated_candidate
from extended_relaxed_cone_screen import cone_point
from joint_collar_fit import kinematics
from radial_continuation import ROOT


def make_field(swirl=0., minimum=0., moment=0.):
    base = load_extended_heated_candidate(
        outer_swirl_bubble_amplitude=swirl)
    if minimum:
        base = BridgePoloidalMode(base, minimum, shape='minimum_curvature')
    if moment:
        base = BridgePoloidalMode(base, moment, shape='moment')
    return base


def response_model(points, tau):
    fields = (make_field(), make_field(swirl=1.),
              make_field(minimum=1.), make_field(moment=1.))
    data = [kinematics(field, points, tau,
                       .001*np.sqrt(field.nu*tau), .00025*tau)
            for field in fields]
    u0, J0, part0 = data[0]
    du = np.stack([row[0]-u0 for row in data[1:]], axis=0)
    dJ = np.stack([row[1]-J0 for row in data[1:]], axis=0)
    dpart = np.stack([row[2]-part0 for row in data[1:]], axis=0)

    def residual(coefficients):
        c = np.asarray(coefficients, float)
        u = u0+np.einsum('k,knj->nj', c, du)
        J = J0+np.einsum('k,knij->nij', c, dJ)
        part = part0+np.einsum('k,knj->nj', c, dpart)
        return part+np.einsum('nij,nj->ni', J, u)

    return residual


def run():
    tau = .5*2**(-5.5)
    X, eta = np.meshgrid([.5859375, .75, 1., 1.25], [.2, .3],
                         indexing='ij')
    X, eta = X.ravel(), eta.ravel()
    reference = make_field()
    points = reference.compact.joined.inner.from_similarity(X, eta, tau)
    residual = response_model(points, tau)
    known = np.array([2., -2.5, 4.5])
    known_max = float(np.max(np.linalg.norm(residual(known), axis=1)))
    candidates = []
    for swirl in np.arange(1., 3.01, .25):
        for minimum in np.arange(-3., -1.99, .25):
            for moment in np.arange(3., 7.01, .5):
                c = [swirl, minimum, moment]
                cost = float(np.max(np.linalg.norm(residual(c), axis=1)))
                candidates.append((cost, c))
    candidates.sort(key=lambda item: item[0])
    cone_checked = []
    for cost, c in candidates[:30]:
        field = make_field(*c)
        cone = [cone_point(field, 1., e, tau, order=8)
                for e in (.2, .3)]
        row = dict(coefficients=c, sampled_momentum_max=cost,
                   cone=[dict(eta=p['eta'], Pc=p.get('Pc'),
                              margin=(p['upper']-p['vs']
                                      if p.get('upper') is not None else None),
                              relaxed_pass=p.get('relaxed_pass'))
                         for p in cone],
                   both_cone_pass=all(p.get('relaxed_pass') for p in cone))
        cone_checked.append(row)
        if row['both_cone_pass']:
            print(json.dumps(row), flush=True)
    report = dict(tau=tau, X=X.tolist(), eta=eta.tolist(),
                  known_coefficients=known.tolist(),
                  known_sampled_momentum_max=known_max,
                  grid_size=len(candidates),
                  unconstrained_best=dict(coefficients=candidates[0][1],
                                          sampled_momentum_max=candidates[0][0]),
                  cone_checked=cone_checked,
                  scope='Quadratic physical-momentum surrogate at eight points; only the 30 lowest sampled-residual grid candidates receive two-point snapshot cone screens. This does not establish continuum/global acceptance.',
                  accepted=False)
    (ROOT/'bridge_quadratic_collocation.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(known_max=known_max,
                          unconstrained_best=report['unconstrained_best'],
                          cone_pass_count=sum(row['both_cone_pass']
                                              for row in cone_checked))))


if __name__ == '__main__':
    run()
