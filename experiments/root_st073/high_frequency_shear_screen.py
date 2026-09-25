"""Numerical scale separation for compact radial shear modulation.

The profile perturbation is B(y) sin(2*pi*N*(y-y0))/(2*pi*N): it has
order-one slope at X=1 and decreasing values/moments as N grows. This is
only a sinusoidal probe; the paper's admissible all-phase shear loop and
nonaxisymmetric wave correction are not implemented here.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from bridge_poloidal_mode import BridgePoloidalMode
from bridge_swirl_moment_balance import cp_tail
from extended_compact_join import load_extended_heated_candidate
from extended_relaxed_cone_screen import cone_point, profile
from joined_field import independent_fd
from radial_continuation import ROOT


SYMMETRIC = -.27311627313724157
AMPLITUDE = 2.


def make_field(cycles, amplitude=AMPLITUDE):
    field = load_extended_heated_candidate(
        swirl_bubble_amplitude=SYMMETRIC,
        outer_swirl_bubble_amplitude=1.,
        shear_swirl_amplitude=amplitude,
        shear_swirl_cycles=cycles)
    field = BridgePoloidalMode(field, -2.75, shape='minimum_curvature')
    return BridgePoloidalMode(field, 4.5, shape='moment')


def moment_probe(base, reference, field, cycles, eta, tau):
    radial = base.compact.joined
    xi, xo = radial.join_X, radial.join_X*radial.outer_ratio**2
    nodes, weights = leggauss(max(128, 4*cycles))
    xs = (xi+xo)/2+(xo-xi)*nodes/2
    weights = (xo-xi)*weights/(4*xs)
    E0 = profile(base, xs, eta, tau)[1]
    Eref = profile(reference, xs, eta, tau)[1]
    E = profile(field, xs, eta, tau)[1]
    return dict(cp_delta_from_unmodified=float(weights@(E**2-E0**2)),
                cp_delta_from_reference=float(weights@(E**2-Eref**2)),
                max_modulation_E=float(np.max(np.abs(E-Eref))),
                min_E=float(np.min(E)))


def run():
    tau = .5*2**(-5.5)
    base = load_extended_heated_candidate()
    reference = make_field(16, 0.)
    rows = []
    for cycles in (8, 16, 32):
        field = make_field(cycles)
        probes = {str(eta): moment_probe(base, reference, field,
                                        cycles, eta, tau)
                  for eta in (.2, .3, .65)}
        rows.append(dict(cycles=cycles, amplitude=AMPLITUDE,
                         probes=probes,
                         downstream_cp_capacity={
                             str(eta): cp_tail(base, 1.5, eta, tau)
                             for eta in (.2, .3, .65)}))
        print(json.dumps(dict(cycles=cycles, probes=probes)), flush=True)
    candidate = make_field(16)
    cones = []
    for label, field in [('reference', reference),
                         ('modulated', candidate)]:
        for eta in (.2, .3):
            point = cone_point(field, 1., eta, tau, order=64)
            cones.append(dict(field=label, eta=eta,
                              Pc=point.get('Pc'),
                              margin=(point['upper']-point['vs']
                                      if point.get('upper') is not None else None),
                              relaxed_pass=point.get('relaxed_pass')))
        print(json.dumps(dict(field=label, cones=cones[-2:])), flush=True)
    x, eta = np.meshgrid([.5859375, .75, 1., 1.25], [.2, .3],
                         indexing='ij')
    points = base.compact.joined.inner.from_similarity(
        x.ravel(), eta.ravel(), tau)
    R, div = independent_fd(
        candidate, points, tau, .001*np.sqrt(candidate.nu*tau),
        .00025*tau)
    report = dict(tau=tau, rows=rows, cones=cones,
                  modulated_eight_point_full_max=float(
                      np.max(np.linalg.norm(R, axis=1))),
                  modulated_eight_point_angular_max=float(
                      np.max(np.abs(R[:, 1]))),
                  modulated_divergence_max=float(np.max(np.abs(div))),
                  scope='Compact sinusoidal shear experiment only. No all-phase admissible cone, exact five-moment restoration, wave stress, or global PDE acceptance.',
                  accepted=False)
    (ROOT/'high_frequency_shear_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(momentum_max=report['modulated_eight_point_full_max'],
                          angular_max=report['modulated_eight_point_angular_max'],
                          divergence_max=report['modulated_divergence_max'])))


if __name__ == '__main__':
    run()
