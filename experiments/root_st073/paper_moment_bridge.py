"""Solve two Appendix-A-style radial moment constraints on joined slices.

The result is a slice-level normalized profile correction, not yet a global
velocity: eta/tau dependence of coefficients must be differentiated in a
streamfunction before claiming incompressibility or momentum improvement.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from joined_field import JoinedField, ROOT
from heat_exterior import tail_moments


XI = 3/64
XB = 4*XI
PATCHES = ((XI, XB), (1., 4.), (4., 16.))


def bump(X, start, end):
    y = (X-start)/(end-start)
    inside = (y > 0) & (y < 1)
    b = np.zeros_like(y)
    db = np.zeros_like(y)
    yy = y[inside]
    b[inside] = 1024*yy**5*(1-yy)**5
    db[inside] = 5120*yy**4*(1-yy)**4*(1-2*yy)/(end-start)
    return b, db


def slice_data(field, eta, tau, order, patch_start, patch_end):
    g, w = leggauss(order)
    xs, ws = [], []
    boundaries = sorted(set((0., XI, XB, patch_start, patch_end)))
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        xs.append((lo+hi)/2+(hi-lo)/2*g)
        ws.append((hi-lo)/2*w)
    X = np.concatenate(xs)
    weights = np.concatenate(ws)
    points = field.inner.from_similarity(X, eta, tau)
    velocity, _ = field.fields(points, tau)
    q = tau/(1-eta**2)
    A = .5+field.inner.h
    U = q**A/np.sqrt(field.nu)*velocity[:, 2]
    E = q**A/np.sqrt(field.nu)*velocity[:, 1]
    b, db = bump(X, patch_start, patch_end)
    _, tail_swirl_square, tail_angular_difference = tail_moments(
        patch_end, eta, c=field.c, h=field.inner.h, n=64)
    tail_swirl_square = float(tail_swirl_square)
    tail_angular_difference = float(tail_angular_difference)
    power_integral = (np.sqrt(2)*field.c*patch_end**(1-field.inner.h)
                      /(1-field.inner.h))
    angular_moment = (float(weights@(np.sqrt(2*X)*E))
                      - power_integral + tail_angular_difference)
    swirl_response = float(weights@(np.sqrt(2*X)*b))
    aE = -angular_moment/swirl_response
    adjusted_E = E+aE*b
    kinetic_baseline = float(weights@(U**2-E**2/2)
                             -tail_swirl_square)
    swirl_linear_response = float(-(weights@(E*b)))
    swirl_quadratic_response = float(-.5*(weights@(b**2)))
    kinetic_moment = float(weights@(U**2-adjusted_E**2/2)
                           -tail_swirl_square)
    linear_U = float(2*weights@(U*db))
    quadratic_U = float(weights@(db**2))
    discriminant = linear_U**2-4*quadratic_U*kinetic_moment
    roots = []
    if discriminant >= 0:
        roots = [(-linear_U+sign*np.sqrt(discriminant))/(2*quadratic_U)
                 for sign in (-1, 1)]
    aU = max(roots) if roots else None
    corrected_kinetic = (kinetic_moment+linear_U*aU+quadratic_U*aU**2
                         if aU is not None else None)
    return {'eta': eta, 'tau': tau, 'order_per_piece': order,
            'patch_start': patch_start, 'patch_end': patch_end,
            'baseline_angular_moment': angular_moment,
            'swirl_bump_response': swirl_response,
            'swirl_amplitude': aE,
            'corrected_angular_moment': angular_moment+aE*swirl_response,
            'kinetic_moment_baseline': kinetic_baseline,
            'kinetic_swirl_linear_response': swirl_linear_response,
            'kinetic_swirl_quadratic_response': swirl_quadratic_response,
            'kinetic_moment_after_swirl': kinetic_moment,
            'meridional_linear_response': linear_U,
            'meridional_quadratic_response': quadratic_U,
            'meridional_discriminant': discriminant,
            'meridional_roots': roots,
            'meridional_amplitude': aU,
            'corrected_kinetic_moment': corrected_kinetic,
            'tail_swirl_square': tail_swirl_square,
            'tail_angular_difference': tail_angular_difference,
            'max_normalized_U_before': float(np.max(np.abs(U))),
            'max_normalized_E_before': float(np.max(np.abs(E)))}


def run():
    field = JoinedField()
    rows = []
    for patch_start, patch_end in PATCHES:
        for tau in (.5/64, .032):
            for eta in (-.3, 0., .3):
                row = slice_data(field, eta, tau, 24,
                                 patch_start, patch_end)
                rows.append(row)
                print(json.dumps({k: row[k] for k in
                                  ('patch_start', 'patch_end', 'eta', 'tau',
                                   'swirl_amplitude', 'meridional_amplitude',
                                   'meridional_discriminant')}), flush=True)
    report = {'radial_X_inner': XI, 'radial_X_heat': XB,
              'patches': PATCHES,
              'heat_amplitude': field.c, 'rows': rows,
              'bump': 'B=1024 y^5(1-y)^5, y=(X-patch_start)/(patch_end-patch_start); delta E=aE B; delta U=aU dB/dX. Both preserve endpoint values and several radial jets. The U bump has zero radial integral. Choose the positive quadratic root for a continuous candidate branch across eta.',
              'scope': 'Appendix-A-style normalized moments of the existing finite-slab JoinedField with a heat exterior tail. Slice-level algebra only; fitted amplitudes are not yet a smooth eta/tau field or physical solenoidal correction. No stress-cone, pointwise momentum, finite-energy, or force acceptance.',
              'accepted': False}
    (ROOT/'paper_moment_bridge.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
