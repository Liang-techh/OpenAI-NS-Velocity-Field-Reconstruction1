"""Audit all five paper radial moments changed by a remote field patch."""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from joined_field import JoinedField, ROOT
from moment_matched_joined_field import MomentMatchedJoinedField


def profiles(field, X, eta, tau):
    points = field.inner.from_similarity(X, eta, tau)
    velocity, _ = field.fields(points, tau)
    q = tau/(1-eta**2)
    factor = q**(.5+field.inner.h)/np.sqrt(field.nu)
    return factor*velocity[:, 2], factor*velocity[:, 1]


def moments(U, E, X, weights):
    H = np.sqrt(2*X)*E
    return np.array([weights@U, weights@H, weights@(U*H),
                     weights@(U**2-E**2/2), weights@(E**2/(2*X))])


def run():
    base = JoinedField()
    variants = (('two_moment', MomentMatchedJoinedField(base=base)),
                ('two_shape', MomentMatchedJoinedField(
                    base=base, meridional_null_amplitude=-5.,
                    meridional_even_amplitude=12.)))
    start, end = variants[0][1].patch_start, variants[0][1].patch_end
    g, w = leggauss(80)
    X = (start+end)/2+(end-start)/2*g
    weights = (end-start)/2*w
    tau = .0084
    rows = []
    for eta in (0., .2, .35):
        U0, E0 = profiles(base, X, eta, tau)
        m0 = moments(U0, E0, X, weights)
        q = tau/(1-eta**2)
        pressure_scale = base.nu*q**(-1-2*base.inner.h)
        for name, field in variants:
            U, E = profiles(field, X, eta, tau)
            delta = moments(U, E, X, weights)-m0
            row = {'variant': name, 'eta': eta, 'tau': tau,
                   'patch_end_moment_difference': dict(zip(
                       ('M', 'I', 'J', 'S', 'Cp'), delta.tolist())),
                   'physical_outer_pressure_increment_required':
                       float(pressure_scale*delta[4])}
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = {'patch_X': [start, end], 'rows': rows,
              'scope': 'Appendix-A/Section-4 five cumulative moment differences across the remote patch, relative to the unmodified joined field. This does not impose the correct exterior target or certify a normalized cone. Nonzero Cp difference means the current unchanged physical pressure cannot obey the radial centrifugal pressure relation after the patch.',
              'accepted': False}
    (ROOT/'remote_five_moment_audit.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
