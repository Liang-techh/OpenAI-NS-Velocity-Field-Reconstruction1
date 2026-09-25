"""Local shear-sign screen for extending the physical cone past X=1.02.

Pressure cannot change lambda squared. This tests shifted compact
solenoidal swirl/curl modes on three eta slices and four radii, before
any expensive full radial stress response is built.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode
from delayed_multimode_cone_fit import BASE_NAME
from delayed_outer_radial_extension import make_extended_modes
from delayed_similarity_curl_screen import CurlPatchedLift, SimilarityCurlMode
from joint_collar_fit import kinematics
from radial_continuation import ROOT


XS = (1.016, 1.02, 1.03, 1.04)
ETAS = (.2, .25, .3)
RADIAL_INTERVALS = ((1.02, 1.06), (1.02, 1.08), (1.02, 1.12),
                    (1.025, 1.09), (1.015, 1.08))
AXIAL_INTERVAL = (.14, .36)


def lambda_squared(u, grad, radius):
    F = u[:, 1]/radius
    shear = np.column_stack((grad[:, 1, 0]-F, grad[:, 2, 0]))
    norm = np.linalg.norm(shear, axis=1)
    Ntheta = shear[:, 0]/norm
    return -2*F*Ntheta*(2*F*Ntheta+norm)


def run():
    source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = CurlPatchedLift(base, make_extended_modes(base),
                           source['coefficients'])
    tau = .5*2**(-5.5)
    X, eta = np.meshgrid(XS, ETAS, indexing='ij')
    points = base.compact.joined.inner.from_similarity(
        X.ravel(), eta.ravel(), tau)
    hs, ht = .0005*np.sqrt(base.nu*tau), .0001*tau
    u0, g0, _ = kinematics(mean, points, tau, hs, ht)
    radius = points[:, 0]
    original = lambda_squared(u0, g0, radius).reshape(len(XS), len(ETAS))
    rows = []
    for kind in ('swirl', 'poloidal'):
        bound = .1 if kind == 'swirl' else .03
        for interval in RADIAL_INTERVALS:
            mode = (SimilaritySwirlMode(base, interval,
                                        eta_interval=AXIAL_INTERVAL)
                    if kind == 'swirl' else
                    SimilarityCurlMode(base, interval,
                                       eta_interval=AXIAL_INTERVAL))
            u1, g1, _ = kinematics(mode, points, tau, hs, ht)
            trials = []
            for a in np.linspace(-bound, bound, 401):
                lam = lambda_squared(u0+a*u1, g0+a*g1,
                                      radius).reshape(len(XS), len(ETAS))
                trials.append(dict(amplitude=float(a),
                                   positive_count=int(np.sum(lam > 0)),
                                   all_X103=bool(np.all(lam[2] > 0)),
                                   all_inner_and_X103=bool(
                                       np.all(lam[:3] > 0)),
                                   min_lambda_squared=float(np.min(lam)),
                                   lambda_X103=lam[2].tolist()))
            valid = [trial for trial in trials
                     if trial['all_inner_and_X103']]
            best = (min(valid, key=lambda trial:
                        abs(trial['amplitude'])) if valid else
                    max(trials, key=lambda trial:
                        (trial['positive_count'],
                         trial['min_lambda_squared'])))
            rows.append(dict(kind=kind, radial_interval=interval,
                             amplitude_bound=bound,
                             valid_amplitude_count=len(valid),
                             selected=best))
            print(json.dumps(dict(kind=kind, interval=interval,
                                  valid_count=len(valid),
                                  selected=best)), flush=True)
    report = dict(source='delayed_outer_admission_search.json',
                  tau=tau, X=XS, eta=ETAS,
                  original_lambda_squared=original.tolist(),
                  axial_interval=AXIAL_INTERVAL,
                  rows=rows,
                  scope='Necessary local shear-sign screen for shifted '
                        'solenoidal velocity modes. No pressure/stress '
                        'primitive, outgoing moment audit, independent '
                        'momentum, continuous cone, wave, or PDE gate.',
                  accepted=False)
    (ROOT/'delayed_outer_shear_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
