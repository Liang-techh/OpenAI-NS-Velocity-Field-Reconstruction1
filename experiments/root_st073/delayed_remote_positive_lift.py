"""Lift independently optimized positive-E slices to a smooth swirl field.

Five compact axial bumps interpolate the slice coefficients exactly.
The resulting axisymmetric swirl perturbation is solenoidal, but
positivity, moment closure, and momentum still require off-slice audits.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import (HOLDOUT_ETAS,
                                          NEAR5_AXIAL_INTERVALS,
                                          current_mean, quadrature)
from delayed_similarity_curl_screen import CurlPatchedLift
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT


def axial_values(etas, intervals):
    return np.array([[bump(np.array([eta]), *interval)[0][0]
                      for interval in intervals] for eta in etas])


def run():
    source = json.loads((ROOT/'delayed_remote_positive_e.json').read_text())
    intervals = source.get('radial_intervals', source.get('intervals'))
    if intervals is None:
        raise KeyError('positive-E source needs radial_intervals')
    intervals = [tuple(map(float, pair)) for pair in intervals]
    rows = sorted(source['rows'], key=lambda row: row['eta'])
    etas = np.array([row['eta'] for row in rows])
    source_eta_interval = tuple(source['mode_eta_interval'])
    # The independent optimizer used a common axial bump in each mode.
    # Convert its physical-mode coefficients to pure radial amplitudes
    # before interpolating with the five new axial bumps.
    source_axial = axial_values(etas, [source_eta_interval])[:, 0]
    slice_coefficients = (np.array([row['coefficients'] for row in rows])
                          *source_axial[:, None])
    if len(etas) != len(NEAR5_AXIAL_INTERVALS):
        raise ValueError('five training slices required for axial lift')
    if slice_coefficients.shape != (len(etas), len(intervals)):
        raise ValueError('slice coefficient dimensions disagree')
    matrix = axial_values(etas, NEAR5_AXIAL_INTERVALS)
    axial_coefficients = np.linalg.solve(matrix, slice_coefficients)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target = make_field(16, 2.)
    modes = [SimilaritySwirlMode(base, radial, eta_interval=axial)
             for axial in NEAR5_AXIAL_INTERVALS for radial in intervals]
    patched = CurlPatchedLift(mean, modes, axial_coefficients.ravel())
    X, weights = quadrature(order=48, radial_intervals=intervals)
    radial_bumps = np.array([bump(X, *interval)[0]
                             for interval in intervals])
    tau = .5*2**(-5.5)
    test_etas = np.sort(np.r_[etas, HOLDOUT_ETAS])
    checks = []
    for eta in test_etas:
        U, E = profile(mean, X, eta, tau)
        U0, E0 = profile(target, X, eta, tau)
        amplitudes = axial_values([eta], NEAR5_AXIAL_INTERVALS)[0]@axial_coefficients
        corrected_E = E+amplitudes@radial_bumps
        moment_defect = (moment_vector(U, corrected_E, X, weights)
                         -moment_vector(U0, E0, X, weights))
        checks.append(dict(eta=float(eta), is_training=bool(np.any(
            np.isclose(eta, etas, atol=1e-12))),
                           moment_defect=moment_defect.tolist(),
                           min_relative_E=float(np.min(corrected_E/E0))))
    # Confirm that the actual Cartesian velocity implements the analytic
    # profile interpolation, rather than just matching its algebra.
    physical_U, physical_E = profile(patched, X, .25, tau)
    mean_U, mean_E = profile(mean, X, .25, tau)
    amp_mid = axial_values([.25], NEAR5_AXIAL_INTERVALS)[0]@axial_coefficients
    expected_E = mean_E+amp_mid@radial_bumps
    cone_points, _, _ = nodes(base, (1.008, 1.016, 1.02, 1.03),
                              (.2, .25, .3), tau)
    old_local, _ = mean.fields(cone_points, tau)
    new_local, _ = patched.fields(cone_points, tau)
    nearby_tau = .5*2**(-5.4)
    remote_points, Xh, etah = nodes(base, (1.1, 1.5, 2., 2.7),
                                   (.225, .275), nearby_tau)
    before, _ = residual(mean, remote_points, nearby_tau)
    after, divergence = residual(patched, remote_points, nearby_tau)
    report = dict(source='delayed_remote_positive_e.json',
                  axial_intervals=NEAR5_AXIAL_INTERVALS,
                  radial_intervals=intervals,
                  axial_interpolation_condition=float(np.linalg.cond(matrix)),
                  axial_coefficients=axial_coefficients.tolist(),
                  checks=checks,
                  max_abs_physical_profile_error=float(max(
                      np.max(np.abs(physical_E-expected_E)),
                      np.max(np.abs(physical_U-mean_U)))),
                  max_abs_local_velocity_change=float(np.max(np.abs(
                      new_local-old_local))),
                  remote_momentum=dict(tau=nearby_tau, X=Xh.tolist(),
                                       eta=etah.tolist(),
                                       before=stats(before), after=stats(after),
                                       max_abs_fd_divergence=float(
                                           np.max(np.abs(divergence)))),
                  scope='Smooth exact-solenoidal swirl interpolation of '
                        'independently fitted E slices. Reports off-slice '
                        'positive-E and moment defects plus sampled full '
                        'Cartesian momentum. No U repair, pressure '
                        'compatibility, continuous bounds, or PDE '
                        'acceptance.', accepted=False)
    (ROOT/'delayed_remote_positive_lift.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(interpolation_condition=(
                        report['axial_interpolation_condition']),
                          max_abs_physical_profile_error=(
                              report['max_abs_physical_profile_error']),
                          max_abs_local_velocity_change=(
                              report['max_abs_local_velocity_change']),
                          train_max_ICp=max(max(abs(check['moment_defect'][i])
                                                for i in (1, 4))
                                            for check in checks
                                            if check['is_training']),
                          holdout_max_ICp=max(max(abs(check['moment_defect'][i])
                                                  for i in (1, 4))
                                              for check in checks
                                              if not check['is_training']),
                          min_relative_E=min(check['min_relative_E']
                                             for check in checks),
                          remote_momentum=report['remote_momentum'])),
          flush=True)


if __name__ == '__main__':
    run()
