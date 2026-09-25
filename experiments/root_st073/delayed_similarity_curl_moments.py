"""Audit fixed-slice moments of the self-similar curl momentum fit."""

import json
import sys

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_similarity_curl_screen import (
    CurlPatchedLift, INTERVALS, SimilarityCurlMode,
)
from delayed_taper_capacity_screen import grid
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


def run(screen_name='delayed_similarity_curl_screen.json'):
    tau = .5*2**(-5.5)
    order = 64 if 'constrained' in screen_name else 48
    X, weights = grid(order=order)
    screen = json.loads((ROOT/screen_name).read_text())
    base = CoupledMomentPhysicalLift(slice_filename=screen['base_slice'])
    if 'eta_intervals' in screen:
        modes = [SimilarityCurlMode(base, radial, axial)
                 for axial in screen['eta_intervals']
                 for radial in screen['radial_intervals']]
        coefficients = screen['coefficients']
    else:
        modes = [SimilarityCurlMode(base, interval) for interval in INTERVALS]
        coefficients = screen['selected']['coefficients']
    patched = CurlPatchedLift(base, modes, coefficients)
    target = make_field(16, 2.)
    rows = []
    for eta in (.2, .25, .3):
        U0, E0 = profile(target, X, eta, tau)
        U_base, E_base = profile(base, X, eta, tau)
        U, E = profile(patched, X, eta, tau)
        before = moment_vector(U_base, E_base, X, weights)
        after = moment_vector(U, E, X, weights)
        target_vector = moment_vector(U0, E0, X, weights)
        row = dict(eta=eta,
                   patch_delta=(after-before).tolist(),
                   total_defect=(after-target_vector).tolist(),
                   max_abs_total_defect=float(np.max(np.abs(
                       after-target_vector))),
                   min_relative_E=float(np.min(E/E0)))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, source_screen=screen_name,
                  quadrature_order=order, rows=rows,
                  scope='Actual physical velocity five moments at two '
                        'reference slices and one intervening slice after '
                        'the fitted exact-curl modes. No continuous '
                        'eta/time moment identity.',
                  accepted=False)
    output_name = screen_name.replace('.json', '_moments.json')
    (ROOT/output_name).write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else
        'delayed_similarity_curl_screen.json')
