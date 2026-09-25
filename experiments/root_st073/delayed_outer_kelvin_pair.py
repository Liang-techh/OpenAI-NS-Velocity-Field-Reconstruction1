"""Frozen Kelvin covariance screen for the pressure-opened ten-mode cone.

This reuses the physical Cartesian velocity jets and the existing affine
pulse integrator. Positive distinct-angular-mode weights are local algebra,
not a compact amplitude/phase solution or a Navier--Stokes cancellation.
"""

import itertools
import json

import numpy as np

from affine_pulse import affine_pulse
from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME, load_cache
from delayed_multimode_pressure_admission import geometry
from delayed_outer_radial_extension import CACHE_NAME, signature
from delayed_pressure_cone_window import pressure_stress_columns
from delayed_similarity_pressure_screen import SimilarityPressurePatch
from radial_continuation import ROOT


ANGULAR_MODES = (1, 2, 4, 8)
AXIAL_RATIOS = (-1, 0, 1)
RADIAL_RATIOS = (-1, 0, 1)


def run():
    source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    pressure_source = json.loads((ROOT/'delayed_outer_pressure_fit.json').read_text())
    c = np.asarray(source['coefficients'])
    pressure_c = np.asarray(pressure_source['pressure_coefficients'])
    rows, _ = load_cache(ROOT/CACHE_NAME, signature())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    patch = SimilarityPressurePatch(base, pressure_c)
    tau = .5*2**(-5.5)
    results = []
    selected_source = None
    for row in rows:
        point = base.compact.joined.inner.from_similarity(
            np.array([row['X']]), np.array([row['eta']]), tau)[0]
        radius = float(point[0])
        u = row['u0']+row['U']@c
        J = row['g0']+np.einsum('abi,i->ab', row['G'], c)
        target, N, K, multiplier, lam2 = geometry(row, c)
        target += pressure_stress_columns(patch, point, tau)@pressure_c
        dot_n, dot_k = float(target@N), float(target@K)
        cone_margin = -dot_n-abs(multiplier*dot_k)
        pulses = []
        for mode in ANGULAR_MODES:
            for axial in AXIAL_RATIOS:
                for radial in RADIAL_RATIOS:
                    tangential = np.array([mode/radius,
                                           axial*mode/radius])
                    n0 = np.r_[radial*np.linalg.norm(tangential),
                               tangential]
                    pulse = affine_pulse(J, base.nu, n0, .1*tau)
                    pulses.append(dict(mode=mode, axial_ratio=axial,
                                       radial_ratio=radial, **pulse))
        options = []
        for i, j in itertools.combinations(range(len(pulses)), 2):
            if pulses[i]['mode'] == pulses[j]['mode']:
                continue
            matrix = np.array([
                pulses[i]['mean_radial_tangential_covariance'],
                pulses[j]['mean_radial_tangential_covariance']]).T
            condition = float(np.linalg.cond(matrix))
            if not np.isfinite(condition) or condition > 1e8:
                continue
            weights = np.linalg.solve(matrix, target)
            if np.min(weights) <= 0:
                continue
            options.append(dict(indices=(i, j),
                                weights=weights.tolist(),
                                modes=(pulses[i]['mode'],
                                       pulses[j]['mode']),
                                weight_sum=float(np.sum(weights)),
                                condition=condition))
        options.sort(key=lambda item: item['weight_sum'])
        result = dict(X=row['X'], eta=row['eta'],
                      point=point.tolist(),
                      lambda_squared=lam2,
                      cone_margin=float(cone_margin),
                      target=target.tolist(),
                      positive_distinct_pair_count=len(options),
                      best_pair=options[0] if options else None)
        results.append(result)
        print(json.dumps(dict(X=row['X'], eta=row['eta'],
                              cone_margin=result['cone_margin'],
                              pair_count=len(options),
                              best_modes=(options[0]['modes']
                                          if options else None))), flush=True)
        if (options and selected_source is None and
                row['X'] == 1.016 and row['eta'] == .25):
            pair = options[0]
            selected_source = dict(point=point.tolist(), tau=tau,
                                   velocity=u.tolist(), gradient=J.tolist(),
                                   local_tangential_stress_primitive=(
                                       target.tolist()),
                                   cone=dict(lambda_squared=lam2,
                                             margin=float(cone_margin)),
                                   kelvin_pulses=[pulses[i]
                                                  for i in pair['indices']],
                                   covariance_nnls=dict(
                                       selected=[dict(index=k,
                                                      weight=pair['weights'][k])
                                                 for k in (0, 1)],
                                       relative_error=0.),
                                   source_mean='delayed_outer_admission_search.json',
                                   source_pressure='delayed_outer_pressure_fit.json',
                                   scope='One-point frozen Kelvin two-mode '
                                         'covariance only. No supported '
                                         'phase/amplitude solution, exact-curl '
                                         'physical wave, or PDE acceptance.',
                                   accepted=False)
    report = dict(source_mean='delayed_outer_admission_search.json',
                  source_pressure='delayed_outer_pressure_fit.json',
                  angular_modes=ANGULAR_MODES,
                  axial_ratios=AXIAL_RATIOS,
                  radial_ratios=RADIAL_RATIOS,
                  duration=.1*tau, rows=results,
                  selected_wave_source=(
                      'delayed_outer_kelvin_source.json'
                      if selected_source is not None else None),
                  scope='Frozen physical-gradient Kelvin ray covariance '
                        'on twelve sampled pressure-opened cone nodes. '
                        'No compact wave, continuous cone, exact moments, '
                        'or full PDE acceptance.', accepted=False)
    (ROOT/'delayed_outer_kelvin_pair.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    if selected_source is not None:
        (ROOT/'delayed_outer_kelvin_source.json').write_bytes(
            (json.dumps(selected_source, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
