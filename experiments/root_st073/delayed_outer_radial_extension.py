"""Add two compact outer-transition velocity controls to the eight-mode mean.

The new swirl and exact-curl modes have radial support shifted toward
X=1.02. A cached ten-mode quadratic response includes their cross-advection
with every existing mode. Only the two new amplitudes are fitted here; the
two incumbent eight-mode candidates are kept as distinct starting fields.
"""

import json

import numpy as np
from scipy.optimize import differential_evolution

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode
from delayed_momentum_tangent_screen import stats
from delayed_multimode_cone_fit import (BASE_NAME, cone, evaluate_holdout,
                                         holdout_response, load_cache,
                                         make_modes, moment_rows,
                                         moment_summary, precompute,
                                         save_cache)
from delayed_multimode_pressure_admission import geometry
from delayed_pressure_cone_geometry import line_interval
from delayed_similarity_curl_screen import SimilarityCurlMode
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT


OUTER_RADIAL = (1.01, 1.05)
OUTER_AXIAL = (.14, .36)
CACHE_NAME = 'delayed_outer_radial_extension_model.npz'
NEARBY_CACHE_NAME = 'delayed_outer_radial_extension_nearby.npz'
SCALES = np.array((.001, .01, .01, .001))


def signature():
    return json.dumps(dict(base=BASE_NAME, outer_radial=OUTER_RADIAL,
                           outer_axial=OUTER_AXIAL, prior_modes=8,
                           full_modes=10, order=48, version=1),
                      sort_keys=True)


def make_extended_modes(base):
    modes, _ = make_modes(base)
    modes += [SimilaritySwirlMode(base, OUTER_RADIAL,
                                  eta_interval=OUTER_AXIAL),
              SimilarityCurlMode(base, OUTER_RADIAL,
                                 eta_interval=OUTER_AXIAL)]
    return modes


def load_or_build(base, modes, tau):
    path = ROOT/CACHE_NAME
    if path.exists():
        rows, moments = load_cache(path, signature())
        print('Loaded ten-mode physical response', flush=True)
    else:
        rows = precompute(base, modes, tau,
                          extra_radial_edges=OUTER_RADIAL)
        moments = moment_rows(base, make_field(16, 2.), modes, tau)
        save_cache(path, rows, moments, signature())
        print('Built ten-mode physical response', flush=True)
    nearby_path = ROOT/NEARBY_CACHE_NAME
    if nearby_path.exists():
        with np.load(nearby_path, allow_pickle=False) as data:
            if str(data['signature']) != signature():
                raise ValueError('Nearby response signature mismatch')
            nearby = {key: data[key] for key in ('tau', 'X', 'eta',
                                                  'R0', 'L', 'Q')}
        print('Loaded ten-mode nearby response', flush=True)
    else:
        nearby = holdout_response(base, modes)
        np.savez_compressed(nearby_path, signature=np.array(signature()),
                            **nearby)
        print('Built ten-mode nearby response', flush=True)
    return rows, moments, nearby


def pressure_admission(rows, c):
    result = []
    for row in rows:
        target, N, K, multiplier, lam2 = geometry(row, c)
        feasible = (lam2 > 0 and multiplier is not None and
                    line_interval([.8*N+sign*multiplier*K
                                   for sign in (-1., 1.)],
                                  target, .1, 1) is not None)
        result.append(dict(X=row['X'], eta=row['eta'],
                           free_axial_080=bool(feasible)))
    return result


def metrics(rows, moments, nearby, c):
    cones = [cone(row, c) for row in rows]
    violation = np.array([row['cone_violation'] for row in cones])
    moment_data = moment_summary(moments, c)
    defects = np.array([row['defect'][1:] for row in moment_data])
    nearby_stats = stats(evaluate_holdout(nearby, c))
    return dict(cone_pass_count=sum(row['strict_pass'] for row in cones),
                cone_max=float(np.max(violation)),
                cone_mean=float(np.mean(violation)),
                outer_cone_max=float(np.max([
                    row['cone_violation'] for row in cones
                    if row['X'] == 1.02])),
                min_margin=float(min(row['margin'] for row in cones)),
                moment_scaled_max=float(np.max(np.abs(defects/SCALES))),
                moment_raw_max=float(np.max(np.abs(defects))),
                min_relative_E=float(min(row['min_relative_E']
                                         for row in moment_data)),
                train_max=float(max(row['residual_norm'] for row in cones)),
                nearby_max=nearby_stats['max'],
                nearby_rms=nearby_stats['rms'])


def run(cache_only=False):
    tau = .5*2**(-5.5)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    modes = make_extended_modes(base)
    rows, moments, nearby = load_or_build(base, modes, tau)
    if cache_only:
        return
    old = json.loads((ROOT/'delayed_multimode_cone_fit.json').read_text())
    tradeoff = json.loads((ROOT/'delayed_multimode_tradeoff_fit.json').read_text())
    baselines = [('momentum_favored', old['coefficients']),
                 ('cone_moment_favored',
                  tradeoff['searches'][0]['coefficients'])]
    base_train = metrics(rows, moments, nearby, np.zeros(10))['train_max']
    base_nearby = stats(nearby['R0'])['max']
    results = []
    for name, incumbent in baselines:
        start = np.r_[incumbent, 0., 0.]
        initial = metrics(rows, moments, nearby, start)

        def score(pair):
            c = np.r_[incumbent, pair]
            m = metrics(rows, moments, nearby, c)
            return (m['cone_max']+.25*m['cone_mean']
                    +.35*m['outer_cone_max']
                    +.08*m['moment_scaled_max']
                    +.03*m['train_max']/base_train
                    +max(0., m['nearby_max']/base_nearby-1.)
                    +.5*max(0., m['nearby_max']/initial['nearby_max']-1.)
                    +100*max(0., .11005-m['min_relative_E']))

        fit = differential_evolution(score, [(-.1, .1), (-.03, .03)],
                                     seed=73077, maxiter=180, popsize=14,
                                     tol=1e-6, polish=True)
        c = np.r_[incumbent, fit.x]
        selected = metrics(rows, moments, nearby, c)
        results.append(dict(incumbent=name,
                            incumbent_coefficients=incumbent,
                            initial=initial, added_coefficients=fit.x.tolist(),
                            coefficients=c.tolist(), selected=selected,
                            initial_pressure=pressure_admission(rows, start),
                            selected_pressure=pressure_admission(rows, c),
                            cone_rows=[cone(row, c) for row in rows],
                            objective=float(fit.fun),
                            optimizer_success=bool(fit.success),
                            optimizer_message=str(fit.message)))
        print(json.dumps(dict(incumbent=name, added=fit.x.tolist(),
                              initial=initial, selected=selected)),
              flush=True)
    report = dict(base_slice=BASE_NAME, tau=tau,
                  outer_radial=OUTER_RADIAL, outer_axial=OUTER_AXIAL,
                  response_cache=CACHE_NAME,
                  nearby_response_cache=NEARBY_CACHE_NAME,
                  results=results,
                  scope='Ten-mode quadratic physical response with two '
                        'new outer-transition amplitudes optimized '
                        'against sampled cone, moments and nearby-time '
                        'momentum. No independent holdout, continuous '
                        'cone, exact moments, wave, or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_outer_radial_extension.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    import sys
    run(cache_only='--cache-only' in sys.argv[1:])
