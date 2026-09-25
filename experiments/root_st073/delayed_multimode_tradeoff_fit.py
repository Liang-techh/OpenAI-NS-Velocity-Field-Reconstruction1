"""Cached eight-mode cone/moment/momentum Pareto screen.

The original selected point overweights sampled cone admission relative to
outgoing moments. These searches expose the tradeoff at three moment weights.
The previously unused disjoint holdout is not used for selection.
"""

import json

import numpy as np
from scipy.optimize import differential_evolution

from delayed_multimode_cone_fit import (CACHE_NAME, HOLDOUT_CACHE_NAME,
                                         cone, evaluate_holdout, load_cache,
                                         moment_summary)
from delayed_momentum_tangent_screen import stats
from delayed_multimode_pressure_admission import geometry
from delayed_pressure_cone_geometry import line_interval
from radial_continuation import ROOT


SCALES = (.001, .01, .01, .001)


def run():
    rows, moments = load_cache(ROOT/CACHE_NAME)
    with np.load(ROOT/HOLDOUT_CACHE_NAME, allow_pickle=False) as data:
        nearby = {key: data[key] for key in ('R0', 'L', 'Q')}
    original = json.loads((ROOT/'delayed_multimode_cone_fit.json').read_text())
    prior = np.asarray(original['coefficients'])
    base_train = max(cone(row, np.zeros(8))['residual_norm'] for row in rows)
    base_nearby = stats(nearby['R0'])['max']

    def metrics(c):
        cones = [cone(row, c) for row in rows]
        violations = np.array([row['cone_violation'] for row in cones])
        moment_data = moment_summary(moments, c)
        defects = np.array([row['defect'][1:] for row in moment_data])
        nearby_stats = stats(evaluate_holdout(nearby, c))
        return dict(cone_max=float(np.max(violations)),
                    cone_mean=float(np.mean(violations)),
                    cone_pass_count=sum(row['strict_pass'] for row in cones),
                    moment_scaled_max=float(np.max(
                        np.abs(defects/np.asarray(SCALES)))),
                    moment_raw_max=float(np.max(np.abs(defects))),
                    min_relative_E=min(row['min_relative_E']
                                       for row in moment_data),
                    train_max=float(max(row['residual_norm'] for row in cones)),
                    nearby_max=nearby_stats['max'],
                    nearby_rms=nearby_stats['rms'])

    def score(m, weight):
        return (m['cone_max']+.25*m['cone_mean']
                +weight*m['moment_scaled_max']
                +.03*m['train_max']/base_train
                +max(0., m['nearby_max']/base_nearby-1.)
                +100*max(0., .11005-m['min_relative_E']))

    def pressure_admission(c):
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

    bounds = [(-.05, .3), (-.05, .3), (-.05, .05), (-.05, .05)]*2
    searches = []
    for weight in (.05, .15, .4):
        fit = differential_evolution(lambda c: score(metrics(c), weight),
                                     bounds, seed=int(73076+weight*100),
                                     maxiter=100, popsize=6, polish=True,
                                     tol=1e-5)
        m = metrics(fit.x)
        searches.append(dict(moment_weight=weight,
                             coefficients=fit.x.tolist(), metrics=m,
                             pressure_admission=pressure_admission(fit.x),
                             objective=float(fit.fun),
                             previous_objective=float(
                                 score(metrics(prior), weight)),
                             optimizer_success=bool(fit.success),
                             optimizer_message=str(fit.message)))
    report = dict(source='delayed_multimode_cone_fit.json',
                  response_cache=CACHE_NAME,
                  nearby_response_cache=HOLDOUT_CACHE_NAME,
                  previous=metrics(prior), searches=searches,
                  scope='Three weighted numerical searches in the '
                        'existing eight-mode response; the separate '
                        'disjoint holdout is excluded from fitting. '
                        'Sampled scores do not establish continuous '
                        'cone, exact moments, or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_multimode_tradeoff_fit.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(previous=report['previous'],
                          searches=[dict(moment_weight=x['moment_weight'],
                                         metrics=x['metrics'],
                                         objective=x['objective'])
                                    for x in searches])), flush=True)


if __name__ == '__main__':
    run()
