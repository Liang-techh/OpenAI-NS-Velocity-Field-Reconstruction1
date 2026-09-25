"""Joint ten-mode search for outer-radius free-axial cone admission.

Free axial stress is a necessary pointwise test for a later compact
pressure correction, not a realizable pressure or oscillatory wave.
The search keeps moment and nearby-time momentum debts in the objective.
"""

import json

import numpy as np
from scipy.optimize import differential_evolution

from delayed_multimode_cone_fit import cone, evaluate_holdout, load_cache
from delayed_momentum_tangent_screen import stats
from delayed_multimode_pressure_admission import geometry
from delayed_outer_radial_extension import (CACHE_NAME, NEARBY_CACHE_NAME,
                                             metrics, pressure_admission,
                                             signature)
from radial_continuation import ROOT


def free_axial_gap(row, c):
    target, N, K, multiplier, lam2 = geometry(row, c)
    if multiplier is None or lam2 <= 0:
        return 1.
    directions = np.array([.8*N-multiplier*K,
                           .8*N+multiplier*K])
    slopes = directions[:, 1]
    intercepts = directions@target+.1
    if slopes[0]*slopes[1] < 0:
        z = (intercepts[1]-intercepts[0])/(slopes[0]-slopes[1])
        violation = max(0., slopes[0]*z+intercepts[0])
    elif np.any(np.abs(slopes) < 1e-14):
        violation = max(0., *intercepts[np.abs(slopes) < 1e-14])
    else:
        violation = 0.
    return float(violation/(1+np.sum(np.abs(intercepts))))


def run():
    rows, moments = load_cache(ROOT/CACHE_NAME, signature())
    with np.load(ROOT/NEARBY_CACHE_NAME, allow_pickle=False) as data:
        nearby = {key: data[key] for key in ('R0', 'L', 'Q')}
    prior = json.loads((ROOT/'delayed_outer_radial_extension.json').read_text())
    incumbent = np.asarray(prior['results'][1]['coefficients'])
    base_train = max(cone(row, np.zeros(10))['residual_norm'] for row in rows)
    base_nearby = stats(nearby['R0'])['max']
    incumbent_nearby = metrics(rows, moments, nearby, incumbent)['nearby_max']
    outer_rows = [row for row in rows if row['X'] == 1.02]

    def score(c):
        m = metrics(rows, moments, nearby, c)
        gap = max(free_axial_gap(row, c) for row in outer_rows)
        return (4*gap+.6*m['cone_max']+.15*m['cone_mean']
                +.08*m['moment_scaled_max']
                +.05*m['train_max']/base_train
                +max(0., m['nearby_max']/base_nearby-1.)
                +.5*max(0., m['nearby_max']/incumbent_nearby-1.)
                +100*max(0., .11005-m['min_relative_E']))

    bounds = ([(-.05, .3), (-.05, .3), (-.05, .05), (-.05, .05)]*2
              +[(-.1, .1), (-.03, .03)])
    fit = differential_evolution(score, bounds, x0=incumbent,
                                 seed=73078, maxiter=160, popsize=8,
                                 tol=1e-6, polish=True)
    c = fit.x
    report = dict(response_cache=CACHE_NAME,
                  nearby_response_cache=NEARBY_CACHE_NAME,
                  incumbent_coefficients=incumbent.tolist(),
                  incumbent_metrics=metrics(rows, moments, nearby,
                                            incumbent),
                  incumbent_outer_gaps=[free_axial_gap(row, incumbent)
                                        for row in outer_rows],
                  coefficients=c.tolist(),
                  selected_metrics=metrics(rows, moments, nearby, c),
                  selected_outer_gaps=[free_axial_gap(row, c)
                                       for row in outer_rows],
                  pressure_admission=pressure_admission(rows, c),
                  cone_rows=[cone(row, c) for row in rows],
                  objective=float(fit.fun),
                  optimizer_success=bool(fit.success),
                  optimizer_message=str(fit.message),
                  scope='Sampled ten-mode pointwise necessary admission '
                        'search with nearby-time momentum and moments. '
                        'No pressure realization, independent holdout, '
                        'continuous cone, wave, or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_outer_admission_search.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(incumbent_gaps=report['incumbent_outer_gaps'],
                          selected_gaps=report['selected_outer_gaps'],
                          selected_metrics=report['selected_metrics'],
                          pressure_count=sum(
                              row['free_axial_080']
                              for row in report['pressure_admission']),
                          objective=report['objective'])), flush=True)


if __name__ == '__main__':
    run()
