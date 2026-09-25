"""Bounded multistart moment feasibility of the cached eight-mode family.

This minimizes the outgoing moment mismatch before attempting further cone
or momentum optimization. Numerical failure is not an impossibility proof.
"""

import json

import numpy as np
from scipy.optimize import least_squares

from delayed_multimode_cone_fit import (CACHE_NAME, cone, load_cache,
                                         moment_summary)
from radial_continuation import ROOT


SCALES = np.array((.001, .01, .01, .001))


def run():
    rows, moments = load_cache(ROOT/CACHE_NAME)
    previous = json.loads((ROOT/'delayed_multimode_cone_fit.json').read_text())
    lower = np.array([-.05]*8)
    upper = np.array([.3, .3, .05, .05]*2)

    def defects(c):
        return np.array([row['defect'][1:]
                         for row in moment_summary(moments, c)])

    def scaled(c):
        return (defects(c)/SCALES).ravel()

    rng = np.random.default_rng(73075)
    starts = [np.zeros(8), np.asarray(previous['coefficients'])]
    starts += [rng.uniform(lower*.3, upper*.3) for _ in range(10)]
    trials = []
    for start in starts:
        fit = least_squares(scaled, start, bounds=(lower, upper),
                            max_nfev=300, xtol=1e-11, ftol=1e-11,
                            gtol=1e-11)
        c = fit.x
        d = defects(c)
        cone_data = [cone(row, c) for row in rows]
        trials.append(dict(coefficients=c.tolist(),
                           scaled_l2=float(np.linalg.norm(scaled(c))),
                           scaled_max=float(np.max(np.abs(d/SCALES))),
                           raw_max=float(np.max(np.abs(d))),
                           max_cone_residual=float(max(
                               row['residual_norm'] for row in cone_data)),
                           cone_pass_count=sum(row['strict_pass']
                                               for row in cone_data),
                           min_relative_E=min(row['min_relative_E']
                                              for row in moment_summary(
                                                  moments, c)),
                           success=bool(fit.success), nfev=int(fit.nfev)))
    trials.sort(key=lambda row: row['scaled_l2'])
    report = dict(response_cache=CACHE_NAME, scales=SCALES.tolist(),
                  bounds=list(zip(lower.tolist(), upper.tolist())),
                  starts=len(starts), best=trials[0],
                  trials=trials,
                  scope='Bounded local multistart least-squares screen '
                        'of three-eta outgoing moments in the eight-mode '
                        'family; neither global infeasibility nor PDE '
                        'acceptance follows from numerical minima.',
                  accepted=False)
    (ROOT/'delayed_multimode_moment_feasibility.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(best=report['best'],
                          scaled_l2_range=(trials[0]['scaled_l2'],
                                           trials[-1]['scaled_l2']))),
          flush=True)


if __name__ == '__main__':
    run()
