"""Choose an axial support width by the complete collar momentum residual."""
import json

import numpy as np

from axial_compact_join import AxiallyCompactField
from joined_field import independent_fd
from radial_continuation import ROOT
from wide_pressure_fit import load_pressure_candidate


def screen(field, k):
    tau = .5*2**(-k)
    fractions = np.array([.25, .5, .75])
    etas = field.eta_flat+(field.eta_outer-field.eta_flat)*fractions
    X, eta = np.meshgrid([.03, .5859375], etas, indexing='ij')
    points = field.joined.inner.from_similarity(X.ravel(), eta.ravel(), tau)
    residual, divergence = independent_fd(
        field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
    norms = np.linalg.norm(residual, axis=1)
    return dict(k=k, tau=tau, eta=etas.tolist(),
                X=X.ravel().tolist(),
                residual_norms=norms.tolist(),
                sampled_max=float(norms.max()),
                vector_rms=float(np.sqrt(np.mean(residual**2))),
                divergence_max=float(np.max(np.abs(divergence))))


def run():
    base = load_pressure_candidate()
    configs = ((.3, .45), (.3, .49), (.25, .49), (.2, .49))
    rows = []
    for flat, outer in configs:
        field = AxiallyCompactField(base, flat, outer)
        result = screen(field, 5.5)
        row = dict(eta_flat=flat, eta_outer=outer, late=result)
        rows.append(row)
        print(json.dumps(dict(eta_flat=flat, eta_outer=outer,
                              late_max=result['sampled_max'],
                              late_rms=result['vector_rms'])), flush=True)
    best = min(rows, key=lambda r: r['late']['sampled_max'])
    chosen = AxiallyCompactField(base, best['eta_flat'], best['eta_outer'])
    mid = screen(chosen, 3.)
    report = dict(rows=rows, selected=dict(
        eta_flat=best['eta_flat'], eta_outer=best['eta_outer'],
        mid_scale=mid),
        scope='Positive axial collar; two radial locations and three relative '
              'axial positions per configuration at k=5.5. One k=3 check '
              'on selected width. All fields are globally callable, '
              'solenoidal and finite-energy at registered positive times; '
              'selection is sampled and does not establish full-domain '
              'momentum acceptance.', accepted=False)
    (ROOT/'axial_cutoff_width_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
