"""Find a wider low-residual radial interface for the frozen ST073-V axis data.

This changes the radial truncation order, not the frozen ST073-V candidate.
Every result is a sampled local diagnostic, not a global PDE certificate.
"""
from dataclasses import asdict, replace
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from radial_continuation import ROOT, FullRadialField


def run():
    frozen = FullRadialField.load(
        ROOT / 'NS_ST073_Full_Local_Recurrence' / 'data' / 'ST073-V.json')
    radii = np.array([1/64, 3/128, 1/32, 3/64, 1/16, 3/32, 1/8, 3/16, 1/4])
    etas = np.linspace(-.4, .4, 9)
    ks = (0, 3, 6)
    rows = []
    for order in (8, 10, 12, 14):
        field = FullRadialField(replace(frozen.p, order=order, X_max=.25))
        for k in ks:
            tau = .5 * 2**(-k)
            for X in radii:
                out = field.evaluate_similarity(
                    np.full(len(etas), X), etas, tau)
                norms = np.linalg.norm(out['residual'], axis=1)
                rows.append(dict(order=order, k=k, X=float(X),
                                 sampled_residual_max=float(norms.max()),
                                 sampled_divergence_max=float(
                                     np.max(np.abs(out['divergence'])))))
            print(json.dumps(dict(order=order, k=k, rows=[
                r for r in rows if r['order'] == order and r['k'] == k])),
                flush=True)
    widest = {}
    for order in (8, 10, 12, 14):
        allowed = [float(X) for X in radii if all(
            next(r['sampled_residual_max'] for r in rows
                 if r['order'] == order and r['k'] == k and r['X'] == X)
            < 1e-3 for k in ks)]
        widest[str(order)] = max(allowed) if allowed else None
    wide = FullRadialField(replace(frozen.p, order=14, X_max=3/32))
    # This is a distinct callable model: never relabel the frozen order-8 field.
    model_path = ROOT / 'NS_ST073_Full_Local_Recurrence' / 'data' / 'ST073-V-wide14.json'
    model = dict(id='ST073-V-wide14', parent='ST073-V',
                 parameters=asdict(wide.p), pde_validated=False,
                 global_field_ready=False,
                 scope='Sampled wider local radial recurrence only; no outer or axial join')
    model_path.write_text(json.dumps(model, indent=2) + '\n')
    gx, wx = leggauss(12)
    ge, we = leggauss(16)
    X = (gx + 1) * wide.p.X_max / 2
    eta = ge * .4
    xx, ee = np.meshgrid(X, eta, indexing='ij')
    quadrature = []
    for k in ks:
        tau = .5 * 2**(-k)
        q = tau / (1 - ee**2)
        volume_weights = (wx[:, None] * wide.p.X_max / 2
                          * we[None, :] * .4
                          * 2*np.pi*wide.nu**1.5*q**(1+wide.D)
                          * (1-2*wide.h*ee**2)/(1-ee**2))
        out = wide.evaluate_similarity(xx, ee, tau)
        norm_sq = np.sum(out['residual']**2, axis=-1)
        quadrature.append(dict(k=k, point_count=int(xx.size),
                               region_volume=float(volume_weights.sum()),
                               sampled_max=float(np.sqrt(norm_sq.max())),
                               physical_volume_L2=float(np.sqrt(
                                   np.sum(volume_weights*norm_sq))),
                               divergence_max=float(np.max(
                                   np.abs(out['divergence'])))))
    report = dict(rows=rows, etas=etas.tolist(), ks=list(ks),
                  wide_candidate=str(model_path.relative_to(ROOT)),
                  wide_quadrature=quadrature,
                  widest_sampled_interface_under_1e_minus3=widest,
                  scope='Frozen ST073-V axis data with new radial truncation orders. '
                        'Nine eta samples at three times. No off-grid or complete '
                        'annulus certificate, outer match, finite energy, or force '
                        'acceptance.', accepted=False)
    (ROOT / 'high_order_interface_screen.json').write_text(
        json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    run()
