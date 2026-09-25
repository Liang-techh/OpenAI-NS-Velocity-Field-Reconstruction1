"""First dynamic check of a widened ST073-V radial join."""
import json

import numpy as np

from radial_continuation import ROOT, FullRadialField
from joined_field import JoinedField, independent_fd


def run():
    inner = FullRadialField.load(
        ROOT / 'NS_ST073_Full_Local_Recurrence' / 'data' / 'ST073-V-wide14.json')
    rows = []
    for outer_ratio, k in ((2., .4), (2., 5.5), (3., 5.5),
                           (4., 5.5), (6., 5.5)):
        field = JoinedField(inner=inner, join_X=3/32,
                            outer_ratio=outer_ratio)
        tau = .5 * 2**(-k)
        fractions = np.array([.2, .5, .8])
        X = (3/32) * (1+(outer_ratio-1)*fractions)**2
        eta = .2
        point = inner.from_similarity(X, np.full(3, eta), tau)
        h = .001 * np.sqrt(inner.nu*tau)
        residual, div = independent_fd(field, point, tau, h, .00025*tau)
        norms = np.linalg.norm(residual, axis=1)
        row = dict(k=k, tau=tau, X=X.tolist(), eta=eta,
                   outer_ratio=outer_ratio,
                   radial_fractions=fractions.tolist(),
                   physical_points=point.tolist(),
                   residual=residual.tolist(),
                   residual_norms=norms.tolist(),
                   sampled_max=float(norms.max()),
                   divergence_max=float(np.max(np.abs(div))))
        rows.append(row)
        print(json.dumps(dict(k=k, outer_ratio=outer_ratio,
                              residual_norms=row['residual_norms'],
                              sampled_max=row['sampled_max'])), flush=True)
    report = dict(join_X=3/32, inner_order=14,
                  heat_amplitude=field.c, rows=rows,
                  scope='C2 radial join from the wider local recurrence to the '
                        'old pure-swirl heat exterior. Fourth-order Cartesian '
                        'finite differences at three radial fractions with several '
                        'outer radius ratios. Ratios are different physical '
                        'points, not identical-domain holdouts. Axial '
                        'closure, five radial moments, finite total energy, '
                        'and global residual remain open.',
                  pde_validated=False, global_field_ready=False)
    (ROOT / 'wide_join_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
