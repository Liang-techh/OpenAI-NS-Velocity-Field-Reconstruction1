"""Check whether implicit axis jets support a wider low-residual ST073 core."""
from dataclasses import asdict
import json

import numpy as np

from extended_axis_jet import ExtendedFullRadialField, ExtendedParameters
from full_radial import FullRadialField
from radial_continuation import ROOT


def run():
    source = FullRadialField.load(
        ROOT/'NS_ST073_Full_Local_Recurrence/data/ST073-V-wide14.json')
    params = asdict(source.p)
    params['eta_max'] = .85
    extended = ExtendedFullRadialField(ExtendedParameters(**params))
    tau = .5*2**(-5.5)
    rows = []
    for eta in (0., .3, .5, .65, .75, .85):
        X = np.array([.01, .03, .06])
        d = extended.evaluate_similarity(X, np.full(len(X), eta), tau)
        row = dict(eta=eta, X=X.tolist(),
                   extended_residual_norms=np.linalg.norm(d['residual'], axis=1).tolist(),
                   extended_speed=np.linalg.norm(d['velocity'], axis=1).tolist())
        if eta <= .5:
            old = source.evaluate_similarity(X, np.full(len(X), eta), tau)
            row['old_new_velocity_gap'] = float(np.max(np.linalg.norm(
                old['velocity']-d['velocity'], axis=1)))
            row['old_new_residual_gap'] = float(np.max(np.linalg.norm(
                old['residual']-d['residual'], axis=1)))
            if (row['old_new_velocity_gap'] > 1e-12
                    or row['old_new_residual_gap'] > 1e-9):
                raise AssertionError('Implicit jets changed the old core')
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, rows=rows,
                  scope='Extended implicit axis jets sampled in the local radial core; no global join or validation.',
                  accepted=False)
    (ROOT/'extended_axis_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
