"""Resolve the radial width of the sampled physical stress-cone window."""
import json

import numpy as np

from extended_physical_cone_map import physical_cone
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def run():
    tau = .5*2**(-5.5)
    field = RadialMomentStep(make_field(16, 2.), -.5)
    radii = (.975, .985, .995, 1., 1.005, 1.015, 1.025)
    rows = []
    for eta in (.2, .3):
        for X in radii:
            row = physical_cone(field, X, eta, tau, order=64)
            rows.append(row)
            print(json.dumps(row), flush=True)

    # The moment step is exactly absent after its restoration radius.
    base = field.base
    outer_points = base.compact.joined.inner.from_similarity(
        np.array([3.05, 3.5, 4.]), np.array([.2, .3, .65]), tau)
    base_u, base_p = base.fields(outer_points, tau)
    step_u, step_p = field.fields(outer_points, tau)
    outer_match = dict(
        max_velocity_difference=float(np.max(np.abs(step_u-base_u))),
        max_pressure_difference=float(np.max(np.abs(step_p-base_p))))
    report = dict(tau=tau, shear=2., moment=-.5, rows=rows,
                  pass_count=sum(row['strict_pass'] for row in rows),
                  outer_match=outer_match,
                  scope='Finite radial sample of physical full-residual cone at two axial slices, not a continuous or all-time certificate. Outer field equality applies to the added step only. No supported wave or PDE acceptance.',
                  accepted=False)
    (ROOT/'moment_shear_band_map.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(pass_count=report['pass_count'],
                          outer_match=outer_match)), flush=True)


if __name__ == '__main__':
    run()
