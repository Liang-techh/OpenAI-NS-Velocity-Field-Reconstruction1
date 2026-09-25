"""Connect the extended ST073 shear background to the existing wave source.

Generate a physical radial stress primitive, a local Kelvin cone, and the
two-direction covariance fit at an actual bridge point. This is a pointwise
preflight only; no wave is promoted into the global candidate.
"""
import json

import numpy as np

from high_frequency_shear_screen import make_field
from radial_peak_cone import run as build_source
from radial_continuation import ROOT


def run():
    tau = .5*2**(-5.5)
    eta = .2
    X = 1.
    rows = []
    for label, amplitude in [('balanced', 0.), ('shear_N16', 2.)]:
        field = make_field(16, amplitude)
        point = field.compact.joined.inner.from_similarity(
            np.array([X]), np.array([eta]), tau)[0]
        source_name = 'st073_'+label+'_wave_source.json'
        build_source(radius=float(point[0]), z=float(point[2]),
                     output_name=source_name, field=field,
                     field_id='extended_'+label, tau=tau,
                     stress_order=64)
        source = json.loads((ROOT/'compact_potential'/source_name).read_text())
        selected = source['covariance_nnls']['selected']
        rows.append(dict(field=label, point=point.tolist(),
                         cone=source['cone'],
                         covariance_relative_error=source['covariance_nnls']['relative_error'],
                         selected_pulse_count=len(selected),
                         selected_modes=[source['kelvin_pulses'][s['index']]['mode']
                                         for s in selected],
                         target=source['local_tangential_stress_primitive']))
    report = dict(tau=tau, X=X, eta=eta, rows=rows,
                  scope='Balanced and high-frequency extended mean fields at one physical bridge point. Kelvin pulses are frozen local ODEs; covariance fit is not a supported exact-curl wave or complete momentum reduction.',
                  accepted=False)
    (ROOT/'extended_wave_source_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run()
