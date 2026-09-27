"""Integrated conservation rows for local state slopes and pressure.

Use the same fixed-radius integral derivative stencil as the independent
moment replay, instead of integrating pointwise second spatial derivatives.
This is still a finite quadrature/stencil operator, not an exact certificate.
"""
import json
from pathlib import Path
import time
import numpy as np
from meridional_state_cache import StatefulMean, value_and_pressure_modes
from midplane_integrated_moment_balance import integral_state, evaluate
from broad_shear_dynamic_control import load_saved_field
from grouped_joined_field import install_in_field


def assemble(base, values, pressures, state, k, breaks, order=96, zfactor=.002):
    tau = .5 * 2.**(-k)
    field = StatefulMean(base, values, pressures, state,
                         np.zeros(len(values)), np.zeros(len(pressures)), k)
    offset = evaluate(field, base.inner, k, order, zfactor=zfactor,
                      radial_breaks=breaks)
    rows = []
    coefficients = np.array([1., -8., 8., -1.])
    for eta in (-.2, .2):
        R, _, z = base.inner.from_similarity(
            [base.inner.p.X_max*16**2], [eta], tau)[0]
        zp = base.inner.from_similarity([base.inner.p.X_max], [eta+.01], tau)[0, 2]
        zm = base.inner.from_similarity([base.inner.p.X_max], [eta-.01], tau)[0, 2]
        hz = zfactor * abs(zp-zm) / .02
        ht = 1e-4*tau
        block = np.zeros((2, len(values)+len(pressures)))
        for i, unit in enumerate(values):
            samples = []
            for j in (-2, -1, 1, 2):
                tt = tau+j*ht
                moment = integral_state(unit, base.inner, R, z, tt, order, breaks)[:2]
                samples.append((-np.log2(2*tt)-k)*moment)
            derivative = coefficients @ np.asarray(samples) / (12*ht)
            block[:, i] = derivative / np.array([R**2, R])
        for i, unit in enumerate(pressures):
            samples = [integral_state(unit, base.inner, R, z+j*hz, tau,
                                      order, breaks)[4] for j in (-2, -1, 1, 2)]
            block[1, len(values)+i] = -(coefficients @ samples)/(12*hz*R)
        rows.append(block)
    return np.concatenate(rows), offset


def run():
    root = Path(__file__).resolve().parent
    saved = json.loads((root/'broad_meridional_constrained.json').read_text())
    replay = json.loads((root/'broad_meridional_constrained_replay.json').read_text())
    dynamic, source = load_saved_field()
    install_in_field(dynamic)
    base, values, pressures, _, _ = value_and_pressure_modes(dynamic, source)
    state = np.zeros(len(values)); state[0] = .01; state[15] = .01
    control = np.asarray(replay['coefficients'])
    breaks = saved['quadrature']['radial_split_breaks']
    started = time.perf_counter()
    E, m = assemble(base, values, pressures, state, source['k'], breaks)
    field = StatefulMean(base, values, pressures, state, control[:len(values)],
                         control[len(values):], source['k'])
    direct = evaluate(field, base.inner, source['k'], 96, radial_breaks=breaks)
    predicted = E @ control + m
    result = dict(accepted=False, pde_validated=False, scale_recursion_established=False,
        scope='Nonzero-state integrated moment linearity check, order96, fixed-radius derivative stencils. Not a trajectory or continuum gate.',
        k=source['k'], state=state.tolist(), order=96,
        predicted=predicted.tolist(), direct=direct.tolist(),
        maximum_absolute_difference=float(np.max(abs(predicted-direct))),
        E=E.tolist(), m=m.tolist(), elapsed_seconds=time.perf_counter()-started)
    root.joinpath('integrated_state_moments.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('E','m','state')}), flush=True)


if __name__ == '__main__':
    run()
