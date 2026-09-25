"""Lagrangian transport inside the low-residual ST073-V wide core only."""
import json
import sys

import numpy as np
from scipy.integrate import solve_ivp

from radial_continuation import ROOT, FullRadialField

sys.path.insert(0, str(ROOT/'NS_ST073_Full_Local_Recurrence/upstream'))
from source_coordinates import coordinates


def similarity(field, radius, z, tau):
    co = coordinates(radius/np.sqrt(field.nu), z/np.sqrt(field.nu), tau,
                     field.h)
    return float(co['X']), float(co['eta'])


def trace(field, X0, eta0, tau0=.128, tau_end=.0084):
    start = field.from_similarity(np.array([X0]), np.array([eta0]), tau0)[0]
    initial = np.array([np.hypot(start[0], start[1]), start[2], 0.])

    def rhs(tau, state):
        radius, z, _ = state
        X, eta = similarity(field, max(radius, 0.), z, tau)
        # The extension is evaluated only to locate the first exit event.
        X = np.clip(X, 0., field.p.X_max)
        eta = np.clip(eta, -field.p.eta_max, field.p.eta_max)
        vel = field.evaluate_similarity(np.array([X]), np.array([eta]), tau)[
            'velocity'][0]
        return [-vel[0], -vel[2], -vel[1]/max(radius, 1e-12)]

    def radial_exit(tau, state):
        return field.p.X_max-similarity(field, max(state[0], 0.), state[1], tau)[0]

    def axial_exit(tau, state):
        return field.p.eta_max-abs(similarity(
            field, max(state[0], 0.), state[1], tau)[1])

    def axis_exit(tau, state):
        return state[0]-1e-10

    for event in (radial_exit, axial_exit, axis_exit):
        event.terminal = True
        event.direction = 0
    times = np.geomspace(tau0, tau_end, 25)
    sol = solve_ivp(rhs, (tau0, tau_end), initial, t_eval=times,
                    events=(radial_exit, axial_exit, axis_exit),
                    rtol=2e-7, atol=2e-9, max_step=.004)
    if not sol.success:
        raise RuntimeError(sol.message)
    reason = ('radial_exit', 'axial_exit', 'axis_exit')
    exit_reason = next((reason[i] for i, ev in enumerate(sol.t_events)
                        if len(ev)), 'reached_final_time')
    rows = []
    for tau, state in zip(sol.t, sol.y.T):
        r, z, theta = state
        X, eta = similarity(field, r, z, tau)
        velocity = field.evaluate_similarity(
            np.array([X]), np.array([eta]), tau)['velocity'][0]
        rows.append(dict(tau=float(tau), radius=float(r), z=float(z),
                         theta=float(theta), turns=float(theta/(2*np.pi)),
                         X=float(X), eta=float(eta),
                         angular_speed=float(velocity[1]/r)))
    return dict(initial_X=X0, initial_eta=eta0, exit_reason=exit_reason,
                sampled_rows=rows, solver_nfev=int(sol.nfev),
                first_to_last_sample_radius_ratio=float(
                    rows[-1]['radius']/rows[0]['radius']),
                first_to_last_sample_turns=float(
                    rows[-1]['turns']-rows[0]['turns']))


def run():
    field = FullRadialField.load(
        ROOT/'NS_ST073_Full_Local_Recurrence/data/ST073-V-wide14.json')
    trajectories = []
    for X0, eta0 in ((.02, 0.), (.05, 0.)):
        row = trace(field, X0, eta0)
        trajectories.append(row)
        print(json.dumps({k: row[k] for k in
                          ('initial_X', 'initial_eta', 'exit_reason',
                           'solver_nfev', 'first_to_last_sample_radius_ratio',
                           'first_to_last_sample_turns')}), flush=True)
    report = dict(model='ST073-V-wide14', trajectories=trajectories,
                  scope='Physical-time particle paths in the registered local '
                        'field, stopped at the first radial/axial exit. '
                        'Accumulated angle is true Lagrangian winding only '
                        'while the particle stays in the local domain. No '
                        'whole-field vortex-core geometry, finite global '
                        'energy, or critical-time limit.', accepted=False)
    (ROOT/'wide_core_trajectories.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
