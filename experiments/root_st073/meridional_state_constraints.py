"""Refresh finite moment/cone constraints from the current nonlinear state.

These sampled constraints are necessary diagnostics, not a continuum proof.
Unlike the initial-time saved problem, shear normals and cone transforms are
recomputed from the evolved velocity. Columns are ordered slopes, pressures.
"""
import json
from pathlib import Path
import time
import numpy as np

from grouped_outer_cache import _group_indices_by_z, _group_quadrature
from meridional_state_cache import StateCache, StatefulMean, value_and_pressure_modes
from affine_momentum import jets, momentum
from broad_shear_dynamic_control import load_saved_field
from grouped_joined_field import install_in_field


class StateConstraints:
    def __init__(self, base, values, pressures, k, locations, breaks, order=24):
        self.k = float(k)
        self.locations = list(locations)
        all_locations = self.locations + [(-.2, 1.), (.2, 1.)]
        tau = .5 * 2.**(-k)
        self.centers = [base.inner.from_similarity(
            [base.inner.p.X_max * (1 + 15*y)**2], [eta], tau)[0]
            for eta, y in all_locations]
        points = list(self.centers)
        panels = {}
        for group in _group_indices_by_z(self.centers):
            panels.update(_group_quadrature(self.centers, group, base.inner,
                                           tau, order, points, breaks))
        self.panels = [panels[i] for i in range(len(self.centers))]
        self.cache = StateCache(base, values, pressures, np.asarray(points), k)

    def assemble(self, state):
        offset, columns = self.cache.derivative_problem(state)
        u = self.cache.base[0] + np.tensordot(state, self.cache.values[0], axes=1)
        grad = self.cache.base[1] + np.tensordot(state, self.cache.values[1], axes=1)
        targets, directions = [], []
        for sl, rr, ww, radius in self.panels:
            targets.append(np.array([-np.dot(ww*rr**2, offset[sl, 1])/radius**2,
                                     -np.dot(ww*rr, offset[sl, 2])/radius]))
            directions.append(np.stack((
                -np.einsum('n,np->p', ww*rr**2, columns[sl, 1])/radius**2,
                -np.einsum('n,np->p', ww*rr, columns[sl, 2])/radius)))
        rows, offsets, diagnostics = [], [], []
        for i, (eta, y) in enumerate(self.locations):
            radius = self.centers[i][0]
            F = u[i, 1] / radius
            shear = np.array([grad[i, 1, 0] - F, grad[i, 2, 0]])
            mag = np.linalg.norm(shear)
            N = shear / max(mag, 1e-30)
            K = np.array([-N[1], N[0]])
            lam = -2*F*N[0]*(2*F*N[0] + mag)
            # -T.N >= .02, and sqrt(lam)*|T.K| <= .95*|2F Ntheta|*(-T.N).
            growth = .95 * abs(2*F*N[0])
            root = np.sqrt(max(lam, 0.))
            H = np.stack((-N, -growth*N-root*K, -growth*N+root*K))
            # Normalize only by current geometry, preserving inequality signs.
            norm = np.maximum(np.linalg.norm(H, axis=1), 1e-30)
            H = H / norm[:, None]
            margin = np.array([.02, 0., 0.]) / norm
            rows.append(H @ directions[i])
            offsets.append(H @ targets[i] - margin)
            diagnostics.append(dict(eta=eta, y=y, lambda_squared=float(lam),
                                    shear_magnitude=float(mag), H=H.tolist()))
        return dict(E=np.concatenate(directions[-2:], axis=0),
                    m=np.concatenate(targets[-2:]),
                    A=np.concatenate(rows, axis=0), b=np.concatenate(offsets),
                    geometry_admissible=all(d['lambda_squared'] > 0 and
                                            d['shear_magnitude'] > 0 for d in diagnostics),
                    diagnostics=diagnostics)


def run():
    root = Path(__file__).resolve().parent
    saved = json.loads((root/'broad_meridional_constrained.json').read_text())
    replay = json.loads((root/'broad_meridional_constrained_replay.json').read_text())
    dynamic, source = load_saved_field()
    install_in_field(dynamic)
    base, values, pressures, _, _ = value_and_pressure_modes(dynamic, source)
    locations = [(d['label']['eta'], d['label']['y'])
                 for d in saved['cone_problem']['diagnostics']]
    started = time.perf_counter()
    engine = StateConstraints(base, values, pressures, source['k'], locations,
                              saved['quadrature']['radial_split_breaks'], order=24)
    initial = engine.assemble(np.zeros(len(values)))
    # Original modes are 25 velocity slopes followed by 19 pressure columns.
    control = np.asarray(replay['coefficients'])
    equality = initial['E'] @ control + initial['m']
    cones = initial['A'] @ control + initial['b']
    shifted = np.zeros(len(values)); shifted[0] = .01; shifted[15] = .01
    changed = engine.assemble(shifted)
    # Independent actual-field replay of the first radial panel at nonzero state.
    sl, rr, ww, radius = engine.panels[0]
    field = StatefulMean(base, values, pressures, shifted, control[:len(values)],
                         control[len(values):], source['k'])
    cache = engine.cache
    direct = momentum(jets(field, cache.points[sl], cache.tau,
                           cache.hspace, cache.htime))
    predicted = cache.residual(shifted, control[:len(values)],
                               control[len(values):])[sl]
    direct_target = np.array([-np.dot(ww*rr**2, direct[:, 1])/radius**2,
                              -np.dot(ww*rr, direct[:, 2])/radius])
    H = np.asarray(changed['diagnostics'][0]['H'])
    # First H row has unit norm, so its physical sign margin is still .02.
    direct_margin = H @ direct_target - np.array([.02, 0., 0.])
    predicted_margin = (changed['A'] @ control + changed['b'])[:3]
    result = dict(accepted=False, pde_validated=False, scale_recursion_established=False,
        scope='Finite state-dependent moment and cone assembly; order 24 diagnostic only. No integrated constrained trajectory.',
        order=24, elapsed_seconds=time.perf_counter()-started,
        initial_geometry_admissible=initial['geometry_admissible'],
        initial_moment_max=float(np.max(abs(equality))),
        initial_minimum_cone_margin=float(np.min(cones)),
        shifted_geometry_admissible=changed['geometry_admissible'],
        changed_cone_matrix_max=float(np.max(abs(changed['A']-initial['A']))),
        changed_moment_offset_max=float(np.max(abs(changed['m']-initial['m']))),
        shifted_first_panel_scaled_residual_error=float(np.max(
            np.linalg.norm(direct-predicted, axis=1) /
            np.maximum(np.linalg.norm(direct, axis=1), 1.))),
        shifted_first_cone_direct_margin=direct_margin.tolist(),
        shifted_first_cone_predicted_margin=predicted_margin.tolist(),
        shifted_first_cone_margin_error_max=float(np.max(abs(direct_margin-predicted_margin))),
        initial_diagnostics=initial['diagnostics'], shifted_diagnostics=changed['diagnostics'])
    root.joinpath('meridional_state_constraints.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if 'diagnostics' not in k}), flush=True)


if __name__ == '__main__':
    run()
