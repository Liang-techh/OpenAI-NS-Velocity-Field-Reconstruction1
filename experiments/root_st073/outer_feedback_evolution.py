"""State-dependent short swirl trajectory; sampled gates, never recursion acceptance.

Physical time is -tau, and k=-log2(2*tau), so a'(k) contributes
+U*a'/(tau*log(2)) to momentum. Constant swirl values and time slopes
are separate: rebuilding at every stage includes the changed advection.
"""
import json
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline
from scipy.optimize import linprog, minimize
from adaptive_bridge_recursive_defect import build_fields
from adaptive_bridge_moment_fit import moment_slices
from affine_momentum import momentum, jets
from grouped_outer_cache import outer_cache
from midplane_resolved_feasibility import ZeroBackground, _integrated_moment_coefficients
from midplane_outer_slope_patch_repair import P_WINDOWS, P_BREAKS
from outer_pressure_modes import OuterPressure, WINDOWS
from outer_swirl_slope import OuterSwirlSlope
from separated_moment_modes import SeparatedMomentModes, flat_bump
from joined_field import coordinates
from midplane_integrated_moment_balance import evaluate
from midplane_outer_residual_source import outer_cones
from radial_continuation import ROOT


class SwirlValue(OuterSwirlSlope):
    """Compact axisymmetric azimuthal velocity; structurally divergence-free."""
    def __init__(self, base, coefficients):
        super().__init__(base, coefficients, 0., P_WINDOWS)

    def fields(self, points, tau):
        points = np.asarray(points)
        u, p = self.base.fields(points, tau)
        r = np.hypot(points[:, 0], points[:, 1])
        safe = np.where(r > 0, r, 1.)
        co = coordinates(r/np.sqrt(self.nu), points[:, 2]/np.sqrt(self.nu), tau, self.inner.h)
        q, eta = np.asarray(co['q']), np.asarray(co['eta'])
        ri = np.sqrt(2*self.nu*q*self.join_X)
        y = (r-ri)/((self.ratio-1)*ri)
        swirl = np.zeros(len(points))
        for j, (lo, hi) in enumerate(self.windows):
            bump, _ = flat_bump(y, lo, hi)
            for power in range(3):
                swirl += self.a[j, power]*np.sqrt(self.nu)*q**(-self.inner.A)*bump*(eta/.3)**power
        u = u.copy()
        u[:, 0] -= swirl*points[:, 1]/safe
        u[:, 1] += swirl*points[:, 0]/safe
        return u, p


def build_current():
    inner, fields = build_fields()
    base = fields['two_sided_cone']
    seed = json.loads((ROOT/'midplane_outer_pressure_staged_k11.json').read_text(encoding='utf-8'))
    velocity = SeparatedMomentModes(base, seed['amplitudes'], windows=WINDOWS,
                                   knots=(11., 15., 19.), axial_powers=(0, 1, 2))
    return inner, base, OuterPressure(velocity, seed['pressure_coefficients'])


def solve_control(inner, base, current, k, state, reference, order=48):
    started = time.perf_counter()
    field = SwirlValue(current, state)
    zero = ZeroBackground(base)
    units = [OuterPressure(zero, e, windows=P_WINDOWS) for e in np.eye(9)]
    units += [OuterSwirlSlope(zero, e, k, P_WINDOWS) for e in np.eye(9)]
    data = moment_slices(inner, base, field, orders=(k,), n=order, unit_fields=units,
                         unit_fields_are_deltas=True, radial_breaks=P_BREAKS)[0]
    m, E, Q = _integrated_moment_coefficients(data)
    # The instantaneous velocity/gradient of every control is exactly zero.
    assert np.max(abs(Q)) < 1e-12
    locations = list(dict.fromkeys(
        [(eta, y) for eta in (-.2, 0., .2) for y in (.5, .75)] +
        [(eta+de, .75+dy) for eta in (-.2, .2)
         for de in (-.005, 0., .005) for dy in (-.005, 0., .005)]))
    cache = outer_cache(field, units, k, order=order, locations=locations)
    residual = momentum(cache['baseline'])
    A, b, lambdas = [], [], []
    for i, (sl, r, w, R) in enumerate(cache['panels']):
        u, J = cache['baseline'][0][i], cache['baseline'][1][i]
        F = u[1]/R
        s = np.array([J[1, 0]-F, J[2, 0]])
        N = s/np.linalg.norm(s)
        K = np.array([-N[1], N[0]])
        lam = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(s))
        lambdas.append(float(lam))
        if lam <= 0:
            raise RuntimeError(f'Lost positive cone geometry at k={k}, node={i}: {lam}')
        L, G = np.sqrt(lam), .95*abs(2*F*N[0])
        T = np.array([-np.dot(w*r*r, residual[sl, 1])/R**2,
                      -np.dot(w*r, residual[sl, 2])/R])
        P = np.stack([-np.einsum('n,pn->p', w*r*r/R**2, cache['modes'][2][:, sl, 1]),
                      -np.einsum('n,pn->p', w*r/R, cache['modes'][2][:, sl, 2])])
        scale = max(abs(T@N), 1.)
        H = np.stack([-N/scale, (-G*N-L*K)/(G*scale), (-G*N+L*K)/(G*scale)])
        A.extend(H@P)
        b.extend(H@T-np.array([.02, 0., 0.]))
    A, b = np.asarray(A), np.asarray(b)
    es = np.maximum(np.max(abs(E), axis=1), 1e-12)
    En, mn = E/es[:, None], m/es
    lp = linprog(np.zeros(18), A_ub=-A, b_ub=b, A_eq=En, b_eq=-mn,
                 bounds=[(None, None)]*18, method='highs')
    if not lp.success:
        raise RuntimeError(f'Feedback infeasible at k={k}: {lp.message}')
    # Fixed reference, not previous solver output: deterministic feedback law.
    fit = minimize(lambda x: (float((x-reference)@(x-reference)), 2*(x-reference)),
                   lp.x, jac=True, method='SLSQP',
                   constraints=[dict(type='eq', fun=lambda x: En@x+mn, jac=lambda x: En),
                                dict(type='ineq', fun=lambda x: A@x+b, jac=lambda x: A)],
                   options=dict(maxiter=250, ftol=1e-10))
    x = fit.x
    if not fit.success or max(abs(En@x+mn)) > 1e-7 or min(A@x+b) < -1e-7:
        raise RuntimeError(f'Feedback optimization failed at k={k}: {fit.message}')
    row = dict(k=float(k), state=state.tolist(), control=x.tolist(),
               moment_max=float(max(abs(E@x+m))), minimum_constraint=float(min(A@x+b)),
               minimum_lambda_squared=min(lambdas), seconds=time.perf_counter()-started)
    print(json.dumps({key: value for key, value in row.items() if key not in ('state', 'control')}), flush=True)
    return x, row


class Trajectory:
    def __init__(self, current, nodes):
        self.base = current
        for name in ('inner', 'nu', 'join_X', 'ratio'):
            setattr(self, name, getattr(current, name))
        k = np.array([row['k'] for row in nodes])
        self.interval = (float(k[0]), float(k[-1]))
        self.state = CubicHermiteSpline(k, [row['state'] for row in nodes],
                                       [row['control'][9:] for row in nodes], extrapolate=False)
        self.pressure = CubicSpline(k, [row['control'][:9] for row in nodes], extrapolate=False)

    def fields(self, points, tau):
        k = -np.log2(2*float(np.asarray(tau).ravel()[0]))
        if not self.interval[0] <= k <= self.interval[1]:
            raise ValueError('Trajectory evaluation outside saved interval')
        return OuterPressure(SwirlValue(self.base, self.state(k)), self.pressure(k),
                             windows=P_WINDOWS).fields(points, tau)


def run():
    inner, base, current = build_current()
    reference = np.array(json.loads((ROOT/'midplane_outer_slope_patch_repair.json').read_text(encoding='utf-8'))['coefficients'])
    report = dict(method='Explicit midpoint, state-dependent constrained pressure and swirl slope; Hermite state replay',
                  accepted=False, scale_recursion_established=False, nodes=[], stages=[], replay=[],
                  scope='Short interval sampled diagnostic; no continuum, energy, full residual or recursive acceptance')
    path = ROOT/'outer_feedback_evolution.json'
    def save():
        path.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    state = np.zeros(9)
    k, step = 11., .0005
    try:
        control, row = solve_control(inner, base, current, k, state, reference)
        report['nodes'].append(row); save()
        for _ in range(2):
            mid, stage = solve_control(inner, base, current, k+step/2, state+step/2*control[9:], reference)
            report['stages'].append(stage); save()
            state = state+step*mid[9:]
            k += step
            control, row = solve_control(inner, base, current, k, state, reference)
            report['nodes'].append(row); save()
        field = Trajectory(current, report['nodes'])
        # Interior holdouts keep the 4th-order time stencil strictly in range.
        for test_k in (11.00035, 11.0005, 11.00065):
            moments = evaluate(field, inner, test_k, 96, .002, radial_breaks=P_BREAKS)
            cone = outer_cones(field, test_k, order=64, radial_breaks=P_BREAKS)
            tau = .5*2.**-test_k
            pts = np.concatenate([inner.from_similarity(inner.p.X_max*(1+15*np.linspace(.01,.99,61))**2,
                                                       np.full(61, eta), tau) for eta in (-.3,-.1,0.,.1,.3)])
            residual = momentum(jets(field, pts, tau, .0005*np.sqrt(inner.nu*tau), .0001*tau))
            row = dict(k=test_k, moments=moments.tolist(), moment_max=float(max(abs(moments))),
                       cones=cone, cone_pass_count=sum(r['cone_pass'] for r in cone),
                       sampled_momentum_peak=float(max(np.linalg.norm(residual, axis=1))))
            report['replay'].append(row); save()
            print(json.dumps({key: value for key, value in row.items() if key not in ('cones','moments')}), flush=True)
    except RuntimeError as error:
        report['failure'] = str(error); save()
        raise


if __name__ == '__main__':
    run()
