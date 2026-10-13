"""Direct Gamma integrals, full tails, native import and ODE history checks.

Finite numerical fixtures use ordinary nonzero parameters. They do not
evaluate the astronomical original radius, or prove final moment closure.
"""
import json
from pathlib import Path
import time
import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import expit, exprel, gamma, roots_genlaguerre
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions as current


def flat_sigma(x):
    x = np.asarray(x)
    s = np.clip(x, 1e-8, 1-1e-8)
    return np.where(x <= 0, 0., np.where(x >= 1, 1., expit(-1/s**2+1/(1-s)**2)))


def flat_phi(x):
    x = np.asarray(x)
    distance = np.where(x < 3, 3-x, 1.)
    return np.where(x < 3, np.exp(-4/distance**2), 0.)


class FixtureEvaluator:
    """Evaluate actual DAG operations; no original field/cap callbacks."""
    def __init__(self, g, delta, Gamma_order=96):
        self.g, self.delta, self.order = g, delta, Gamma_order
        self.variables = {row['name']: i for i, row in enumerate(g.nodes)
            if row['operation'] == 'bound_variable'}
        self.Gamma_rules = {}

    def __call__(self, ref, env=None):
        return float(self.value(ref.node, env or {}))

    def value(self, index, env):
        cache = {}
        def at(i):
            if i in env:
                return env[i]
            if i in cache:
                return cache[i]
            row = self.g.nodes[i]
            op = row['operation']
            if op == 'exact_rational':
                result = row['numerator']/row['denominator']
            elif op == 'sum':
                result = sum(at(j) for j in row['arguments'])
            elif op == 'product':
                result = 1.
                for j in row['arguments']:
                    result = result*at(j)
            elif op == 'negative':
                result = -at(row['argument'])
            elif op == 'positive_quotient':
                result = at(row['numerator'])/at(row['denominator'])
            elif op == 'analytic_unary':
                x = at(row['argument'])
                name = row['name']
                if name == 'exprel':
                    result = exprel(x)
                elif name == 'log1prel':
                    xx = np.asarray(x)
                    result = np.divide(np.log1p(xx), xx, out=np.ones_like(xx, dtype=float), where=xx != 0)
                elif name == 'gamma':
                    result = gamma(x)
                else:
                    result = getattr(np, name)(x)
            elif op == 'exact_original_source_flat_sigma':
                result = flat_sigma(at(row['argument']))
            elif op == 'exact_original_source_heat_phi':
                result = flat_phi(at(row['argument']))
            elif op == 'mathematical_positive_infinity':
                result = np.inf
            elif op == 'function_substitution':
                result = self.value(row['expression'], {**env, row['variable']: at(row['value'])})
            elif op == 'definite_integral':
                variable = self.variables[row['variable']]
                lo, hi = at(row['lower']), at(row['upper'])
                body = lambda x: self.value(row['integrand'], {**env, variable: x})
                if row.get('exact_positive_Gamma_expectation'):
                    alpha = float(at(self.g.nodes[row['Gamma_normalization']]['argument']))
                    if alpha not in self.Gamma_rules:
                        self.Gamma_rules[alpha] = roots_genlaguerre(self.order, alpha)
                    points, weights = self.Gamma_rules[alpha]
                    result = np.dot(weights, body(points)/(np.exp(-points)*points**alpha))
                elif np.isinf(hi):
                    # An exact substitution on the ENTIRE [lo,infinity).
                    rate = self.delta/2 if row['full_exact_Gamma_future_defect'] == 'theta' else 1.
                    result = quad(lambda u: body(lo-np.log(u)/rate)/(rate*u),
                        0, 1, epsabs=2e-10, epsrel=2e-10, limit=180)[0]
                else:
                    points = [lo, 1., hi] if lo < 1 < hi else [lo, hi]
                    result = sum(quad(body, left, right, epsabs=2e-10,
                        epsrel=2e-10, limit=120)[0] for left, right in zip(points, points[1:]))
            else:
                raise ValueError('No mathematical fixture oracle for ' + op)
            cache[i] = result
            return result
        return at(index)


def exact_identities():
    x, a, S, D, Dz, pre = sy.symbols('x a S D Dz pre', real=True)
    assert sy.simplify(1-sy.exp(-a*sy.log(1+x))-
        a*sy.log(1+x)*(sy.exp(-a*sy.log(1+x))-1)/(-a*sy.log(1+x))) == 0
    square = D*(2*pre-a*S*D)
    assert sy.simplify(sy.diff(square, D)*Dz-2*(pre-a*S*D)*Dz) == 0
    return dict(passed=True,exact_divided_Gamma_deficit_identity=True,
        exact_square_derivative_uses_two_actual_brackets=True,
        endpoint_continuations_log1prel_and_exprel_are_one=True)


def Gamma_and_tail_fixture():
    delta, eps, S = .3, .02, .12
    aa = delta/2
    g = current.source.FunctionTransportGraph()
    functions, pieces = current.Gamma_future_functions(g, g.constant('3/10'),
        g.constant('1/50'), g.constant('3/25'))
    evaluate = FixtureEvaluator(g, delta)
    z = g.symbol('Z')
    t = pieces['offset_variable']
    # Independent direct expectation of 1-(1+xi*v)^(-a), in the
    # original v^a measure. It does not evaluate the DAG's divided body.
    v, w = roots_genlaguerre(128, aa)
    def direct_D(Z, offset):
        xx = 2*(1-Z*Z)*S*np.exp(-offset)*v
        return np.dot(w, -np.expm1(-aa*np.log1p(xx)))/(aa*S*gamma(1+aa))
    def direct_D_mp(Z, offset):
        with mp.workdps(35):
            Z, offset, s, a = map(mp.mpf, (str(Z), str(offset), '.12', '.15'))
            xi = 2*(1-Z*Z)*s*mp.exp(-offset)
            return float(mp.quad(lambda v: mp.exp(-v)*v**a*
                (-mp.expm1(-a*mp.log1p(xi*v))), [0, 1, mp.inf])/(a*s*mp.gamma(1+a)))
    errors = []
    for Z, offset in ((-.65, .7), (0., 0.), (.65, 4.), (1., .7)):
        env = {z.node: Z, t.node: offset}
        value = evaluate(functions['Dhat'].value, env)
        want = direct_D_mp(Z, offset)
        errors.append(abs(value-want))
        h = 1e-4
        derivative = ((3*want-4*direct_D_mp(Z-h, offset)+direct_D_mp(Z-2*h, offset))/(2*h)
            if Z == 1. else (direct_D_mp(Z+h, offset)-direct_D_mp(Z-h, offset))/(2*h))
        assert abs(evaluate(functions['Dhat'].Z, env)-derivative) < 4e-8
        assert abs(value-want) < 2e-10
    def direct_piece(name, Z, lo, hi):
        def body(offset):
            sig, ph = float(flat_sigma(offset)), float(flat_phi(offset))
            C, pre = 1-eps*ph, 1-eps*(1-sig+sig*ph)
            D = sig*C*direct_D(Z, offset)
            if name == 'theta':
                return aa*np.exp((1-aa)*offset)*D
            square = D*(2*pre-aa*S*D)
            return (np.exp(-(1+delta)*offset)/2 if name == 'pressure'
                    else np.exp(-delta*offset))*square
        if np.isinf(hi):
            rate = aa if name == 'theta' else 1.
            return quad(lambda u: body(lo-np.log(u)/rate)/(rate*u), 0, 1,
                epsabs=1e-10, epsrel=1e-10, limit=180)[0]
        return quad(body, lo, 1, epsabs=1e-10, epsrel=1e-10)[0] + \
            quad(body, 1, hi, epsabs=1e-10, epsrel=1e-10)[0]
    values = {}
    for name in ('theta', 'pressure', 'energy'):
        for Z in (.65, 1.):
            env = {z.node: Z}
            collar = evaluate(pieces['collar'][name].value, env)
            exterior = evaluate(pieces['exterior'][name].value, env)
            expected = direct_piece(name, Z, 0, 3)+direct_piece(name, Z, 3, np.inf)
            value = evaluate(functions[name].value, env)
            assert abs(value-collar-exterior) < 2e-10
            assert abs(value-expected) < 3e-9
            h = 1e-4
            full = lambda zz: direct_piece(name, zz, 0, 3)+direct_piece(name, zz, 3, np.inf)
            expected_Z = ((3*expected-4*full(Z-h)+full(Z-2*h))/(2*h) if Z == 1.
                else (full(Z+h)-full(Z-h))/(2*h))
            assert abs(evaluate(functions[name].Z, env)-expected_Z) < 2e-7
            if Z == 1.:
                assert value == 0 and evaluate(functions[name].Z, env) < 0
            else:
                assert value > 0 and exterior > 0 and evaluate(functions[name].Z, env) < 0
                values[name] = dict(value=value,collar=collar,entire_exterior=exterior)
        assert abs(evaluate(functions[name].Z, {z.node: 0.})) == 0
    return dict(passed=True,mathematical_fixture_only=True,fixture_delta=delta,
        fixture_nonzero_S=S,fixture_epsilon=eps,original_parameters_not_replaced=True,
        direct_original_Gamma_expectation_max_error=max(errors),
        all_three_complete_value_and_first_Z_functions_checked=True,
        exact_end_values_zero_but_endpoint_derivatives_nonzero=True,
        full_infinite_exterior_evaluated_without_finite_cutoff=True,values=values)


def history_fixture():
    mu, delta, eps, logRp = .6, .3, .0003, 1.4
    r, k = 1-mu, 1-delta/2
    g = current.source.FunctionTransportGraph()
    z = g.symbol('Z')
    Xp = current.source.C1Function(g.add(g.constant(3),g.mul(g.constant('1/7'),z),
        g.mul(g.constant('1/9'),z,z)),g.add(g.constant('1/7'),g.mul(g.constant('2/9'),z)))
    schedule, histories = current.angular_history_functions(g, g.constant('3/5'),
        g.constant('3/10'), g.constant('3/10000'), Xp, g.constant('7/5'))
    evaluate = FixtureEvaluator(g,delta)
    def ode_endpoint(Z):
        X0, D0 = 3+Z/7+Z*Z/9, 1/7+2*Z/9
        memory = np.exp(-13*r/mu)
        Xv, Dv = 1/r+(X0-1/r)*memory, D0*memory
        def flatten(t,y):
            s = t/100
            sig = float(flat_sigma(s))
            ds = sig*(1-sig)*(2/s**3+2/(1-s)**3)/100 if 0<s<1 else 0.
            coefficient = r+np.log((1+Z*Z)/2)*ds
            return [1-coefficient*y[0],-coefficient*y[1]-2*Z/(1+Z*Z)*ds*y[0]]
        Xf, Df = solve_ivp(flatten,[0,100],[Xv,Dv],rtol=2e-12,atol=2e-13,max_step=.5).y[:,-1]
        Lrel=-30*np.log(mu)
        Xrel=1/r+(Xf-1/r)*np.exp(-r*Lrel);Drel=Df*np.exp(-r*Lrel)
        def unit(y,coeff):
            return solve_ivp(lambda t,x:[1-coeff(t)*x[0],-coeff(t)*x[1]],
                [0,1],y,rtol=2e-12,atol=2e-13,max_step=.01).y[:,-1]
        Xin,Din=unit([Xrel,Drel],lambda t:r*(1-float(flat_sigma(t))))
        Ts=4*np.log(2/delta)
        Xt,Dt=unit([Xin+Ts,Din],lambda t:k*float(flat_sigma(t)))
        return dict(Xv=(Xv,Dv),Xf=(Xf,Df),Xrel=(Xrel,Drel),Xt=(Xt,Dt)),Ts
    for Z in (-.4,.5):
        expected,_=ode_endpoint(Z)
        for name,want in expected.items():
            pair=histories[name]
            assert abs(evaluate(pair.value,{z.node:Z})-want[0])<2e-9,name
            assert abs(evaluate(pair.Z,{z.node:Z})-want[1])<2e-9,name
    expected,Ts=ode_endpoint(0.)
    J=quad(lambda t:np.exp(k*t)*(1-float(flat_sigma(t))+float(flat_sigma(t))*float(flat_phi(t))),
        0,3,points=[1],epsabs=1e-12)[0]
    wait=(np.log(expected['Xt'][0]-1/k)+np.log1p(-eps)-np.log(eps)-np.log(1/k+J))/k
    assert wait>0
    assert abs(evaluate(schedule['waiting'])-wait)<2e-9
    want_logRtail=logRp+13/mu+102-30*np.log(mu)+Ts+wait
    assert abs(evaluate(schedule['logRtail'])-want_logRtail)<2e-9
    assert evaluate(schedule['S'])>0
    assert abs(np.log(evaluate(schedule['S']))+want_logRtail)<2e-9
    return dict(passed=True,mathematical_fixture_only=True,nonconstant_native_X_fixture=True,
        original_native_X_import_checked_separately=True,
        nonzero_pulse_memory_and_first_Z_rows_checked_against_independent_ODE=True,
        flatten_steep_waiting_and_actual_Rtail_composition_checked=True,
        original_waiting_scalar_uses_Xt_at_Z0=True)


def run(owner=None):
    began=time.monotonic()
    owner=owner or current.CurrentOriginalN1GammaFutureFunctions(require_checked=False)
    owner.assert_graph()
    raw=current.parent.linear.read(current.NAME)
    assert raw['exact_graph_nodes']==owner.g.nodes
    assert raw['function_bindings']==owner.snapshot
    assert not any(raw[key] for key in current.OPEN)
    for pair in owner.functions.values():
        assert type(pair) is current.source.C1Function
        owner.g.ids([pair.value,pair.Z])
    frame=current.parent.linear.read(current.parent.quiet.FRAME_NAME)
    original=current.parent.linear.realization.restore_graph(frame['exact_graph_nodes'])
    importer=current.parent.GraphImport(owner.g,original)
    native=frame['actual_Rp_native_frame_functions']
    expected=[importer(current.source.FunctionRef(original,i)).node
        for i in native['actual_identified_Rp_native_frame_C3']['X'][:2]]
    assert expected==current.parent.encode(owner.histories['Xp'])
    assert importer(current.source.FunctionRef(original,native['original_logRp'])).node==owner.schedule['logRp'].node
    assert owner.g.nodes[owner.schedule['S'].node]['operation']=='analytic_unary'
    assert owner.g.nodes[owner.schedule['S'].node]['name']=='exp'
    assert owner.g.nodes[owner.pieces['infinity'].node]['operation']=='mathematical_positive_infinity'
    fixtures=dict(exact_identities=exact_identities(),Gamma_and_infinite_tails=Gamma_and_tail_fixture(),
        native_history_and_waiting=history_fixture(),source_native_import=dict(passed=True,
            actual_X_C1_and_logRp_bodies_relocated_from_checked_frame=True,
            all_ancestor_graph_prefixes_unchanged=True))
    owner.assert_graph()
    for name,digest in owner.hashes.items():
        assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=fixtures,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_EXACT_GAMMA_FUTURE_C1_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':
    run()
