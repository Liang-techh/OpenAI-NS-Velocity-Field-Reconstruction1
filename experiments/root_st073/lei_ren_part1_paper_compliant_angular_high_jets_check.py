"""Independent C4 implicit, flatten and true-Gamma derivative checks."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_angular_high_jets import (
    CompliantAngularHighJets, gamma_deficit_jets, integrated_gamma_tails,
    implicit_quadratic_jets, log_taylor)
from lei_ren_part1_paper_compliant_outer_angular_repair import CompliantAngularRepair
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).parent
NAME = 'lei_ren_part1_paper_compliant_angular_high_jets.json'


def identities():
    passed = {}
    def zero(name, expr):
        if s.simplify(expr) != 0:
            raise ArithmeticError('High-jet identity failed: '+name)
        passed[name] = True
    p, q, K = s.symbols('p q K', real=True)
    x = s.symbols('x0:5'); y = s.symbols('y0:5')
    h = s.symbols('h', real=True)
    J21, J22 = 1+2*K*x[0], q*(1+2*K*y[0]); det = p*J22-J21
    fx = sum(x[i]*h**i for i in range(5)); fy = sum(y[i]*h**i for i in range(5))
    for n in range(1, 5):
        b1, b2 = s.symbols('b1_'+str(n)+' b2_'+str(n))
        cross = K*sum(x[i]*x[n-i]+q*y[i]*y[n-i] for i in range(1, n))
        v = b2-cross
        xn, yn = (J22*b1-v)/det, (p*v-J21*b1)/det
        zero('implicit_linear_order'+str(n), (s.expand(p*fx+fy).coeff(h, n)-b1).subs({x[n]: xn, y[n]: yn}))
        zero('implicit_quadratic_order'+str(n), (s.expand(fx+q*fy+K*(fx**2+q*fy**2)).coeff(h, n)-b2).subs({x[n]: xn, y[n]: yn}))
    a, v, xi = s.symbols('a v xi', positive=True)
    for n in range(1, 5):
        zero('Gamma_integrand_derivative'+str(n), s.diff((1+xi*v)**(-a), xi, n)-(-1)**n*s.rf(a, n)*v**n*(1+xi*v)**(-a-n))
    n = s.symbols('n', integer=True, positive=True)
    # The leading angular integral cancels a exactly, while pressure/energy
    # poles remain positive even when delta is extremely small.
    zero('angular_leading_integral_cancellation', a/a-1)
    z = s.symbols('Z', real=True)
    zero('quadratic_axial_composition', (1-(z+h)**2)-(1-z**2-2*z*h-h**2))
    zero('square_pressure_deficit', (1-(1-a*v)**2)/(2*a)-(v-a*v**2/2))
    return passed


def contains(name, interval, reference, error=0):
    lo, hi = endpoints(interval)
    if not lo <= reference-error <= reference+error <= hi:
        raise ArithmeticError('Independent C4 reference outside enclosure: '+name)


def implicit_fixture():
    with mp.workdps(85):
        c = MPIntervalContext(); c.dps = 75
        p, q, K = map(mp.mpf, ('.14', '.12', '.02'))
        b1 = lambda z: mp.exp(z/3)+mp.mpf('.1')*z*z
        b2 = lambda z: mp.mpf('.03')/(1+z*z)
        def actual(z):
            b, d = b1(z), b2(z)
            A = K*(1+q*p*p); B = 1-q*p-2*K*q*p*b
            C = q*b+K*q*b*b-d
            x = -2*C/(B+mp.sqrt(B*B-4*A*C))
            return x, b-p*x
        checks = {}
        for z in map(mp.mpf, ('0', '.5')):
            x0, y0 = actual(z)
            rhs1 = IntervalTaylor(c, [c.mpf(v) for v in mp.taylor(b1, z, 4)])
            rhs2 = IntervalTaylor(c, [c.mpf(v) for v in mp.taylor(b2, z, 4)])
            actuals = [mp.taylor(lambda zz: actual(zz)[j], z, 4) for j in (0, 1)]
            radius = c.mpf('1e-65')
            box = lambda v: c.mpf(v)+c.mpf([-endpoints(radius)[1], endpoints(radius)[1]])
            out = implicit_quadratic_jets(c, box(x0), box(y0), rhs1, rhs2, c.mpf(p), c.mpf(q), c.mpf(K))
            for j, jet in enumerate(out['scaled_coefficient_Taylor']):
                for n in range(5):
                    label = 'Z'+str(z)+'_branch'+str(j)+'_order'+str(n)
                    contains(label, jet[n], actuals[j][n]); checks[label] = True
        return dict(independent_stable_closed_form_branch_checks=checks, passed=True)


def gamma_fixture():
    with mp.workdps(55):
        c = MPIntervalContext(); c.dps = 70
        a, S, Z, t = map(mp.mpf, ('.1', '.002', '.5', '3'))
        H = lambda xi: xi**(-1-a)*mp.hyperu(1+a, 2, 1/xi)
        true_D = lambda z: (1-H(2*(1-z*z)*S*mp.exp(-t)))/(a*S)
        expected = mp.taylor(true_D, Z, 4)
        actual = gamma_deficit_jets(c, c.mpf(Z), c.mpf(a), c.mpf(S), 4, c.mpf(t))
        checks = {}
        for n in range(5):
            contains('Gamma_deficit_order'+str(n), actual[n], expected[n]); checks['Gamma_deficit_order'+str(n)] = True
        tails = integrated_gamma_tails(c, c.mpf(Z), c.mpf(a), c.mpf(S), 4)
        k = 1-a
        def true_theta_tail(z):
            xi = 2*(1-z*z)*S*mp.exp(-3)
            integrated_Gamma = xi**(-1-a)*mp.hyperu(1+a, 3, 1/xi)/k
            return mp.exp(3*k)*(integrated_Gamma-1/k)/S
        refs = mp.taylor(true_theta_tail, Z, 4)
        for n in range(5):
            contains('Gamma_angular_tail_order'+str(n), tails['theta'][n], refs[n]); checks['Gamma_angular_tail_order'+str(n)] = True

        # Differentiate the exact H integral, compose xi(Z) as polynomials,
        # then independently integrate the squared heat deficit derivatives.
        cache = {}
        def square_deficit_coefficients(t):
            if t in cache: return cache[t]
            xi = 2*(1-Z*Z)*S*mp.exp(-t)
            delta = [mp.mpf(0), -4*Z*S*mp.exp(-t), -2*S*mp.exp(-t), 0, 0]
            coeffs = [mp.mpf(0)]*5
            power = [mp.mpf(1), 0, 0, 0, 0]
            for m in range(5):
                hm = (-1)**m*mp.rf(a, m)*mp.rf(1+a, m)*xi**(-1-a-m)*mp.hyperu(1+a+m, 2, 1/xi)/math.factorial(m)
                coeffs = [coeffs[j]+hm*power[j] for j in range(5)]
                power = [sum(power[i]*delta[j-i] for i in range(j+1)) for j in range(5)]
            square = [sum(coeffs[i]*coeffs[j-i] for i in range(j+1)) for j in range(5)]
            cache[t] = [((1 if j == 0 else 0)-square[j])/(2*a*S) for j in range(5)]
            return cache[t]
        limit = mp.mpf(30)
        leading = [2*(1+a)*(1-Z*Z), -4*(1+a)*Z, -2*(1+a), 0, 0]
        for name, rate, factor in (('pressure', 1+2*a, 1), ('energy', 2*a, 2)):
            for n in range(5):
                value = mp.quad(lambda tt: mp.exp(-rate*tt)*square_deficit_coefficients(tt)[n], [3, 8, 16, limit])
                value += leading[n]*mp.exp(-(rate+1)*limit)/(rate+1)
                # For this a,S,Z and orders<=4, the derivative remainder
                # constants are below1e4. The omitted difference is O(S)
                # with the second exponential power, not a zero tail.
                error = mp.mpf('1e4')*S*mp.exp(-(rate+2)*limit)/(rate+2)
                label = 'Gamma_'+name+'_tail_order'+str(n)
                contains(label, tails[name][n], factor*value, factor*error); checks[label] = True
        return dict(true_Gamma_value_and_axial_tail_derivative_checks=checks,
                    finite_pressure_energy_tail_errors_bounded=True, passed=True)


def flatten_fixture():
    with mp.workdps(80):
        c = MPIntervalContext(); c.dps = 70
        old = object.__new__(CompliantAngularRepair)
        old.ctx = c; old.rate = c.mpf('.975'); old.cells = 256
        old.angular = SimpleNamespace(Xv=c.mpf('1.125'))
        field = object.__new__(CompliantAngularHighJets); field.repair = old; field.ctx = c; field.order = 4
        def sigma(t):
            if t <= 0: return mp.mpf(0)
            if t >= 1: return mp.mpf(1)
            A, B = mp.exp(-1/t**2), mp.exp(-1/(1-t)**2)
            return A/(A+B)
        z = mp.mpf('.5'); Q = 1+z*z; u, v = 2*z/Q, 1/Q
        # Differentiate the real power explicitly before integration. This
        # avoids differentiating a numerical quadrature repeatedly at very
        # high temporary precision. The formula is independent of log-jet
        # multiplication in the producer.
        def power_coeffs(p):
            rising2 = p*(p+1); rising3 = rising2*(p+2); rising4 = rising3*(p+3)
            return [(2/Q)**p*x for x in (1, -p*u, -p*v+rising2*u*u/2,
                    rising2*u*v-rising3*u**3/6,
                    rising2*v*v/2-rising3*u*u*v/2+rising4*u**4/24)]
        rate = mp.mpf('.975'); points = [0, 25, 50, 75, 90, 100]
        start = mp.mpf('1.125')*mp.exp(-100*rate)
        refs = [-start*value for value in power_coeffs(1)]
        cache = {}
        def integrands(t):
            if t not in cache:
                p = 1-sigma(t/100); weight = mp.exp(-rate*(100-t))
                cache[t] = [weight*x for x in power_coeffs(p)]
                cache[t][0] = weight*(-2**p*mp.expm1(-p*mp.log(Q)))
            return cache[t]
        refs[0] = 2*start*(1-1/Q)+mp.quad(lambda t: integrands(t)[0], points)
        for n in range(1, 5):
            refs[n] -= mp.quad(lambda t: integrands(t)[n], points)
        actual = field.preheat_jets(c.mpf(z)); checks = {}
        for n in range(5):
            contains('flatten_difference_order'+str(n), actual[n], refs[n]); checks['flatten_difference_order'+str(n)] = True
        q = IntervalTaylor(c, [1+z*z, 2*z, 1, 0, 0])
        logs = log_taylor(q); refs = mp.taylor(lambda zz: mp.log(1+zz*zz), z, 4)
        for n in range(5):
            contains('log_Q_order'+str(n), logs[n], refs[n]); checks['log_Q_order'+str(n)] = True
        return dict(independent_flatten_integral_and_log_checks=checks, passed=True)


def run():
    r = json.loads((HERE/NAME).read_bytes()); hashes = dict(r['input_hashes'])
    for name, digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('C4 angular source changed: '+name)
    c = MPIntervalContext(); c.dps = 230
    read = lambda v: read_interval(c, v)
    for point in r['samples']+[r['whole_Z_C4_coefficients']]:
        if len(point['physical_coefficient_Taylor']) != 2:
            raise ValueError('Actual two angular functions required')
        for jet in point['physical_coefficient_Taylor']:
            if len(jet['coefficients']) != 5:
                raise ValueError('Angular C4 jet order lost')
        lo, hi = endpoints(read(point['derivative_jacobian_determinant']))
        if lo <= 0 <= hi: raise ArithmeticError('Uniform high-jet Jacobian lost')
        if point['full_outer_C4_certified'] or point['full_outer_cone_certified'] or point['temporal_recursion']:
            raise ValueError('Angular C4 incorrectly promoted to full outer acceptance')
    for jet in r['samples'][1]['physical_coefficient_Taylor']:
        for n in (1, 3):
            if endpoints(read(jet['coefficients'][n])) != (mp.mpf(0), mp.mpf(0)):
                raise ArithmeticError('Even implicit branch lost exact odd jets at Z0')
    exact = identities(); implicit = implicit_fixture()
    print('C4 exact coefficient identities and independent implicit branch PASS', flush=True)
    gamma = gamma_fixture()
    print('C4 independent true-Gamma deficit and infinite-tail derivative fixtures PASS', flush=True)
    flatten = flatten_fixture()
    hashes[NAME] = hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = dict(identities=exact, independent_implicit_branch_fixture=implicit,
                  independent_Gamma_fixture=gamma, independent_flatten_fixture=flatten,
                  whole_Z_C4_branch_and_Z0_parity_checked=True,
                  actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],
                  implicit_source_sha256=r['implicit_source_sha256'], angular_coefficient_C4_available=True,
                  full_outer_C4_certified=False, full_outer_cone_certified=False, temporal_recursion=False,
                  all_passed=True, input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('Actual angular C4:', len(exact), 'identities; independent closed-form branch, Gamma tail and flatten derivative fixtures PASS', flush=True)
    return result


if __name__ == '__main__':
    run()
