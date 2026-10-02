"""Higher axial jets of the actual implicit two-bump angular repair.

All jets are ordinary Taylor coefficients (derivative divided by factorial).
The same C0 branch and the same exact right sides are retained. Gamma bounds
are finite derivative bounds for the true integral, not a truncated definition.
Positive formal S factors are enclosed by caps, never replaced by their caps.
This layer does not certify the full outer C4 or stress cone.
"""
import hashlib
import json
import math
import operator
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_outer_angular_repair import CompliantAngularRepair, intersect, magnitude
from lei_ren_part1_paper_compliant_outer_initial import stable_sigma
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def log_taylor(jet):
    """Finite log jets from (log f)'=f'/f, retaining every requested order."""
    c = jet.ctx
    if endpoints(jet[0])[0] <= 0:
        raise ValueError('Positive Taylor value required for logarithm')
    reciprocal = jet.reciprocal()
    out = [c.ln(jet[0])]
    for n in range(1, jet.order+1):
        out.append(sum((i*jet[i]*reciprocal[n-i] for i in range(1, n+1)), c.mpf(0))/n)
    return IntervalTaylor(c, out)


def pochhammer(c, value, n):
    result = c.mpf(1)
    for j in range(n):
        result *= value+j
    return result


def gamma_leading_error_weights(c, Z, a, S_cap, order):
    """Dhat_j(t)=A_j exp(-t)+error, |error|<=sum_m W_jm exp(-m t).

    Differentiation in d=1-Z^2 uses the exact positive Gamma expectation.
    W coefficients enclose the actual S powers; no infinite Taylor expansion
    in S is evaluated. Composition uses delta_d=-2Z*h-h^2.
    """
    Z = c.mpf(Z); order = operator.index(order)
    if order < 1 or endpoints(Z)[0] < -1 or endpoints(Z)[1] > 1:
        raise ValueError('Positive order and real Z in [-1,1] required')
    d = 1-Z**2; dmax = c.mpf(endpoints(d)[1]); Smax = c.mpf(endpoints(S_cap)[1])
    if endpoints(a)[0] <= 0 or endpoints(S_cap)[0] < 0:
        raise ValueError('Positive a and nonnegative inverse-radius bound required')
    djet = IntervalTaylor(c, [d, -2*Z]+([-c.mpf(1)] if order >= 2 else [])+[0]*max(0, order-2))
    delta = IntervalTaylor(c, [0]+list(djet.coefficients[1:]))
    leading = djet*(2*(1+a)); B = (1+a)**2*(2+a)
    weights = [{2: c.mpf(0)} for _ in range(order+1)]
    weights[0][2] = 2*B*dmax**2*Smax
    for j in range(1, order+1):
        weights[j][2] += 4*B*dmax*Smax*c.mpf(magnitude(delta[j]))
    for n in range(2, order+1):
        # (a)_n/a=(1+a)_(n-1), avoiding division of two tiny values.
        bound = pochhammer(c, 1+a, n-1)*pochhammer(c, 1+a, n)*2**n*Smax**(n-1)/math.factorial(n)
        power = delta**n
        for j in range(n, order+1):
            weights[j][n] = weights[j].get(n, c.mpf(0))+bound*c.mpf(magnitude(power[j]))
    return leading, weights


def symmetric(c, radius):
    upper = max(mp.mpf(0), endpoints(radius)[1])
    return c.mpf([-upper, upper])


def gamma_deficit_jets(c, Z, a, S_cap, order, t):
    leading, weights = gamma_leading_error_weights(c, Z, a, S_cap, order)
    t = c.mpf(t)
    if endpoints(t)[0] < 0:
        raise ValueError('Nonnegative heat log offset required')
    coeffs = []
    for j in range(order+1):
        radius = sum((w*c.exp(-m*t) for m, w in weights[j].items()), c.mpf(0))
        coeffs.append(leading[j]*c.exp(-t)+symmetric(c, radius))
    # The exact Gamma deficit is nonnegative. This narrows an enclosure,
    # and is not a substitution of a zero or an upper cap for the function.
    lo, hi = endpoints(coeffs[0])
    if hi < 0:
        raise ArithmeticError('Gamma deficit enclosure is negative')
    coeffs[0] = c.mpf([max(mp.mpf(0), lo), hi])
    return IntervalTaylor(c, coeffs)


def integrated_gamma_tails(c, Z, a, S_cap, order, start=3):
    """Analytic infinite-tail derivative enclosures; no finite radial cutoff."""
    A, errors = gamma_leading_error_weights(c, Z, a, S_cap, order)
    start = c.mpf(start); Smax = c.mpf(endpoints(S_cap)[1])
    # Full Dhat absolute coefficient envelopes include their leading power.
    full = [dict(w) for w in errors]
    for j in range(order+1):
        full[j][1] = c.mpf(magnitude(A[j]))
    square = [{} for _ in range(order+1)]
    for j in range(order+1):
        for i in range(j+1):
            for m, wm in full[i].items():
                for n, wn in full[j-i].items():
                    square[j][m+n] = square[j].get(m+n, c.mpf(0))+wm*wn

    angular = []
    for j in range(order+1):
        err = sum((w*a*c.exp(-(m-1+a)*start)/(m-1+a) for m, w in errors[j].items()), c.mpf(0))
        # a/(1-k)=1 is exact; do not divide interval representations of a.
        angular.append(A[j]*c.exp(-a*start)+symmetric(c, err))

    def square_tail(rate):
        out = []
        for j in range(order+1):
            err = sum((w*c.exp(-(rate+m)*start)/(rate+m) for m, w in errors[j].items()), c.mpf(0))
            err += sum((w*(a*Smax/2)*c.exp(-(rate+m)*start)/(rate+m) for m, w in square[j].items()), c.mpf(0))
            out.append(A[j]*c.exp(-(rate+1)*start)/(rate+1)+symmetric(c, err))
        return IntervalTaylor(c, out)
    return dict(theta=IntervalTaylor(c, angular), pressure=square_tail(1+2*a), energy=square_tail(2*a)*2)


def heat_future_jets(c, Z, a, eps, S_cap, order=4, cells=256):
    """True scaled collar/Gamma angular, pressure and energy defect jets."""
    cells = operator.index(cells)
    if cells <= 0:
        raise ValueError('Positive collar cell count required')
    zero = IntervalTaylor.constant(c, 0, order)
    theta, pressure, energy = zero, zero, zero
    Sbox = c.mpf([0, endpoints(S_cap)[1]])
    for i in range(cells):
        left = c.mpf(3)*i/cells; right = c.mpf(3)*(i+1)/cells
        t = c.mpf([endpoints(left)[0], endpoints(right)[1]])
        sig = stable_sigma(c, t)[0]
        def phi(v):
            return c.exp(-4/(3-c.mpf(v))**2) if v < 3 else c.mpf(0)
        phi_box = c.mpf([endpoints(phi(endpoints(t)[1]))[0], endpoints(phi(endpoints(t)[0]))[1]])
        C = 1-eps*phi_box; pre = 1-eps*(1-sig+sig*phi_box)
        D = gamma_deficit_jets(c, Z, a, S_cap, order, t)*(sig*C)
        squaredifference = D*(D*(-a*Sbox)+2*pre)
        dt = c.mpf(3)/cells
        theta += D*(dt*a*c.exp((1-a)*t))
        pressure += squaredifference*(dt*c.exp(-(1+2*a)*t)/2)
        energy += squaredifference*(dt*c.exp(-2*a*t))
    tails = integrated_gamma_tails(c, Z, a, S_cap, order)
    return dict(theta=theta+tails['theta'], pressure=pressure+tails['pressure'], energy=energy+tails['energy'])


def implicit_quadratic_jets(c, x0, y0, b1, b2, p, q, nonlinear, prior_C1=None):
    """Differentiate exact equations, using one C0 Jacobian at every order."""
    if b1.order != b2.order:
        raise ValueError('RHS Taylor orders must agree')
    J21 = 1+2*nonlinear*x0; J22 = q*(1+2*nonlinear*y0)
    det = p*J22-J21
    if endpoints(det)[0] <= 0 <= endpoints(det)[1]:
        raise ArithmeticError('Implicit high-jet Jacobian is singular')
    xs, ys = [x0], [y0]
    for n in range(1, b1.order+1):
        known = nonlinear*sum((xs[i]*xs[n-i]+q*ys[i]*ys[n-i] for i in range(1, n)), c.mpf(0))
        v = b2[n]-known
        xn, yn = (J22*b1[n]-v)/det, (p*v-J21*b1[n])/det
        if n == 1 and prior_C1 is not None:
            xn = intersect(c, xn, prior_C1[0][1]); yn = intersect(c, yn, prior_C1[1][1])
        xs.append(xn); ys.append(yn)
    return dict(scaled_coefficient_Taylor=[IntervalTaylor(c, xs), IntervalTaylor(c, ys)],
                derivative_jacobian_determinant=det,
                recurrence='J*(x_n,y_n)=(b1_n,b2_n-K*scale*sum_i=1..n-1(x_i*x_(n-i)+q*y_i*y_(n-i)))')


class CompliantAngularHighJets:
    def __init__(self, order=4, cells=256):
        self.order = operator.index(order)
        if self.order < 2 or self.order > 4:
            raise ValueError('This admission covers orders2..4')
        self.repair = CompliantAngularRepair(); self.ctx = self.repair.ctx
        self.cells = operator.index(cells)
        if self.cells <= 0:
            raise ValueError('Positive collar cell count required')
        self.cache = {}; self.hashes = dict(self.repair.hashes)
        for stem, gate in (('compliant_outer_angular_repair_check', 'actual_implicit_functional_angular_pressure_repair_independently_checked'),
                           ('compliant_absolute_moment_closure_check', 'all_five_terminal_moment_identities_certified')):
            name = PREFIX+stem+'.json'; r = json.loads((HERE/name).read_bytes())
            if not r.get(gate) or r['actual_five_defect_family_sha256'] != self.repair.angular.initial.family:
                raise ValueError('High-jet prerequisite/family mismatch: '+name)
            for source, digest in r['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest() != digest:
                    raise ValueError('High-jet source changed: '+source)
                self.hashes[source] = digest
            self.hashes[name] = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def preheat_jets(self, Z):
        r = self.repair; c = self.ctx; Z = c.mpf(Z); N = self.order
        q = IntervalTaylor(c, [1+Z**2, 2*Z, 1]+[0]*(N-2)); lq = log_taylor(q)
        Xf = q.reciprocal()*(2*r.angular.Xv*c.exp(-100*r.rate))
        for i in range(r.cells):
            a = c.mpf(100)*i/r.cells; b = c.mpf(100)*(i+1)/r.cells
            sig = stable_sigma(c, c.mpf([endpoints(a)[0], endpoints(b)[1]])/100)[0]; power = 1-sig
            weight = (c.exp(-r.rate*(100-b))-c.exp(-r.rate*(100-a)))/r.rate
            Xf += (lq*(-power)+power*c.ln(2)).exp()*weight
        low = r.preheat_difference(Z)
        # Keep correlated value and slope instead of subtracting two
        # independent near-equal Xf enclosures at Z close to zero.
        return IntervalTaylor(c, list(low.coefficients)+[-Xf[n] for n in range(2, N+1)])

    def coefficients(self, Z):
        c = self.ctx; Z = c.mpf(Z); key = tuple(endpoints(Z))
        if key in self.cache:
            return self.cache[key]
        with mp.workdps(210):
            r = self.repair; pre = self.preheat_jets(Z)
            heat = heat_future_jets(c, Z, r.delta/2, r.heat.epsilon, r.strong_S_cap, self.order, self.cells)
            old = r.coefficients(Z); oldrhs = old['defects']
            def scaled_heat(jet, cap, oldjet):
                value = jet*c.mpf([0, endpoints(cap)[1]])
                coeffs = list(value.coefficients)
                # These providers enclose the same exact Gamma integrals.
                for j in (0, 1):
                    coeffs[j] = intersect(c, coeffs[j], oldjet[j])
                return IntervalTaylor(c, coeffs)
            rh = scaled_heat(heat['theta'], r.theta_heat_over_scale_cap, oldrhs['rheat_scaled_enclosure'])
            sh = scaled_heat(heat['pressure'], r.pressure_heat_over_scale_cap, oldrhs['s_scaled'])
            rhs = pre+rh
            b1 = rhs*(c.exp(r.rate)/r.weights['A']); b2 = sh*(c.exp(-3*r.prate)/r.weights['B'])
            prior = old['scaled_coefficient_Taylor']
            out = implicit_quadratic_jets(c, prior[0][0], prior[1][0], b1, b2, r.p, r.q, r.k*r.scale, prior)
            physical = [v*r.scale for v in out['scaled_coefficient_Taylor']]
            bounds = [sum((math.factorial(j)*c.mpf(magnitude(v[j])) for v in physical), c.mpf(0)) for j in range(self.order+1)]
            if any(endpoints(v-r.mu/100)[1] >= 0 for v in bounds):
                raise ArithmeticError('Uniform angular coefficient derivative smallness failed')
            out.update(Z=Z, physical_coefficient_Taylor=physical, preheat_difference_scaled_Taylor=pre,
                       scaled_Gamma_future_defect_Taylor=heat, r_scaled_Taylor=rhs, s_scaled_Taylor=sh,
                       physical_derivative_sum_bounds_by_order=bounds,
                       true_implicit_branch_retained=True, ordinary_Taylor_order=self.order,
                       true_Gamma_integral_derivatives_enclosed=True, exact_positive_S_not_replaced_by_cap=True,
                       angular_coefficient_C4_available=self.order == 4,
                       fifth_derivative_Taylor_remainder_available=False,
                       full_outer_C4_certified=False, full_outer_cone_certified=False, temporal_recursion=False)
            self.cache[key] = out
            return out

    def report(self):
        with mp.workdps(210):
            samples = [self.coefficients(z) for z in ('-1', '0', '.5', '1')]
            whole = self.coefficients([-1, 1])
            return dict(samples=samples, whole_Z_C4_coefficients=whole,
                ordinary_Taylor_order=self.order, angular_coefficient_C4_available=self.order == 4,
                implicit_source_sha256=self.repair.angular.initial.datum.source_sha,
                actual_five_defect_family_sha256=self.repair.angular.initial.family,
                exact_heat_definition=self.repair.heat.heat_definition,
                Gamma_derivative_bound='|d_d^m Dhat/m!| <= (1+a)_(m-1)*(1+a)_m*2^m*S^(m-1)*exp(-m*t)/m!, m>=2',
                finite_derivative_bounds_not_infinite_S_series=True,
                entire_Gamma_tail_integrals_included=True,
                angular_coefficient_derivatives_through4_below_mu_over100=True,
                fifth_derivative_Taylor_remainder_available=False,
                full_outer_C4_certified=False, full_outer_cone_certified=False,
                temporal_recursion=False, input_hashes=self.hashes)


def run():
    result = CompliantAngularHighJets().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)), indent=2)+'\n', encoding='utf-8')
    print('Actual implicit angular repair C4 jets and full Gamma-defect derivatives generated; full outer C4/cone/temporal pending', flush=True)
    return result


if __name__ == '__main__':
    run()
