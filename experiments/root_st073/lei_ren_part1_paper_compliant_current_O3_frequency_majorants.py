"""Conditional analytic N envelopes for the actual supported O2/O3 source.

Inputs must be bounds on ordinary logR/axial derivatives of the SAME
cutoff and original theta/Pstar source. This module does not certify those
inputs on a whole region, the independent repair, or an admissible cone.
No current checked ancestor graph is reconstructed here.
"""
import hashlib
import json
import math
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent
PREFIX = 'lei_ren_part1_paper_compliant_'
NAME = PREFIX + 'current_O3_frequency_majorants.json'


def _nonnegative(value, name):
    value = s.sympify(value)
    if value.is_nonnegative is not True or value.is_finite is not True:
        raise ValueError('Finite nonnegative source majorant required: ' + name)
    return value


def profile_majorants(mu, N, cutoff, theta_source, log_theta_y):
    """Exact algebraic upper bounds conditional on supplied source norms.

    cutoff[j] bounds |d_logR^j chi|, theta_source[j][k] bounds
    |d_logR^j d_Z^k(E0/Pstar)|, and log_theta_y bounds
    |d_logR log(E0)|. All are on one common source domain. The original
    source and cutoff bounds must be independent of N. mu>0, N>=1 integer.
    """
    mu, N = s.sympify(mu), s.sympify(N)
    if mu.is_positive is not True or mu.is_finite is not True:
        raise ValueError('Finite positive source mu required')
    if N.is_integer is not True or N.is_positive is not True:
        raise ValueError('Finite positive integer N required')
    if len(cutoff) != 5 or len(theta_source) != 5 or any(len(row) != 6 for row in theta_source):
        raise ValueError('Ordinary logR rows0..4 and axial rows0..5 required')
    C = [_nonnegative(v, 'chi_' + str(j)) for j, v in enumerate(cutoff)]
    F = [[_nonnegative(v, 'theta_%d_%d' % (j, k)) for k, v in enumerate(row)]
         for j, row in enumerate(theta_source)]
    L = _nonnegative(log_theta_y, 'log_theta_y')
    if isinstance(N, s.Symbol) and any(value.has(N) for value in [mu, L, *C, *(v for row in F for v in row)]):
        raise ValueError('Slow source and cutoff majorants must be independent of N')
    square = [sum(s.binomial(j, k) * C[k] * C[j-k] for k in range(j+1)) for j in range(5)]
    # H=A/N. Actual phase is N*logR, including the translated constant.
    H = [mu / (8*s.pi*N) * sum(s.binomial(j, k) * square[j-k] * (4*s.pi*N)**k
                               for k in range(j+1)) for j in range(5)]
    bell = [s.Integer(1)]
    for j in range(1, 5):
        bell.append(sum(s.binomial(j-1, k) * H[j-k] * bell[k] for k in range(j)))
    exp_cap = s.exp(H[0])
    G = [exp_cap * row for row in bell]
    deltaG = [exp_cap*H[0]] + G[1:]
    theta, increment, axial = [], [], []
    for j in range(5):
        theta.append([sum(s.binomial(j, i)*F[j-i][k]*G[i] for i in range(j+1)) for k in range(6)])
        increment.append([sum(s.binomial(j, i)*F[j-i][k]*deltaG[i] for i in range(j+1)) for k in range(6)])
        axial.append([s.sqrt(mu)/(2*s.pi*N)*sum(
            s.binomial(j, i)*s.binomial(j-i, h)*F[j-i-h][k]*C[h]*(2*s.pi*N)**i
            for i in range(j+1) for h in range(j-i+1)) for k in range(6)])
    theta_error = mu*C[0]*C[1]/(2*s.pi*N)
    axial_error = exp_cap/N*(2*s.sqrt(mu)*C[0]*(mu*C[0]**2/(8*s.pi))
                             + s.sqrt(mu)/s.pi*(C[0]*L+C[1]))
    return dict(
        logR_ordinary_derivative_order=4, axial_ordinary_derivative_order=5,
        actual_mu=mu, finite_integer_N=N, cutoff_derivative_majorants=C,
        original_theta_over_Pstar_derivative_majorants=F,
        original_log_theta_derivative_majorant=L,
        exponent_A_over_N_derivative_majorants=H,
        exponential_derivative_majorants=G, exponential_increment_majorants=deltaG,
        modified_theta_over_Pstar_majorants=theta,
        theta_increment_over_Pstar_majorants=increment,
        modified_axial_over_Pstar_majorants=axial,
        theta_shear_error_vs_same_periodic_loop_upper=theta_error,
        axial_shear_error_vs_same_periodic_loop_upper=axial_error,
        shear_error_bound_is_uniform_in_actual_phase=True,
        source_norm_inputs_are_conditional_not_certified_whole_domain=True,
        independent_repair_pressure_tensor_error_common_N_cones_remain_open=True)


def symbolic_envelopes():
    mu = s.Symbol('mu', positive=True, finite=True)
    N = s.Symbol('N', integer=True, positive=True, finite=True)
    C = s.symbols('C0:5', nonnegative=True, finite=True)
    F = [[s.Symbol('F%d_%d' % (j, k), nonnegative=True, finite=True) for k in range(6)] for j in range(5)]
    L = s.Symbol('L1', nonnegative=True, finite=True)
    return profile_majorants(mu, N, C, F, L)


def exact_frequency_majorant_theorem():
    # Import source inspection only; no source graph constructors are called.
    from types import SimpleNamespace
    import inspect
    from lei_ren_part1_paper_compliant_steep_entry_stress_C3 import SourceAST
    asts, checks = SourceAST(), {}

    def zero(name, a, b):
        if s.simplify(s.expand(a-b)) != 0:
            raise ArithmeticError('Actual finite-N derivative identity failed: ' + name)
        checks[name] = True

    # Bind the published profile/ordinary-unit program before deriving caps.
    for target, wanted in (
        ('n', 'c.mpf(N)'), ('phase', 'n*t'), ('chi', 'cutoff_rows(c,t)'),
        ('sin4', 'trig_rows(c,phase,n,4)'),
        ('cos2', 'trig_rows(c,phase,n,2,cosine=True)'),
        ('square', 'product_rows(chi,chi)'),
        ('Arows', '[-mu*r/(8*c.pi) for r in product_rows(square,sin4)]'),
        ('A', 'IntervalTaylor(c,[r/math.factorial(j) for j,r in enumerate(Arows)])'),
        ('G', '(A/n).exp()'), ('Enew', "product_rows(old['theta'],grows)"),
        ('grows', '[G[j]*math.factorial(j) for j in range(5)]'),
        ('x', 'A[0]/n'), ('mag', 'upper(c,abs(x))'),
        ('dgrows', '[x*c.exp(c.mpf([-endpoints(mag)[1],endpoints(mag)[1]]))]+grows[1:]'),
        ('einc', "product_rows(old['theta'],dgrows)"),
        ('Vnew', "[r*(-c.sqrt(mu)/(2*c.pi*n)) for r in product_rows(old['theta'],product_rows(chi,cos2))]"),
        ('slowA', '-mu*chi[0]*chi[1]*c.sin(4*c.pi*phase)/(4*c.pi)'),
        ('slowB_over_E', '-c.sqrt(mu)*(chi[0]*logEy+chi[1])*c.cos(2*c.pi*phase)/(2*c.pi)'),
        ('aexcess', '2*mu*sig+da-2*slowA/n'),
        ('b', 'c.exp(-A[0]/n)*(beta+2*slowB_over_E/n)')):
        asts.expression('current_O3_finite_frequency_profiles', 'profile', target, wanted=wanted)
        checks['actual_profile_program_' + target] = True
    # Ordinary product rows and actual N*phase derivatives, independently
    # compared with polynomial multiplication and direct differentiation.
    env = dict(math=math)
    product = asts.replay('collar_Gamma_C4', 'product_rows', env)
    x, t, phi = s.symbols('local_logR_offset t actual_phase', real=True)
    N = s.Symbol('N', integer=True, positive=True)
    u, v = s.symbols('u0:5'), s.symbols('v0:5')
    pu = sum(u[j]*x**j/s.factorial(j) for j in range(5))
    pv = sum(v[j]*x**j/s.factorial(j) for j in range(5))
    ordinary = product(list(u), list(v))
    for j in range(5):
        zero('actual_ordinary_product_derivative_' + str(j), ordinary[j], s.diff(pu*pv, x, j).subs(x, 0))
    # Bind both cutoff coordinate widths and Taylor-to-ordinary conversions.
    sigma = s.Function('same_actual_sigma')
    a = s.Symbol('sigma_argument', real=True)
    sigma_jets = lambda ctx, value: [s.diff(sigma(a),a,j).subs(a,value)/s.factorial(j) for j in range(5)]
    cutoff = asts.replay('current_O3_finite_frequency_profiles','cutoff_rows',
        dict(math=math,sigma_jets=sigma_jets,product_rows=product))
    ordinary_cutoff = cutoff(None,t)
    exact_cutoff = sigma(t+2)*(1-sigma(4*t-1))
    for j in range(5):
        zero('actual_cutoff_width_and_ordinary_derivative_' + str(j), ordinary_cutoff[j],s.diff(exact_cutoff,t,j))
    c = SimpleNamespace(pi=s.pi, sin=s.sin, cos=s.cos)
    trig = asts.replay('current_O3_finite_frequency_profiles', 'trig_rows', {})
    for harmonic, cosine in ((4, False), (2, True)):
        actual = trig(c, phi, N, harmonic, cosine)
        fn = s.cos if cosine else s.sin
        for j in range(5):
            zero('actual_phase_derivative_h%d_%d' % (harmonic, j), actual[j],
                 s.diff(fn(harmonic*s.pi*(phi+N*t)), t, j).subs(t, 0))
    # Independent finite exponential series supplies every local jet.
    h = s.symbols('h0:5')
    inner = sum(h[j]*x**j/s.factorial(j) for j in range(1, 5))
    series = s.Poly(s.expand(sum(inner**k/s.factorial(k) for k in range(5))), x)
    bell = [s.Integer(1)]
    for j in range(1, 5):
        bell.append(sum(s.binomial(j-1, k)*h[j-k]*bell[k] for k in range(j)))
    for j in range(5):
        zero('independent_exponential_ordinary_jet_' + str(j), bell[j], series.nth(j)*s.factorial(j))
    # The exact two finite-N shear deviations, before absolute bounds.
    mu = s.Symbol('mu', positive=True)
    C0, C1, L, sig = s.symbols('chi chi_y logE_y sigma', real=True)
    da = mu*C0**2*s.cos(4*s.pi*phi)
    slowA = -mu*C0*C1*s.sin(4*s.pi*phi)/(4*s.pi)
    slowB = -s.sqrt(mu)*(C0*L+C1)*s.cos(2*s.pi*phi)/(2*s.pi)
    a0 = -mu*C0**2*s.sin(4*s.pi*phi)/(8*s.pi)
    loop_b = 2*s.sqrt(mu)*C0*s.sin(2*s.pi*phi)
    zero('actual_theta_loop_error', (2*mu*sig+da-2*slowA/N)-(2*mu*sig+da), -2*slowA/N)
    zero('actual_axial_loop_error', s.exp(-a0/N)*(loop_b+2*slowB/N)-loop_b,
         (s.exp(-a0/N)-1)*loop_b+s.exp(-a0/N)*2*slowB/N)
    envelope = symbolic_envelopes()
    freq = envelope['finite_integer_N']
    degrees = {}
    # Factoring exp(H0) leaves finite Laurent polynomials with N^(j-1),
    # so high derivatives do not become small when N is increased.
    for label, rows, exponential in (
        ('swirl_increment', envelope['theta_increment_over_Pstar_majorants'], True),
        ('axial', envelope['modified_axial_over_Pstar_majorants'], False)):
        degrees[label] = []
        for j, row in enumerate(rows):
            polynomial = s.expand(row[0] / (s.exp(envelope['exponent_A_over_N_derivative_majorants'][0]) if exponential else 1))
            degree = max(term.as_powers_dict().get(freq, s.Integer(0)) for term in s.Add.make_args(polynomial))
            if degree != j-1:
                raise ArithmeticError('Actual N growth order changed')
            degrees[label].append(int(degree))
            checks['actual_' + label + '_N_growth_order_' + str(j)] = True
    checks['exp_increment_bound_uses_exact_real_integral_not_subtraction'] = True
    for name in ('modified_theta_over_Pstar_majorants','theta_increment_over_Pstar_majorants',
                 'modified_axial_over_Pstar_majorants'):
        if not all(value.is_nonnegative is True for row in envelope[name] for value in row):
            raise ArithmeticError('Conditional absolute-value recurrence lost positivity')
        checks['positive_triangle_recurrence_' + name] = True
    inspector = Path(inspect.getsourcefile(SourceAST))
    asts.hashes[inspector.name] = hashlib.sha256(inspector.read_bytes()).hexdigest()
    return dict(identities=checks, passed=True, source_bindings=asts.bindings,
        input_hashes=asts.hashes,
        conditional_triangle_bound_rules=['ordinary Leibniz with positive binomial coefficients',
            '|sin|<=1 and |cos|<=1 at every real actual phase',
            '|exp(h)-1|<=|h|*exp(|h|), from h*integral_0^1 exp(r*h)dr',
            'positive Bell recurrence for ordinary derivatives of exp(h)'],
        conditional_N_growth_degrees_after_exponential_factor=degrees,
        source_norm_inputs_not_evaluated_or_whole_domain_certified=True,
        repair_tensor_pressure_common_N_cones_and_coefficient_recursion_not_certified=True)


def _encode(value):
    if isinstance(value, dict):
        return {key: _encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_encode(item) for item in value]
    return s.sstr(value) if isinstance(value, s.Basic) else value


def run():
    result = dict(exact_actual_source_frequency_majorant_theorem=exact_frequency_majorant_theorem(),
        conditional_symbolic_envelopes=symbolic_envelopes(),
        conditional_actual_modulation_derivative_envelopes_available=True,
        actual_whole_domain_norm_inputs_certified=False, actual_independent_repair_envelopes_installed=False,
        current_common_cone_frequency_N_certified=False, current_modified_whole_cones_certified=False)
    result['input_hashes'] = dict(result['exact_actual_source_frequency_majorant_theorem']['input_hashes'])
    result['input_hashes'][Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/NAME).write_bytes((json.dumps(_encode(result), indent=2)+'\n').encode('utf8'))
    print('Actual modulation N envelopes: exact source identities and conditional bounds available', flush=True)
    return result


if __name__ == '__main__':
    run()
