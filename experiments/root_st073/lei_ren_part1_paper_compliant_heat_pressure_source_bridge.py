"""Defining-function and unit bridge from absolute closure to Gamma C4.

Returned interval boxes need not be equal. The two evaluators enclose one
full Gamma function, one amplitude and one original pressure history. This
bridge binds their production origins and proves the exact unit identities.
"""
import ast
from pathlib import Path

import sympy as s
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def source_bridge():
    proofs = {}

    def zero(name, expression):
        if s.simplify(s.expand_power_exp(expression)) != 0:
            raise ArithmeticError('Heat pressure source bridge failed: '+name)
        proofs[name] = True

    def syntax(stem, method, target, expected):
        tree = ast.parse((HERE/(PREFIX+stem+'.py')).read_text(encoding='utf8'))
        function = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == method)
        nodes = [n.value for n in ast.walk(function) if isinstance(n, ast.Assign)
                 and any(ast.unparse(t) == target for t in n.targets)]
        expressions = expected if isinstance(expected, list) else [expected]
        if len(nodes) != len(expressions) or any(
                ast.dump(node) != ast.dump(ast.parse(expression, mode='eval').body)
                for node, expression in zip(nodes, expressions)):
            raise ValueError('Heat pressure defining source changed: '+stem+':'+target)
        proofs['production_'+stem+'_'+target] = True

    # Both actual inlet amplitudes are evaluations of the same original O3
    # endpoint at Z=0. Different integration cells provide enclosures of
    # this function, not alternative amplitudes.
    syntax('compliant_axial_high_jets', '_incoming_constants', 'p0', 'buffer.power(0, 1, cells=128)')
    syntax('compliant_axial_high_jets', '_incoming_constants', 'U', "box(p0['Utheta_over_Pstar'][0])")
    syntax('compliant_outer_angular_candidate', '__init__', 'inlet', "self.buffer.power('0', 1, cells)")
    syntax('compliant_outer_angular_candidate', '__init__', 'self.loguRp0', "c.ln(inlet['Utheta_over_Pstar'][0])")
    # The exact inverse radius originates in the same actual heat component;
    # strong_S_cap is a bound, not the definition of S.
    syntax('compliant_exact_heat_component', '__init__', 'self.logradius_terms',
           'dict(selected_reference=self.logRref, pulse_term=13/self.params.mu, finite_offset=self.tail_finite)')
    syntax('compliant_outer_angular_repair', '__init__', 'logRtail', 'self.heat.logRref+13/self.mu+self.heat.tail_finite')
    syntax('compliant_collar_Gamma_C4', '__init__', 'self.Scap', 'self.steep.S')
    syntax('compliant_steep_waiting_C4', '__init__', 'self.S',
           "read_interval(c, json.loads((HERE/(PREFIX+'compliant_outer_angular_repair.json')).read_bytes())['strong_inverse_radius_positive_cap'])")
    syntax('compliant_corrected_outer_field', 'heat_local', 'Sbox',
           ['self.repair.strong_S_cap*self.pulse.factor(-t)', 'c.mpf([0,endpoints(Sbox)[1]])'])
    syntax('compliant_collar_Gamma_C4', 'local_Gamma', 'Sc', 'self.S*c.exp(-t)')
    syntax('compliant_collar_Gamma_C4', 'gamma_deficit_mixed', 'xi', 'base*S')
    syntax('compliant_collar_Gamma_C4', 'gamma_deficit_mixed', 'base', 'd*(2*decay)')
    syntax('compliant_collar_Gamma_C4', 'gamma_deficit_mixed', 'd', 'IntervalTaylor(c, [1-Z**2, -2*Z, -1, 0, 0, 0])')
    syntax('compliant_exact_heat_component', '__init__', 'self.heat_definition',
           repr('H_delta(xi)=Gamma(1+a)^-1 integral_0^infinity exp(-v)*v^a*(1+xi*v)^-a dv, a=delta/2'))

    mu, a, length, Ts, tau, lone, t = s.symbols('mu a Lrel Ts tau lone t', real=True)
    bp, bh, rate, k, p = s.Rational(1, 2)+mu, s.Rational(1, 2)+a, 1-mu, 1-a, 1+2*a
    env = {'self.bp': bp, 'self.bh': bh, 'self.rate': rate, 'self.k': k,
           'self.outer.Lrel': length, 'self.Lrel': length, 'self.Ts': Ts,
           'self.params.Ts': Ts, 'self.mu': mu, 'self.a': a,
           'self.phrate': p, 'self.prate': 1+2*mu, 'self.angular.waiting': tau,
           'self.angular.waiting_logone': lone}
    theta = {}
    for target in ('self.thetaR', 'self.thetaS', 'self.thetaQ', 'self.thetaT'):
        theta[target] = assignment('compliant_steep_waiting_C4', '__init__', target, {**env, **theta})
    logs = {}
    for target in ('self.Lf', 'self.LR', 'self.Ls', 'self.Lq', 'self.Lt', 'self.Ltail'):
        logs[target] = assignment('compliant_corrected_outer_field', '__init__', target, {**env, **logs})
    for theta_name, log_name in (('self.thetaR', 'self.LR'), ('self.thetaS', 'self.Ls'),
                                 ('self.thetaQ', 'self.Lq'), ('self.thetaT', 'self.Lt')):
        zero('same_amplitude_'+theta_name[5:], theta[theta_name]-s.exp(logs[log_name]))
    theta_base = assignment('compliant_collar_Gamma_C4', '__init__', 'self.theta_base',
                            {'self.steep.thetaT': theta['self.thetaT'], 'self.bh': bh,
                             'self.steep.wait': tau, 'self.steep.logone': lone})
    zero('same_actual_heat_velocity_amplitude', theta_base-s.exp(logs['self.Ltail']-lone))
    weights = {}
    for target in ('self.Pf', 'self.Prel', 'self.Ps', 'self.Pq', 'self.Pt', 'self.Ptail', 'self.tailmult'):
        weights[target] = assignment('compliant_corrected_outer_field', '__init__', target, {**env, **weights})
    zero('Rtail_pressure_units_equal_theta_base_squared',
         weights['self.Prel']*weights['self.Ptail']*weights['self.tailmult']-theta_base**2)
    common_scale = s.symbols('Ev0_squared_over_Pstar_squared', positive=True)
    zero('production_C4_pressure_scale',
         assignment('compliant_collar_Gamma_C4', '__init__', 'self.pressure_scale',
                    {'self.Ev2': common_scale, 'self.theta_base': theta_base})-common_scale*theta_base**2)
    syntax('compliant_flatten_mixed_C4', '__init__', 'self.logEv2_parts',
           'dict(inlet_log=2*c.ln(self.U), inverse_mu_term=-13/self.mu, finite_offset=c.mpf(-26))')
    syntax('compliant_corrected_outer_field', '__init__', 'self.logEv0Pstar2_parts',
           'dict(inlet_log=2*self.angular.loguRp0, inverse_mu_term=-13/self.mu, finite_offset=c.mpf(-26))')

    # Shared defining Gamma expectation; the accepted low-order and C4
    # moment bounds are different enclosures of these very same integrals.
    z, S, H = s.symbols('Z S H', real=True)
    d = 1-z*z; xi = 2*d*S*s.exp(-t)
    deficit_Rtail = (1-H)/(a*S)
    deficit_current = (1-H)/(a*S*s.exp(-t))
    zero('same_full_Gamma_K_from_Rtail_or_current_deficit',
         1-a*S*deficit_Rtail-(1-a*S*s.exp(-t)*deficit_current))
    zero('same_full_Gamma_argument', xi-2*d*(S*s.exp(-t)))
    zero('full_Gamma_pressure_deficit_integrand',
         H**2-(1-a*S*s.exp(-t)*deficit_current*(2-a*S*s.exp(-t)*deficit_current)))
    A, C = s.symbols('A_pressure C_pressure', real=True)
    zero('pressure_tail_normalization_cancels_K_squared', (A/H**2)*H**2-A)
    zero('same_exterior_tail_pressure',
         -C*theta_base**2*s.exp(-p*t)*A
         +C*weights['self.Prel']*weights['self.Ptail']*weights['self.tailmult']*s.exp(-p*t)*A)

    # Retained forward history: the common inlet P0+Mp(Rv) plus the same
    # original swirl integral determines both pressures uniquely. The
    # earlier closure proof fixes the infinity constant; no P3=-tail
    # substitution is used to manufacture this bridge.
    P0, Mp, Prv, A3 = s.symbols('P0 Mp_Rv Prv A3', real=True)
    P3 = P0+Mp+C*(Prv-theta_base**2*s.exp(-3*p)*A3)
    zero('retained_P3_plus_heat_tail_equals_original_infinity_offset',
         P3+C*theta_base**2*s.exp(-3*p)*A3-(P0+Mp+C*Prv))
    Pforward = P3+C*theta_base**2*(A3*s.exp(-3*p)-A*s.exp(-p*t))
    zero('same_forward_pressure_for_entire_exterior',
         Pforward-(-C*theta_base**2*s.exp(-p*t)*A)-(P0+Mp+C*Prv))
    syntax('compliant_collar_Gamma_C4', 'data', 'Ptail', "terminal['pressure_over_Pstar_squared_Taylor']")
    syntax('compliant_collar_Gamma_C4', 'data', 'pressure3', 'self.forward_pressure(Z, 3, Ptail)')
    syntax('compliant_collar_Gamma_C4', 'forward_pressure', 'K', "self.shape(Z, v, False)['K_rows'][0]")
    return dict(identities=proofs,
                complete_defining_function_history_bridge_verified=False,
                unresolved_bindings=[
                    'Bind both inlet buffer objects to the same SharedOuterBuffer callable',
                    'Bind exact S=exp(-logRtail), independently of all numerical caps',
                    'Bind actual C4 Ptail/forward_pressure integrals to retained closure P0/Mp/Prv histories',
                    'Bind canonical full Gamma defining expectation to the C4 derivative evaluator'],
                shared_defining_heat_function='Full positive Gamma expectation with a=delta/2, xi=2*(1-Z^2)*S*exp(-t)',
                shared_pressure_tail='A_p=integral_0^infinity exp(-(1+delta)*v)*H(xi*exp(-v))^2/2 dv',
                pressure_infinity_constant_source='Original P0+Mp(Rv)+(Ev0^2/Pstar^2)*Prv; zero by the independently admitted absolute closure',
                pressure_equivalence_uses_interval_overlap=False,
                interval_boxes_or_caps_treated_as_exact_functions=False)
