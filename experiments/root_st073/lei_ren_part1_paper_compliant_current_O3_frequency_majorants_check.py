"""Focused independent differentiation check of conditional N envelopes.

Synthetic analytic cutoff data exercise the general conditional inequalities;
they are not a certification of the current graph's whole-domain source norms.
"""
import hashlib
import json
from pathlib import Path
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_frequency_majorants as source

RECEIPT = source.PREFIX + 'current_O3_frequency_majorants_check.json'


def run():
    data = json.loads((source.HERE/source.NAME).read_bytes())
    for name, digest in data['input_hashes'].items():
        if hashlib.sha256((source.HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Conditional envelope source differs: ' + name)
    theorem = source.exact_frequency_majorant_theorem()
    if source._encode(theorem) != data['exact_actual_source_frequency_majorant_theorem']:
        raise ValueError('Exact source-bound finite-N derivation differs')
    if source._encode(source.symbolic_envelopes()) != data['conditional_symbolic_envelopes']:
        raise ValueError('Published symbolic envelope differs')
    # Independently differentiate the defining functions, rather than the
    # production recurrence. These are synthetic jets for conditional tests.
    t, z = s.symbols('t z', real=True)
    mu = s.Symbol('mu', positive=True)
    N = s.Symbol('N', integer=True, positive=True)
    chi = (1-t*t)**5
    et, cz = s.exp(-t/2), 1/(1+z*z)
    h = -mu*chi**2*s.sin(4*s.pi*N*t)/(8*s.pi*N)
    profiles = dict(theta=et*s.exp(h), increment=et*(s.exp(h)-1),
                    axial=-s.sqrt(mu)*et*chi*s.cos(2*s.pi*N*t)/(2*s.pi*N))
    radial = {name: [s.diff(value, t, j) for j in range(5)] for name, value in profiles.items()}
    axial_shape = [s.diff(cz, z, k) for k in range(6)]
    tested = 0
    cases = ((s.Rational(-371,1000), s.Rational(137,1000), 3),
             (s.Rational(127,1000), s.Rational(-237,1000), 13),
             (s.Rational(631,1000), s.Rational(337,1000), 73))
    for tv, zv, frequency in cases:
        m = s.Rational(1,37)
        env = {t:tv, z:zv, mu:m, N:frequency}
        cutoff = [abs(s.diff(chi,t,j).subs(t,tv)) for j in range(5)]
        F = [[abs((s.diff(et,t,j)*axial_shape[k]).subs({t:tv,z:zv})) for k in range(6)] for j in range(5)]
        caps = source.profile_majorants(m, s.Integer(frequency), cutoff, F, s.Rational(1,2))
        for name, key in (('theta','modified_theta_over_Pstar_majorants'),
                          ('increment','theta_increment_over_Pstar_majorants'),
                          ('axial','modified_axial_over_Pstar_majorants')):
            radial_values = [value.subs(env).evalf(80) for value in radial[name]]
            z_values = [value.subs(z,zv).evalf(80) for value in axial_shape]
            for j in range(5):
                for k in range(6):
                    actual = abs(radial_values[j]*z_values[k])
                    upper = caps[key][j][k].evalf(80)
                    if actual > upper + s.Float('1e-65',80)*max(s.Integer(1),upper):
                        raise ArithmeticError('Independent conditional derivative exceeds envelope: %s/%d/%d' % (name,j,k))
                    tested += 1
        # Both exact shear deviations against the same supported loop.
        slow_a = -m*chi*s.diff(chi,t)*s.sin(4*s.pi*frequency*t)/(4*s.pi)
        loop_b = 2*s.sqrt(m)*chi*s.sin(2*s.pi*frequency*t)
        slow_b = -s.sqrt(m)*(chi*(-s.Rational(1,2))+s.diff(chi,t))*s.cos(2*s.pi*frequency*t)/(2*s.pi)
        H = h.subs({mu:m,N:frequency})
        for actual, key in ((-2*slow_a/frequency,'theta_shear_error_vs_same_periodic_loop_upper'),
                           ((s.exp(-H)-1)*loop_b+s.exp(-H)*2*slow_b/frequency,'axial_shear_error_vs_same_periodic_loop_upper')):
            value = abs(actual.subs(t,tv).evalf(80))
            upper = caps[key].evalf(80)
            if value > upper + s.Float('1e-65',80)*max(s.Integer(1),upper):
                raise ArithmeticError('Independent finite-N shear error exceeds envelope')
            tested += 1
    args = source.symbolic_envelopes()
    rejected = 0
    for mu_bad, n_bad in ((-1,3),(s.oo,3),(1,0),(1,s.Rational(3,2)),(1,True)):
        try:
            source.profile_majorants(mu_bad,n_bad,[1]*5,[[1]*6 for _ in range(5)],1)
        except ValueError:
            rejected += 1
        else:
            raise ArithmeticError('Invalid parameter accepted')
    try:
        source.profile_majorants(s.Rational(1,37),args['finite_integer_N'],
            [args['finite_integer_N']]*5,[[1]*6 for _ in range(5)],1)
    except ValueError:
        rejected += 1
    else:
        raise ArithmeticError('N-dependent slow norm accepted')
    hashes = dict(data['input_hashes'])
    for name in (source.NAME,Path(__file__).name):
        hashes[name] = hashlib.sha256((source.HERE/name).read_bytes()).hexdigest()
    receipt = dict(all_passed=True,exact_source_frequency_identities=len(theorem['identities']),
        independent_synthetic_derivative_and_shear_comparisons=tested,
        synthetic_cases_are_not_actual_whole_domain_norm_certification=True,
        rejected_invalid_or_N_dependent_inputs=rejected,
        conditional_N_growth_degrees=theorem['conditional_N_growth_degrees_after_exponential_factor'],
        conditional_actual_modulation_derivative_envelopes_available=True,
        actual_whole_domain_norm_inputs_certified=False,actual_independent_repair_envelopes_installed=False,
        current_common_cone_frequency_N_certified=False,current_modified_whole_cones_certified=False,
        input_hashes=hashes)
    (source.HERE/RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode('utf8'))
    print('Conditional finite-N envelopes PASS:',len(theorem['identities']),'source identities;',tested,
          'independent derivative/shear comparisons;',rejected,'invalid inputs rejected',flush=True)
    return receipt


if __name__ == '__main__':
    run()
