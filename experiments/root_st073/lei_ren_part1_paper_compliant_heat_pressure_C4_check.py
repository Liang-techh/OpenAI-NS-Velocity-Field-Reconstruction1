"""Independent integral/differentiation checks of the absolute heat pressure."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_heat_pressure_C4 import pressure_y_rows
from lei_ren_part1_paper_compliant_absolute_moment_closure import source_identities
from lei_ren_part1_paper_compliant_heat_pressure_source_bridge import source_bridge
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import gamma_deficit_mixed
from lei_ren_part1_paper_compliant_angular_high_jets import integrated_gamma_tails
from lei_ren_part1_paper_compliant_pulse_physical_bounds import P
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).parent
NAME = 'lei_ren_part1_paper_compliant_heat_pressure_C4.json'


def independent_integral_fixture():
    """Full hyperu heat function, direct pressure-tail quadrature and jets.

    Finite parameters exercise the genuine Gamma source without fitting a
    pressure constant. This fixture does not certify the actual parameter
    choices. Actual enclosures come from admitted same-source receipts.
    """
    with mp.workdps(45):
        c = MPIntervalContext(); c.dps = 100
        a, S, Z, t = map(mp.mpf, ('.15', '.04', '.4', '3'))
        rate = 1+2*a

        def H(xi, n=0):
            if not xi:
                return (-1)**n*mp.rf(a, n)*mp.rf(1+a, n)
            return (-1)**n*mp.rf(a, n)*mp.rf(1+a, n)*xi**(-1-a-n)*mp.hyperu(1+a+n, 2, 1/xi)

        xi = 2*(1-Z**2)*S*mp.exp(-t)
        direct_H = mp.quad(lambda v: mp.exp(-v)*v**a*(1+xi*v)**(-a), [0, 1, mp.inf])/mp.gamma(1+a)
        if abs(H(xi)-direct_H) > mp.mpf('1e-38'):
            raise ArithmeticError('Independent full Gamma representation disagrees')
        D = gamma_deficit_mixed(c, c.mpf(Z), c.mpf(a), c.mpf(S), c.mpf(t))
        one = IntervalTaylor.constant(c, 1, 5)
        K = [one-D[0]*(a*S)]+[-row*(a*S) for row in D[1:]]
        tails = integrated_gamma_tails(c, c.mpf(Z), c.mpf(a), c.mpf(S*mp.exp(-t)), 5, 0)
        numerator = one/(2*rate)-tails['pressure']*(a*S*mp.exp(-t))
        rows = pressure_y_rows(K, numerator, c.mpf(rate), c.mpf(1), c.mpf(t))

        # s=exp(-v) maps the entire radial tail to s in[0,1].
        target = -mp.exp(-rate*t)*mp.quad(lambda u: u**(rate-1)*H(xi*u)**2, [0, '.25', 1])/2
        if not endpoints(rows[0][0])[0] <= target <= endpoints(rows[0][0])[1]:
            raise ArithmeticError('Direct infinite pressure tail outside enclosure')
        # Differentiate under the true infinite pressure integral in Z.
        target_Z = -mp.exp(-rate*t)*mp.quad(
            lambda u: u**(rate-1)*H(xi*u)*H(xi*u, 1)*(-4*Z*S*mp.exp(-t)*u),
            [0, '.25', 1])
        if not endpoints(rows[0][1])[0] <= target_Z <= endpoints(rows[0][1])[1]:
            raise ArithmeticError('Direct pressure axial derivative outside enclosure')

        prime = lambda y, z: mp.exp(-rate*y)*H(2*(1-z*z)*S*mp.exp(-y))**2/2
        checks = 0
        for k, n in ((1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (3, 0), (4, 0)):
            value = mp.diff(lambda z: mp.diff(lambda y: prime(y, z), t, k-1), Z, n)/math.factorial(n)
            lo, hi = endpoints(rows[k][n])
            if not lo <= value <= hi:
                raise ArithmeticError('Independent heat pressure mixed derivative failed: '+str((k, n)))
            checks += 1
        return dict(full_Gamma_quadrature_checked=True, infinite_pressure_integral_checks=2,
                    independent_mixed_derivative_checks=checks, moderate_fixture_only=True, passed=True)


def pressure_integral_identities():
    y, rate = s.symbols('y rate', positive=True)
    h = s.Function('H')(y)
    source = s.exp(-rate*y)*h**2/2
    identities = {}
    for order in range(1, 5):
        leibniz = s.exp(-rate*y)*sum(
            s.binomial(order-1, j)*(-rate)**(order-1-j)*s.diff(h**2, y, j)
            for j in range(order))/2
        if s.simplify(s.diff(source, y, order-1)-leibniz) != 0:
            raise ArithmeticError('Pressure FTC/Leibniz identity failed')
        identities['pressure_derivative_'+str(order)+'_from_defining_integral'] = True
    return identities


def run():
    record = json.loads((HERE/NAME).read_bytes())
    hashes = dict(record['input_hashes'])
    for filename, digest in hashes.items():
        if hashlib.sha256((HERE/filename).read_bytes()).hexdigest() != digest:
            raise ValueError('Heat pressure dependency changed: '+filename)
    closure = json.loads((HERE/'lei_ren_part1_paper_compliant_absolute_moment_closure.json').read_bytes())
    if closure['actual_five_defect_family_sha256'] != record['actual_five_defect_family_sha256'] or closure['source_sha256'] != record['implicit_source_sha256']:
        raise ValueError('Absolute closure and heat C4 use different sources')
    absolute_proof = source_identities()
    if absolute_proof != closure['proof']:
        raise ValueError('Admitted absolute pressure source proof changed')
    bridge = source_bridge()
    if bridge != record['defining_source_bridge']:
        raise ValueError('Heat defining source/scale bridge changed')
    if not bridge['complete_defining_function_history_bridge_verified'] or bridge['unresolved_bindings']:
        raise ValueError('Complete actual defining-function/history transfer required')
    c = MPIntervalContext(); c.dps = 270
    read = lambda q: read_interval(c, q)
    overlaps = 0
    for point in record['samples']+[record['whole_Z_inlet'], record['whole_Z_unbounded_exterior']]:
        if not point['original_pressure_datum_and_forward_history_retained']:
            raise ValueError('Original pressure history lost')
        if (not point['pressure_infinity_offset_exactly_zero_from_source_closure']
                or not point['defining_source_and_pressure_scale_bridge_verified']
                or not point['absolute_pressure_same_source_mixed4_available']
                or point['source_history_transfer_conditional']):
            raise ValueError('Admitted pressure source/history transfer lost')
        if point['heat_exterior_stress_identity_certified'] or point['global_admissible_stress_lift_constructed'] or point['temporal_recursion']:
            raise ValueError('Pressure adapter overclaims stress or recursion')
        ratio = point['pressure_over_Utheta_squared_Taylor']['coefficients'][0]
        if endpoints(read(ratio))[1] >= 0:
            raise ArithmeticError('Remaining pressure tail is not strictly positive in normalized units')
        for key, new in point['physical_mixed_derivatives_total_order_le4'][P].items():
            nl, nh = endpoints(read(new)); ol, oh = endpoints(read(point['original_forward_pressure_mixed_bounds'][key]))
            if max(nl, ol) > min(nh, oh):
                raise ArithmeticError('Same-source forward and tail pressure bounds disagree: '+key)
            overlaps += 1
    # Independent pre-existing C1 absolute closure at the full-Z heat inlet.
    reference = closure['samples'][-1]
    if reference['coordinate']['offset']['lower'] != '3.0' or reference['Z']['lower'] != '-1.0' or reference['Z']['upper'] != '1.0':
        raise ValueError('Unexpected absolute closure reference domain')
    for new, old in zip(record['whole_Z_inlet']['pressure_over_Pstar_squared_Taylor']['coefficients'], reference['P_over_Pstar_squared']['coefficients']):
        nl, nh = endpoints(read(new)); ol, oh = endpoints(read(old))
        if max(nl, ol) > min(nh, oh):
            raise ArithmeticError('C4 heat pressure and admitted absolute C1 pressure disagree')
    identities = pressure_integral_identities()
    fixture = independent_integral_fixture()
    hashes[NAME] = hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = dict(actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],
                  implicit_source_sha256=record['implicit_source_sha256'],
                  exact_pressure_source_identities=identities, independent_integral_fixture=fixture,
                  admitted_absolute_pressure_source_proof_recomputed=True,
                  defining_source_bridge=bridge,
                  actual_forward_tail_overlap_diagnostics=overlaps,
                  independent_admitted_C1_absolute_pressure_consistent=True,
                  absolute_pressure_same_source_mixed4_available=True,
                  candidate_pressure_tail_mixed4_available=True,
                  source_history_transfer_conditional=False,
                  heat_exterior_stress_identity_certified=False,
                  global_admissible_stress_lift_constructed=False, temporal_recursion=False,
                  candidate_checks_passed=True, all_passed=True, input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf8')
    print('Source-bound absolute Gamma pressure integral and mixed derivatives PASS', flush=True)
    return result


if __name__ == '__main__':
    run()
