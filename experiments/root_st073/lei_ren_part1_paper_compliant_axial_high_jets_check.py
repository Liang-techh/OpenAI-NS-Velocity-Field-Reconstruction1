"""Source shape propagation and independent positive-root C4 checks."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_axial_high_jets import positive_quadratic_jets
from lei_ren_part1_paper_compliant_outer_pulse_map import SharedOuterPulseMap
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).parent
NAME = 'lei_ren_part1_paper_compliant_axial_high_jets.json'


def source_identities():
    checks = {}
    def zero(name, expression):
        if s.simplify(expression) != 0: raise ArithmeticError('Incoming source identity failed: '+name)
        checks[name] = True
    z, h, U, EZ, EQ, t, invP2, KB2, energy = s.symbols('z h U EZ EQ t invP2 KB2 energy', real=True)
    raw = EZ*z*z+EQ/(1+z*z)**2
    env = {'z': z, 'self.invP2': invP2, 'qi': 1/(1+z*z), 'yy': t,
           'integrals[2]': energy}
    zero('slope_separated_positive_axial_energy', assignment('compliant_outer_initial', 'slope', 'e', env)-
        (16*invP2*z*z-(s.Rational(5,12)+energy/2)*s.exp(-t)/(1+z*z)**2))
    env = {'z': z, "get('Mztheta_over_R_Pstar_squared')": raw,
           'decay': s.exp(-t), 'self.invP2': invP2, "K['B_squared_mass']": KB2,
           'u1': U/(1+z*z), 't': t}
    zero('axial_separated_energy_transport', assignment('compliant_outer_initial', 'axial', 'e', env)-
        ((EZ*s.exp(-t)+16*invP2*KB2)*z*z+(EQ-U*U*t/2)*s.exp(-t)/(1+z*z)**2))
    env.update({"kernels['energy']": energy, 'decay_integral(c, 2 * mu, t)': energy})
    for method in ('slope_mu', 'power'):
        zero(method+'_separated_energy_transport', assignment('compliant_outer_buffer', method, 'e', env)-
            (EZ*s.exp(-t)*z*z+(EQ-U*U*energy/2)*s.exp(-t)/(1+z*z)**2))
    M, K, invP = s.symbols('M K invP', real=True)
    zero('incoming_m1_shape', M*z/(U/(1+z*z))*invP-invP*M/U*(z+z**3))
    zero('incoming_m2_shape', K*z/(1+z*z)/(U/(1+z*z))**2*invP-invP*K/U**2*(z+z**3))
    zero('incoming_energy_shape_without_extra_invP', raw/(U/(1+z*z))**2-(EQ/U**2+EZ/U**2*(z*z+2*z**4+z**6)))
    g = [z+z**3, 1+3*z*z, 3*z, 1, 0]
    hp = [z*z+2*z**4+z**6, 2*z+8*z**3+6*z**5, 1+12*z*z+15*z**4, 8*z+20*z**3, 2+15*z*z]
    zero('ordinary_moment_Taylor_coefficients', s.expand((z+h)+(z+h)**3)-sum(g[n]*h**n for n in range(5)))
    zero('ordinary_energy_Taylor_coefficients_through4', s.Poly(s.expand((z+h)**2+2*(z+h)**4+(z+h)**6-sum(hp[n]*h**n for n in range(5))), h).as_expr().series(h, 0, 5).removeO())
    return checks


def independent_root_fixture():
    with mp.workdps(85):
        c = MPIntervalContext(); c.dps = 100; radius = mp.mpf('1e-75'); Z = mp.mpf('.3')
        jet = lambda f: IntervalTaylor(c, [c.mpf([v-radius,v+radius]) for v in mp.taylor(f, Z, 4)])
        u = [lambda z: z+z**3, lambda z: -.2*z+.1*z*z]
        v = [mp.mpf('-.7'), mp.mpf('.9')]; nu = [mp.mpf('.012'), mp.mpf('.009')]
        K = mp.mpf('1.3'); B = lambda z: 2+mp.exp(-z*z)+z**4/10
        A2 = K+sum(nu[j]*v[j]**2 for j in range(2))
        A1f = lambda z: 2*sum(nu[j]*u[j](z)*v[j] for j in range(2))
        A0f = lambda z: sum(nu[j]*u[j](z)**2 for j in range(2))-B(z)
        root = lambda z: (-A1f(z)+mp.sqrt(A1f(z)**2-4*A2*A0f(z)))/(2*A2)
        a0 = root(Z); ap, D = positive_quadratic_jets(c, c.mpf([a0-radius,a0+radius]), c.mpf(A2), jet(A1f), jet(A0f))
        checks = {}
        for n, value in enumerate(mp.taylor(root, Z, 4)):
            lo, hi = endpoints(ap[n])
            if not lo <= value <= hi: raise ArithmeticError('Independent scalar root C4 failed: '+str(n))
            checks['direct_positive_root_order'+str(n)] = True
        for j in range(2):
            end = jet(u[j])+ap*c.mpf(v[j])
            for n, value in enumerate(mp.taylor(lambda z: u[j](z)+v[j]*root(z), Z, 4)):
                lo, hi = endpoints(end[n])
                if not lo <= value <= hi: raise ArithmeticError('Independent end coefficient C4 failed')
                checks['end'+str(j+1)+'_order'+str(n)] = True
        if endpoints(D)[0] <= 0: raise ArithmeticError('Fixture positive inverse failed')
        return dict(checks=checks, actual_Md40_source_admission=False, passed=True)


def run():
    r = json.loads((HERE/NAME).read_bytes()); hashes = dict(r['input_hashes'])
    for name, digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Axial C4 dependency changed: '+name)
    c = MPIntervalContext(); c.dps = 230; read = lambda v: read_interval(c, v)
    k = r['incoming_constants']
    for name in ('U', 'M', 'K', 'E_Z', 'C1', 'C2', 'C_E'):
        if endpoints(read(k[name]))[0] <= 0: raise ArithmeticError('Positive incoming constant lost: '+name)
    for name in ('E_Q', 'C0'):
        if endpoints(read(k[name]))[1] >= 0: raise ArithmeticError('Negative swirl constant lost: '+name)
    pulse = SharedOuterPulseMap()
    with mp.workdps(210):
        for Z, point in zip(('-1', '0', '.5', '1', [-1,1]), r['samples']+[r['whole_Z_C4']]):
            if endpoints(read(point['positive_root_derivative_denominator']))[0] <= 0:
                raise ArithmeticError('Actual selected inverse failed')
            trial = pulse.coefficients(Z, 1)
            for source, actual in zip(point['incoming']['moment_Taylor'], trial['actual_incoming_moment_Taylor_enclosures']):
                for n in (0,1):
                    lo, hi = endpoints(read(source['coefficients'][n])); al, ah = endpoints(actual[n])
                    if hi < al or ah < lo: raise ArithmeticError('Actual incoming source C1 disagrees')
            for n in (0,1):
                lo, hi = endpoints(read(point['incoming']['energy_Taylor']['coefficients'][n]))
                al, ah = endpoints(trial['actual_normalized_incoming_energy'][n])
                if hi < al or ah < lo: raise ArithmeticError('Actual incoming energy C1 disagrees')
            for jet in [point['selected_ap_Taylor']]+point['selected_scaled_end_coefficient_Taylor']:
                if len(jet['coefficients']) != 5: raise ValueError('Actual selected C4 lost')
            controls = point['selected_scaled_end_coefficient_Taylor']
            if not endpoints(read(controls[0]['coefficients'][0]))[1] < 0 < endpoints(read(controls[1]['coefficients'][0]))[0]:
                raise ArithmeticError('Actual selected coefficient signs lost')
            for key in ('full_pulse_C4_installed', 'full_outer_C4_certified', 'whole_outer_cone_certified', 'temporal_recursion'):
                if point[key]: raise ValueError('Axial C4 scope overclaimed')
    proof = source_identities(); fixture = independent_root_fixture()
    hashes[NAME] = hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out = dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],
        implicit_source_sha256=r['implicit_source_sha256'], incoming_function_source_identities=proof,
        independent_root_C4_fixture=fixture, actual_incoming_C1_source_agreement_checked=True,
        actual_selected_ap_c1_c2_C4_available=True, all_passed=True,
        full_pulse_C4_installed=False, full_outer_C4_certified=False,
        whole_outer_cone_certified=False, temporal_recursion=False, input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
    print('Actual incoming shapes, independent root/end C4 fixture and whole-Z positive inverse/sign gates PASS', flush=True)
    return out


if __name__ == '__main__': run()
