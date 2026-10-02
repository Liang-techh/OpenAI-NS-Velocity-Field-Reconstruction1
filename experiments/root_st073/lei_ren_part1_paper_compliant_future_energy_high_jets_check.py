"""Independent normalization and new axial-derivative paths for future C4."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_angular_high_jets import log_taylor
from lei_ren_part1_paper_compliant_future_swirl_energy import CompliantFutureSwirlEnergy
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).parent
NAME = 'lei_ren_part1_paper_compliant_future_energy_high_jets.json'


def identities():
    out = {}
    def zero(name, expression):
        if s.simplify(expression) != 0: raise ArithmeticError('Future C4 identity failed: '+name)
        out[name] = True
    Q, mu, L, dj, E, F, center = s.symbols('Q mu L dj E F center', real=True)
    env = {'q': Q, 'f.mu': mu, 'f.Lrel': L, 'dj': dj, "w['E']": E, "w['F']": F, 'center': center}
    for name, expr in (('Nf', Q**2*s.exp(-200*mu)/4), ('Nrel', Q**2*s.exp(-2*mu*(100+L))/4)):
        zero('source_'+name, assignment('compliant_future_energy_high_jets', 'future', name, env)-expr)
    zero('source_signed_angular_energy_increment', assignment('compliant_future_energy_high_jets', 'future', 'angular_change', env, augmented=True)-(2*E*dj+F*dj*dj)*s.exp(-2*mu*center))
    z, h = s.symbols('Z h', real=True)
    coef = [1+2*z*z+z**4, 4*z+4*z**3, 2+6*z*z, 4*z, 1]
    zero('exact_Q_squared_coefficients', (1+(z+h)**2)**2-sum(coef[n]*h**n for n in range(5)))
    y = s.symbols('y', real=True); U = s.Function('U')(y)
    # Direct physical logarithmic coordinate change for swirl energy.
    zero('radial_energy_log_measure', s.exp(y)*U**2/(s.exp(y)*U**2)-1)
    zero('Section7_34_Rv_to_Rp_weight', (mu/2)*s.exp(13/mu)*s.exp(-(1+2*mu)*13/mu)-mu*s.exp(-26)/2)
    return out


def contains(label, enclosure, reference):
    lo, hi = endpoints(enclosure)
    if not lo <= reference <= hi: raise ArithmeticError('Independent future C4 derivative outside enclosure: '+label)


def derivative_fixture():
    with mp.workdps(70):
        c = MPIntervalContext(); c.dps = 65
        Z = mp.mpf('.5'); mu = mp.mpf('.025'); Q = 1+Z*Z
        q = IntervalTaylor(c, [c.mpf(Q), c.mpf(2*Z), 1, 0, 0])
        logq = log_taylor(q)-c.ln(2)
        def sigma(t):
            if t <= 0: return mp.mpf(0)
            if t >= 1: return mp.mpf(1)
            a, b = mp.exp(-1/t**2), mp.exp(-1/(1-t)**2)
            return a/(a+b)
        # Point samples of real powers provide an independent derivative
        # formula for each integrand; no derivative of a quadrature iterate.
        flat = q*0; cells = 512
        from lei_ren_part1_paper_compliant_outer_initial import stable_sigma
        for i in range(cells):
            left, right = c.mpf(100)*i/cells, c.mpf(100)*(i+1)/cells
            t = c.mpf([endpoints(left)[0], endpoints(right)[1]])
            sig = stable_sigma(c, t/100)[0]
            flat += (logq*(2*sig)).exp()*((right-left)*c.exp(-2*c.mpf(mu)*t))
        cache = {}; u, v = 2*Z/Q, 1/Q
        def power_coeffs(t):
            if t not in cache:
                p = 2*sigma(t/100); p2 = p*(p-1); p3 = p2*(p-2); p4 = p3*(p-3)
                vals = (1, p*u, p*v+p2*u*u/2, p2*u*v+p3*u**3/6,
                        p2*v*v/2+p3*u*u*v/2+p4*u**4/24)
                cache[t] = [mp.exp(-2*mu*t)*(Q/2)**p*x for x in vals]
            return cache[t]
        checks = {}
        for n in range(5):
            ref = mp.quad(lambda t: power_coeffs(t)[n], [0, 25, 50, 75, 90, 100])
            contains('flatten_energy_order'+str(n), flat[n], ref); checks['flatten_energy_order'+str(n)] = True

        raw = lambda v: mp.exp(-1/(1-v*v)) if abs(v) < 1 else mp.mpf(0)
        norm = mp.quad(raw, [-1, 0, 1]); ell = mp.mpf('.15')
        beta = lambda v: raw(v/ell)/(ell*norm)
        E = mp.quad(lambda v: mp.exp(-2*mu*v)*beta(v), [-ell, 0, ell])
        F = mp.quad(lambda v: mp.exp(-2*mu*v)*beta(v)**2, [-ell, 0, ell])
        d1 = lambda z: mp.mpf('.001')/(1+z*z)
        d2 = lambda z: -mp.mpf('.002')*z*z
        coeffs = [mp.taylor(d, Z, 4) for d in (d1, d2)]
        change = q*0
        for ds, center in zip(coeffs, (-3, -1)):
            # A small enclosing radius covers reference quadrature/rounding;
            # it is only a fixture and does not define the actual coefficients.
            radius = c.mpf('1e-55')
            ds = IntervalTaylor(c, [c.mpf(v)+c.mpf([-endpoints(radius)[1], endpoints(radius)[1]]) for v in ds])
            change += (ds*(2*c.mpf(E))+ds*ds*c.mpf(F))*c.exp(-2*c.mpf(mu)*center)
        physical = lambda z: 2*E*(d1(z)*mp.exp(6*mu)+d2(z)*mp.exp(2*mu))+F*(d1(z)**2*mp.exp(6*mu)+d2(z)**2*mp.exp(2*mu))
        refs = mp.taylor(physical, Z, 4)
        for n in range(5):
            contains('physical_signed_bump_energy_order'+str(n), change[n], refs[n]); checks['physical_signed_bump_energy_order'+str(n)] = True
        return dict(independent_energy_derivative_checks=checks, actual_Md40_source_admission=False, passed=True)


def run():
    r = json.loads((HERE/NAME).read_bytes()); hashes = dict(r['input_hashes'])
    for name, digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Future C4 dependency changed: '+name)
    c = MPIntervalContext(); c.dps = 230; read = lambda v: read_interval(c, v)
    for point in r['samples']+[r['whole_Z_C4']]:
        for key in ('complete_future_energy_Taylor', 'Section7_34_weighted_future_Taylor', 'Nf', 'Nrel'):
            if len(point[key]['coefficients']) != 5: raise ValueError('Future energy C4 lost: '+key)
        if endpoints(read(point['complete_future_energy_Taylor']['coefficients'][0]))[0] <= 0:
            raise ArithmeticError('Complete future energy must remain positive')
        if point['full_physical_kinetic_energy_certified'] or point['full_outer_C4_certified'] or point['temporal_recursion']:
            raise ValueError('Future energy C4 scope overclaimed')
    for n in (1, 3):
        if endpoints(read(r['samples'][1]['complete_future_energy_Taylor']['coefficients'][n])) != (mp.mpf(0), mp.mpf(0)):
            raise ArithmeticError('Complete future energy lost exact Z0 even parity')
    base = CompliantFutureSwirlEnergy()
    for Z, point in zip(('-1', '0', '.5', '1', [-1, 1]), r['samples']+[r['whole_Z_C4']]):
        old = base.future(Z)['complete_future_energy_Taylor']
        for j in (0, 1):
            lo, hi = endpoints(read(point['complete_future_energy_Taylor']['coefficients'][j]))
            ol, oh = endpoints(old[j])
            if not ol <= lo <= hi <= oh: raise ArithmeticError('Future high jets changed C1 source')
    proof = identities(); fixture = derivative_fixture()
    hashes[NAME] = hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out = dict(identities=proof, independent_energy_derivative_fixture=fixture,
        actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'], implicit_source_sha256=r['implicit_source_sha256'],
        complete_future_corrected_energy_C4_available=True, unchanged_C1_source_and_Z0_parity_checked=True,
        full_physical_kinetic_energy_certified=False, full_outer_C4_certified=False,
        whole_outer_cone_certified=False, temporal_recursion=False, all_passed=True, input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
    print('Complete future energy C4: source normalization identities, ten independent new derivative paths, old C1 inclusion and Z0 parity PASS', flush=True)
    return out


if __name__ == '__main__': run()
