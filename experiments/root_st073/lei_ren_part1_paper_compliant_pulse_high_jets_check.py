"""Independent physical radial derivatives and same-source C4 chart checks."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_pulse_high_jets import CompliantPulseHighJets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).parent
NAME = 'lei_ren_part1_paper_compliant_pulse_high_jets.json'


def source_identities():
    checks = {}
    def zero(name, expression):
        if s.simplify(expression) != 0: raise ArithmeticError('Pulse high-jet identity failed: '+name)
        checks[name] = True
    Z, b, m, mz, delta = s.symbols('Z b m mz delta', real=True)
    q, d, L = 1+Z*Z, 1-Z*Z, 1-delta*Z*Z
    env = dict(z=Z, b=b, m=m, mz=mz, q=q, d=d, **{'self.delta':delta})
    numerator = 2*Z*b-(1-delta)*Z*m-d*(mz-2*Z*m/q)
    zero('production_paper_3_9_numerator', assignment('compliant_pulse_high_jets', 'radial', 'numerator', env)-numerator)
    t, mu = s.symbols('t mu', real=True)
    B, M = s.Function('B')(t,Z), s.Function('M')(t,Z)
    N = numerator.subs({b:B, m:M, mz:s.diff(M,Z)})/L
    # Differentiating normalized Mz transport and the true physical factor
    # sqrt(R/2)*Utheta is equivalent to paper (3.8), not a new velocity fit.
    my = B-(s.Rational(1,2)-mu)*M
    Ny = s.diff(N,t).subs(s.diff(M,t,Z), s.diff(my,Z)).subs(s.diff(M,t),my)-mu*N
    thetaZ = -2*Z/q
    source_3_8 = ((1+delta)*Z*B-d*(s.diff(B,Z)+thetaZ*B)+2*Z*(s.diff(B,t)-(s.Rational(1,2)+mu)*B))/L-N/2
    zero('radial_log_derivative_from_paper_3_8', Ny-source_3_8)
    theta, A = s.Function('theta')(Z), s.Function('A')(Z)
    zero('physical_radial_Z_product_rule', s.diff(theta*A,Z)/theta-(s.diff(A,Z)+s.diff(theta,Z)/theta*A))
    # Arbitrary smooth histories can be differentiated at every chart
    # interface once the exact selected equations identify those functions.
    ap, e0, ev, end, K = s.symbols('ap e0 ev end K', real=True)
    chosen = (1-s.exp(-26))/4-mu*e0+mu*s.exp(-26)*(ev-end)
    forward = s.exp(22)*(e0+ap*ap*K/mu-(1-s.exp(-22))/(4*mu))
    backward = s.exp(-4)*ev+(1-s.exp(-4))/(4*mu)-s.exp(-4)*end
    zero('main_gap_energy_arbitrary_functions', (forward-backward).subs(ap*ap*K, chosen))
    a = 1+Z*Z
    coeff = [a**-2, -4*Z*a**-3, 2*(5*Z*Z-1)*a**-4,
             4*Z*(3-5*Z*Z)*a**-5, (35*Z**4-42*Z*Z+3)*a**-6]
    for n in range(5):
        zero('pressure_shape_order'+str(n), s.diff(a**-2,Z,n)/s.factorial(n)-coeff[n])
    return checks


def independent_physical_radial_fixture():
    with mp.workdps(85):
        c = MPIntervalContext(); c.dps = 75; radius = mp.mpf('1e-65')
        delta, mu = mp.mpf('.03'), mp.mpf('.04')
        ap = lambda z: mp.sqrt(1+z*z)+mp.mpf('.02')*z
        b = lambda z: ap(z)*mp.mpf('1.1')
        by = lambda z: ap(z)*mp.mpf('.6')
        m = lambda z: mp.mpf('.03')*(z+z**3)+mp.mpf('.3')*ap(z)
        theta = lambda z: 1/(1+z*z)
        my = lambda z: b(z)-(mp.mpf('.5')-mu)*m(z)
        # Direct integrated paper formula in physical factor units, with
        # an independently differentiated Mz=theta*m, not its normalized jet.
        ur = lambda z: (2*z*theta(z)*b(z)-(1-delta)*z*theta(z)*m(z)
            -(1-z*z)*mp.diff(lambda x: theta(x)*m(x),z))/(1-delta*z*z)
        ury = lambda z: ((1+delta)*z*theta(z)*b(z)
            -(1-z*z)*mp.diff(lambda x: theta(x)*b(x),z)
            +2*z*theta(z)*(by(z)-(mp.mpf('.5')+mu)*b(z)))/(1-delta*z*z)-ur(z)/2
        field = CompliantPulseHighJets.__new__(CompliantPulseHighJets)
        field.ctx=c; field.delta=c.mpf(delta); field.mu=c.mpf(mu)
        checks = {}
        for Z in (mp.mpf(-1),mp.mpf(0),mp.mpf('.3'),mp.mpf(1)):
            jet = lambda f: IntervalTaylor(c,[c.mpf([v-radius,v+radius]) for v in mp.taylor(f,Z,4)])
            u = jet(theta); B = jet(b); M = jet(m)
            field.data = lambda z: (None,None,u,None,None)
            point = dict(Z=c.mpf(Z), Uz_over_Utheta=B,
                Uz_y_over_Utheta=jet(by)-B*(c.mpf('.5')+field.mu),
                Mz_over_R_Utheta=M, Mtheta_over_sqrt2_R_3half_Utheta=c.mpf(1),
                Ur_over_sqrt_R_over_2_Utheta=field.radial(c.mpf(Z),B,M))
            point = field._high_packet(point)
            tests = [('physical_Ur',point['factored_Ur_over_sqrt_R_over_2_Pstar_Taylor'],ur),
                ('physical_Ur_y',point['Ur_y_over_sqrt_R_over_2_Utheta']*u.truncate(3),ury),
                ('physical_Ur_Z',point['Ur_Z_over_sqrt_R_over_2_Utheta']*u.truncate(2),lambda z:mp.diff(ur,z)),
                ('physical_Ur_yZ',point['Ur_yZ_over_sqrt_R_over_2_Utheta']*u.truncate(2),lambda z:mp.diff(ury,z))]
            for label, actual, function in tests:
                for n, value in enumerate(mp.taylor(function,Z,actual.order)):
                    lo, hi = endpoints(actual[n])
                    if not lo <= value <= hi: raise ArithmeticError('Independent radial derivative outside enclosure: '+label)
                    checks[label+'_Z'+str(Z)+'_order'+str(n)] = True
        return dict(checks=checks, actual_Md40_source_admission=False, passed=True)


def run():
    r = json.loads((HERE/NAME).read_bytes()); hashes = dict(r['input_hashes'])
    for name, digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Pulse C4 source changed: '+name)
    c = MPIntervalContext(); c.dps = 240; read = lambda v:read_interval(c,v)
    keys = ('Uz_over_Utheta','Uz_y_over_Utheta','Mz_over_R_Utheta',
        'Mtheta_z_over_sqrt2_R_3half_Utheta_squared','Mtheta_over_sqrt2_R_3half_Utheta',
        'Mztheta_over_R_Utheta_squared')
    radial_keys = ('Ur_over_sqrt_R_over_2_Utheta','Ur_y_over_sqrt_R_over_2_Utheta')
    groups = ('main_samples','entrance_samples','gap_samples','gap_end_samples','end_samples')
    samples = [p for name in groups for p in r[name]]+[r[name] for name in
        ('whole_Z_terminal','whole_Z_main','whole_Z_end_support','whole_Z_support_crossing')]
    old = json.loads((HERE/'lei_ren_part1_paper_compliant_axial_pulse_field.json').read_bytes())
    for point in samples:
        for key in keys:
            if len(point[key]['coefficients']) != 5: raise ArithmeticError('Partial primitive C4 lost: '+key)
        for key in radial_keys:
            if len(point[key]['coefficients']) != 4: raise ArithmeticError('Physical radial C3 lost: '+key)
        for key in ('Ur_Z_over_sqrt_R_over_2_Utheta','Ur_yZ_over_sqrt_R_over_2_Utheta'):
            if len(point[key]['coefficients']) != 3: raise ArithmeticError('Radial derivative order lost')
        pressure = point['pressure']
        if not pressure['exact_pressure_source_unchanged'] or not point['divergence_preserved_by_Mz_recovery']:
            raise ValueError('Source restoration/radial recovery lost')
        for key in ('Mp_over_Pstar_squared','P0_over_Pstar_squared','P_over_Pstar_squared','P_y_over_Pstar_squared'):
            if len(pressure[key]['coefficients']) != 5: raise ValueError('Pressure C4 lost')
        for mv,p0v,pv in zip(pressure['Mp_over_Pstar_squared']['coefficients'],pressure['P0_over_Pstar_squared']['coefficients'],pressure['P_over_Pstar_squared']['coefficients']):
            expected=read(mv)+read(p0v); actual=read(pv)
            if not endpoints(actual)[0] <= endpoints(expected)[0] <= endpoints(expected)[1] <= endpoints(actual)[1]:
                raise ArithmeticError('Absolute pressure no longer retains P0+Mp')
        if not point['Ur_Z_available']: raise ValueError('Radial derivative omitted')
        for key in ('radial_velocity_C4_available','full_pulse_C4_installed','full_outer_C4_certified','whole_outer_cone_certified','temporal_recursion','fifth_derivative_Taylor_remainder_available'):
            if point[key]: raise ValueError('Pulse high-jet scope overclaimed: '+key)
    lower_checks = 0
    for name in groups:
        for actual, prior in zip(r[name],old[name]):
            for key in keys:
                if key == 'Mtheta_over_sqrt2_R_3half_Utheta': continue
                for n in (0,1):
                    a,b = endpoints(read(actual[key]['coefficients'][n])); d,e = endpoints(read(prior[key]['coefficients'][n]))
                    if max(a,d) > min(b,e): raise ArithmeticError('Higher pulse jets changed the source: '+key)
                    lower_checks += 1
    whole = r['whole_Z_terminal']
    for key in ('Uz_over_Utheta','Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared')+radial_keys+('Ur_Z_over_sqrt_R_over_2_Utheta','Ur_yZ_over_sqrt_R_over_2_Utheta'):
        for value in whole[key]['coefficients']:
            if endpoints(read(value)) != (mp.mpf(0),mp.mpf(0)):
                raise ArithmeticError('Exact whole-Z terminal derivative closure lost: '+key)
    if endpoints(read(whole['Mztheta_over_R_Utheta_squared']['coefficients'][0]))[0] <= 0:
        raise ArithmeticError('Positive future-energy target lost')
    for actual, expected in zip(whole['Mztheta_over_R_Utheta_squared']['coefficients'],whole['positive_terminal_future_energy_Taylor']['coefficients']):
        a,b=endpoints(read(actual)); d,e=endpoints(read(expected))
        if not a <= d <= e <= b: raise ArithmeticError('Terminal future-energy derivative identity lost')
    interfaces = {}
    for label,left,right in (('main_gap_xi11',r['main_samples'][-1],r['gap_samples'][0]),
        ('gap_coordinate_overlap',r['gap_samples'][2],r['gap_end_samples'][0]),
        ('gap_end_s_minus4',r['gap_end_samples'][-1],r['end_samples'][0])):
        for key in keys+radial_keys+('Ur_Z_over_sqrt_R_over_2_Utheta','Ur_yZ_over_sqrt_R_over_2_Utheta'):
            for av,bv in zip(left[key]['coefficients'],right[key]['coefficients']):
                a,b=endpoints(read(av)); d,e=endpoints(read(bv))
                if max(a,d)>min(b,e): raise ArithmeticError('High-jet interface disjoint: '+label+' '+key)
        for key in ('Mp_over_Pstar_squared','P0_over_Pstar_squared','P_over_Pstar_squared','P_y_over_Pstar_squared'):
            for av,bv in zip(left['pressure'][key]['coefficients'],right['pressure'][key]['coefficients']):
                a,b=endpoints(read(av)); d,e=endpoints(read(bv))
                if max(a,d)>min(b,e): raise ArithmeticError('Pressure C4 interface disjoint: '+label+' '+key)
        interfaces[label]=True
    proof=source_identities(); fixture=independent_physical_radial_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],
        implicit_source_sha256=r['implicit_source_sha256'],source_and_physical_derivative_identities=proof,
        independent_physical_radial_fixture=fixture,directed_high_jet_interfaces=interfaces,
        unchanged_C1_partial_field_comparisons=lower_checks,whole_Z_terminal_derivatives_checked=True,
        pulse_primitives_axial_C4_installed=True,pulse_radial_velocity_axial_C3_installed=True,
        radial_velocity_Z_and_yZ_available=True,all_passed=True,
        radial_velocity_C4_available=False,full_pulse_C4_installed=False,
        full_outer_C4_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Pulse primitive C4/radial C3: source identities, independent physical derivatives, three chart overlaps and whole-Z terminal gates PASS',flush=True)
    return result


if __name__ == '__main__': run()
